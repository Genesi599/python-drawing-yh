
from create_para import *


def trim_outer_border(dataframe):
    # Remove the first and last rows, and the first and last columns
    trimmed_df = dataframe.iloc[1:-1, 1:-1]
    return trimmed_df

def read_file(file_path, para_df, drugnames):
    # Open and process the text file
    data = []
    current_iteration = None
    with open(file_path, "r") as file:
        for line in file:
            # Strip whitespace from the line
            line = line.strip()
            # Extract metadata (Plate and Wavelength)
            # Detect new iteration
            if line.startswith("Iteration:"):
                current_iteration = int(line.split(":")[1].strip())
                # Parse numerical grid data
            elif line and current_iteration is not None:
                # Split the line into individual values
                values = list(map(float, line.split()))
                # Append the iteration and parsed values as a row
                data.append([current_iteration] + values)

            # Create a DataFrame with the parsed data
    columns = ["Iteration"] + [f"Col{i + 1}" for i in range(len(data[0]) - 1)]  # Dynamic column names
    df = pd.DataFrame(data, columns=columns)
    # 1. Split DataFrames based on "Iteration"
    df1 = df[df["Iteration"] == 1].reset_index(drop=True)  # DataFrame for Iteration 1
    df2 = df[df["Iteration"] == 2].reset_index(drop=True)  # DataFrame for Iteration 2
    # 2. Drop the "Iteration" column for numerical operations
    df1_vals = df1.drop(columns=["Iteration"])
    df2_vals = df2.drop(columns=["Iteration"])
    # 3. Calculate cell-wise average
    df_mean = (df1_vals + df2_vals) / 2
    # 4. Add "Iteration" information to the averaged DataFrame if needed
    df_final = df_mean.copy()
    df_final = trim_outer_border(df_final)
    df_final.columns = range(df_final.shape[1])




    def replace_with_photometry_dict(df):
        # Replace each value in the DataFrame with a dictionary
        # where the key is 'photometry' and the value is the original value
        df = df.map(lambda x: {'photometry': x})
        return df

    # Replace values in the DataFrame with dictionaries
    df_final = replace_with_photometry_dict(df_final)
    df_final = df_final.reset_index(drop=True)
    df_final.columns = range(len(df_final.columns))


    def unify_indexes_and_columns(df):
        # 完整重置
        df = df.reset_index(drop=True)
        df.columns = range(len(df.columns))
        return df


    dfs = [df_final] + list(para_df)


    merged_list = []
    merged = {}
    for i in range(len(dfs[0])):
        for j in range(len(dfs[0].columns)):
            for df in dfs:
                merged = {**merged, **df.iloc[i,j]}
            merged_list.append(merged)

    return merged_list

def get_viablity(merged_list, cell_name, drug_names, date):
    new_list = []
    #  计算背景值并减去
    backround = []
    for i in merged_list:
        if i['cell'] is None:
            backround.append(i['photometry'])
    # 计算总和
    total_sum = sum(backround)
    # 计算元素数量
    num_elements = len(backround)
    # 计算平均值
    average = total_sum / num_elements

    for i in range(len(merged_list)):
        if merged_list[i]['cell'] is not None:
            merged_list[i]['photometry'] -= average
            new_list.append(merged_list[i])
    full_live = []
    for i in new_list:
        if i[drug_names[0]] == 0 and i[drug_names[1]] == 0:
            full_live.append(i['photometry'])
    # 计算总和
    total_sum = sum(full_live)
    # 计算元素数量
    num_elements = len(full_live)
    # 计算平均值
    full_live = total_sum / num_elements

    for i in range(len(new_list)):
        new_list[i]['viability'] = new_list[i]['photometry'] / full_live
    df = pd.DataFrame(new_list)
    df.to_csv(f'data/viability/{cell_name}-{drug_names[0]}+{drug_names[1]}-{date}.csv', index=False)






