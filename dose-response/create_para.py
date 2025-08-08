import pandas as pd


def gradient_dilution(start_concentration, dilution_factor, num_gradients):
    """
    计算梯度稀释的浓度列表。

    参数:
        start_concentration (float): 起始浓度。
        dilution_factor (float): 稀释倍数。
        num_gradients (int): 总的梯度数量。

    返回:
        list: 包含浓度梯度的列表。
    """
    if dilution_factor <= 1:
        raise ValueError("稀释倍数必须大于1")
    if num_gradients <= 0:
        raise ValueError("梯度数量必须是正整数")

    # 初始化浓度列表
    concentrations = []
    current_concentration = start_concentration

    # 计算每个梯度的浓度
    for _ in range(num_gradients):
        concentrations.append(current_concentration)
        current_concentration /= dilution_factor

    return concentrations


def create_cell_df(cell_name):
    # Initialize the cell DataFrame with 'MOLM-13'
    cell_df = pd.DataFrame([[cell_name] * 10] * 6)

    # Set the last column to null (empty values)
    cell_df.iloc[:, 0] = None

    # Convert each cell value into a dictionary with key 'cell'
    cell_df = cell_df.map(lambda x: {'cell': x})
    return cell_df


def create_date_df(date):
    # Define the dictionary value
    date_dict = {'date': date}

    # Use a nested list comprehension to create a 6x10 grid of the dictionary
    data = [[date_dict for _ in range(10)] for _ in range(6)]

    # Create the DataFrame
    date_df = pd.DataFrame(data)

    return date_df


def create_drug_df(drugnames, file_order):
    # 药物-浓度 (nM)
    drug_name_1 = drugnames[0]
    concentration1 = None
    # Create a 6x10 DataFrame
    drug_df = pd.DataFrame([[{drug_name_1: concentration1} for _ in range(10)] for _ in range(6)])
    # Define the new drug to add
    drug_name_2 = drugnames[1]
    concentration_2 = 0
    gradient = gradient_dilution(256000, 4, 8)
    for row_index in range(len(drug_df)):  # For each row
        for col_index in range(1, drug_df.shape[1] - 1):  # For all columns except the last two
            # Calculate the dilution for 'dBET1'
            drug_df.iloc[row_index, col_index][drug_name_1] = gradient[col_index-1]
            # Set the second-to-last column ('dBET1') concentration to 0
            drug_df.iloc[row_index, 0][drug_name_1] = 0
            drug_df.iloc[row_index, -1][drug_name_1] = 0


    # Function to add the new drug to the dictionary
    def add_new_drug(cell):
        cell[drug_name_2] = concentration_2  # Add a new key-value pair
        return cell
    # Apply the function to every cell in the DataFrame
    drug_df = drug_df.map(add_new_drug)

    # Update concentrations for CC-90009
    gradient = gradient_dilution(256000, 4, 7)
    gradient = gradient + [0]
    for col_index in range(drug_df.shape[1]):  # Iterate through each column

        for row_index in range(3):  # Iterate through each row except the last
            drug_df.iloc[row_index, col_index][drug_name_2] = gradient[2*file_order]
        for row_index in range(3,6):  # Iterate through each row except the last
            drug_df.iloc[row_index, col_index][drug_name_2] = gradient[2*file_order+1]
    return drug_df


def create_para_df(cell_name, drugnames, file_order, date):
    cell_df = create_cell_df(cell_name)
    cell_df.to_csv(f'data/parameters/cell.csv', index=False)
    date_df = create_date_df(date)
    date_df.to_csv(f'data/parameters/date.csv', index=False)
    drug_df = create_drug_df(drugnames, file_order)
    drug_df.to_csv(f'data/parameters/drug.csv', index=False)
    return cell_df, date_df, drug_df

