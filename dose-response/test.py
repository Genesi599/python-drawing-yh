# 指定文件夹路径
import os

folder_path = "data/raw"

# 获取文件夹中的所有文件和文件夹名称
entries = os.listdir(folder_path)

# 打印结果
for entry in entries:
    print(entry)