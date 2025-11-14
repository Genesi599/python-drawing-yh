import gseapy as gp

# 获取所有可用的基因集
available_gene_sets = gp.get_library_name()
print(available_gene_sets)