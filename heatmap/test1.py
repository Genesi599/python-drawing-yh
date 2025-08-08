import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pandas import Series,DataFrame
import seaborn as sns
import palettable

plt.rcParams['font.sans-serif']=['SimHei']  # 用于显示中文
plt.rcParams['axes.unicode_minus'] = False  # 用于显示中文

from sklearn import datasets
iris=datasets.load_iris()
x, y = iris.data, iris.target
pd_iris = pd.DataFrame(np.hstack((x, y.reshape(150, 1))),columns=['sepal length(cm)','sepal width(cm)','petal length(cm)','petal width(cm)','class'] )

plt.figure(dpi=200)
sns.clustermap(data=pd_iris,#仅仅需传入绘图数据集
              )