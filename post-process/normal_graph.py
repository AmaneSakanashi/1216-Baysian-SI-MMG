import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import norm

# plot params
plt.rcParams["font.family"] = "Times New Roman"       # 使用するフォント
# plt.rcParams["font.serif"] = "Times New Roman"
plt.rcParams["font.size"] = 16              # 基本となるフォントの大きさ
plt.rcParams["mathtext.cal"] = "serif"      # TeX表記に関するフォント設定
plt.rcParams["mathtext.rm"] = "serif"       # TeX表記に関するフォント設定
plt.rcParams["mathtext.it"] = "serif:italic"# TeX表記に関するフォント設定
plt.rcParams["mathtext.bf"] = "serif:bold"  # TeX表記に関するフォント設定
plt.rcParams["mathtext.fontset"] = "cm"     # TeX表記に関するフォント設定

mean_result = pd.read_csv("opt_result/mean_result.csv", header=None)
var_result = pd.read_csv("opt_result/var_result.csv", header=None)

m_params = np.array(mean_result.values.flatten())
v_params = np.array(var_result.values.flatten())  

fig = plt.figure( figsize = (12,24))

x = np.arange(-10,10,0.1)
for n in range(len(m_params)):
    y = norm.pdf( x,m_params[n], v_params[n] )
    plt.subplot(len(m_params),1,n+1)
    plt.title(n)
    plt.plot(x,y)

plt.tight_layout()
savepath = "output/"
fig.savefig( savepath + "params_normalgraph.pdf")

plt.show()


        