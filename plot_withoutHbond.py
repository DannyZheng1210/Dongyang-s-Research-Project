

import pandas as pd
import matplotlib.pyplot as plt

# 读取数据
file_path = r"E:\new_HPMCAS\script\10ASD_new_RMSE.csv"
df = pd.read_csv(file_path)

# ========== 参数设置 ==========
percent_to_keep = 100   # 比如 10 表示取最小10%的点，100 表示全部

# ========== 筛选数据 ==========
df = df.sort_values(by="Distance_to_origin", ascending=True).reset_index(drop=True)
n_keep = int(len(df) * percent_to_keep / 100)
if n_keep < 1:
    n_keep = 1
df = df.iloc[:n_keep]

# ========== 绘图 ==========
plt.figure(figsize=(8, 6))
point_size = 5

plt.scatter(
    df["API_RMSE_all"],
    df["HPMCAS_RMSE_all"],
    c="gray",
    alpha=0.7,
    s=point_size,
    edgecolors="none"
)

plt.title("10%ASD", fontsize=18, fontweight="bold")

# 坐标轴
plt.xlabel("RMSE $_{API}$", fontsize=18, fontweight="bold")
plt.ylabel("RMSE $_{HPMCAS}$", fontsize=18, fontweight="bold")

plt.grid(True)
plt.tight_layout()
plt.show()
