import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

HARTREE_TO_KJMOL = 2625.5  # 单位换算常数

def calc_cv(series):
    """计算变异系数 (CV, %)"""
    mu = series.mean()
    sigma = series.std()
    return abs(sigma / mu) * 100

def calc_sem(series):
    """计算均值的标准误差 (SEM)"""
    sigma = series.std()
    n = series.count()
    return sigma / np.sqrt(n) if n > 0 else np.nan

def load_and_clean(filepath):
    """读取 + 单位换算 + 数据清洗"""
    df = pd.read_csv(filepath)
    df["E_HPMCAS"] = df["E_HPMCAS"] * HARTREE_TO_KJMOL
    df["E_Para"]   = df["E_Para"] * HARTREE_TO_KJMOL

    # 数据清洗：只保留三个能量都 < 0 的点
    df = df[(df["Ebinding(kJ/mol)"] < 0) &
            (df["E_Para"] < 0) &
            (df["E_HPMCAS"] < 0)]
    return df

def plot_scatter(ax, filepath, title):
    """原来的散点图"""
    df = pd.read_csv(filepath)
    x = df["Distance_to_origin"]
    y = df["Ebinding(kJ/mol)"]

    mask_neg = y < 0
    mask_pos = y >= 0

    total = len(y)
    pct_neg = 100 * mask_neg.sum() / total
    pct_pos = 100 * mask_pos.sum() / total

    ax.scatter(x[mask_neg], y[mask_neg], c="blue", s=5,
               label=f"Ebinding < 0 ({pct_neg:.1f}%)")
    ax.scatter(x[mask_pos], y[mask_pos], c="red", s=5,
               label=f"Ebinding ≥ 0 ({pct_pos:.1f}%)")

    ax.set_title(title, fontsize=18, fontweight="bold")
    ax.set_xlabel(r"$d$", fontsize=18, fontweight="bold")
    ax.set_ylabel(r"$E_{binding}$ (kJ/mol)", fontsize=18, fontweight="bold")

    ax.set_ylim(-100, 500)
    ax.legend(fontsize=12, markerscale=2)

# ========== 第一幅图：原始散点 ==========
fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharey=True)
plot_scatter(axes[0], r"E:\new_HPMCAS\script\10ASD_new_RMSE.csv", "10%ASD")
plot_scatter(axes[1], r"E:\new_HPMCAS\script\20ASD_new_RMSE.csv", "20%ASD")
axes[1].set_ylabel("")
plt.suptitle("Scatter plots", fontsize=20, fontweight="bold")
plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.show()

# ========== 第二幅图：平均值 + SEM ==========
# df10 = load_and_clean(r"E:\new_HPMCAS\script\10ASD_new_RMSE.csv")  # 先不画10%
df20 = load_and_clean(r"E:\new_HPMCAS\script\20ASD_new_RMSE.csv")

labels = ["E_Para", "E_HPMCAS"]

# 只计算20%ASD
means20 = [df20["E_Para"].mean(), df20["E_HPMCAS"].mean()]
sems20  = [calc_sem(df20["E_Para"]), calc_sem(df20["E_HPMCAS"])]

fig, ax = plt.subplots(figsize=(6, 6))  # 只要一个子图

bars = ax.bar(labels, means20, yerr=sems20, capsize=8,
              color=["skyblue", "lightcoral"], alpha=0.8)
ax.set_title("20%ASD", fontsize=18, fontweight="bold")
ax.set_ylabel("Energy (kJ/mol)", fontsize=14, fontweight="bold")
ax.grid(axis="y", linestyle="--", alpha=0.6)

plt.suptitle("Average energies (after filtering) ± SEM", fontsize=20, fontweight="bold")
plt.tight_layout(rect=[0, 0, 1, 0.95])
plt.show()