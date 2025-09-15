import pandas as pd
import matplotlib.pyplot as plt

# 读取数据
file_path = r"E:\new_HPMCAS\script\10ASD_new_RMSE.csv"
df = pd.read_csv(file_path)

# ========== 先筛选 E_Binding < 0 ==========
df = df[df["E_Binding(eV)"] < 0].copy()

# ========== 参数设置 ==========
percent_to_keep = 0.1   # ⚠️ 修改这里，比如10表示取最小10%的点，100表示全部

# ========== 筛选数据 ==========
df = df.sort_values(by="Distance_to_origin", ascending=True).reset_index(drop=True)
n_keep = int(len(df) * percent_to_keep / 100)
if n_keep < 1:
    n_keep = 1
df = df.iloc[:n_keep]

# ========== 归类规则（允许多个类别） ==========
def has_substituent(subst, key):
    if pd.isna(subst) or subst.strip() == "":
        return False
    return key in subst.split(",")

df["is_S"] = df["Substituent"].apply(lambda x: has_substituent(x, "S") or has_substituent(x, "O6S"))
df["is_M"] = df["Substituent"].apply(lambda x: has_substituent(x, "M") or has_substituent(x, "O6M"))
df["is_A"] = df["Substituent"].apply(lambda x: has_substituent(x, "A") or has_substituent(x, "O6A"))
df["is_P"] = df["Substituent"].apply(lambda x: has_substituent(x, "P") or has_substituent(x, "O6P"))

# Cyclic H-bond（S 的子集）
df["is_Cyclic"] = df["is_S"] & (df["whether cyclic Hbond"].astype(str).str.lower() == "true")

# ========== 颜色定义 ==========
colors = {
    "S": "red",
    "Cyclic H-bond": "green",
    "M": "blue",
    "A": "purple",
    "P": "orange",
    "others": "gray"
}

# ========== 百分比计算（允许重复统计） ==========
total = len(df)
counts = {
    "S": df["is_S"].sum() / total * 100,
    "Cyclic H-bond": df["is_Cyclic"].sum() / total * 100,
    "M": df["is_M"].sum() / total * 100,
    "A": df["is_A"].sum() / total * 100,
    "P": df["is_P"].sum() / total * 100,
    "others": (  # others = 不属于 S/M/A/P 的点
        (~(df["is_S"] | df["is_M"] | df["is_A"] | df["is_P"]))
        .sum() / total * 100
    )
}

# ========== 绘制顺序 ==========
plot_order = ["others", "A", "P", "M", "S", "Cyclic H-bond"]

plt.figure(figsize=(8, 6))
point_size = 5

for cls in plot_order:
    if cls == "S":
        subset = df[df["is_S"]]
        label = f"S ({counts['S']:.1f}%)"
    elif cls == "Cyclic H-bond":
        subset = df[df["is_Cyclic"]]
        label = f"   Cyclic H-bond ({counts['Cyclic H-bond']:.1f}%)"
    elif cls == "M":
        subset = df[df["is_M"]]
        label = f"M ({counts['M']:.1f}%)"
    elif cls == "A":
        subset = df[df["is_A"]]
        label = f"A ({counts['A']:.1f}%)"
    elif cls == "P":
        subset = df[df["is_P"]]
        label = f"P ({counts['P']:.1f}%)"
    else:
        subset = df[~(df["is_S"] | df["is_M"] | df["is_A"] | df["is_P"])]
        label = f"others ({counts['others']:.1f}%)"

    plt.scatter(
        subset["API_RMSE_all"],
        subset["HPMCAS_RMSE_all"],
        c=colors[cls],
        label=label,
        alpha=0.7,
        s=point_size,
        edgecolors="none",
        zorder=plot_order.index(cls)
    )

# ========== Legend 设置 ==========
order_for_legend = ["S", "Cyclic H-bond", "M", "A", "P", "others"]
handles, labels = plt.gca().get_legend_handles_labels()
new_handles, new_labels = [], []
for name in order_for_legend:
    for h, l in zip(handles, labels):
        if l.startswith(name) or l.strip().startswith(name):
            new_handles.append(h)
            new_labels.append(l)

legend = plt.legend(
    new_handles, new_labels,
    title="Substituents",
    scatterpoints=1,
    markerscale=3,
    loc="upper right"
)

# 缩小 Cyclic H-bond 的字体
for text in legend.get_texts():
    if "Cyclic H-bond" in text.get_text():
        text.set_fontsize(8)

plt.xlabel("RMSE $_{API}$", fontsize=18, fontweight="bold")
plt.ylabel("RMSE $_{HPMCAS}$", fontsize=18, fontweight="bold")
# plt.title("10%ASD H-bond Analysis", fontsize=18, fontweight="bold")
plt.title(
          f"10%ASD H-bond Analysis(Top {percent_to_keep} %)", 
        #   f"10%ASD H-bond Analysis",
          fontsize=18, fontweight="bold"
          )
plt.grid(True)
plt.tight_layout()
plt.show()
