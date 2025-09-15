import pandas as pd
import matplotlib.pyplot as plt

# ========== 读取数据 ==========
file_path = r"E:\new_HPMCAS\script\20ASD_new_RMSE.csv"
df = pd.read_csv(file_path)

# ========== 筛选 Ebinding(kJ/mol) < 0 ==========
df = df[df["Ebinding(kJ/mol)"] < 0].copy()

# ========== 参数设置 ==========
percent_to_keep = 1   # ⚠️ 修改这里，比如10表示取最小10%的点，100表示全部

# ========== 筛选 Top X% ==========
df = df.sort_values(by="Distance_to_origin", ascending=True).reset_index(drop=True)
n_keep = int(len(df) * percent_to_keep / 100)
if n_keep < 1:
    n_keep = 1
df = df.iloc[:n_keep]

# ========== 分类规则 ==========
def classify_substituent(subst, row):
    if pd.isna(subst) or str(subst).strip() == "":
        # without H-bond
        if (pd.isna(row["124N"]) or str(row["124N"]).strip() == "") \
           and (pd.isna(row["125O"]) or str(row["125O"]).strip() == "") \
           and (pd.isna(row["110H"]) or str(row["110H"]).strip() == ""):
            return "without H-bond"
        else:
            return "unknown"

    subst_str = str(subst)
    if "S" in subst_str:
        return "S"
    if "M" in subst_str:
        return "M"
    if "A" in subst_str:
        return "A"
    if "P" in subst_str:
        return "P"

    return "unknown"

# 主分类（用于统计）
df["Class"] = df.apply(lambda r: classify_substituent(r["Substituent"], r), axis=1)

# 绘图分类（在 S 里单独分出 cyclic）
df["PlotClass"] = df["Class"].copy()
df.loc[(df["Class"] == "S") & (df["whether cyclic Hbond"].astype(str).str.lower() == "true"),
       "PlotClass"] = "Cyclic H-bond"

# ========== 颜色定义 ==========
colors = {
    "S": "red",
    "Cyclic H-bond": "green",
    "M": "blue",
    "A": "purple",
    "P": "orange",
    "without H-bond": "brown",   # 改成棕色
    "unknown": "gray"            # unknown 用灰色
}

# ========== 百分比计算 ==========
total = len(df)
counts = {
    "S": len(df[df["Class"] == "S"]) / total * 100 if total > 0 else 0,
    "M": len(df[df["Class"] == "M"]) / total * 100 if total > 0 else 0,
    "A": len(df[df["Class"] == "A"]) / total * 100 if total > 0 else 0,
    "P": len(df[df["Class"] == "P"]) / total * 100 if total > 0 else 0,
    "without H-bond": len(df[df["Class"] == "without H-bond"]) / total * 100 if total > 0 else 0,
    # Cyclic H-bond 单独统计
    "Cyclic H-bond": len(df[df["whether cyclic Hbond"].astype(str).str.lower() == "true"]) / total * 100 if total > 0 else 0
}

# ========== 绘制顺序（底层到顶层） ==========
plot_order = ["without H-bond", "unknown", "A", "P", "M", "S", "Cyclic H-bond"]

plt.figure(figsize=(8, 6))
point_size = 5

for cls in plot_order:
    subset = df[df["PlotClass"] == cls]
    if len(subset) == 0:
        continue
    if cls == "Cyclic H-bond":
        label = f"   Cyclic H-bond ({counts['Cyclic H-bond']:.1f}%)"
    elif cls == "unknown":
        label = None  # 不在 legend 显示
    else:
        label = f"{cls} ({counts[cls]:.1f}%)"
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
order_for_legend = ["S", "Cyclic H-bond", "M", "A", "P", "without H-bond"]
handles, labels = plt.gca().get_legend_handles_labels()
new_handles, new_labels = [], []
for name in order_for_legend:
    for h, l in zip(handles, labels):
        if l and l.strip().startswith(name):  # 排除 None
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
plt.title(f"20%ASD H-bond Analysis (Top {percent_to_keep}%)",
          fontsize=18, fontweight="bold")
plt.grid(True)
plt.tight_layout()
plt.show()
