import pandas as pd
import matplotlib.pyplot as plt

def process_file(file_path):
    df = pd.read_csv(file_path)
    df_sub = df[["Substituent", "E_Binding(eV)", "whether cyclic Hbond"]].copy()
    df_sub = df_sub[df_sub["E_Binding(eV)"] < 0]  # 只保留结合能为负值

    # 分类函数：只保留单一 M, P, S, A
    def simplify_substituent(x):
        if pd.isna(x):
            return None
        if "," in str(x):  # 多取代基直接丢弃
            return None
        if str(x) in ["M", "P", "S", "A"]:
            return str(x)
        return None

    df_sub["Substituent_grouped"] = df_sub["Substituent"].apply(simplify_substituent)
    df_sub = df_sub.dropna(subset=["Substituent_grouped"])

    # 单位换算
    df_sub["E_Binding(kJ/mol)"] = df_sub["E_Binding(eV)"] * 96.485

    # 规则：cyclic Hbond == True 仅允许 S
    mask_invalid = (df_sub["whether cyclic Hbond"] == True) & (df_sub["Substituent_grouped"] != "S")
    df_sub = df_sub[~mask_invalid]

    # 按 Substituent 和 cyclic Hbond 分组计算平均结合能
    avg_binding_grouped = (
        df_sub.groupby(["Substituent_grouped", "whether cyclic Hbond"])["E_Binding(kJ/mol)"]
        .mean()
        .reset_index()
    )

    # 排序
    order = ["M", "P", "S", "A"]
    avg_binding_grouped["Substituent_grouped"] = pd.Categorical(
        avg_binding_grouped["Substituent_grouped"], categories=order, ordered=True
    )
    avg_binding_grouped = avg_binding_grouped.sort_values("Substituent_grouped")

    return avg_binding_grouped

# 文件路径
file_10 = r"E:\new_HPMCAS\script\10ASD_new_RMSE.csv"
file_20 = r"E:\new_HPMCAS\script\20ASD_new_RMSE.csv"

# 分别处理
avg_10 = process_file(file_10)
avg_20 = process_file(file_20)

# 作图：左右两个子图
fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharey=True)

colors = {True: "royalblue", False: "skyblue"}

# 左图：10% ASD
for hbond in [True, False]:
    data = avg_10[avg_10["whether cyclic Hbond"] == hbond]
    axes[0].bar(
        data["Substituent_grouped"],
        data["E_Binding(kJ/mol)"],
        width=0.4,
        align="center" if hbond else "edge",
        color=colors[hbond],
        label=f"Cyclic H-bond = {hbond}"
    )
axes[0].set_title("10% ASD", fontsize=18, fontweight="bold")
axes[0].set_xlabel("Substituent", fontsize=16, fontweight="bold")
axes[0].set_ylabel("E$_{binding}$ (kJ/mol)", fontsize=16, fontweight="bold")
axes[0].legend(loc="lower right")
axes[0].grid(axis="y", linestyle="--", alpha=0.6)

# 右图：20% ASD
for hbond in [True, False]:
    data = avg_20[avg_20["whether cyclic Hbond"] == hbond]
    axes[1].bar(
        data["Substituent_grouped"],
        data["E_Binding(kJ/mol)"],
        width=0.4,
        align="center" if hbond else "edge",
        color=colors[hbond],
        label=f"Cyclic H-bond = {hbond}"
    )
axes[1].set_title("20% ASD", fontsize=18, fontweight="bold")
axes[1].set_xlabel("Substituent", fontsize=16, fontweight="bold")
axes[1].legend(loc="lower right")
axes[1].grid(axis="y", linestyle="--", alpha=0.6)

plt.tight_layout()
plt.show()
