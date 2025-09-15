from ase.io import read
from shiftml.ase import ShiftML
import numpy as np
import csv
import os

# 初始化 ShiftML 模型
calculator = ShiftML("ShiftML3")

# 输入/输出路径
input_folder = "20ASD"
output_folder = "20ASD_shiftml"
os.makedirs(output_folder, exist_ok=True)

# 找到所有 .cif 文件
cif_files = [f for f in os.listdir(input_folder) if f.endswith(".cif")]

for idx, input_file in enumerate(cif_files, 1):
    input_path = os.path.join(input_folder, input_file)
    base_name = os.path.splitext(input_file)[0]
    output_csv = os.path.join(output_folder, f"{base_name}_ShiftML_results.csv")

    # 读取结构
    frame = read(input_path)

    # 预测屏蔽值和不确定性
    cs_iso = np.array(calculator.get_cs_iso(frame))
    cs_committee_iso = calculator.get_cs_iso_ensemble(frame)
    uncertainty = np.std(cs_committee_iso, axis=1)
    atom_types = frame.get_chemical_symbols()

    # 收集当前文件结果
    csv_rows = []
    for i, (atom_type, sigma, unc) in enumerate(zip(atom_types, cs_iso, uncertainty)):
        if atom_type == "C":
            shift = -0.9732 * sigma + 166.23
        elif atom_type == "H":
            shift = -0.9024 * sigma + 28.05
        elif atom_type == "N":
            shift = -1.0250 * sigma + 183.34
        else:
            shift = np.nan
        atom_index = i + 1
        csv_rows.append([atom_index, atom_type, sigma, shift, unc])

    # 写入 CSV 文件
    with open(output_csv, mode="w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "Atom Index", "Atom Type", 
            "Shielding (ppm)", "Chemical Shift (ppm)", "Uncertainty (ppm)"
        ])
        for row in csv_rows:
            writer.writerow([
                row[0],
                row[1],
                f"{row[2]:.6f}",   # 保留 6 位小数，避免精度丢失
                f"{row[3]:.6f}" if not np.isnan(row[3]) else "",
                f"{row[4]:.6f}"
            ])

    # 打印进度提示
    print(f"{input_file} 计算完成 ({idx}/{len(cif_files)})")

print(f"\n全部完成！结果保存在文件夹: {output_folder}")
