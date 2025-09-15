
import os
import numpy as np
import pandas as pd

# ========== filepath ==========
folder = "/mnt/fastscratch/users/sgdzheng/20ASD_xyz/"   
output_csv = "/mnt/fastscratch/users/sgdzheng/20ASD_Hbonds_results_2.csv"

# ========== function of angel ==========
def distance(a, b):
    return np.linalg.norm(a - b)

def angle(a, b, c):
    v1 = a - b
    v2 = c - b
    cosang = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
    return np.degrees(np.arccos(np.clip(cosang, -1.0, 1.0)))

def parse_xyz(filepath):
    with open(filepath, "r") as f:
        lines = f.readlines()
    atoms = []
    for i, line in enumerate(lines[2:], start=1):  # skip the first two lines
        parts = line.split()
        if len(parts) >= 4:
            atom = parts[0]
            x, y, z = map(float, parts[1:4])
            atoms.append((i, atom, x, y, z))
    return atoms

# ========== atom index ==========
donor_N = 124
donor_H_for_N = 127
acceptor_O = 125
donor_H = 110
donor_H_parentO = 116

# ========== Substituent mapping rules ==========
substituent_map = {
    "M": ["O23", "O35", "O48", "O50", "O37", "O52"],
    "P": ["O47"],
    "A": ["O42"],
    "S": ["O28", "O32", "O33"],
    "O6P": ["O43"],
    "O6A": ["O39"],
    "O6S": ["O26"]
}

def get_substituent(hb_124N):
    found = []
    for label, targets in substituent_map.items():
        for t in targets:
            if t in hb_124N:
                found.append(label)
                break
    return ",".join(found)

# ========== processing ==========
results = []
files = [f for f in os.listdir(folder) if f.endswith(".xyz")]
total = len(files)

for idx, fname in enumerate(files, start=1):
    atoms = parse_xyz(os.path.join(folder, fname))
    coords = np.array([(x, y, z) for (_, _, x, y, z) in atoms])
    elements = [a for (_, a, _, _, _) in atoms]

    hb_124N, hb_125O, hb_110H = [], [], []

    # --- 124 N donor（H127）
    idx_N = donor_N - 1
    idx_HN = donor_H_for_N - 1
    for j, elem_j in enumerate(elements):
        if j not in [idx_N, idx_HN] and elem_j in ["O", "N"]:
            d_HA = distance(coords[idx_HN], coords[j])
            if d_HA <= 2.5:
                ang = angle(coords[idx_N], coords[idx_HN], coords[j])
                if ang >= 130:
                    hb_124N.append(f"{elem_j}{j+1}")

    # --- 125 O acceptor
    idx_O = acceptor_O - 1
    for i, elem_i in enumerate(elements):
        if elem_i == "H":
            for j, elem_j in enumerate(elements):
                if j != i and elem_j in ["O", "N"]:
                    if distance(coords[j], coords[i]) < 1.2:
                        d_HA = distance(coords[i], coords[idx_O])
                        if d_HA <= 2.5:
                            ang = angle(coords[j], coords[i], coords[idx_O])
                            if ang >= 130:
                                hb_125O.append(f"{elem_j}{j+1}")

    # --- 110 H donor（O 116）
    idx_H = donor_H - 1
    idx_OH = donor_H_parentO - 1
    for j, elem_j in enumerate(elements):
        if j not in [idx_H, idx_OH] and elem_j in ["O", "N"]:
            d_HA = distance(coords[idx_H], coords[j])
            if d_HA <= 2.5:
                ang = angle(coords[idx_OH], coords[idx_H], coords[j])
                if ang >= 130:
                    hb_110H.append(f"{elem_j}{j+1}")

    # --- cyclic hydrogen bond condition
    cyclic_hbond = False
    if "O28" in hb_124N and "O33" in hb_125O:
        cyclic_hbond = True

    # --- Substituent
    substituent = get_substituent(hb_124N)

    results.append({
        "file": fname,
        "Substituent": substituent,
        "124N": ", ".join(hb_124N) if hb_124N else "",
        "125O": ", ".join(hb_125O) if hb_125O else "",
        "110H": ", ".join(hb_110H) if hb_110H else "",
        "whether cyclic Hbond": "true" if cyclic_hbond else "false"
    })

    print(f"Processed {idx}/{total} : {fname}")

# ========== 保存结果 ==========
df = pd.DataFrame(results)
df.to_csv(output_csv, index=False, encoding="utf-8")
print(f"save result to {output_csv}")
