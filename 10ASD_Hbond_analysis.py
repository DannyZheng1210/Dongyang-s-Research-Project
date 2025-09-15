import os
import numpy as np
import pandas as pd

# ========== Paths ==========
folder = "/mnt/fastscratch/users/sgdzheng/10ASD_xyz/"   # <-- change to your xyz directory
output_csv = "/mnt/fastscratch/users/sgdzheng/hbonds_results_10ASD.csv"

# ========== Utility functions ==========
def distance(a, b):
    return np.linalg.norm(a - b)

def angle(a, b, c):
    """angle a-b-c, with b as vertex"""
    v1 = a - b
    v2 = c - b
    cosang = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
    return np.degrees(np.arccos(np.clip(cosang, -1.0, 1.0)))

def parse_xyz(filepath):
    with open(filepath, "r") as f:
        lines = f.readlines()
    atoms = []
    for i, line in enumerate(lines[2:], start=1):  # skip first two lines
        parts = line.split()
        if len(parts) >= 4:
            atom = parts[0]
            x, y, z = map(float, parts[1:4])
            atoms.append((i, atom, x, y, z))
    return atoms

# ========== Specific atom indices (10ASD = 20ASD + 109) ==========
OFFSET = 109
donor_N = 124 + OFFSET        # 233 (amide N)
donor_H_for_N = 127 + OFFSET  # 236 (amide N-H)
acceptor_O = 125 + OFFSET     # 234 (amide C=O O)
donor_H = 110 + OFFSET        # 219 (phenolic OH H)
donor_H_parentO = 116 + OFFSET # 225 (phenolic OH O)

# ========== Substituent mapping (both HPMCAS chains) ==========
substituent_map = {
    "M": ["O23", "O35", "O48", "O50", "O37", "O52",
          "O132", "O144", "O157", "O159", "O146", "O161"],
    "P": ["O47", "O156"],
    "A": ["O42", "O151"],
    "S": ["O28", "O32", "O33", "O137", "O141", "O142"],
    "O6P": ["O43", "O152"],
    "O6A": ["O39", "O148"],
    "O6S": ["O26", "O135"]
}

def get_substituent(hb_124N):
    found = []
    for label, targets in substituent_map.items():
        for t in targets:
            if t in hb_124N:
                found.append(label)
                break
    return ",".join(found)

# ========== Main processing ==========
results = []
files = [f for f in os.listdir(folder) if f.endswith(".xyz")]
total = len(files)

for idx, fname in enumerate(files, start=1):
    atoms = parse_xyz(os.path.join(folder, fname))
    coords = np.array([(x, y, z) for (_, _, x, y, z) in atoms])
    elements = [a for (_, a, _, _, _) in atoms]

    hb_124N, hb_125O, hb_110H = [], [], []

    # --- API amide N donor (233N-H236)
    idx_N = donor_N - 1
    idx_HN = donor_H_for_N - 1
    for j, elem_j in enumerate(elements):
        if j not in [idx_N, idx_HN] and elem_j in ["O", "N"]:
            d_HA = distance(coords[idx_HN], coords[j])
            if d_HA <= 2.5:
                ang = angle(coords[idx_N], coords[idx_HN], coords[j])
                if ang >= 130:
                    hb_124N.append(f"{elem_j}{j+1}")

    # --- API amide O acceptor (234)
    idx_O = acceptor_O - 1
    for i, elem_i in enumerate(elements):
        if elem_i == "H":
            for j, elem_j in enumerate(elements):
                if j != i and elem_j in ["O", "N"]:
                    if distance(coords[j], coords[i]) < 1.2:  # covalent bond
                        d_HA = distance(coords[i], coords[idx_O])
                        if d_HA <= 2.5:
                            ang = angle(coords[j], coords[i], coords[idx_O])
                            if ang >= 130:
                                hb_125O.append(f"{elem_j}{j+1}")

    # --- API phenolic OH donor (219H, O225)
    idx_H = donor_H - 1
    idx_OH = donor_H_parentO - 1
    for j, elem_j in enumerate(elements):
        if j not in [idx_H, idx_OH] and elem_j in ["O", "N"]:
            d_HA = distance(coords[idx_H], coords[j])
            if d_HA <= 2.5:
                ang = angle(coords[idx_OH], coords[idx_H], coords[j])
                if ang >= 130:
                    hb_110H.append(f"{elem_j}{j+1}")

    # --- Cyclic hydrogen bond condition ---
    cyclic_hbond = False

    # HPMCAS1: API 233N-H236 → O32 AND H81 → API O234
    if ("O32" in hb_124N) and ("H81" in hb_125O):
        cyclic_hbond = True

    # HPMCAS2: API 233N-H236 → O141 AND H190 → API O234
    if ("O141" in hb_124N) and ("H190" in hb_125O):
        cyclic_hbond = True

    # --- Substituent mapping
    substituent = get_substituent(hb_124N)

    results.append({
        "file": fname,
        "Substituent": substituent,
        "124N_donor": ", ".join(hb_124N) if hb_124N else "",
        "125O_acceptor": ", ".join(hb_125O) if hb_125O else "",
        "110H_donor": ", ".join(hb_110H) if hb_110H else "",
        "cyclic_Hbond": "true" if cyclic_hbond else "false"
    })

    print(f"Processed {idx}/{total} : {fname}")

# ========== Save results ==========
df = pd.DataFrame(results)
df.to_csv(output_csv, index=False, encoding="utf-8")
print(f"✅ Results saved to {output_csv}")
