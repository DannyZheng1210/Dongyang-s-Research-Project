import os
import pandas as pd
import numpy as np
from collections import defaultdict
from sklearn.metrics import mean_squared_error
from tqdm import tqdm   # NEW: progress bar

# === 10ASD experimental values -> Atom Index ===

experiment_pairs_C_10asd = [
    # HPMCAS1 (1–109)
    (58 , [34,36,38,49,51,53]),
    (70 , [44]),
    (61 , [45]),
    (17 , [46]),
    (173, [40]),
    (21 , [41]),
    (177, [27]),
    (29 , [29]),
    (29 , [30]),
    (177, [31]),
    (102, [13,4]),
    (75 , [19,18,12,11,3,2]),
    (84 , [17,10]),
    (61 , [22,15,6]),
    (70 , [24,16,9]),
    # HPMCAS2 (110–218)
    (58 , [143,145,147,158,160,162]),
    (70 , [153]),
    (61 , [154]),
    (17 , [155]),
    (173, [149]),
    (21 , [150]),
    (177, [136]),
    (29 , [138]),
    (29 , [139]),
    (177, [140]),
    (102, [122,113]),
    (75 , [128,127,121,120,112,111]),
    (84 , [126,119]),
    (61 , [131,124,115]),
    (70 , [133,125,118]),
    # API (219–238)
    (24.1, [232]),
    (171 , [231]),
    (131 , [230]),
    (121 , [223,238]),
    (116 , [224,227]),
    (154 , [222]),
]

experiment_pairs_H_10asd = [
    # HPMCAS1 (1–109)
    (3.0, [82,83,84,85,86,87,104,105,106,101,102,103,88,89,90,107,108,109]),
    (3.9, [94,95]),
    (3.4, [96]),
    (1.1, [91,92,93]),
    (2.0, [77,78,79,80]),
    (4.8, [72,73,65,57]),
    (3.4, [71,70,74,64,63,66,56,55,58]),
    (2.5, [69,62,54]),
    (3.9, [75,76,67,68,60,61]),
    # HPMCAS2 (110–218)
    (3.0, [191,192,193,194,195,196,213,214,215,210,211,212,197,198,199,216,217,218]),
    (3.9, [203,204]),
    (3.4, [205]),
    (1.1, [200,201,202]),
    (2.0, [186,187,188,189]),
    (4.8, [181,182,174,166]),
    (3.4, [180,179,183,173,172,175,165,164,167]),
    (2.5, [178,171,163]),
    (3.9, [184,185,176,177,169,170]),
    # API (219–238)
    (1.1, [228,229,237]),
    (7.4, [220,235]),
    (6.8, [221,226]),
    (8.5, [219,236]),
]

def merge_pairs_to_mapping(pairs):
    merged = defaultdict(list)
    for exp_val, idx_list in pairs:
        merged[float(exp_val)].extend(int(i) for i in idx_list)
    out = {}
    for exp_val, lst in merged.items():
        seen = set()
        uniq = []
        for x in lst:
            if x not in seen:
                seen.add(x)
                uniq.append(x)
        out[exp_val] = uniq
    return out

experiment_mapping_C = merge_pairs_to_mapping(experiment_pairs_C_10asd)
experiment_mapping_H = merge_pairs_to_mapping(experiment_pairs_H_10asd)

def mapping_to_one2one(experiment_mapping: dict) -> dict:
    out = {}
    for exp_val, atom_list in experiment_mapping.items():
        for idx in atom_list:
            out[int(idx)] = float(exp_val)
    return out

exp_map_C = mapping_to_one2one(experiment_mapping_C)
exp_map_H = mapping_to_one2one(experiment_mapping_H)

# API atom sets (10ASD)
api_atoms_C = {232,231,230,223,238,224,227,222}
api_atoms_H = {228,229,237,220,235,221,226,219,236}

def calculate_rmse(exp_map, pred_map, include=None, exclude=None):
    if include is not None:
        indices = [i for i in include if i in pred_map]
    elif exclude is not None:
        indices = [i for i in exp_map if i not in exclude and i in pred_map]
    else:
        indices = [i for i in exp_map if i in pred_map]

    if not indices:
        return np.nan

    exp_vals = np.array([exp_map[i] for i in indices], dtype=float)
    pred_vals = np.array([pred_map[i] for i in indices], dtype=float)
    return np.sqrt(mean_squared_error(exp_vals, pred_vals))

# === Main program ===
folder_path = "/mnt/fastscratch/users/sgdzheng/10ASD_shiftml_result/"
csv_files = sorted(f for f in os.listdir(folder_path) if f.endswith('.csv'))

rows = []
for csv_file in tqdm(csv_files, desc="Processing files"):   # NEW: progress bar
    file_path = os.path.join(folder_path, csv_file)

    # read once
    df = pd.read_csv(file_path)
    df['Atom Index'] = pd.to_numeric(df['Atom Index'], errors='coerce').astype('Int64')
    df['Chemical Shift (ppm)'] = pd.to_numeric(df['Chemical Shift (ppm)'], errors='coerce')

    # build lookup dicts
    pred_C = df[df['Atom Type']=='C'].set_index('Atom Index')['Chemical Shift (ppm)'].to_dict()
    pred_H = df[df['Atom Type']=='H'].set_index('Atom Index')['Chemical Shift (ppm)'].to_dict()

    # API
    rmse_C_api = calculate_rmse(exp_map_C, pred_C, include=api_atoms_C)
    rmse_H_api = calculate_rmse(exp_map_H, pred_H, include=api_atoms_H)
    rmse_api_all = calculate_rmse({**exp_map_C, **exp_map_H}, {**pred_C, **pred_H},
                                  include=api_atoms_C | api_atoms_H)

    # HPMCAS (everything except API)
    rmse_C_nonapi = calculate_rmse(exp_map_C, pred_C, exclude=api_atoms_C)
    rmse_H_nonapi = calculate_rmse(exp_map_H, pred_H, exclude=api_atoms_H)
    rmse_hpmcas_all = calculate_rmse({**exp_map_C, **exp_map_H}, {**pred_C, **pred_H},
                                     exclude=api_atoms_C | api_atoms_H)

    # distance to origin
    if not np.isnan(rmse_api_all) and not np.isnan(rmse_hpmcas_all):
        distance = np.sqrt(rmse_api_all**2 + rmse_hpmcas_all**2)
    else:
        distance = np.nan

    rows.append([
        csv_file,
        rmse_C_api, rmse_H_api, rmse_api_all,
        rmse_C_nonapi, rmse_H_nonapi, rmse_hpmcas_all,
        distance
    ])

results_df = pd.DataFrame(
    rows,
    columns=[
        'File Name',
        'API_13C_RMSE', 'API_1H_RMSE', 'API_RMSE_all',
        'HPMCAS_13C_RMSE', 'HPMCAS_1H_RMSE', 'HPMCAS_RMSE_all',
        'Distance_to_origin'
    ]
)

results_df.to_csv('10ASD_new_RMSE.csv', index=False)
print("✅ All files processed, results saved.")
