import os 
import pandas as pd
import numpy as np
from collections import defaultdict
from sklearn.metrics import mean_squared_error

# === 13C 实验值 -> Atom Index
experiment_pairs_C = [
    #atoms in API
    (	24.1	,	[	123	]	)	,
    (	171	,	[	122	]	)	,
    (	131	,	[	121	]	)	,
    (	121	,	[	114,129	]	)	,
    (	116	,	[	115,118	]	)	,
    (	154	,	[	113	]	)	,

    #following atoms is the atoms in substitution
    (	58	,	[	34,36,38,49,51,53	]	)	,
    (	70	,	[	44	]	)	,
    (	61	,	[	45	]	)	,
    (	17	,	[	46	]	)	,
    (	173	,	[	40	]	)	,
    (	21	,	[	41	]	)	,
    (	177	,	[	27	]	)	,
    (	29	,	[	29	]	)	,
    (	29	,	[	30	]	)	,
    (	177	,	[	31	]	)	,    
    
    # #following atoms is the atom in celluloses' rings
    (	102	,	[	13,4]	)	,
    (	75	,	[	19,18,12,11,3,2]	)	,
    (	84	,	[	17,	10]	)	,
    (	61	,	[	22,	15,	6]	)	,
    (	70	,	[	24	,16	, 9]	)	,
]

# === 1H 实验值 -> Atom Index
experiment_pairs_H = [
    #API
    (1.1, [119,120,128]),
    (7.4, [111,126]),
    (6.8, [112,117]),
    (8.5, [110,127]),
    #Atoms in substitution
    (3.0, [82,83,84, 85,86,87,	104,105,106, 101,102,103, 88,89,90,	107,108,109]),
    (3.9, [94,95]),	
    (3.4, [96]),	
    (1.1, [91,92,93]),	
    (2.0, [77,78,79,80]),
    #atoms in rings
    (4.8, [72,73, 65, 57]),
    (3.4, [71,70,74, 64,63,66, 56,55,58]),
    (2.5, [69, 62, 54]),
    (3.9, [75,76, 67,68, 60,61]),

]

def merge_pairs_to_mapping(pairs):
    """将 [(exp, [idxs]), ...] 合并为 {exp: [idxs...]}，避免重复键覆盖。"""
    from collections import defaultdict
    merged = defaultdict(list)
    for exp_val, idx_list in pairs:
        merged[float(exp_val)].extend(int(i) for i in idx_list)
    # 去重但保持顺序
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

# 合并后的映射
experiment_mapping_C = merge_pairs_to_mapping(experiment_pairs_C)
experiment_mapping_H = merge_pairs_to_mapping(experiment_pairs_H)

def mapping_to_one2one(experiment_mapping: dict) -> dict:
    """把 {exp: [idx, idx,...]} 展开为 {idx: exp} 的一对一映射"""
    out = {}
    for exp_val, atom_list in experiment_mapping.items():
        for idx in atom_list:
            out[int(idx)] = float(exp_val)
    return out

def _get_aligned_arrays(
    file_path: str,
    element: str,
    experiment_mapping: dict,
    include_indices: set | None = None,
    exclude_indices: set | None = None,
):
    """返回按索引对齐后的 (experimental, predicted) 数组。"""
    try:
        df = pd.read_csv(file_path)
        df = df[df['Atom Type'] == element].copy()

        df['Atom Index'] = pd.to_numeric(df['Atom Index'], errors='coerce').astype('Int64')
        df['Chemical Shift (ppm)'] = pd.to_numeric(df['Chemical Shift (ppm)'], errors='coerce')

        exp_one2one = mapping_to_one2one(experiment_mapping)

        atoms_to_extract = list(exp_one2one.keys())
        if include_indices is not None:
            atoms_to_extract = [a for a in atoms_to_extract if a in include_indices]
        if exclude_indices is not None:
            atoms_to_extract = [a for a in atoms_to_extract if a not in exclude_indices]

        if not atoms_to_extract:
            return None, None

        sub = df[df['Atom Index'].isin(atoms_to_extract)].copy()

        have = set(sub['Atom Index'].dropna().astype(int).tolist())
        need = set(atoms_to_extract)
        missing = sorted(list(need - have))
        if missing:
            print(f"[⚠️] {os.path.basename(file_path)} 缺少 {element} 原子索引: {missing}")
            return None, None

        sub = sub.set_index('Atom Index').loc[atoms_to_extract]

        if sub['Chemical Shift (ppm)'].isnull().any():
            print(f"[⚠️] 存在空值: {file_path} ({element})")
            return None, None

        predicted = sub['Chemical Shift (ppm)'].astype(float).values
        experimental = np.array([exp_one2one[idx] for idx in atoms_to_extract], dtype=float)
        return experimental, predicted
    except Exception as e:
        print(f"[❌] 处理失败 {file_path} ({element}): {e}")
        return None, None

def calculate_rmse_for_element(file_path, element, experiment_mapping, include_indices=None, exclude_indices=None):
    """元素级别 RMSE（API 或 HPMCAS 分别统计）"""
    experimental, predicted = _get_aligned_arrays(
        file_path, element, experiment_mapping,
        include_indices=include_indices, exclude_indices=exclude_indices
    )
    if experimental is None:
        return np.nan
    return np.sqrt(mean_squared_error(experimental, predicted))

def calculate_rmse_api_all(file_path):
    exp_C, pred_C = _get_aligned_arrays(file_path, 'C', experiment_mapping_C, include_indices=api_atoms_C)
    exp_H, pred_H = _get_aligned_arrays(file_path, 'H', experiment_mapping_H, include_indices=api_atoms_H)
    if exp_C is None or exp_H is None:
        return np.nan
    return np.sqrt(mean_squared_error(np.concatenate([exp_C, exp_H]),
                                      np.concatenate([pred_C, pred_H])))

def calculate_rmse_hpmcas_all(file_path):
    exp_C, pred_C = _get_aligned_arrays(file_path, 'C', experiment_mapping_C, exclude_indices=api_atoms_C)
    exp_H, pred_H = _get_aligned_arrays(file_path, 'H', experiment_mapping_H, exclude_indices=api_atoms_H)
    if exp_C is None or exp_H is None:
        return np.nan
    return np.sqrt(mean_squared_error(np.concatenate([exp_C, exp_H]),
                                      np.concatenate([pred_C, pred_H])))

# === 定义 API 原子集合（顺序变了的话也要改这里） ===
api_atoms_C = {123,122,121,114,129,115,118,113}
api_atoms_H = {119,120,128,111,126,112,117,110,127}

# === 主程序 ===
folder_path = '20ASD_shiftml_result'  # 结果文件夹
csv_files = sorted(f for f in os.listdir(folder_path) if f.endswith('.csv'))

rows = []
for csv_file in csv_files:
    file_path = os.path.join(folder_path, csv_file)

    # API 部分
    rmse_C_api = calculate_rmse_for_element(file_path, 'C', experiment_mapping_C, include_indices=api_atoms_C)
    rmse_H_api = calculate_rmse_for_element(file_path, 'H', experiment_mapping_H, include_indices=api_atoms_H)
    rmse_api_all = calculate_rmse_api_all(file_path)

    # HPMCAS 部分
    rmse_C_nonapi = calculate_rmse_for_element(file_path, 'C', experiment_mapping_C, exclude_indices=api_atoms_C)
    rmse_H_nonapi = calculate_rmse_for_element(file_path, 'H', experiment_mapping_H, exclude_indices=api_atoms_H)
    rmse_hpmcas_all = calculate_rmse_hpmcas_all(file_path)

    # 到原点的距离
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

results_df.to_csv('20ASD_new_RMSE.csv', index=False)
print("✅ 所有文件处理完成，结果已保存。")
