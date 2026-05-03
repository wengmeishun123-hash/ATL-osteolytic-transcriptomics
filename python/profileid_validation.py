import csv
import os

CCLE_DIR = r"D:\ATL research\analysis\data\ccle"
OUT_DIR = r"D:\ATL research\analysis\results\celline"

known_profiles = {
    'HUT102': 'PR-Myvrcc', 'MJ': 'PR-wNql0Y', 'ATN1': 'PR-xbYnX2',
    'HUT78': 'PR-xQbaOc', 'JURKAT': 'PR-tYD2aE', 'MOLT4': 'PR-ns3q2x', 'CCRFCEM': 'PR-DDO590',
}

profile_files = sorted([f for f in os.listdir(CCLE_DIR) if f.startswith("OmicsProfiles")])
print(f"Found {len(profile_files)} OmicsProfiles files\n")

for bn in profile_files:
    fpath = os.path.join(CCLE_DIR, bn)
    with open(fpath, 'r', encoding='utf-8', errors='replace') as f:
        reader = csv.DictReader(f)
        header = reader.fieldnames
        rows = list(reader)
    
    prof_col = 'ProfileID' if 'ProfileID' in header else None
    model_col = 'ModelID' if 'ModelID' in header else None
    datatype_col = 'DataType' if 'DataType' in header else ('Datatype' if 'Datatype' in header else None)
    
    all_pids = set(r.get(prof_col, '') for r in rows) if prof_col else set()
    
    print(f"=== {bn} ({len(rows)} rows) ===")
    print(f"  Columns: {header}")
    
    for cl, pid in known_profiles.items():
        found = pid in all_pids if prof_col else False
        if found:
            matches = [r for r in rows if r.get(prof_col, '') == pid]
            m = matches[0]
            model_id = m.get(model_col, 'N/A') if model_col else 'N/A'
            dtype = m.get(datatype_col, 'N/A') if datatype_col else 'N/A'
            print(f"  {cl}: {pid} -> FOUND (ModelID={model_id}, DataType={dtype})")
        else:
            print(f"  {cl}: {pid} -> NOT FOUND")
    
    if model_col:
        ach_to_pids = {}
        for r in rows:
            mid = r.get(model_col, '')
            pid = r.get(prof_col, '')
            dt = r.get(datatype_col, '')
            if mid and pid and dt == 'rna':
                if mid not in ach_to_pids:
                    ach_to_pids[mid] = []
                ach_to_pids[mid].append(pid)
        
        cellline_ach = {}
        for cl in known_profiles:
            for r in rows:
                mid = r.get(model_col, '')
                if cl.lower() in mid.lower():
                    cellline_ach[cl] = mid
                    break
        
        if cellline_ach:
            print(f"  ACH mapping found: {cellline_ach}")
    print()

# Also check expression file
expr_file = os.path.join(CCLE_DIR, "OmicsExpressionGenesExpectedCountProfile.csv")
if os.path.exists(expr_file):
    with open(expr_file, 'r', encoding='utf-8', errors='replace') as f:
        reader = csv.DictReader(f)
        header = reader.fieldnames
        first_col = header[0]
        expr_ids = set()
        for row in reader:
            expr_ids.add(row[first_col])
    
    print(f"=== Expression file ({len(expr_ids)} profiles) ===")
    print(f"  First column: {first_col}")
    for cl, pid in known_profiles.items():
        found = pid in expr_ids
        print(f"  {cl}: {pid} -> {'FOUND' if found else 'NOT FOUND'} in expression data")
    
    # Also check ACH- prefixed IDs
    ach_ids = {k: v for k, v in known_profiles.items()}
    for cl in known_profiles:
        for eid in expr_ids:
            if cl.lower() in eid.lower():
                print(f"  {cl}: found via name match -> {eid}")

print("\n=== TASK_D verification complete ===")
