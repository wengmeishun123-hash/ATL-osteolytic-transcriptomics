import csv
import os
import re
from collections import defaultdict

OUT_DIR = r"D:\ATL research\analysis\results\cmap"

def normalize_drug_name(name):
    if not name or name == "-666":
        return None
    name = name.strip().lower()
    name = re.sub(r'[^a-z0-9]', '', name)
    return name if name else None

def parse_enrichr_term(term):
    if not term:
        return None, None, None
    parts = term.split('-')
    compound = None
    dose = None
    cell_line = None
    if len(parts) >= 2:
        header = parts[0].strip()
        header_parts = header.split()
        if len(header_parts) >= 2:
            cell_line = header_parts[1]
        last = parts[-1].strip()
        try:
            dose = float(last)
            if len(parts) >= 3:
                compound = parts[-2].strip()
        except ValueError:
            compound = last
    if compound and compound == "-666":
        compound = None
    return compound, cell_line, dose

def load_enrichr_v2(filepath, direction='down'):
    results = []
    with open(filepath, 'r', encoding='utf-8', errors='replace') as f:
        reader = csv.reader(f)
        header = next(reader)
        for row in reader:
            if len(row) < 5:
                continue
            term = row[1]
            try:
                pval = float(row[2])
                zscore = float(row[3])
                combined = float(row[4])
            except:
                continue
            compound, cell_line, dose = parse_enrichr_term(term)
            if compound and compound != "-666":
                results.append({
                    'compound': compound,
                    'pval': pval,
                    'zscore': zscore,
                    'combined_score': combined,
                    'source': f'Enrichr_LINCS_{direction}',
                    'cell_line': cell_line,
                    'raw_term': term
                })
    return results

def load_l1000cds2(filepath):
    results = []
    with open(filepath, 'r', encoding='utf-8', errors='replace') as f:
        reader = csv.DictReader(f)
        for row in reader:
            name = row.get('pert_desc', '').strip()
            if name == "-666" or not name:
                continue
            score = row.get('score', '0')
            try:
                score = float(score)
            except:
                score = 0.0
            results.append({
                'compound': name,
                'score': score,
                'source': 'L1000CDS2',
                'pert_id': row.get('pert_id', ''),
                'cell_id': row.get('cell_id', ''),
            })
    return results

def load_clue_new(filepath, reversal_only=True, min_sigs=2):
    results = []
    with open(filepath, 'r', encoding='utf-8', errors='replace') as f:
        reader = csv.DictReader(f)
        for row in reader:
            compound = row.get('compound', '').strip()
            if compound == "-666" or not compound:
                continue
            if compound.startswith('BRD-'):
                continue
            try:
                avg_norm_cs = float(row.get('avg_norm_cs', '0'))
                min_norm_cs = float(row.get('min_norm_cs', '0'))
                max_fdr = float(row.get('max_fdr_q_nlog10', '0'))
                n_sigs = int(row.get('n_signatures', '0'))
                n_hiq = int(row.get('n_hiq', '0'))
            except:
                continue
            if n_sigs < min_sigs:
                continue
            if reversal_only and avg_norm_cs >= 0:
                continue
            results.append({
                'compound': compound,
                'avg_norm_cs': avg_norm_cs,
                'min_norm_cs': min_norm_cs,
                'max_fdr_q_nlog10': max_fdr,
                'n_signatures': n_sigs,
                'n_hiq': n_hiq,
                'moa': row.get('moa', ''),
                'target': row.get('target', ''),
                'source': 'clue.io',
            })
    return results

def main():
    print("=" * 80)
    print("CMap Convergence Analysis v4 - ALL L1000-filtered sources")
    print("=" * 80)
    
    all_results = defaultdict(list)
    
    l1000cds2_file = os.path.join(OUT_DIR, "l1000cds2_results.csv")
    if os.path.exists(l1000cds2_file):
        print(f"\n>>> Loading ORIGINAL L1000CDS2 results (server down, using cached)...")
        l1000cds2 = load_l1000cds2(l1000cds2_file)
        print(f"    Loaded {len(l1000cds2)} results")
        for r in l1000cds2:
            norm = normalize_drug_name(r['compound'])
            if norm:
                all_results[norm].append(r)
    
    enrichr_down_v2 = os.path.join(OUT_DIR, "lincs_l1000_chem_pert_down_v2.csv")
    if os.path.exists(enrichr_down_v2):
        print(f"\n>>> Loading NEW Enrichr DOWN results (L1000-filtered)...")
        ed = load_enrichr_v2(enrichr_down_v2, 'down')
        print(f"    Loaded {len(ed)} results")
        for r in ed:
            norm = normalize_drug_name(r['compound'])
            if norm:
                all_results[norm].append(r)
    
    enrichr_up_v2 = os.path.join(OUT_DIR, "lincs_l1000_chem_pert_up_v2.csv")
    if os.path.exists(enrichr_up_v2):
        print(f"\n>>> Loading NEW Enrichr UP results (L1000-filtered)...")
        eu = load_enrichr_v2(enrichr_up_v2, 'up')
        print(f"    Loaded {len(eu)} results")
        for r in eu:
            norm = normalize_drug_name(r['compound'])
            if norm:
                all_results[norm].append(r)
    
    enrichr_dn_from_dn = os.path.join(OUT_DIR, "lincs_l1000_chem_pert_down_from_DN_v2.csv")
    if os.path.exists(enrichr_dn_from_dn):
        print(f"\n>>> Loading NEW Enrichr DOWN-from-DN results (L1000-filtered)...")
        edn = load_enrichr_v2(enrichr_dn_from_dn, 'down_from_DN')
        print(f"    Loaded {len(edn)} results")
        for r in edn:
            norm = normalize_drug_name(r['compound'])
            if norm:
                all_results[norm].append(r)
    
    enrichr_up_from_dn = os.path.join(OUT_DIR, "lincs_l1000_chem_pert_up_from_DN_v2.csv")
    if os.path.exists(enrichr_up_from_dn):
        print(f"\n>>> Loading NEW Enrichr UP-from-DN results (L1000-filtered)...")
        eud = load_enrichr_v2(enrichr_up_from_dn, 'up_from_DN')
        print(f"    Loaded {len(eud)} results")
        for r in eud:
            norm = normalize_drug_name(r['compound'])
            if norm:
                all_results[norm].append(r)
    
    clue_new_file = os.path.join(OUT_DIR, "clue_aggregated_results_v2_new_query.csv")
    if os.path.exists(clue_new_file):
        print(f"\n>>> Loading NEW clue.io results (L1000-filtered query, reversal only)...")
        clue = load_clue_new(clue_new_file, reversal_only=True, min_sigs=2)
        print(f"    Loaded {len(clue)} reversal compounds")
        for r in clue:
            norm = normalize_drug_name(r['compound'])
            if norm:
                all_results[norm].append(r)
    
    print(f"\n>>> Total unique compounds: {len(all_results)}")
    
    convergence = []
    for norm_name, entries in all_results.items():
        raw_sources = set(e['source'] for e in entries)
        canonical_sources = set()
        for s in raw_sources:
            if s.startswith('Enrichr_LINCS'):
                canonical_sources.add('Enrichr')
            else:
                canonical_sources.add(s)
        n_sources = len(canonical_sources)
        
        l1000_scores = []
        enrichr_scores = []
        clue_scores = []
        clue_min_norm_cs = None
        clue_max_fdr = 0
        clue_n_sigs = 0
        best_pval = 1.0
        best_l1000_score = 0.0
        display_name = entries[0]['compound']
        moa = ''
        target = ''
        
        for e in entries:
            if e['source'] == 'L1000CDS2':
                s = e.get('score', 0)
                l1000_scores.append(s)
                if s > best_l1000_score:
                    best_l1002_score = s
            elif e['source'].startswith('Enrichr_LINCS'):
                enrichr_scores.append(e.get('combined_score', 0))
                pval = e.get('pval', 1)
                if pval < best_pval:
                    best_pval = pval
            elif e['source'] == 'clue.io':
                clue_scores.append(e.get('avg_norm_cs', 0))
                if clue_min_norm_cs is None or e.get('min_norm_cs', 0) < clue_min_norm_cs:
                    clue_min_norm_cs = e.get('min_norm_cs', 0)
                if e.get('max_fdr_q_nlog10', 0) > clue_max_fdr:
                    clue_max_fdr = e.get('max_fdr_q_nlog10', 0)
                clue_n_sigs = max(clue_n_sigs, e.get('n_signatures', 0))
                if e.get('moa', ''):
                    moa = e['moa']
                if e.get('target', ''):
                    target = e['target']
        
        avg_l1000 = sum(l1000_scores) / len(l1000_scores) if l1000_scores else None
        best_enrichr = max(enrichr_scores) if enrichr_scores else None
        avg_clue = sum(clue_scores) / len(clue_scores) if clue_scores else None
        clue_reversal = avg_clue is not None and avg_clue < 0
        
        convergence.append({
            'compound': display_name,
            'normalized': norm_name,
            'n_sources': n_sources,
            'sources': '; '.join(sorted(canonical_sources)),
            'l1000cds2_score': best_l1000_score if l1000_scores else None,
            'l1000cds2_n': len(l1000_scores),
            'enrichr_best_score': best_enrichr,
            'enrichr_n': len(enrichr_scores),
            'clue_avg_norm_cs': avg_clue,
            'clue_min_norm_cs': clue_min_norm_cs,
            'clue_max_fdr': clue_max_fdr,
            'clue_n': len(clue_scores),
            'clue_n_sigs': clue_n_sigs,
            'clue_reversal': clue_reversal,
            'best_pval': best_pval,
            'moa': moa,
            'target': target,
            'n_entries': len(entries)
        })
    
    convergence.sort(key=lambda x: (
        -x['n_sources'],
        x['clue_avg_norm_cs'] if x['clue_avg_norm_cs'] is not None else 0,
        x['best_pval'],
    ))
    
    three_src = [c for c in convergence if c['n_sources'] >= 3]
    two_src = [c for c in convergence if c['n_sources'] == 2]
    two_clue_rev = [c for c in two_src if c['clue_reversal']]
    
    print(f"\n{'='*80}")
    print("CONVERGENCE RESULTS (v4 - ALL L1000-filtered)")
    print(f"{'='*80}")
    
    print(f"\n--- 3-source convergence ({len(three_src)} compounds) ---")
    print(f"{'Compound':<25} {'L1000':<8} {'Enrichr':<10} {'ClueAvg':<10} {'ClueMin':<10} {'MOA'}")
    print("-" * 100)
    for c in three_src[:20]:
        l1k = f"{c['l1000cds2_score']:.3f}" if c['l1000cds2_score'] else "-"
        enr = f"{c['enrichr_best_score']:.1f}" if c['enrichr_best_score'] else "-"
        ca = f"{c['clue_avg_norm_cs']:.4f}" if c['clue_avg_norm_cs'] else "-"
        cm = f"{c['clue_min_norm_cs']:.4f}" if c['clue_min_norm_cs'] else "-"
        print(f"{c['compound']:<25} {l1k:<8} {enr:<10} {ca:<10} {cm:<10} {c['moa'][:35]}")
    
    print(f"\n--- 2-source with clue.io reversal (top 30) ---")
    print(f"{'Compound':<25} {'Sources':<20} {'Enrichr':<10} {'ClueAvg':<10} {'MOA'}")
    print("-" * 100)
    for c in two_clue_rev[:30]:
        enr = f"{c['enrichr_best_score']:.1f}" if c['enrichr_best_score'] else "-"
        ca = f"{c['clue_avg_norm_cs']:.4f}" if c['clue_avg_norm_cs'] else "-"
        print(f"{c['compound']:<25} {c['sources']:<20} {enr:<10} {ca:<10} {c['moa'][:35]}")
    
    out_csv = os.path.join(OUT_DIR, "cmap_convergence_results_v4.csv")
    fields = ['compound', 'normalized', 'n_sources', 'sources',
              'l1000cds2_score', 'l1000cds2_n',
              'enrichr_best_score', 'enrichr_n',
              'clue_avg_norm_cs', 'clue_min_norm_cs', 'clue_max_fdr', 'clue_n',
              'clue_n_sigs', 'clue_reversal',
              'best_pval', 'moa', 'target', 'n_entries']
    with open(out_csv, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(convergence)
    print(f"\n>>> Saved to: {out_csv}")
    
    print(f"\n>>> ATL-RELEVANT CANDIDATES:")
    atl_keywords = ['hdac', 'bcr-abl', 'flt3', 'pdgfr', 'pi3k', 'mtor', 'jak', 'stat',
                     'nfkb', 'cdk', 'egfr', 'vegfr', 'src', 'syk', 'btk', 'hsp90', 'fgfr',
                     'parp', 'aurora', 'proteasome', 'bcl', 'mdm', 'wnt', 'gsk', 'mek',
                     'kit', 'topoisomerase', 'dna inhibitor', 'p53', 'pkc', 'raf']
    
    atl_candidates = []
    for c in convergence:
        if c['n_sources'] < 2:
            continue
        compound_lower = c['compound'].lower()
        moa_lower = c.get('moa', '').lower()
        target_lower = c.get('target', '').lower()
        all_text = f"{compound_lower} {moa_lower} {target_lower}"
        if any(kw in all_text for kw in atl_keywords):
            atl_candidates.append(c)
    
    atl_candidates.sort(key=lambda x: (-x['n_sources'], x.get('clue_avg_norm_cs', 0) or 0))
    print(f"{'Compound':<25} {'N_Src':<6} {'L1000':<8} {'Enrichr':<10} {'ClueAvg':<10} {'Rev':<5} {'MOA'}")
    print("-" * 110)
    for c in atl_candidates[:30]:
        l1k = f"{c['l1000cds2_score']:.3f}" if c['l1000cds2_score'] else "-"
        enr = f"{c['enrichr_best_score']:.1f}" if c['enrichr_best_score'] else "-"
        ca = f"{c['clue_avg_norm_cs']:.4f}" if c['clue_avg_norm_cs'] else "-"
        rev = "REV" if c['clue_reversal'] else "SIM"
        print(f"{c['compound']:<25} {c['n_sources']:<6} {l1k:<8} {enr:<10} {ca:<10} {rev:<5} {c['moa'][:40]}")
    
    print(f"\n>>> SUMMARY:")
    print(f"    Total compounds: {len(convergence)}")
    print(f"    3-source: {len(three_src)}")
    print(f"    2-source: {len(two_src)}")
    print(f"    2-source with clue.io reversal: {len(two_clue_rev)}")
    print(f"    ATL-relevant: {len(atl_candidates)}")
    print(f"    NOTE: L1000CDS2 results are from original query (server currently down)")
    print(f"    When L1000CDS2 is available again, re-submit for improved 3-source coverage")

if __name__ == "__main__":
    main()
