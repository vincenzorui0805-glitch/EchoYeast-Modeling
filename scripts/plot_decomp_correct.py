#!/usr/bin/env python3
"""
Per-residue MM-GBSA decomposition - only first 51 frames (correct).
Complex/Ligand 各取第一组 51 帧，按帧号配对相减。
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

CSV = 'FINAL_DECOMP_100ns.csv'
lines = open(CSV).read().splitlines()

def find_idx(marker):
    for i, l in enumerate(lines):
        if l.strip() == marker:
            return i
    return None

c_s = find_idx('Complex:') + 1
r_s = find_idx('Receptor:')
l_s = find_idx('Ligand:') + 1

def parse_first_group(block):
    """每个残基只取前 51 行（第一组），返回 {residue: {frame: total}}"""
    d = {}
    for line in block:
        p = [x for x in line.split(',') if x.strip()]
        if len(p) >= 8:
            try:
                frame = int(p[0])
                res   = p[1]
                total = float(p[-1])
                if len(d.setdefault(res, {})) < 51:
                    d[res][frame] = total
            except Exception:
                pass
    return d

cplx = parse_first_group(lines[c_s:r_s])
recp = parse_first_group(lines[r_s:l_s])
lig  = parse_first_group(lines[l_s:])

def mean_diff(a, b):
    """按帧号配对相减，取平均"""
    out = {}
    for res in a:
        if res not in b:
            continue
        frames = sorted(set(a[res]) & set(b[res]))
        if not frames:
            continue
        out[res] = float(np.mean([a[res][f] - b[res][f] for f in frames]))
    return out

rec_contrib = {r: v for r, v in mean_diff(cplx, recp).items() if r.startswith('R:')}
lig_contrib = {r: v for r, v in mean_diff(cplx, lig).items()  if r.startswith('L:')}

rec_items = sorted(rec_contrib.items(), key=lambda kv: kv[1])
lig_items = sorted(lig_contrib.items(), key=lambda kv: kv[1])
rec_labels, rec_vals = zip(*rec_items)
lig_labels, lig_vals = zip(*lig_items)

# ================== 画图 ==================
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 8),
                               gridspec_kw={'width_ratios': [1.3, 1]})

# ----- (a) 受体 HA -----
x1 = np.arange(len(rec_vals))
c1 = ['#c0392b' if v < -0.5 else '#e67e22' if v < -0.2 else '#3498db'
      for v in rec_vals]
ax1.bar(x1, rec_vals, color=c1, edgecolor='black', linewidth=0.4)
ax1.axhline(0, color='black', linewidth=0.8)
ax1.set_xticks(x1)
ax1.set_xticklabels([l.replace('R:A:', '') for l in rec_labels],
                    rotation=90, fontsize=8)
ax1.set_ylabel('MM-GBSA contribution (kcal/mol)', fontsize=12)
ax1.set_title('(a) Receptor (HA) per-residue contribution',
              fontsize=13, fontweight='bold')
ax1.grid(axis='y', alpha=0.3)
for i in range(min(5, len(rec_vals))):
    ax1.text(i, rec_vals[i] - 0.02, f'{rec_vals[i]:.2f}',
             ha='center', va='top', fontsize=8, fontweight='bold')

# ----- (b) 配体 纳米抗体 -----
x2 = np.arange(len(lig_vals))
c2 = []
for lab, v in zip(lig_labels, lig_vals):
    if 'ARG:31' in lab:
        c2.append('#27ae60')       # 热点绿色
    elif v < -3:
        c2.append('#c0392b')
    elif v < -1.5:
        c2.append('#e67e22')
    else:
        c2.append('#3498db')
ax2.bar(x2, lig_vals, color=c2, edgecolor='black', linewidth=0.4)
ax2.axhline(0, color='black', linewidth=0.8)
ax2.set_xticks(x2)
ax2.set_xticklabels([l.replace('L:B:', '') for l in lig_labels],
                    rotation=90, fontsize=8)
ax2.set_ylabel('MM-GBSA contribution (kcal/mol)', fontsize=12)
ax2.set_title('(b) Ligand (nanobody) per-residue contribution',
              fontsize=13, fontweight='bold')
ax2.grid(axis='y', alpha=0.3)
for i, lab in enumerate(lig_labels):
    if 'ARG:31' in lab:
        ax2.text(i, lig_vals[i] - 0.08, f'ARG-B:31\n{lig_vals[i]:.2f}',
                 ha='center', va='top', fontsize=9,
                 fontweight='bold', color='#27ae60')
        break
for i in range(min(5, len(lig_vals))):
    if 'ARG:31' in lig_labels[i]:
        continue
    ax2.text(i, lig_vals[i] - 0.08, f'{lig_vals[i]:.2f}',
             ha='center', va='top', fontsize=8)

plt.tight_layout()
plt.savefig('per_residue_decomp_100ns_correct.png', dpi=300, bbox_inches='tight')
print('Saved per_residue_decomp_100ns_correct.png')

print('\n=== Receptor (HA) Top 10 ===')
for lab, v in zip(rec_labels[:10], rec_vals[:10]):
    print(f'  {v:9.3f}  {lab}')

print('\n=== Ligand (nanobody) Top 10 ===')
for lab, v in zip(lig_labels[:10], lig_vals[:10]):
    print(f'  {v:9.3f}  {lab}')

print('\n=== ARG-B:31 ===')
for lab, v in zip(lig_labels, lig_vals):
    if 'ARG:31' in lab:
        print(f'  {v:9.3f}  {lab}')
