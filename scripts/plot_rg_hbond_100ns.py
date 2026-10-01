#!/usr/bin/env python3
"""Plot Rg and H-bond count over 100 ns MD (no titles, for iGEM Wiki)."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def read_xvg(fname):
    x, y = [], []
    with open(fname) as f:
        for line in f:
            if line.startswith(('#', '@')):
                continue
            p = line.split()
            if len(p) >= 2:
                try:
                    x.append(float(p[0]) / 1000.0)  # ps → ns
                    y.append(float(p[1]))
                except ValueError:
                    pass
    return x, y

# --- Rg ---
t_rg, rg = read_xvg('gyrate_100ns.xvg')

# --- H-bonds ---
t_hb, hb = read_xvg('hbonds_100ns.xvg')

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

# Left: Rg
ax1.plot(t_rg, rg, color='steelblue', linewidth=0.8)
ax1.set_xlabel('Time (ns)', fontsize=12)
ax1.set_ylabel('Radius of gyration (nm)', fontsize=12)
ax1.grid(alpha=0.3)

# Right: H-bonds
ax2.plot(t_hb, hb, color='seagreen', linewidth=0.8)
ax2.set_xlabel('Time (ns)', fontsize=12)
ax2.set_ylabel('Number of H-bonds', fontsize=12)
ax2.grid(alpha=0.3)

plt.tight_layout()
plt.savefig('rg_hbond_100ns.png', dpi=300, bbox_inches='tight')
print('Saved rg_hbond_100ns.png')

# Print stats
import statistics
print(f"\nRg 90-100 ns: mean = {statistics.mean([v for t,v in zip(t_rg,rg) if t>=90]):.4f} nm")
print(f"H-bonds full: mean = {statistics.mean(hb):.2f}, SD = {statistics.stdev(hb):.2f}")
