import matplotlib
matplotlib.use('Agg')  # 无图形界面，直接保存文件
import matplotlib.pyplot as plt

# 读取 RMSD 数据
time_rmsd, rmsd = [], []
with open('rmsd_100ns.xvg') as f:
    for line in f:
        if line.startswith(('#', '@')):
            continue
        parts = line.split()
        time_rmsd.append(float(parts[0]))
        rmsd.append(float(parts[1]))

# 读取 RMSF 数据
resid, rmsf = [], []
with open('rmsf_100ns.xvg') as f:
    for line in f:
        if line.startswith(('#', '@')):
            continue
        parts = line.split()
        resid.append(int(parts[0]))
        rmsf.append(float(parts[1]))

# 创建两个子图
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))

# 左图：RMSD
ax1.plot(time_rmsd, rmsd, color='steelblue', linewidth=1)
ax1.set_xlabel('Time (ns)')
ax1.set_ylabel('RMSD (nm)')
ax1.set_title('RMSD of R1aB6–H1N1 complex')
ax1.grid(alpha=0.3)

# 右图：RMSF
ax2.plot(resid, rmsf, color='coral', linewidth=1)
ax2.set_xlabel('Residue number')
ax2.set_ylabel('RMSF (nm)')
ax2.set_title('RMSF per residue')
ax2.grid(alpha=0.3)

plt.tight_layout()
plt.savefig('rmsd_rmsf_100ns.png', dpi=300)
print('Saved rmsd_rmsf_100ns.png')
