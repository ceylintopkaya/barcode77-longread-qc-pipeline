# -*- coding: utf-8 -*-
import sys
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Komut satiri argumanlarini al
input_file = sys.argv[1]
output_file = sys.argv[2]

# Sample adini cikti dosyasinin ISMINDEN cikar
SAMPLE = os.path.basename(output_file).split('_qc_plots_python.png')[0]
if SAMPLE == "":
    SAMPLE = os.path.basename(os.path.dirname(output_file))

df = pd.read_csv(input_file)

# N50 hesapla
sorted_lengths = sorted(df['length'], reverse=True)
total_bases = sum(sorted_lengths)
cumsum = 0
n50 = 0
for l in sorted_lengths:
    cumsum += l
    if cumsum >= total_bases / 2:
        n50 = l
        break

# Kalite kategorileri
total_reads = len(df)
q10_pct = (df['mean_quality'] >= 10).sum() / total_reads * 100
q20_pct = (df['mean_quality'] >= 20).sum() / total_reads * 100
q30_pct = (df['mean_quality'] >= 30).sum() / total_reads * 100
total_yield_mb = total_bases / 1_000_000

print("=== SUMMARY STATISTICS ===")
print(f"GC Content   - Mean: {df['gc_content'].mean():.2f} | Median: {df['gc_content'].median():.2f}")
print(f"Read Length  - Mean: {df['length'].mean():.2f} | Median: {df['length'].median():.2f} | N50: {n50}")
print(f"Quality Score- Mean: {df['mean_quality'].mean():.2f} | Median: {df['mean_quality'].median():.2f}")
print(f"\nQ10 Above: %{q10_pct:.1f} | Q20 Above: %{q20_pct:.1f} | Q30 Above: %{q30_pct:.1f}")
print(f"Total Yield: {total_yield_mb:.2f} Mb")

sns.set_theme(style="whitegrid")
fig, axes = plt.subplots(2, 2, figsize=(16, 11))
fig.suptitle(f"{SAMPLE.capitalize()} - Long-read QC Report (Python)", fontsize=15, fontweight='bold')

# -----------------------------------------------------------------------------
# GRAFIK 1: READ LENGTH 
# -----------------------------------------------------------------------------
axes[0,0].hist(df['length'], bins=50, color='steelblue', edgecolor='black', alpha=0.7)
axes[0,0].set_xscale('log')
axes[0,0].axvline(df['length'].median(), color='orange', linestyle='dashed', linewidth=1.5, label='Median')
axes[0,0].axvline(n50, color='red', linestyle='dashed', linewidth=1.5, label='N50')
axes[0,0].set_title('Read Length Distribution (Log Scale)')
axes[0,0].set_xlabel('Length (bp) [log scale]')
axes[0,0].set_ylabel('Read Count')
axes[0,0].legend()
stats_text = f"Mean: {df['length'].mean():.1f}\nMedian: {df['length'].median():.1f}\nN50: {n50}\nTotal: {len(df)}"
axes[0,0].text(0.97, 0.97, stats_text, transform=axes[0,0].transAxes, fontsize=8,
               verticalalignment='top', horizontalalignment='right',
               bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

# -----------------------------------------------------------------------------
# GRAFIK 2: QUALITY SCORE 
# -----------------------------------------------------------------------------
axes[0,1].hist(df['mean_quality'], bins=50, color='coral', edgecolor='black', alpha=0.7)
axes[0,1].axvline(df['mean_quality'].median(), color='orange', linestyle='dashed', linewidth=1.5, label='Median')
axes[0,1].axvline(10, color='darkred', linestyle='dotted', linewidth=2, label='Q10')
axes[0,1].axvline(20, color='purple', linestyle='dotted', linewidth=2, label='Q20')
axes[0,1].axvline(30, color='blue', linestyle='dotted', linewidth=2, label='Q30')
axes[0,1].set_title('Mean Quality Score Distribution')
axes[0,1].set_xlabel('Mean Quality Score')
axes[0,1].set_ylabel('Read Count')
axes[0,1].legend()
stats_text = f"Mean: {df['mean_quality'].mean():.1f}\nMedian: {df['mean_quality'].median():.1f}\nQ10+: %{q10_pct:.1f}\nQ20+: %{q20_pct:.1f}\nQ30+: %{q30_pct:.1f}"
axes[0,1].text(0.03, 0.97, stats_text, transform=axes[0,1].transAxes, fontsize=8,
               verticalalignment='top', horizontalalignment='left',
               bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

# -----------------------------------------------------------------------------
# GRAFIK 3: GC CONTENT 
# -----------------------------------------------------------------------------
axes[1,0].hist(df['gc_content'], bins=50, color='seagreen', edgecolor='black', alpha=0.7)
axes[1,0].axvline(df['gc_content'].median(), color='orange', linestyle='dashed', linewidth=1.5, label='Median')
axes[1,0].axvline(52, color='purple', linestyle='dashed', linewidth=1.5, label='Expected GC (52%)')
axes[1,0].set_title('GC Content Distribution')
axes[1,0].set_xlabel('GC Content (%)')
axes[1,0].set_ylabel('Read Count')
axes[1,0].legend()
stats_text = f"Mean: {df['gc_content'].mean():.1f}%\nMedian: {df['gc_content'].median():.1f}%\nTotal: {len(df)}"
axes[1,0].text(0.03, 0.97, stats_text, transform=axes[1,0].transAxes, fontsize=8,
               verticalalignment='top', horizontalalignment='left',
               bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

# -----------------------------------------------------------------------------
# GRAFIK 4: LENGTH VS QUALITY - HEXBIN (YlGnBu skalasi)
# -----------------------------------------------------------------------------
hb = axes[1,1].hexbin(df['length'], df['mean_quality'],
                      gridsize=40, cmap='YlGnBu',
                      xscale='log', bins='log', mincnt=1, alpha=0.8)
plt.colorbar(hb, ax=axes[1,1], label='Read Count (log scale)')
axes[1,1].axhline(10, color='darkred', linestyle='dotted', linewidth=2, label='Q10')
axes[1,1].axhline(20, color='purple', linestyle='dotted', linewidth=2, label='Q20')
axes[1,1].axhline(30, color='blue', linestyle='dotted', linewidth=2, label='Q30')
axes[1,1].set_title('Read Length vs Quality (Density)')
axes[1,1].set_xlabel('Read Length (bp) [log scale]')
axes[1,1].set_ylabel('Mean Quality Score')
axes[1,1].legend(loc='lower right')

# Ultra-uzun okumayi isaretle (varsa)
max_idx = df['length'].idxmax()
max_length = df.loc[max_idx, 'length']
if max_length > 100000:  # 100kb'den buyukse
    axes[1,1].plot(max_length, df.loc[max_idx, 'mean_quality'],
                   'r*', markersize=15, label=f'Max: {max_length:,.0f} bp')
    axes[1,1].legend()

plt.tight_layout()
plt.savefig(output_file, dpi=150, bbox_inches='tight')
print(f"Plot saved as {output_file}!")