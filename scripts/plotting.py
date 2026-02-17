from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.formula.api as smf
import seaborn as sns
import ptitprince as pt

basedir = Path(__file__).parent.parent
pathin = basedir / "data"
pathout = basedir / "nsp"

# Load data for distribution plots
df_balance_simple = pd.read_csv(pathout / "nsp_balance_zscores.csv")

# Get all FC columns from df_balance_simple
fc_columns = df_balance_simple.columns.tolist()
n_cols = len(fc_columns)

# Create subplots for distributions
fig, axes = plt.subplots(2, 3, figsize=(15, 10))
axes = axes.flatten()

# Plot each FC type
for idx, col in enumerate(fc_columns):
    axes[idx].hist(df_balance_simple[col], bins=20, color='skyblue', edgecolor='black', alpha=0.7)
    axes[idx].set_title(col, fontsize=12, fontweight='bold')
    axes[idx].set_xlabel('Balance Coefficient (z-scores)')
    axes[idx].set_ylabel('Frequency')
    axes[idx].grid(axis='y', alpha=0.3)

# Hide any unused subplots if you have fewer than 6 FC types
for idx in range(n_cols, len(axes)):
    axes[idx].axis('off')

plt.tight_layout()
plt.show()

# ===== Twin Pair Comparison Plot =====

# Load corrected wide format data with complete twin pairs
df_wide = pd.read_excel(pathin / "balance_rs.xlsx")

# Reshape from wide to long format for plotting
twin_pairs_list = []
for idx, row in df_wide.iterrows():
    twin_pairs_list.append({
        'Family_ID': row['Family_ID'],
        'zyg': row['zyg'],
        'Subject': f"{row['Family_ID']}_T1",
        'FC_rs': row['FC_rs_T1']
    })
    twin_pairs_list.append({
        'Family_ID': row['Family_ID'],
        'zyg': row['zyg'],
        'Subject': f"{row['Family_ID']}_T2",
        'FC_rs': row['FC_rs_T2']
    })

twin_pairs = pd.DataFrame(twin_pairs_list).groupby('Family_ID')

# Create figure
fig, ax = plt.subplots(figsize=(14, 7))

# Track x position for plotting
x_pos = 0
pair_labels = []
pair_positions = []

# Define colors for MZ and DZ
colors = {'MZ': 'blue', 'DZ': 'red'}

# Plot each twin pair
for family_id, group in twin_pairs:
    zygosity = group['zyg'].iloc[0]
    color = colors.get(zygosity, 'gray')
    
    # Get FC_rs values for both twins
    hbal_values = group['FC_rs'].values
    
    # Plot both twins at the same x position
    ax.scatter([x_pos, x_pos], hbal_values, s=100, color=color, alpha=0.7, edgecolors='black', linewidth=1.5)
    
    # Draw a line connecting the twins
    ax.plot([x_pos, x_pos], hbal_values, color=color, alpha=0.4, linewidth=2)
    
    pair_labels.append(f"{family_id}\n({zygosity})")
    pair_positions.append(x_pos)
    x_pos += 1

# Customize plot
ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5, linewidth=1)
ax.set_xlabel('Twin Pairs', fontsize=12, fontweight='bold')
ax.set_ylabel('FC_rs (z-scores)', fontsize=12, fontweight='bold')
ax.set_title('Twin Pair Comparison of FC_rs (124 Complete Pairs)', fontsize=14, fontweight='bold')
ax.set_xticks(pair_positions[::5])  # Show every 5th label to avoid crowding
ax.set_xticklabels(pair_labels[::5], rotation=45, ha='right', fontsize=9)
ax.grid(axis='y', alpha=0.3)

# Add legend
from matplotlib.patches import Patch
legend_elements = [
    Patch(facecolor='blue', edgecolor='black', alpha=0.7, label='MZ (Monozygotic)'),
    Patch(facecolor='red', edgecolor='black', alpha=0.7, label='DZ (Dizygotic)')
]
ax.legend(handles=legend_elements, loc='best', fontsize=11)

plt.tight_layout()
plt.show()

# ===== Raincloud Plot for MZ vs DZ =====

# Create data for raincloud plot from corrected wide format
fc_type = 'FC_rs'

# Extract MZ and DZ data from the corrected wide format
mz_indices = df_wide[df_wide['zyg'] == 'MZ'].index
dz_indices = df_wide[df_wide['zyg'] == 'DZ'].index

# Combine both twin values for each group
mz_data = np.concatenate([df_wide.loc[mz_indices, 'FC_rs_T1'].values, 
                          df_wide.loc[mz_indices, 'FC_rs_T2'].values])
dz_data = np.concatenate([df_wide.loc[dz_indices, 'FC_rs_T1'].values, 
                          df_wide.loc[dz_indices, 'FC_rs_T2'].values])

# Create raincloud plot
fig, ax = plt.subplots(figsize=(12, 6))

# Define positions and colors
positions = [1, 0]
colors_dict = {0: 'blue', 1: 'red'}
data_list = [mz_data, dz_data]
labels = ['MZ (Monozygotic)', 'DZ (Dizygotic)']

# Plot each group
for i, data in enumerate(data_list):
    # Add individual points (rain) with jitter at the bottom
    jitter = np.random.normal(positions[i], 0.04, size=len(data))
    ax.scatter(data, jitter, alpha=0.5, s=60, color=colors_dict[i], edgecolors='black', linewidth=0.5)
    
    # Add density plot (cloud) using KDE extending upward
    from scipy.stats import gaussian_kde
    kde = gaussian_kde(data)
    x_range = np.linspace(data.min() - 1, data.max() + 1, 200)
    density = kde(x_range)
    
    # Scale density for visualization
    density_scaled = density / density.max() * 0.3
    
    # Plot cloud (density curve) extending upward from the point row
    ax.fill_between(x_range, positions[i], positions[i] + density_scaled, 
                     alpha=0.4, color=colors_dict[i], label=labels[i])
    
    # Add median line
    median = np.median(data)
    ax.vlines(median, positions[i], positions[i] + 0.15, colors=colors_dict[i], 
              linewidth=3, alpha=0.8)

# Customize plot
ax.axvline(x=0, color='gray', linestyle='--', alpha=0.5, linewidth=1)
ax.set_xlabel(f'hbal Coefficient (z-scores) - {fc_type}', fontsize=12, fontweight='bold')
ax.set_ylabel('Zygosity', fontsize=12, fontweight='bold')
ax.set_title(f'Raincloud Plot: hbal Distribution by Zygosity ({fc_type})', fontsize=14, fontweight='bold')
ax.set_yticks([1, 0])
ax.set_yticklabels(labels, fontsize=11)
ax.grid(axis='x', alpha=0.3)
ax.legend(loc='best', fontsize=11)

# Add sample sizes
ax.text(ax.get_xlim()[0] + 0.2, positions[0], f'n={len(mz_data)}', 
        va='center', fontsize=10, fontweight='bold')
ax.text(ax.get_xlim()[0] + 0.2, positions[1], f'n={len(dz_data)}', 
        va='center', fontsize=10, fontweight='bold')

plt.tight_layout()
plt.show()


# ===== Intra-Pair Difference Boxplot (MZ vs DZ) =====

# Compute absolute differences from corrected wide format
pair_diffs = df_wide.copy()
pair_diffs['abs_diff'] = np.abs(pair_diffs['FC_rs_T1'] - pair_diffs['FC_rs_T2'])
pair_diffs['Zygosity'] = pair_diffs['zyg']

print("Number of valid twin pairs:", len(pair_diffs))


plt.figure(figsize=(7,6))

sns.boxplot(
    data=pair_diffs,
    x='Zygosity',
    y='abs_diff',
    palette={'MZ': 'blue', 'DZ': 'red'},
    width=0.5
)

sns.swarmplot(
    data=pair_diffs,
    x='Zygosity',
    y='abs_diff',
    color='black',
    alpha=0.6
)

plt.title(f'Intra-Pair Absolute Differences: {fc_type}', fontsize=14, fontweight='bold')
plt.xlabel('Twin Type', fontsize=12)
plt.ylabel('Absolute Difference in hbal (z-score)', fontsize=12)
plt.grid(axis='y', alpha=0.3)

plt.tight_layout()
plt.show()

# ===== "multilevel plot" =====
df_multilevel = df_wide.copy()
sns.lmplot(x='FC_rs_T1', y='FC_rs_T2', hue='zyg', data=df_multilevel)
plt.title("pairwise linear regression by zygosity (124 Complete Pairs)", fontsize=14, fontweight='bold')

# ===== Correlation Plot by Zygosity (Side-by-side) =====
from scipy.stats import pearsonr

# Separate by zygosity
mz_data_corr = df_wide[df_wide['zyg'] == 'MZ']
dz_data_corr = df_wide[df_wide['zyg'] == 'DZ']

# Calculate correlations and fit lines
mz_corr, mz_pval = pearsonr(mz_data_corr['FC_rs_T1'], mz_data_corr['FC_rs_T2'])
dz_corr, dz_pval = pearsonr(dz_data_corr['FC_rs_T1'], dz_data_corr['FC_rs_T2'])

# Fit lines for regression
mz_z = np.polyfit(mz_data_corr['FC_rs_T1'], mz_data_corr['FC_rs_T2'], 1)
mz_p = np.poly1d(mz_z)

dz_z = np.polyfit(dz_data_corr['FC_rs_T1'], dz_data_corr['FC_rs_T2'], 1)
dz_p = np.poly1d(dz_z)

# Create figure with two subplots
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# MZ plot
axes[0].scatter(mz_data_corr['FC_rs_T1'], mz_data_corr['FC_rs_T2'], alpha=0.6, s=80, color='#1f77b4', edgecolors='black', linewidth=0.5)
x_range = np.linspace(mz_data_corr['FC_rs_T1'].min(), mz_data_corr['FC_rs_T1'].max(), 100)
axes[0].plot(x_range, mz_p(x_range), 'r-', linewidth=2, label='Regression line')
axes[0].set_xlabel('Balance coefficient Twin pair 1 (z-score)', fontsize=12, fontweight='bold')
axes[0].set_ylabel('Balance coefficient Twin pair 2 (z-score)', fontsize=12, fontweight='bold')
axes[0].set_title(f'Monozygotes (MZ)\nN = {len(mz_data_corr)}, r = {mz_corr:.4f}', fontsize=13, fontweight='bold')
axes[0].grid(True, alpha=0.3)
axes[0].legend()

# DZ plot
axes[1].scatter(dz_data_corr['FC_rs_T1'], dz_data_corr['FC_rs_T2'], alpha=0.6, s=80, color='#ff7f0e', edgecolors='black', linewidth=0.5)
x_range = np.linspace(dz_data_corr['FC_rs_T1'].min(), dz_data_corr['FC_rs_T1'].max(), 100)
axes[1].plot(x_range, dz_p(x_range), 'r-', linewidth=2, label='Regression line')
axes[1].set_xlabel('Balance coefficient Twin pair 1 (z-score)', fontsize=12, fontweight='bold')
axes[1].set_ylabel('Balance coefficient Twin pair 2 (z-score)', fontsize=12, fontweight='bold')
axes[1].set_title(f'Dizygotes (DZ)\nN = {len(dz_data_corr)}, r = {dz_corr:.4f}', fontsize=13, fontweight='bold')
axes[1].grid(True, alpha=0.3)
axes[1].legend()

plt.tight_layout()
plt.show()