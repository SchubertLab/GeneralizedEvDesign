# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.4.2
#   kernelspec:
#     display_name: Python [conda env:genev]
#     language: python
#     name: conda-env-genev-py
# ---

# %load_ext autoreload
# %autoreload 2

# +
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib as mpl
import pandas as pd
from scipy.stats import pearsonr, spearmanr
import seaborn as sns
from scipy.cluster import hierarchy
from scipy.spatial.distance import pdist
from collections import defaultdict
from matplotlib import gridspec
from scipy.stats import pearsonr, spearmanr
from scipy import stats
from tqdm.auto import tqdm

# use LaTeX fonts in the plot
# https://ercanozturk.org/2017/12/16/python-matplotlib-plots-in-latex/
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
# -

FIG_FORMAT = 'pdf'
FIG_ORDER = {f: i + 3 for i, f in enumerate([
    'pareto', 'cocktail', 'advantage', 'entropy', 'bound', 'immunogens', 'shuffled', 'bound_igs'
])}

# # comparison with fischer

xx = np.random.random((10, 2))
yy = np.array([1, -1])

# +
base_dir = 'experiments/results/'
data_frames = []
for dd in os.listdir(base_dir):
    if not dd.startswith('nef-300-'):
        continue
    
    for df_name in ['fischer', 'fischer-val', 'genev', 'genev-val']:
        df = pd.read_csv('%s/%s/mosaic-%s-vaccine-evaluation.csv' % (base_dir, dd, df_name))
        df['name'] = df_name
        data_frames.append(df)

df = pd.concat(data_frames).reset_index()
df.head()


# -

def compare(df, feature, name1, name2):
    data1 = df[df.name == name1][feature].values
    data2 = df[df.name == name2][feature].values
    diff = data2 - data1
    tt = np.sqrt(len(diff)) * diff.mean() / diff.std(ddof=1)
    if np.isfinite(tt):
        pp = stats.t(df=len(diff) - 1).cdf(tt)
        if pp > 0.5:
            pp = 1 - pp
    else:
        tt, pp = 0, 0.5
    
    return np.mean(data1), tt, pp, np.mean(diff), np.std(diff)


features = ['norm_prot_coverage', 'conservation', 'rel_pop_coverage', 'immunogen']

for val in ['', '-val']:
    print('\ncomparing fischer and genev on %s data' % ('training' if not val else 'validation'))
    
    for feature in features:
        vv, tt, pp, dd, ss = compare(df, feature, 'genev', 'fischer')
        print('    feature: %20s | value: %.3f difference: %+.3f (%+.3f) - t value: %+.3f - p value: %.3f' % (feature, vv, dd, ss, tt, pp))

for method in ['fischer', 'genev']:
    print('\ncomparing %s on training and validation data' % method)
    
    for feature in features:
        vv, tt, pp, dd, ss = compare(df, feature, method, method + '-val')
        print('    feature: %20s | value: %.3f difference: %+.3f (%+.3f) - t value: %+.3f - p value: %.3f' % (feature, vv, dd, ss, tt, pp))

# # plots

sns.set()
sns.set_style('whitegrid')

colorblindfriendly = [ '#1b9e77', '#d95f02', '#7570b3']
sns.palplot(colorblindfriendly)


# +
def to_hex(r, g, b):
    return '#%02x%02x%02x' % (r,g,b)

to_hex(255,150, 9)
# -

theme_dark = ['#3344a7', '#ff9609', '#d22d2f']
sns.palplot(theme_dark)

theme_light = [ '#a8d6ff', '#ffca68', '#e35e7e']
sns.palplot(theme_light)


def set_font_size(font_size):
    plt.rc('font', size=font_size)          # controls default text sizes
    plt.rc('axes', titlesize=font_size)     # fontsize of the axes title
    plt.rc('axes', labelsize=font_size)     # fontsize of the x and y labels
    plt.rc('xtick', labelsize=font_size)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=font_size)    # fontsize of the tick labels
    plt.rc('legend', fontsize=font_size)    # legend fontsize
    plt.rc('figure', titlesize=font_size)   # fontsize of the figure title


plt.rc('axes', lw=0.8)
plt.rc('grid', lw=0.8)
plt.rc('lines', lw=0.8, markersize=4)
plt.rc('ytick.major', size=3, width=0.6)


# # bound plot (fig. 7)

# +
def get_aggregate_mean_std(immunogen):
    df = pd.read_csv(f'experiments/results/ig-bound-{immunogen}-evaluation-aggregate.csv')
    cols = df.columns.tolist()
    cols[0] = 'metric'
    df.columns = cols

    vaxmap = {
        'coverage': 'Pathogen\ncoverage',
        'norm_prot_coverage': 'Pathogen\ncoverage',
        'immunogen': 'Immunogenicity',
        'conservation': 'Conservation',
        'rel_pop_coverage': 'Population\ncoverage',
    }

    pivot_mean = df[
        df.metric.isin(['conservation', 'immunogen', 'norm_prot_coverage', 'rel_pop_coverage'])
    ].pivot_table('mean', 'aa', ['metric', 'vax']).rename(columns=vaxmap)

    pivot_std = df[
        df.metric.isin(['conservation', 'immunogen', 'norm_prot_coverage', 'rel_pop_coverage'])
    ].pivot_table('std', 'aa', ['metric', 'vax']).rename(columns=vaxmap)

    max_ig = pivot_mean['Immunogenicity'].max().max()
    pivot_mean['Immunogenicity'] /= max_ig
    pivot_std['Immunogenicity'] /= max_ig
    
    return max_ig, pivot_mean, pivot_std


max_ig, pivot_mean, pivot_std = get_aggregate_mean_std('netmhcpan')
pivot_mean.T

# +
set_font_size(6)
plt.rc('legend', title_fontsize=8)

fig = plt.figure(figsize=(5, 5 / 4.), dpi=300)

axes = fig.subplots(1, 4, gridspec_kw={
    'height_ratios': [1],
    'width_ratios': [1, 1, 1, 1],
    'hspace': 0.05,
    'wspace': 0.1,
    'left': 0,
    'bottom': 0,
    'right': 1,
    'top': 1,
})


letters = ['(a) ', '(b) ', '(c) ', ' (d) ']
for i, (ax, col) in enumerate(zip(axes, pivot_mean.columns.levels[0])):
    pivot_mean[col][
        # ensure same ordering
        pivot_mean.columns.levels[1]
    ].plot.line(
        ax=ax, marker='.', yerr=pivot_std[col], color=theme_dark
    )

    ax.set_ylim(0, 1.05)
    ax.set_xlim(40, 140)

    # title
    if col == 'Immunogenicity':
        ax.set_title(letters[i] + 'Immunogenicity\n100\%% = %.3f' % max_ig)
    else:
        ax.set_title(letters[i] + col)

    # y axis
    ax.set_yticks([0., 0.25, 0.5, 0.75, 1.0])
    ax.set_yticklabels([] if i > 0 else [f'{i}\%' for i in range(0, 101, 25)])

    # x axis
    ax.set_xlabel('')
    ax.set_xticks([45, 72, 90, 135])

    # legend
    if i == 3:
        ax.legend(loc='lower right', title='Optimize')
    else:
        ax.get_legend().remove()

    ax.grid(False, axis='x')
    sns.despine(fig, ax)

fig.suptitle('Length in amino acids', y=-0.25)
plt.savefig(f'plots/Fig{FIG_ORDER["bound"]}.{FIG_FORMAT}', bbox_inches='tight')
# -

# # bound plot with different immunogenicities (fig. 10)

# +
pretty_ig = dict(zip([
    'netmhcpan', 'netmhcpan-rank', 'pickpocket', 'mhcflurry'
], ['NetMHCpan (IC$_{50}$)', 'NetMHCpan (\%rank)', 'PickPocket (IC$_{50}$)', 'MHCflurry (IC$_{50}$)']))

bound_mean_data, bound_std_data = [], []
for ig in ['netmhcpan', 'netmhcpan-rank', 'pickpocket', 'mhcflurry']:
    _, pivot_mean, pivot_std = get_aggregate_mean_std(ig)
    
    # here we add a new level to the index
    # https://stackoverflow.com/a/42094658
    bound_mean_data.append(pd.concat([pivot_mean.T], keys=[pretty_ig[ig]], names=['Immunogenicity']))
    bound_std_data.append(pd.concat([pivot_std.T], keys=[pretty_ig[ig]], names=['Immunogenicity']))

bound_mean_data = pd.concat(bound_mean_data).T
bound_std_data = pd.concat(bound_std_data).T

bound_mean_data

# +
differences = pd.DataFrame([
    pd.Series(np.abs((mean1.values - mean2.values).ravel()))
    for i1 in range(len(bound_mean_data.columns.levels[0]) - 1)
    for i2 in range(i1 + 1, len(bound_mean_data.columns.levels[0]))
    for m in bound_mean_data.columns.levels[2]
    for mean1, mean2 in [[
        bound_mean_data.loc[:, (bound_mean_data.columns.levels[0][i1], slice(None), m)],
        bound_mean_data.loc[:, (bound_mean_data.columns.levels[0][i2], slice(None), m)],
    ]]
])

differences.index = pd.MultiIndex.from_tuples([
    ('{}\n{}'.format(
        bound_mean_data.columns.levels[0][i1],
        bound_mean_data.columns.levels[0][i2],
    ), m)
    for i1 in range(len(bound_mean_data.columns.levels[0]) - 1)
    for i2 in range(i1 + 1, len(bound_mean_data.columns.levels[0]))
    for m in bound_mean_data.columns.levels[2]
])

#differences = differences.T
#differences.head()
differences
# -

pd.DataFrame(differences.values.ravel()).describe().T

# +
set_font_size(6)
plt.rc('legend', title_fontsize=6, fontsize=5)

fig = plt.figure(figsize=(5, 3), dpi=300)

root_gs = mpl.gridspec.GridSpec(
    1, 2, figure=fig, width_ratios=[2, 1],
)

left_gs = mpl.gridspec.GridSpecFromSubplotSpec(
    2, 2, subplot_spec=root_gs[0], hspace=0.3, wspace=0.1
)

ax = fig.add_subplot(root_gs[1])
for i, idx in enumerate(differences.index):
    box = ax.boxplot(differences.loc[idx], vert=False, patch_artist=True, positions=[i], sym='.')
    for key in box:
        for item in box[key]:
            item.set_color(theme_dark[i % 3])

ax.yaxis.tick_right()
ax.set_yticks([
    y - 0.5 for y in range(0, len(differences), len(differences.index.levels[1]))
], minor=True)
ax.set_yticks(range(1, len(differences), len(differences.index.levels[1])), minor=False)
ax.set_yticklabels(differences.index.levels[0])
ax.grid(True, axis='y', which='minor')
ax.grid(False, axis='y', which='major')
ax.set_xticklabels([f'{x * 100:.0f}\%' for x in ax.get_xticks()])
ax.set_title('(e) Absolute difference')
ax.set_xlim(ax.get_xlim()[::-1])
sns.despine(ax=ax, left=True, right=False, top=True, bottom=False)


axes = [fig.add_subplot(left_gs[i, j]) for i in range(2) for j in range(2)]

letters = ['(a) ', '(b) ', '(c) ', ' (d) ']
linestyles = ['-', '--', ':', '-.']
markerstyles = ['o', '^', 's', 'v']
for ls, ms, ig in zip(linestyles, markerstyles, bound_mean_data.columns.levels[0]):
    for i, (ax, col) in enumerate(zip(axes, bound_mean_data.columns.levels[1])):
        bound_mean_data[ig][col][
            # keep same ordering
            bound_mean_data.columns.levels[2]
        ].plot.line(ax=ax, marker=ms, yerr=pivot_std[col], color=theme_dark, linestyle=ls, markersize=2)

        ax.set_ylim(0, 1.05)
        ax.set_xlim(40, 140)
        
        ax.set_title(letters[i] + col)

        # y axis
        ax.set_yticks([0., 0.25, 0.5, 0.75, 1.0])
        ax.set_yticklabels([] if i == 1 or i == 3 else [f'{i}\%' for i in range(0, 101, 25)])

        # x axis
        ax.set_xlabel('Length in amino acids' if i > 1 else '')
        ax.set_xticks([45, 72, 90, 135] if i > 1 else [])

        ax.get_legend().remove()

        ax.grid(False, axis='x')
        sns.despine(ax=ax)


axes[-1].legend([
    mpl.patches.Patch(alpha=0)  # fake patch to hold the title
] + [
    plt.Line2D((-1, -2), (-1, -2), c='k', marker=ms, markersize=2, linestyle=ls)
    for ms, ls in zip(markerstyles, linestyles)
] + [
    mpl.patches.Patch(alpha=0)  # fake patch to hold the title
] + [
    plt.Line2D((-1, -2), (-1, -2), c=c) for c in theme_dark
], [
    r'\footnotesize Immunogenicity'
] + [
    ig for ig in bound_mean_data.columns.levels[0]
] + [
    r'\footnotesize Optimize'
] + [
    opt.replace('\n', ' ') for opt in bound_mean_data.columns.levels[1]
], loc='lower right', ncol=2)

fig.tight_layout()
plt.savefig(f'plots/Fig{FIG_ORDER["bound_igs"]}.{FIG_FORMAT}', bbox_inches='tight')
# -

# # pareto plot (fig. 3)

dfs = []
for i in range(1, 6):
    df = pd.read_csv('experiments/results/nef-300-%d/made-tradeoff.csv' % i)
    #df = pd.read_csv('/mnt/edo/home/edo/phd/extern/repos/GeneralizedEvDesign/experiments/results/nef-300-%d/made-tradeoff.csv' % i)
    df['rep'] = i
    dfs.append(df)
df = pd.concat(dfs)
df.head()

# +
set_font_size(8)
fig, ax = plt.subplots(figsize=(3, 2), dpi=300)
cm = plt.get_cmap('RdBu_r')
vmin, vmax = 0.1, 1.2
for rep in df.rep.unique():
    data = df[df.rep == rep].values[:, :2]
    for i in range(len(data) - 1):
        (y1, x1), (y2, x2) = data[i], data[i + 1]
        slope = ((y2 - y1) / 3) / ((x2 - x1) / 40)
        color = 255 * (np.arctan(slope) - vmin) / (vmax - vmin)
        ax.plot([-x1, -x2], [y1, y2], '.-', color=theme_dark[0])


plt.xlabel('Cleavage Score')
plt.ylabel('Immunogenicity')
plt.subplots_adjust(hspace=0, wspace=0, left=0, bottom=0, right=1, top=1)
#sns.despine()

plt.savefig(f'plots/Fig{FIG_ORDER["pareto"]}.{FIG_FORMAT}', bbox_inches='tight')
# -

# # advantage (fig. 5)

df = pd.read_csv('experiments/results/advantage-evaluation-summary.csv')
df.head()

aaa = df.columns.tolist()
aaa[0] = 'metric'
df.columns = aaa

df['vax'] = df.apply(lambda x: x['vax'] if '-' not in x['size'] else x['vax'] + '-' + x['size'].split('-')[1], axis=1)
df['size'] = df['size'].apply(lambda x: int(x.split('-')[0]))

df.head()

df[(df.vax == 'OptiTope') & (df.metric == 'norm_prot_coverage')].sort_values('size')

df = df.replace({
    'optitope': 'EpiMix',
    'mosaic-o4': '4-mosaic',
    'mosaic-o8': '8-mosaic',
})

# +
set_font_size(8)

fig = plt.figure(figsize=(5, 2.5), dpi=300)
#gs = gridspec.GridSpec(2, 2, hspace=0.025, wspace=0.025, left=0, bottom=0, right=1, top=1)
gs = gridspec.GridSpec(2, 2, hspace=0.15, wspace=0.05, left=0, bottom=0, right=1, top=1)
ax1 = plt.subplot(gs[0])
ax2 = plt.subplot(gs[1])
ax3 = plt.subplot(gs[2])
ax4 = plt.subplot(gs[3])

#fig, (ax1, ax2, ax3) = plt.subplots(3, 1)

for i, vax in enumerate(df.vax.unique()):
    data = df[(df.vax == vax) & (df['size'] < 1000)].sort_values('size')
    
    data_ig = data[data.metric == 'immunogen']
    ax1.errorbar(data_ig['size'], data_ig['mean'], yerr=data_ig['std'], fmt='.-', label=vax, color=theme_dark[i])
    
    data_eps = data[data.metric == 'rel_pop_coverage']
    ax2.errorbar(data_eps['size'], data_eps['mean'], yerr=data_eps['std'], fmt='.-', label=vax, color=theme_dark[i])
    
    data_cov = data[data.metric == 'norm_prot_coverage']
    ax3.errorbar(data_cov['size'], data_cov['mean'], yerr=data_cov['std'], fmt='.-', label=vax, color=theme_dark[i])
    
    data_cons = data[data.metric == 'conservation']
    ax4.errorbar(data_cons['size'], data_cons['mean'], yerr=data_cons['std'], fmt='.-', label=vax, color=theme_dark[i])


xticks = range(0, 800, 90)
ax1.set_xlim(0, 750)
ax1.set_xticks(xticks)
ax2.set_xlim(0, 750)
ax2.set_xticks(xticks)
ax3.set_xlim(0, 750)
ax3.set_xticks(xticks)
ax4.set_xlim(0, 750)
ax4.set_xticks(xticks)

ax1.tick_params(axis='x', labelbottom=False)
ax2.tick_params(axis='x', labelbottom=False)

ax1.yaxis.set_ticks([0, 10, 20, 30, 40, 50])
ax2.yaxis.tick_right()
ax2.yaxis.set_ticks([0.9, 0.95, 1.0])
ax2.yaxis.set_ticklabels(['90\%', '95\%', '100\%'])
ax2.yaxis.set_label_coords(1.22, 0.5)
ax3.yaxis.set_ticks([0, 0.25, 0.50, 0.75, 1.0])
ax3.yaxis.set_ticklabels(['0\%', '25\%', '50\%', '75\%', '100\%'])
ax3.set_ylim(-0.05, 1.1)
ax4.yaxis.tick_right()
ax4.set_ylim(0, 0.09)
ax4.yaxis.set_ticklabels(['0.0\%', '2.5\%', '5.0\%', '7.5\%', '10.0\%'])
ax4.yaxis.set_label_coords(1.22, 0.5)

ax3.set_xlabel('Length in Aminoacids')
ax4.set_xlabel('Length in Aminoacids')

ax1.set_ylabel('Immunogenicity')
ax2.set_ylabel('Pop. Coverage')
ax3.set_ylabel('Pathogen Cov.')
ax4.set_ylabel('Conservation')

tt1 = ax1.set_title('(a)')
tt2 = ax2.set_title('(b)')
tt3 = ax3.set_title('(c)')
tt4 = ax4.set_title('(d)')

tt1.set_position((0.53, 0.82))
tt2.set_position((0.53, 0.82))
tt3.set_position((0.53, 0.0))
tt4.set_position((0.53, 0.0))

ax1.grid(False, axis='x')
ax2.grid(False, axis='x')
ax3.grid(False, axis='x')
ax4.grid(False, axis='x')

sns.despine(fig, ax1)
sns.despine(fig, ax2)
sns.despine(fig, ax3)
sns.despine(fig, ax4)

ax1.legend(loc='upper left')

plt.savefig(f'plots/Fig{FIG_ORDER["advantage"]}.{FIG_FORMAT}', bbox_inches='tight')
# -

# # cocktail (fig. 4)

# We parse the log file and look for this piece:
#
# ```
# 2019-09-24 10:51:58,687 INFO: The polypeptide has 184 epitopes 
# 2019-09-24 10:51:58,687 INFO: The epitopes have immunogenicity 2.316
# 2019-09-24 10:51:58,687 INFO: The epitopes cover 27 alleles
# 2019-09-24 10:51:58,687 INFO: The maximum population coverage is 91.29%
# 2019-09-24 10:51:58,687 INFO: The epitopes cover 91.29% of the population (100.00% of the maximum)
# 2019-09-24 10:51:58,688 INFO: The epitopes cover 1736 proteins (90.56% of the total)
# 2019-09-24 10:51:58,688 INFO: The average epitope conservation is 13.17%
# ```
#
# Every piece is about a different chain, exept for the last one which is about the vaccine as a whole.

# +
import re

epis, immunogs, popc, protc, cons = [], [], [], [], []
with open('experiments/results/hiv1bc-full/mosaic-4cocktail-evaluation.log') as f:
    for row in f:
        match = re.search(r'The polypeptide has (\d+) epitopes', row)
        if match:
            epis.append(int(match.group(1)))
            continue
        
        match = re.search(r'The epitopes have immunogenicity ([0-9.]+)', row)
        if match:
            immunogs.append(float(match.group(1)))
            continue
        
        match = re.search(r'The epitopes cover [0-9.]+% of the population \(([0-9.]+)% of the maximum\)', row)
        if match:
            popc.append(float(match.group(1)) / 100)
            continue
        
        match = re.search(r'The epitopes cover \d+ proteins \(([0-9.]+)% of the total\)', row)
        if match:
            protc.append(float(match.group(1)) / 100)
            continue
        
        match = re.search(r'The average epitope conservation is ([0-9.]+)%', row)
        if match:
            cons.append(float(match.group(1)) / 100)
# -

mos_df = pd.DataFrame({
    'num_epitopes': epis,
    'immunogen': immunogs,
    'rel_pop_coverage': popc,
    'norm_prot_coverage': protc,
    'conservation': cons,
    'vax': ['chain-%d' % i for i in range(len(epis) - 1)] + ['mosaic'],
})
mos_df

ot_df = pd.read_csv('experiments/results/hiv1bc-full/optitope-p99-evaluation.csv')
ot_df = ot_df.drop(['prot_coverage', 'alleles', 'pop_coverage', 'max_pop_coverage'], axis=1)
ot_df['vax'] = 'optitope'
ot_df

highig_df = pd.read_csv('experiments/results/hiv1bc-full/mosaic-highig-evaluation.csv')
highig_df = highig_df.drop(['prot_coverage', 'alleles', 'pop_coverage', 'max_pop_coverage'], axis=1)
highig_df['vax'] = 'highig'
highig_df

df = pd.concat([mos_df, ot_df, highig_df], sort=True)
df


# +
def num_epis_to_num_aas(x):
    if x['vax'].startswith('chain') or x['vax'] == 'highig':
        return x['num_epitopes'] + 8
    elif x['vax'] == 'mosaic':
        return 0.0
    elif x['vax'] == 'optitope':
        return x['num_epitopes'] * 9
    
df['num_aas'] = df.apply(num_epis_to_num_aas, axis=1)
df.loc[(df.vax == 'mosaic'), 'num_aas'] = df[df.vax.str.startswith('chain-')].num_aas.sum()

max_short, max_aas = df.immunogen.max(), df.num_aas.max()
df.immunogen = df.immunogen / max_short
df.num_aas = df.num_aas / max_aas
# -

res = df.drop(['vax', 'num_epitopes'], axis=1).T
res.columns = df.vax
res

res = res.rename(index={
    'conservation': 'Epitope\nConservation',
    'immunogen': 'Immunogen.\n100\\%%=%.3f' % max_short,
    'norm_prot_coverage': 'Pathogen\nCoverage',
    'rel_pop_coverage': 'Population\nCoverage',
    'num_aas': 'Amino Acids\n100\\%%=%d' % max_aas,
})

# +
set_font_size(8)
fig, ax = plt.subplots(figsize=(5, 5 * 1 / 2.), dpi=300)

#color = [theme_light[0]] * 4 + [theme_light[1], theme_light[2]]
#ax.bar(x=[x - 0.25 for x in range(5)], width=0.32, height=res.mosaic, align='edge', edgecolor='k', linewidth=2, color='white')
#res.drop('mosaic', axis=1).plot.bar(color=color, rot=0, legend=False, ax=ax, edgecolor='black')

color = [theme_light[0]] * 4 + [theme_dark[0], theme_light[1], theme_light[2]]
res.plot.bar(color=color, rot=0, legend=False, ax=ax, edgecolor='black')

plt.yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0], ['0\%', '20\%', '40\%', '60\%', '80\%', '100\%'])

plt.legend(handles=[
    mpatches.Patch(facecolor=color[0], edgecolor='k', label='Polypep.'),
    #mpatches.Patch(facecolor='#ffffff', edgecolor='k', label='Cocktail'),
    mpatches.Patch(facecolor=color[4], edgecolor='k', label='Cocktail'),
    mpatches.Patch(facecolor=color[5], edgecolor='k', label='EpiMix'),
    mpatches.Patch(facecolor=color[6], edgecolor='k', label='High imm.'),
])
plt.grid(False, axis='x')

plt.subplots_adjust(hspace=0, wspace=0, left=0, bottom=0, right=1, top=1)
#sns.despine()

plt.savefig(f'plots/Fig{FIG_ORDER["cocktail"]}.{FIG_FORMAT}', bbox_inches='tight')


# -

# # correlate coverage with variations and immunogenicity (fig. 6)

# ## find covered positions

# +
def read_fasta(fname):
    with open(fname) as f:
        prots = {}
        for row in f:
            row = row.strip()
            if row.startswith('>'):
                pid = row[1:]
                prots[pid] = []
            else:
                prots[pid].append(row)
    return {pid: ''.join(ps) for pid, ps in prots.items()}

prots = read_fasta('experiments/resources/hiv1-bc-nef.fasta')


# +
def read_vax(vax_file):
    with open(vax_file) as f:
        _ = next(f)
        peps = [r.split(',')[-1].strip() for r in f]
    return peps


vax_short = read_vax('experiments/results/hiv1bc-full/mosaic-entropy-short.csv')
vax_sob = read_vax('experiments/results/hiv1bc-full/optitope-entropy.csv')
vax_mid = read_vax('experiments/results/hiv1bc-full/mosaic-entropy-mid.csv')
vax_long = read_vax('experiments/results/hiv1bc-full/mosaic-entropy-long.csv')


# +
def find_coverage_by_protein(prots, vax, sort=False):
    # step 1: find all occurrences of all epitopes
    vax = set(vax)
    peps = defaultdict(list)
    for i, pp in enumerate(prots.values()):
        for j in range(len(pp) - 8):
            seq = pp[j:j+9]
            if seq in vax:
                peps[seq].append((i, j))
    
    # step 2: aggregate the positions by protein
    mat = np.zeros((len(prots), max(len(p) for p in prots.values())))
    for poss in peps.values():
        for i, j in poss:
            mat[i, j:j+9] += 1
    
    # step 3: sort by similarity
    if sort:
        link = hierarchy.linkage(mat, 'ward')
        idx = hierarchy.leaves_list(link)
        mat = mat[idx]
    
    return mat


mat_short = find_coverage_by_protein(prots, vax_short)
mat_mid = find_coverage_by_protein(prots, vax_mid)
mat_long = find_coverage_by_protein(prots, vax_long)
mat_sob = find_coverage_by_protein(prots, vax_sob)
# -

plt.matshow(mat_sob.T)

idx = np.random.choice(len(prots), len(prots))
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 6))
ax1.matshow(mat_short[idx].T, cmap='Reds')
ax1.grid(False)
ax1.set_xlabel('Protein')
ax1.xaxis.set_label_position('top') 
ax1.set_ylabel('Position')
ax2.set_ylabel('Position')
ax2.matshow(mat_long[idx].T, cmap='Reds')
ax2.grid(False)
ax2.set_title('Top: Short mosaic - Bottom: Long mosaic')
ax2.tick_params(axis='x', labeltop=False)
fig.set_tight_layout(True)
fig.subplots_adjust(left=0, bottom=-0.35, right=1, top=1, wspace=0, hspace=0)
#fig.savefig('plots/epitopes-positions.tiff', bbox_inches='tight')

# ## find immunogenicity by position

# +
with open('experiments/results/hiv1bc-full/made-epitopes.csv') as f:
    _ = next(f)
    epitope_immunogen = {row.split(',')[-1].strip(): float(row.split(',')[0]) for row in f}

peps = defaultdict(list)
for i, pp in enumerate(prots.values()):
    for j in range(len(pp) - 8):
        seq = pp[j:j+9]
        if seq in epitope_immunogen:
            peps[seq].append((i, j))

# step 2: aggregate the positions by protein
mat = np.zeros((len(prots), max(len(p) for p in prots.values())))
counts = np.zeros((len(prots), max(len(p) for p in prots.values())))
for pep, positions in peps.items():
    for i, j in positions:
        mat[i, j:j+9] += epitope_immunogen[pep]
        counts[i, j:j+9] += 1

#avg = mat.sum(axis=0) / counts.sum(axis=0)
#ig_val = np.where(np.isfinite(avg), avg, 0.0)
ig_mat = mat
# -

plt.matshow(ig_mat.T)

# ## find entropy by position

# from https://mafft.cbrc.jp/alignment/server/
aligned = read_fasta('plots/aligned.fasta')

aminos = sorted(set(a for p in prots.values() for a in p)) + ['-']
counts = np.zeros((max(len(a) for a in aligned.values()), len(aminos) + 1))
for p in aligned.values():
    for i, a in enumerate(p):
        counts[i, aminos.index(a)] += 1

# remove gaps from aligned counts
#aminocounts = np.array([c for c in counts if np.argmax(c) < len(aminos) - 1])
aminocounts = counts

plt.matshow(aminocounts.T)
plt.grid(False)

# +
# aggregate counts for same position
# idea 1: just sum how many different aminoacids
ss = np.sum(aminocounts > 0, axis=1)

# idea 2: entropy
ps = aminocounts / np.sum(aminocounts, axis=1, keepdims=True)
es = -ps * np.log(ps)
ss = np.sum(np.where(np.isfinite(es), es, 0.0), axis=1)

# smooth
variants = pd.Series(ss).rolling(9, center=True, win_type='triang').mean()
# -

plt.plot(variants)

# ## align immunogenicity and coverage

# +
# step 3: use alignments
ordered_alignments = [aligned[k] for k in prots]


def align_matrix(mat):
    # mat is (proteins,positions)
    # expand the positions to account for the alignment
    res = np.zeros((len(mat), max(len(a) for a in ordered_alignments)))
    for i, values in enumerate(mat):
        cursor = 0
        for j, a in enumerate(ordered_alignments[i]):
            if a != '-':
                res[i, j] = mat[i, cursor]
                cursor += 1
    return res


ig_ali = align_matrix(ig_mat)
short_ali = align_matrix(mat_short)
mid_ali = align_matrix(mat_mid)
long_ali = align_matrix(mat_long)
sob_ali = align_matrix(mat_sob)
# -

plt.matshow(sob_ali.T)

# +
consensus_positions = np.argmax(counts, axis=1) < len(aminos) - 1

cov_long = long_ali.sum(axis=0)[consensus_positions]
cov_long /= cov_long.max()

cov_mid = mid_ali.sum(axis=0)[consensus_positions]
cov_mid /= cov_mid.max()

cov_short = short_ali.sum(axis=0)[consensus_positions]
cov_short /= cov_short.max()

cov_sob = sob_ali.sum(axis=0)[consensus_positions]
cov_sob /= cov_sob.max()

cov_ig = ig_ali.sum(axis=0)[consensus_positions] / len(prots)
ent_smo = variants[consensus_positions]
ent = ss[consensus_positions]
# -

# ## correlations

agg = pd.DataFrame({
    'Immuno,': cov_ig,
    'Entropy': ent_smo,
    'Long Mosaic': cov_long,
    'Short Mosaic': cov_short,
    'OptiTope': cov_sob,
})

corrs = agg.corr()
corrs

# ## plot everything together

# +
set_font_size(6)

fig = plt.figure(figsize=(5.2, 5.2 / 4), dpi=300)

gs = gridspec.GridSpec(2, 2, height_ratios=[1, 1], width_ratios=[3, 1],
                       hspace=0.05, wspace=0.1, left=0, bottom=0, right=1, top=1,
                       figure=fig)
ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[1, 0])

ax12 = ax1.twinx()
p1 = ax1.plot(ent_smo.values, color='#1f78b4', label='Entropy')
ax1.plot(ent, color='#1f78b4', alpha=0.2)
p2 = ax12.plot(cov_ig, color='#33a02c', label='Immuno.')

p3 = ax2.plot(cov_mid, color=theme_dark[0], label='Long Mosaic')
p5 = ax2.plot(cov_sob, color=theme_dark[1], label='OptiTope')
p4 = ax2.plot(cov_short, color=theme_dark[2], label='Short Mosaic')

# adjust chart
ax1.set_title('(a)')
ax2.set_title('(b)', pad=-8)
ax1.set_ylabel('Entropy')
ax1.set_yticks([0, 1, 2])
ax12.set_ylabel('Immuno.')
ax12.set_yticks([0, 0.5, 1])

ax12.yaxis.set_label_coords(0.97, 0.5)
ax1.tick_params(axis='x', labelbottom=False)
ax1.grid(True, axis='x')
ax12.grid(False)

ax1.legend(p1 + p2, [p.get_label() for p in p1 + p2], ncol=2, loc='upper left', bbox_to_anchor=(-0.012, 1.06))
ax2.legend(ncol=1, loc='upper left',
           bbox_to_anchor=(-0.012, 1.07))
ax2.set_xlabel('Position')
ax2.set_ylabel('Rel. Coverage')

ax1.set_xticks(range(0, 210, 18))
ax2.set_xticks(range(0, 210, 18))
ax2.set_yticks([0, 0.5, 1.0])
ax1.set_ylim(0, 2.15)
ax12.set_ylim(0, 1.1)
ax2.set_ylim(0, 1.1)

# scatter matrix
cm = plt.get_cmap('bwr')
gs00 = gridspec.GridSpecFromSubplotSpec(5, 5, subplot_spec=gs[:,1], hspace=0, wspace=0)
for i in range(len(corrs)):
    for j in range(len(corrs)):
        plot = fig.add_subplot(gs00[i, j])
        
        plot.set_xticks([])
        plot.set_yticks([])

        # write correlation coeff and pval in the lower triangular half
        xs, ys = agg.iloc[4:-4, j], agg.iloc[4:-4, i]
        if i > j:
            corr = spearmanr(xs, ys)
            plot.set_facecolor(cm(int(128 + corr[0] * 128)))
            
            # bold if significant (with bonferroni correction)
            if corr[1] < 0.05 / 10:
                fmt = '\\textbf{{{:.3f}}}\n\\textbf{{{:.0e}}}'
            else:
                fmt = '{:.3f}\n{:.0e}'

            plot.text(0.5, 0.5, fmt.format(*corr),
                      fontdict={'size': 'small'},
                      horizontalalignment='center',
                      verticalalignment='center',
                      weight='bold',
                      transform=plot.transAxes)

        # add regression line to upper triangular half
        elif i < j:
            plot.scatter(xs, ys, alpha=0.2, marker='.')
            xs = np.vstack([xs, np.ones(len(xs))])
            ws = np.linalg.inv(xs.dot(xs.T) + 1e-6 * np.eye(2)).dot(xs).dot(ys)
            xlims, ylims = plot.get_xlim(), plot.get_ylim()
            reg = xlims[0] * ws[0] + ws[1], xlims[1] * ws[0] + ws[1]
            plot.plot(xlims, reg, 'C3')
            plot.set_ylim(ylims)
            plot.set_xlim(xlims)

        # histogram on the diagonal
        else:
            plot.hist(xs)
            plot.set_yscale('log')
            plot.grid(False)
            plot.tick_params(axis='both', which='both', left=False,
                             bottom=False, labelleft=False)

        # put labels on the bottom / right
        if j == len(corrs) - 1:
            plot.set_ylabel(corrs.columns[i], rotation=-45, ha='left')
            plot.yaxis.set_label_coords(1.4, 0.5)
        if i == len(corrs) - 1:
            plot.set_xlabel(corrs.columns[j], rotation=-45, ha='left')
        if i == 0 and j == 2:
            plot.set_title('(c)')

#sns.despine(fig, ax2, top=True, right=False)

plt.savefig(f'plots/Fig{FIG_ORDER["entropy"]}.{FIG_FORMAT}', bbox_inches='tight')


# -

# # compare rank and ic50 immunogenicities (fig. 8)

# +
def get_epitope_immunogen(bootstrap, method):
    if method == 'netmhcpan':
        mm = ''
    else:
        mm = method + '-'

    igs = {}
    with open(f'experiments/results/nef-300-{bootstrap}/made-{mm}epitopes.csv') as f:
        header = next(f).strip().split(',')
        for row in f:
            row = dict(zip(header, row.strip().split(',')))
            igs[row['epitope']] = float(row['immunogen'])
    return igs


def get_tradeoff_ranks(bootstrap, method, immunogens):
    if method == 'netmhcpan':
        mm = ''
    else:
        mm = method + '-'

    epi_igs = []
    with open(f'experiments/results/nef-300-{bootstrap}/made-{mm}tradeoff.csv') as f:
        header = next(f).strip().split(',')
        for row in f:
            row = dict(zip(header, row.strip().strip().split(',')))
            epi_igs.append({
                method: sum(ig[e] for e in row['epitopes'].split(';'))
                for method, ig in immunogens.items()
            })
            epi_igs[-1]['optimize'] = method
            epi_igs[-1]['epitopes'] = row['epitopes']
            epi_igs[-1]['cleavage'] = float(row['cleavage'])
    return epi_igs


# -

all_igs = []
for i in range(1, 6):
    immunogens = {
        method: get_epitope_immunogen(i, method)
        for method in ['netmhcpan', 'netmhcpan-rank', 'pickpocket', 'mhcflurry']
    }
    for method in ['netmhcpan', 'netmhcpan-rank', 'pickpocket', 'mhcflurry']:
        all_igs.extend(get_tradeoff_ranks(i, method, immunogens))

df = pd.DataFrame([
    {
        'method': m,
        'epitope': e,
        'immunog': i
    }
    for m, es in immunogens.items()
    for e, i in es.items()
]).pivot(index='epitope', columns='method', values='immunog')
df

df.corr()

# within-method immunogenicity variation
for col in df:
    vals = df[col].sort_values()
    windows = []
    j = k = 0
    for i in range(len(vals)):
        while j < i and (vals[i] - vals[j]) / vals[i] > 0.005:
            j += 1
        while k < len(vals) and (vals[k] - vals[i]) / vals[i] < 0.005:
            k += 1
        windows.append(k - j + 1)

    print(col)
    print(pd.Series(windows).describe())

df = pd.DataFrame(all_igs)
df

# +
set_font_size(6)
fig = plt.figure(figsize=(6.7261, 6.7261 * 3 / 4), dpi=300)
axes = fig.subplots(4, 4)

renames = {
    'netmhcpan': 'NetMHCpan (IC$_{50}$)',
    'netmhcpan-rank': 'NetMHCpan (\%rank)',
    'pickpocket': 'PickPocket (IC$_{50}$)',
    'mhcflurry': 'MHCflurry (IC$_{50}$)',
}

for i, ig_method in enumerate(['netmhcpan', 'netmhcpan-rank', 'pickpocket', 'mhcflurry']):
    for j, optimize in enumerate(['netmhcpan', 'netmhcpan-rank', 'pickpocket', 'mhcflurry']):
        ax = axes[i, j]
        data = df[df['optimize'] == optimize]
        
        if i == j:
            sns.distplot(data[optimize].values, ax=ax, rug=True, kde=False, norm_hist=False, bins=8)
        else:
            ax.scatter(data[optimize], data[ig_method], c=data['cleavage'], cmap='viridis', marker='.')
            ax.set_ylim((df[ig_method].min(), df[ig_method].max()))
        
            e1 = [set(e.split(';')) for e in df[df['optimize'] == optimize].sort_values('cleavage')['epitopes'].values]
            e2 = [set(e.split(';')) for e in df[df['optimize'] == ig_method].sort_values('cleavage')['epitopes'].values]
            ious = [len(a & b) / len(a | b) for a, b in zip(e1, e2)]
            shareds = [len(a & b) for a, b in zip(e1, e2)]
            iou = '{:.2f} ({:.2f} - {:.2f})'.format(*np.percentile(ious, [50, 25, 75]))
            shared = '{:.0f} ({:.0f} - {:.0f})'.format(*np.percentile(shareds, [50, 25, 75]))
            
            if 'rank' not in optimize:
                cname, (corr, pval) = 'Pearson', pearsonr(data[optimize], data[ig_method])
            else:
                cname, (corr, pval) = 'Spearman', spearmanr(data[optimize], data[ig_method])
            
            if i == 1:
                xy = 0.975, 0.025
                ha, va = 'right', 'bottom'
            else:
                xy = 0.025, 0.975
                ha, va = 'left', 'top'
                
            ax.annotate(f'IoU: {iou}\nShared epitopes: {shared}\n{cname}: {corr:.2f} ($p$={pval:.0e})',
                        xy=xy, xycoords='axes fraction', ha=ha, va=va)
            
        ax.tick_params(axis='y', rotation=90)
        ax.set_xlim((data[optimize].min(), data[optimize].max()))
        sns.despine(fig, ax)
        ax.grid(False)
        
        if i == 3:
            ax.set_xlabel('Optimize\n' + renames.get(optimize, optimize))
        if j == 0:
            ax.set_ylabel('Immunogenicity\n' + renames.get(ig_method, ig_method))
            

fig.tight_layout()
plt.savefig(f'plots/Fig{FIG_ORDER["immunogens"]}.{FIG_FORMAT}', bbox_inches='tight')
# -

# # compare optimized cleavage with shuffled epitopes (fig. 9)

from data_preparation import get_cleavage_score_process

records = []
for i in tqdm(range(1, 6)):
    with open(f'experiments/results/nef-300-{i}/made-tradeoff.csv') as f:
        header = next(f).strip().split(',')
        for j, row in tqdm(enumerate(f), leave=False, total=10):
            row = dict(zip(header, row.strip().split(',')))
            epis = row['epitopes'].split(';')

            scores = get_cleavage_score_process(
                penalty=0.1, cleavage_model='pcm', window_size=5,
                epitopes=list(zip(epis[:-1], epis[1:]))
            )
            assert np.allclose(sum(x[2] for x in scores), float(row['cleavage']))
            
            records.append({
                'immunogen': float(row['immunogenicity']),
                'cleavage': float(row['cleavage']),
                'positive': sum(x[2] > 0 for x in scores),
                'shuffled': False,
                'bootstrap': i,
                'index': j,
                'epis': row['epitopes'],
            })
            
            for k in range(50):
                np.random.shuffle(epis)
                scores = get_cleavage_score_process(
                    # always make sure these params are equal to
                    # what was used for the optimization
                    penalty=0.1, cleavage_model='pcm', window_size=5,
                    epitopes=list(zip(epis[:-1], epis[1:]))
                )
                
                records.append({
                    'immunogen': float(row['immunogenicity']),
                    'cleavage': sum(x[2] for x in scores),
                    'positive': sum(x[2] > 0 for x in scores),
                    'shuffled': True,
                    'bootstrap': i,
                    'index': j,
                    'rep': k,
                    'epis': ';'.join(epis)
                })

df = pd.DataFrame(records)
df

for i, g in df.groupby(['bootstrap', 'index']):
    opt = g[g['rep'].isna()].iloc[0]
    best = g[~g['rep'].isna()].loc[g[~g['rep'].isna()].cleavage.idxmin()]
    
    print(opt['cleavage'], opt['epis'])
    print(best['cleavage'], best['epis'])
    print('---')

# +
# number of cleavage sites with positive score
# False = shuffled, True = optimized

df.groupby(df['rep'].isna()).apply(lambda g: g['positive'].describe())


# -

def legend_without_duplicate_labels(ax, **kwargs):
    # https://stackoverflow.com/a/56253636/521776
    handles, labels = ax.get_legend_handles_labels()
    unique = [(h, l) for i, (h, l) in enumerate(zip(handles, labels)) if l not in labels[:i]]
    ax.legend(*zip(*unique), **kwargs)


# +
fig = plt.figure(figsize=(5, 3), dpi=300)
ax2, ax1 = fig.subplots(1, 2)
ax11 = ax1.twiny()

order = df.groupby([
    'bootstrap', 'index'
]).agg({
    'immunogen': 'min',
    'cleavage': 'mean'
}).sort_values([
    'bootstrap',
    'immunogen'
]).reset_index()


for i, row in order.iterrows():
    g = df[(df['bootstrap'] == row['bootstrap']) & (df['index'] == row['index'])]
    cleavs = g['cleavage'].values
    
    assert g['rep'].isna().sum() == 1
    reduction = g[~g['rep'].isna()]['cleavage'] / g[g['rep'].isna()]['cleavage'].values[0]
    reduction = 100 * (1 - reduction)
    
    #ax2.boxplot([reduction], positions=[i])
    ax2.scatter(reduction, [i] * len(reduction), marker='.',
                c=f'C{int(row["bootstrap"]) + 1}')
    ax2.scatter([np.median(reduction)], [i], marker='x', c='k', label='Median')

pos = 0
for i, g in order.groupby(['bootstrap']):
    ax1.plot(g['immunogen'], [pos + j for j in range(len(g))], 'C0')
    ax11.plot(-g['cleavage'], [pos + j for j in range(len(g))], 'C1')
    pos += len(g)

ax1.grid(False)
ax11.grid(False)
ax2.grid(False)
ax1.set_yticks([])
ax2.set_yticks([])
ax1.set_title('(b)')
ax2.set_title('(a)')
ax11.set_xlabel('Cleavage Score', color='C1')
ax1.set_xlabel('Vaccine immunogenicity', color='C0')

fig.subplots_adjust(wspace=0.)
ax1.set_ylim(ax2.get_ylim())
ax1.grid(False, axis='x')
ax11.grid(False, axis='x')
ax2.grid(False, axis='x')

group_sizes = [0] + order.groupby(['bootstrap']).agg('size').cumsum().tolist()
ax2.set_yticks([(a + b) / 2 for a, b in zip(group_sizes[:-1], group_sizes[1:])])
ax2.set_yticklabels(range(1, len(group_sizes)))
ax2.set_ylabel('Bootstrap')
ax2.set_xlabel('Percent decrease in cleavage score')
sns.despine(fig, ax1)
sns.despine(fig, ax2)

legend_without_duplicate_labels(ax2)

fig.tight_layout()
plt.savefig(f'plots/Fig{FIG_ORDER["shuffled"]}.{FIG_FORMAT}', bbox_inches='tight')


# -

# # show peptides overlap

# +
def load_epis(fname):
    with open(fname) as f:
        _ = next(f)
        epis = []
        for row in f:
            epis.append(row[:-1].split(',')[-1])
    return {e: i for i, e in enumerate(sorted(epis))}


epis = [
    load_epis('experiments/results/nef-300-%d/made-epitopes.csv' % i)
    for i in range(1, 6)
]

all_epis = load_epis('experiments/results/hiv1bc-full/made-epitopes.csv')
# -

# pairwise overlaps
ovs = np.zeros((len(epis), len(epis)))
for i in range(len(epis)):
    for j in range(len(epis)):
        ovs[i, j] = len(epis[i].keys() & epis[j].keys()) / max(len(epis[i].keys()), len(epis[j].keys()))
ovs

# +
coverage = np.zeros((len(epis), len(all_epis)))
for e, i in all_epis.items():
    for j, b in enumerate(epis):
        if e in b:
            coverage[j, i] = 1

np.mean(coverage, axis=1)  # portion of covered epitopes 
# -

# try to sort stuff decently
mat = np.array(sorted(coverage.T, key=lambda x: sum(x))).T
link = hierarchy.linkage(mat, 'ward')
idx = hierarchy.leaves_list(link)
mat = mat[idx]

# +
fig = plt.figure(figsize=(10, 5))
ax1, ax2 = fig.subplots(1, 2,)

ax1.imshow(mat, aspect='auto')
ax1.grid(False)

im = ax2.imshow(ovs, vmin=0.4, vmax=0.5)
fig.colorbar(im, ax=ax2)
# -




