"""Plot two completed solver cases without rerunning or changing their data."""
from pathlib import Path
import hashlib
import json
import sys

ROOT = Path(__file__).resolve().parents[1]
RUN = ROOT

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scientific_plot_style import configure_matplotlib, style_axes, save_png_pdf

OUT = ROOT / 'figures'
CASES = ['domain30_dx1_np64_r8_20261016', 'domain30_dx0p5_np64_r32_20261016']
configure_matplotlib()
fig, axes = plt.subplots(1, 2, figsize=(18, 9.5), sharey=True)
fig.subplots_adjust(left=.075, right=.98, bottom=.12, top=.68, wspace=.13)
records = []
for index, (ax, name) in enumerate(zip(axes, CASES)):
    folder = RUN / 'outputs' / name
    data = np.load(folder / 'results.npz')
    summary = json.loads((folder / 'summary.json').read_text())
    p = summary['parameters']
    for k, expected in {'start':-15, 'end':15, 'u':10, 'D':.5, 'alpha':-1,
                        'kappa':.1, 'T':1, 'K':2, 'Lp':16, 'Np':64, 'mu':-10}.items():
        assert p[k] == expected
    error = np.linalg.norm(data['phi_trotter']-data['phi_nonlinear']) / np.linalg.norm(data['phi_nonlinear'])
    assert np.isclose(error, summary['errors']['phit_vs_finite_difference'])
    for key, label, color, marker, ls in [
        ('phi_nonlinear', 'Finite difference', '#222222', 'o', '-'),
        ('phi_trotter', 'phit', '#CC0077', 'D', '--')]:
        ax.plot(data['x'], data[key].real, label=label, color=color, marker=marker,
                markersize=6.5, markerfacecolor='white', linestyle=ls, linewidth=2)
    style_axes(ax, xlabel=r'$x$', ylabel=r'$\phi$' if index == 0 else '')
    ax.set_xlim(-15,15); ax.set_ylim(-.08,.23)
    ax.set_xticks(np.arange(-15,16,5)); ax.set_yticks(np.arange(-.05,.201,.05))
    ax.set_title(
        rf"({'ab'[index]}) $\Delta x={p['deltax']:g},\ N_x={len(data['x'])},\ r={p['r']}$"
        + '\n' + rf'$E={100*error:.2f}\%$', fontsize=21, fontweight='bold', pad=14)
    records.append(dict(case=name, parameters=p, relative_l2_error=float(error),
                        results_sha256=hashlib.sha256((folder/'results.npz').read_bytes()).hexdigest()))
fig.text(.5,.96,r'$x\in[-15,15),\ \mathrm{periodic},\quad \phi(0,x)=e^{-(x+10)^2}$',
         ha='center', va='center', fontsize=20)
fig.text(.5,.91,r'$u=10,\quad D=0.5,\quad \alpha=-1,\quad \kappa=0.1$',
         ha='center', va='center', fontsize=20)
fig.text(.5,.86,r'$K=2,\quad L_p=16,\quad N_p=64,\quad T=1$',
         ha='center', va='center', fontsize=20)
handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc='center', bbox_to_anchor=(.5,.80), ncol=2,
           frameon=False, prop={'family':'Times New Roman','weight':'bold','size':20})
assert all(line.get_linewidth()==2 for ax in axes for line in ax.lines)
save_png_pdf(fig, OUT/'selected_grid_results')
plt.close(fig)
(OUT/'figure_sources.json').write_text(json.dumps(records,indent=2),encoding='utf-8')
print(json.dumps(records,indent=2))
