import matplotlib.pyplot as plt
import numpy as np

import matplotlib.ticker as tkr
import matplotlib.dates as mdates
import datetime


import seaborn as sns
from matplotlib.colors import ListedColormap
sns.set_style("whitegrid")
sns.set_style("ticks")

import scipy.io as sio

from matplotlib import rc
rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})
rc('text', usetex=True)


def load_data():
    mat_contents = sio.loadmat('../IO/MyLakeResults.mat')
    MyLake_results = mat_contents['MyLake_results']
    Sediment_results = mat_contents['Sediment_results']
    return MyLake_results, Sediment_results


def contour_plot(results, elem, lbl=False, years_ago=1, cmap=ListedColormap(sns.color_palette("Blues", 51))):
    plt.figure(figsize=(6, 4), dpi=192)
    start = -365 * years_ago
    end = -365 * (years_ago - 1) - 1
    X, Y = np.meshgrid(results['days'][0, 0][0][start:end] - 366, -results['z'][0, 0])
    z = 0
    for e in elem:
        z += results[e][0, 0][:, start:end]
    CS = plt.contourf(X, Y, z, 51, cmap=cmap, origin='lower')
#     plt.clabel(CS, inline=1, fontsize=10, colors='w')
    cbar = plt.colorbar(CS)

    if not results['O2zt'][0, 0][:, 0].shape[0] == 256:
        ice_thickness = results['His'][0, 0][0, start:end]
        plt.fill_between(results['days'][0, 0][0][start:end] - 366, 0, -ice_thickness, where=-ice_thickness <= 0, facecolor='red', interpolate=True)

    plt.ylabel('Depth, [m]')
    ax = plt.gca()
    ax.ticklabel_format(useOffset=False)
    ax.xaxis.set_major_locator(mdates.MonthLocator(interval=2))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b\n%Y'))
    if lbl:
        cbar.ax.set_ylabel(lbl)
    else:
        cbar.ax.set_ylabel(elem)
    plt.show()


def plot_profile(results, elem):
    plt.figure(figsize=(6, 4), dpi=192)
    for e in elem:
        plt.plot(results[e][0, 0][1:-1, -1], -results['z'][0, 0][1:-1, -1], lw=3, label=e)
    plt.xlabel('mmol/L')
    plt.ylabel('Depth, m')
    ax = plt.gca()
    ax.ticklabel_format(useOffset=False)
    ax.grid(linestyle='-', linewidth=0.2)
    plt.legend()
    plt.ylim([-results['z'][0, 0][-1], -results['z'][0, 0][0]])
    plt.tight_layout()
    plt.show()


def plot_fit(MyLake_results):
    fig, axes = plt.subplots(4, 1, sharex='col', figsize=(8, 5), dpi=192)

    inx = sum(MyLake_results['z'][0, 0] < 4)[0]
    TOTP = np.mean(MyLake_results['Pzt'][0, 0][0:inx, :], axis=0) + \
        np.mean(MyLake_results['PPzt'][0, 0][0:inx, :], axis=0) + \
        np.mean(MyLake_results['DOPzt'][0, 0][0:inx, :], axis=0) + \
        np.mean(MyLake_results['DOCzt'][0, 0][0:inx, :], axis=0) + \
        np.mean(MyLake_results['Chlzt'][0, 0][0:inx, :], axis=0) + \
        np.mean(MyLake_results['Czt'][0, 0][0:inx, :], axis=0)  + \
        np.mean(MyLake_results['POCzt'][0, 0][0:inx, :], axis=0)
    Chl = np.mean(MyLake_results['Czt'][0, 0][0:inx, :], axis=0) + np.mean(MyLake_results['Chlzt']
                                                                           [0, 0][0:inx, :], axis=0)
    PO4 = np.mean(MyLake_results['Pzt'][0, 0][0:inx, :], axis=0)
    Part = np.mean(MyLake_results['PPzt'][0, 0][0:inx, :], axis=0) + np.mean(MyLake_results['POCzt'][0, 0][0:inx, :], axis=0)

    axes[0].plot(-366 + MyLake_results['days'][0, 0][0], TOTP, c=sns.xkcd_rgb["denim blue"], lw=3, label='Total P')
    axes[1].plot(-366 + MyLake_results['days'][0, 0][0], Chl, c=sns.xkcd_rgb["denim blue"], lw=3, label='Chl-a')
    axes[2].plot(-366 + MyLake_results['days'][0, 0][0], PO4, c=sns.xkcd_rgb["denim blue"], lw=3, label='PO_4')
    axes[3].plot(-366 + MyLake_results['days'][0, 0][0], Part, c=sns.xkcd_rgb["denim blue"], lw=3, label='Solid P')

    TOTP = np.loadtxt('../obs/store_obs/TOTP.dat', delimiter=',')
    Chl = np.loadtxt('../obs/store_obs/Cha_aquaM_march_2017.dat', delimiter=',')
    PO4 = np.loadtxt('../obs/store_obs/PO4.dat', delimiter=',')
    Part = np.loadtxt('../obs/store_obs/Part.dat', delimiter=',')
    axes[0].plot(-366 + TOTP[:, 0], TOTP[:, 1], 'bo', c=sns.xkcd_rgb["pale red"], markersize=4)
    axes[1].plot(-366 + Chl[:, 0], Chl[:, 1], 'bo', c=sns.xkcd_rgb["pale red"], markersize=4)
    axes[2].plot(-366 + PO4[:, 0], PO4[:, 1], 'bo', c=sns.xkcd_rgb["pale red"], markersize=4)
    axes[3].plot(-366 + Part[:, 0], Part[:, 1], 'bo', c=sns.xkcd_rgb["pale red"], markersize=4)

    axes[3].xaxis.set_major_locator(mdates.MonthLocator(interval=12))
    axes[3].xaxis.set_major_formatter(mdates.DateFormatter('%b\n%Y'))
    axes[1].set_xlim([732313 - 366, 735234 - 366 * 2])
    for ax in axes:
        ax.grid(linestyle='-', linewidth=0.2)
        ax.set_ylim([0, 50])
        ax.set_ylabel(r'$mg / m^3$')
        ax.legend(loc=1)


def plot_flux(results, elem, lbl, years_ago=1):
    plt.figure(figsize=(6, 4), dpi=192)
    start = -365 * years_ago
    end = -365 * (years_ago - 1) - 1
    try:
        plt.plot(results['days'][0, 0][0][start:end] - 366, results['Bioirrigation_fx_zt'][0, 0][elem]
                 [0, 0][0][start:end], sns.xkcd_rgb["medium green"], lw=3, label='Bioirrigation')
        plt.plot(results['days'][0, 0][0][start:end] - 366, results['sediment_SWI_fluxes'][0, 0][elem][0, 0][0][start:end], sns.xkcd_rgb["denim blue"], lw=3, label='Diffusive')
    except:
        plt.plot(results['days'][0, 0][0][start:end] - 366, results['sediment_SWI_fluxes'][0, 0][elem][0, 0][0][start:end], sns.xkcd_rgb["denim blue"], lw=3, label=elem)
    ax = plt.gca()
    ax.set_ylabel(lbl)
    ax.ticklabel_format(useOffset=False)
    ax.xaxis.set_major_locator(mdates.MonthLocator(interval=2))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b\n%Y'))
    ax.set_xlim([results['days'][0, 0][0][start:end][0] - 366, results['days'][0, 0][0][start:end][-1] - 366])
    ax.grid(linestyle='-', linewidth=0.2)
    ax.legend(loc=1)
    legend = plt.legend(frameon=1)
    frame = legend.get_frame()
    frame.set_facecolor('white')
    plt.tight_layout()
    plt.show()
