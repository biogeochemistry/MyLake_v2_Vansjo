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


def load_data():
    mat_contents = sio.loadmat('../IO/MyLakeResults.mat')
    MyLake_results = mat_contents['MyLake_results']
    Sediment_results = mat_contents['Sediment_results']
    return MyLake_results, Sediment_results


def contour_plot(results, elem):
    plt.figure(figsize=(6, 4), dpi=192)
    X, Y = np.meshgrid(results['days'][0, 0][0][-365:], -results['z'][0, 0])
    z = results[elem][0, 0][:, -365:]
    CS = plt.contourf(X, Y, z, 51, cmap=ListedColormap(sns.color_palette("Blues", 51)), origin='lower')
#     plt.clabel(CS, inline=1, fontsize=10, colors='w')
    cbar = plt.colorbar(CS)
    plt.ylabel('Depth, [m]')
    ax = plt.gca()
    ax.ticklabel_format(useOffset=False)
    ax.xaxis.set_major_locator(mdates.MonthLocator(interval=2))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b'))
    cbar.ax.set_ylabel(elem)
    plt.show()


def plot_profile(results, elem):
    plt.figure(figsize=(6, 4), dpi=192)
    plt.plot(results[elem][0, 0][:, -1], -results['z'][0, 0], sns.xkcd_rgb["denim blue"], lw=3, label=elem)
    plt.xlabel(elem)
    plt.ylabel('Depth, m')
    ax = plt.gca()
    ax.ticklabel_format(useOffset=False)
    ax.grid(linestyle='-', linewidth=0.2)
    plt.legend()
    plt.tight_layout()
    plt.show()


def plot_fit(MyLake_results):
    fig, axes = plt.subplots(4, 1, sharex='col', figsize=(8, 5), dpi=192)
    TOTP = np.loadtxt('../obs/store_obs/TOTP.dat', delimiter=',')
    Chl = np.loadtxt('../obs/store_obs/Cha_aquaM_march_2017.dat', delimiter=',')
    PO4 = np.loadtxt('../obs/store_obs/PO4.dat', delimiter=',')
    Part = np.loadtxt('../obs/store_obs/Part.dat', delimiter=',')
    axes[0].plot(-366 + TOTP[:, 0], TOTP[:, 1], 'bo', c=sns.xkcd_rgb["denim blue"], markersize=4)
    axes[1].plot(-366 + Chl[:, 0], Chl[:, 1], 'bo', c=sns.xkcd_rgb["denim blue"], markersize=4)
    axes[2].plot(-366 + PO4[:, 0], PO4[:, 1], 'bo', c=sns.xkcd_rgb["denim blue"], markersize=4)
    axes[3].plot(-366 + Part[:, 0], Part[:, 1], 'bo', c=sns.xkcd_rgb["denim blue"], markersize=4)
    axes[3].xaxis.set_major_locator(mdates.MonthLocator(interval=12))
    axes[3].xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
    axes[1].set_xlim([732313 - 366, 735234 - 366 * 2])
    for ax in axes:
        ax.grid(linestyle='-', linewidth=0.2)
    axes[0].set_ylim([0, 50])
    axes[1].set_ylim([0, 50])
    axes[2].set_ylim([0, 50])
    axes[3].set_ylim([0, 50])

    inx = sum(MyLake_results['z'][0, 0] < 4)[0]
    TOTP = np.sum(MyLake_results['Pzt'][0:inx, :][0]) + np.sum(MyLake_results['PPzt'][0:inx, :][0]) + np.sum(MyLake_results['DOPzt'][0:inx, :][0])
    axes[0].plot(-366 + MyLake_results['days'][0, 0][0], TOTP, c=sns.xkcd_rgb["denim blue"])
