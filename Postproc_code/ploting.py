import matplotlib.pyplot as plt
import numpy as np
import matplotlib.dates as mdates
import seaborn as sns
from matplotlib.colors import ListedColormap
import scipy.io as sio
from matplotlib import rc

sns.set_style("whitegrid")
sns.set_style("ticks")

rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})
rc('text', usetex=True)


solid_species = ['OMzt', 'OMbzt', 'FeOOHzt', 'FeSzt', 'S8zt', 'FeS2zt', 'AlOH3zt', 'Ca3PO42zt', 'PO4adsazt', 'PO4adsbzt', 'OMSzt', 'FeOH3zt']
solute_species = ['O2zt', 'NO3zt', 'SO4zt', 'NH4zt', 'Fe2zt', 'H2Szt', 'S0zt', 'PO4zt', 'Ca2zt',
                  'HSzt', 'Hzt', 'OHzt', 'CO2zt', 'CO3zt', 'HCO3zt', 'NH3zt', 'H2CO3zt', 'DOM1zt', 'DOM2zt']
molar_masses = {
    'O2': 31998.8,
    'OM': 30973.762,
    'OM': 30973.762,
    'DOM1': 12010.7,
    'DOM2': 12010.7,
    'NO3': 62004,
    'FeOH3': 106867.0,
    'SO4': 96062,
    'NH4': 18038,
    'Fe2': 55845,
    'H2S': 34080.9,
    'HS': 33072.9,
    'P': 30973.762,
    'Al3': 78003.6,
    'PP': 30973.762,
    'Ca2': 80156.0,
    'CO2': 44009.5,
    'POC': 12010.7,
    'OMb': 12010.7}


def load_data():
    mat_contents = sio.loadmat('../IO/MyLakeResults.mat')
    MyLake_results = mat_contents['MyLake_results']
    Sediment_results = mat_contents['Sediment_results']
    return MyLake_results, Sediment_results


def plot_intime(results, elem):
    plt.figure(figsize=(6, 4), dpi=192)
    for e in elem:
        plt.plot(-366 + results['days'][0, 0][0], results[e][0, 0].T, lw=3, label=e)
    ax = plt.gca()
    ax.ticklabel_format(useOffset=False)
    ax.grid(linestyle='-', linewidth=0.2)
    plt.legend()
    plt.tight_layout()
    ax = plt.gca()
    ax.ticklabel_format(useOffset=False)
    ax.xaxis.set_major_locator(mdates.MonthLocator(interval=2))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b\n%Y'))
    plt.show()


class ResultsPlotter:
    """docstring for ResultsPlotter"""

    def __init__(self, years_ago=1):
        MyLake_results, Sediment_results = load_data()
        self.myLake_results = MyLake_results
        self.sediment_results = Sediment_results
        self.years_ago = years_ago

    def env_getter(self, env):
        if env == 'sediment':
            results = self.sediment_results
        elif env == 'water-column':
            results = self.myLake_results
        return results

    def plot_flux(self, elem, lbl):
        results = self.sediment_results
        plt.figure(figsize=(6, 4), dpi=192)
        start = -365 * self.years_ago
        end = -365 * (self.years_ago - 1) - 1
        try:
            plt.plot(results['days'][0, 0][0][start:end] - 366, results['Bioirrigation_fx_zt'][0, 0][elem]
                     [0, 0][0][start:end], sns.xkcd_rgb["medium green"], lw=3, label='Bioirrigation')
            plt.plot(results['days'][0, 0][0][start:end] - 366, results['sediment_SWI_fluxes'][0, 0]
                     [elem][0, 0][0][start:end], sns.xkcd_rgb["denim blue"], lw=3, label='Diffusive')
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

    def plot_fit(self):
        fig, axes = plt.subplots(4, 1, sharex='col', figsize=(8, 5), dpi=192)

        inx = sum(self.myLake_results['z'][0, 0] < 4)[0]
        TOTP = np.mean(self.myLake_results['Pzt'][0, 0][0:inx, :], axis=0) + \
            np.mean(self.myLake_results['PPzt'][0, 0][0:inx, :], axis=0) + \
            np.mean(self.myLake_results['DOPzt'][0, 0][0:inx, :], axis=0) + \
            np.mean(self.myLake_results['Chlzt'][0, 0][0:inx, :], axis=0) + \
            np.mean(self.myLake_results['Czt'][0, 0][0:inx, :], axis=0)
        # np.mean(self.myLake_results['DOCzt'][0, 0][0:inx, :], axis=0) + \
        # np.mean(self.myLake_results['POCzt'][0, 0][0:inx, :], axis=0) + \
        Chl = np.mean(self.myLake_results['Czt'][0, 0][0:inx, :], axis=0) + np.mean(self.myLake_results['Chlzt']
                                                                                    [0, 0][0:inx, :], axis=0)
        PO4 = np.mean(self.myLake_results['Pzt'][0, 0][0:inx, :], axis=0)
        Part = np.mean(self.myLake_results['PPzt'][0, 0][0:inx, :], axis=0)
        # + np.mean(self.myLake_results['POCzt'][0, 0][0:inx, :], axis=0)

        axes[0].plot(-366 + self.myLake_results['days'][0, 0][0], TOTP, c=sns.xkcd_rgb["denim blue"], lw=3, label='Total P')
        axes[1].plot(-366 + self.myLake_results['days'][0, 0][0], Chl, c=sns.xkcd_rgb["denim blue"], lw=3, label='Chl-a')
        axes[2].plot(-366 + self.myLake_results['days'][0, 0][0], PO4, c=sns.xkcd_rgb["denim blue"], lw=3, label='PO_4')
        axes[3].plot(-366 + self.myLake_results['days'][0, 0][0], Part, c=sns.xkcd_rgb["denim blue"], lw=3, label='Solid P')
        # axes[4].plot(-366 + self.myLake_results['days'][0, 0][0], np.mean(self.myLake_results['DOPzt'][0, 0][0:inx, :], axis=0), lw=3, label='DOPzt')
        # axes[4].plot(-366 + self.myLake_results['days'][0, 0][0], np.mean(self.myLake_results['Pzt'][0, 0][0:inx, :], axis=0), lw=3, label='Pzt')
        # axes[4].plot(-366 + self.myLake_results['days'][0, 0][0], np.mean(self.myLake_results['POCzt'][0, 0][0:inx, :], axis=0), lw=3, label='POCzt')

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

        # axes[-1].grid(linestyle='-', linewidth=0.2)
        # axes[-1].set_ylabel(r'$mg / m^3$')
        # axes[-1].legend(loc=1)

    def unit_converter(self, convert_units, env, e):
        if convert_units:
            if env == 'sediment':
                coef = molar_masses[e[:-2]]
                units = '[$mg/L$]'
            elif env == 'water-column':
                coef = 1 / molar_masses[e[:-2]]
                units = '[$mmol/L$]'
        else:
            coef = 1
            if env == 'sediment':
                units = '[$mmol/L$]'
            elif env == 'water-column':
                units = '[$mg/L$]'
        return coef, units

    def plot_profile(self, env, elem, convert_units=False):
        results = self.env_getter(env)
        plt.figure(figsize=(6, 4), dpi=192)
        for e in elem:
            coef, units = self.unit_converter(convert_units, env, e)
            plt.plot(results[e][0, 0][:, -1] * coef, -results['z'][0, 0][:, -1], lw=3, label=e[:-2])
        plt.xlabel(units)
        plt.ylabel('Depth, m')
        ax = plt.gca()
        ax.ticklabel_format(useOffset=False)
        ax.grid(linestyle='-', linewidth=0.2)
        plt.legend()
        plt.ylim([-results['z'][0, 0][-1], -results['z'][0, 0][0]])
        plt.tight_layout()
        plt.show()

    def contour_plot(self, env, elem, convert_units=False, cmap=ListedColormap(sns.color_palette("Blues", 51))):
        results = self.env_getter(env)
        plt.figure(figsize=(6, 4), dpi=192)
        start = -365 * self.years_ago
        end = -365 * (self.years_ago - 1) - 1
        X, Y = np.meshgrid(results['days'][0, 0][0][start:end] - 366, -results['z'][0, 0])
        z = 0
        for e in elem:
            coef, units = self.unit_converter(convert_units, env, e)
            z += results[e][0, 0][:, start:end] * coef
        CS = plt.contourf(X, Y, z, 51, cmap=cmap, origin='lower')
    #     plt.clabel(CS, inline=1, fontsize=10, colors='w')
        cbar = plt.colorbar(CS)

        if env == 'water-column':
            ice_thickness = results['His'][0, 0][0, start:end]
            plt.fill_between(results['days'][0, 0][0][start:end] - 366, 0, -ice_thickness, where=-ice_thickness <= 0, facecolor='red', interpolate=True)

        plt.ylabel('Depth, [m]')
        ax = plt.gca()
        ax.ticklabel_format(useOffset=False)
        ax.xaxis.set_major_locator(mdates.MonthLocator(interval=2))
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%b\n%Y'))
        lbl = ''
        for e in elem:
            lbl += e[:-2] + ', '
        cbar.ax.set_ylabel(lbl + units)
        if elem[0] == 'Tzt':
            cbar.ax.set_ylabel('Temperature, [C]')
        if elem[0] == 'H_sw_zt':
            cbar.ax.set_ylabel('Light limiting function (group 1), [-]')
        if elem[0] == 'H_sw_zt_2':
            cbar.ax.set_ylabel('Light limiting function(group 2), [-]')
        plt.show()
