import matplotlib.pyplot as plt
import numpy as np
import matplotlib.dates as mdates
import pandas as pd
import seaborn as sns
from matplotlib.colors import ListedColormap
import scipy.io as sio
from matplotlib import rc
from scipy.interpolate import UnivariateSpline

sns.set_style("whitegrid")
sns.set_style("ticks")

rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})
rc('text', usetex=True)


solid_species = ['OMzt', 'OMbzt', 'FeOOHzt', 'FeSzt', 'S8zt', 'FeS2zt', 'AlOH3zt', 'Ca3PO42zt', 'PO4adsazt', 'PO4adsbzt', 'OMSzt', 'FeOH3zt']
solute_species = ['O2zt', 'NO3zt', 'SO4zt', 'NH4zt', 'Fe2zt', 'H2Szt', 'S0zt', 'PO4zt', 'Ca2zt',
                  'HSzt', 'Hzt', 'OHzt', 'CO2zt', 'CO3zt', 'HCO3zt', 'NH3zt', 'H2CO3zt', 'DOM1zt', 'DOM2zt']
molar_masses = {
    'O2': 31998.8,
    'Ox': 31998.8,
    'OM1': 30973.762,
    'DOP': 30973.762,
    'OM': 30973.762,
    'OM2': 12010.7,
    'OMb': 12010.7,
    'DOM1': 30973.762,
    'DOM2': 12010.7,
    'NO3': 62004,
    'FeOH3': 106867.0,
    'Fe3': 106867.0,
    'SO4': 96062,
    'NH4': 18038,
    'Fe2': 55845,
    'H2S': 34080.9,
    'HS': 33072.9,
    'P': 30973.762,
    'PO4adsa': 30973.762,
    'PO4': 30973.762,
    'Al3': 78003.6,
    'PP': 30973.762,
    'Ca2': 80156.0,
    'CO2': 44009.5,
    'POC': 12010.7,
    'OMb': 12010.7}


solid = ['OM', 'OMb', 'FeOH3', 'PO4adsa', 'OMb']
disolved = ['O2', 'DOM1', 'DOM2', 'NO3', 'SO4', 'NH4', 'Fe2', 'H2S', 'HS', 'PO4', 'Al3', 'Ca2', 'CO2']


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

    def __init__(self, years_ago=0.):
        MyLake_results, Sediment_results = load_data()
        self.myLake_results = MyLake_results
        self.sediment_results = Sediment_results
        self.years_ago = years_ago + 1

    def env_getter(self, env, basin=1):
        if env == 'sediment':
            results = self.sediment_results['basin' + str(basin)][0, 0]
        elif env == 'water-column':
            results = self.myLake_results['basin' + str(basin)][0, 0]
        return results

    def unit_converter(self, convert_units, env, e):
        if convert_units:
            if env == 'sediment':
                coef = molar_masses[e[:-2]]
                units = '[$mg/m^3$]'
            elif env == 'water-column':
                coef = 1 / molar_masses[e[:-2]]
                units = '[$mmol/L$]'
        else:
            coef = 1
            if env == 'sediment':
                units = '[$mmol/L$]'
            elif env == 'water-column':
                units = '[$mg/m^3$]'
        return coef, units

    def phosphorus_fit(self):
        fig, axes = plt.subplots(4, 1, sharex='col', figsize=(8, 5), dpi=192)

        results = self.env_getter('water-column', basin=1)

        inx = sum(results['z'][0, 0] < 4)[0]
        TOTP = np.mean(results['Pzt'][0, 0][0:inx, :], axis=0) + \
            np.mean(results['PPzt'][0, 0][0:inx, :], axis=0) + \
            np.mean(results['DOPzt'][0, 0][0:inx, :], axis=0) + \
            np.mean(results['Chlzt'][0, 0][0:inx, :], axis=0) + \
            np.mean(results['Czt'][0, 0][0:inx, :], axis=0)
        # np.mean(results['DOCzt'][0, 0][0:inx, :], axis=0) + \
        # np.mean(results['POCzt'][0, 0][0:inx, :], axis=0) + \
        Chl = np.mean(results['Czt'][0, 0][0:inx, :], axis=0) + np.mean(results['Chlzt']
                                                                        [0, 0][0:inx, :], axis=0)
        PO4 = np.mean(results['Pzt'][0, 0][0:inx, :], axis=0)
        Part = np.mean(results['PPzt'][0, 0][0:inx, :], axis=0)
        # + np.mean(results['POCzt'][0, 0][0:inx, :], axis=0)

        axes[0].plot(-366 + results['days'][0, 0][0], TOTP, c=sns.xkcd_rgb["denim blue"], lw=3, label='Total P')
        axes[1].plot(-366 + results['days'][0, 0][0], Chl, c=sns.xkcd_rgb["denim blue"], lw=3, label='Chl-a')
        axes[2].plot(-366 + results['days'][0, 0][0], PO4, c=sns.xkcd_rgb["denim blue"], lw=3, label='PO_4')
        axes[3].plot(-366 + results['days'][0, 0][0], Part, c=sns.xkcd_rgb["denim blue"], lw=3, label='Solid P')
        # axes[4].plot(-366 + results['days'][0, 0][0], np.mean(results['DOPzt'][0, 0][0:inx, :], axis=0), lw=3, label='DOPzt')
        # axes[4].plot(-366 + results['days'][0, 0][0], np.mean(results['Pzt'][0, 0][0:inx, :], axis=0), lw=3, label='Pzt')
        # axes[4].plot(-366 + results['days'][0, 0][0], np.mean(results['POCzt'][0, 0][0:inx, :], axis=0), lw=3, label='POCzt')

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
        axes[1].set_xlim([732313 - 366, 735234 - 366 * 1])

        for ax in axes:
            ax.grid(linestyle='-', linewidth=0.2)
            ax.set_ylim([0, 50])
            ax.set_ylabel(r'$[mg / m^3]$')
            ax.legend(loc=1)

    def plot_flux(self, elem, convert_units=False, smoothing_factor=False):
        results = self.env_getter('sediment', basin=1)

        plt.figure(figsize=(6, 4), dpi=192)
        start = int(-365 * self.years_ago)
        end = int(-365 * (self.years_ago - 1) - 1)
        x = results['days'][0, 0][0][start:end] - 366
        y = results['sediment_transport_fluxes'][0, 0][elem][0, 0][0][start:end]
        total = {}
        lines = {}
        if convert_units:
            y = y / (molar_masses[elem] * 10**4 / 365 / 10**6)
            lbl = elem + ' flux, $[umol/cm^{2}/y]$'
            total['D'] = np.trapz(y, x / 365)
        else:
            lbl = elem + ' flux, $[mg/m^{2}/d]$'
            total['D'] = np.trapz(y, x)
        if smoothing_factor:
            spl = UnivariateSpline(x, y)
            spl.set_smoothing_factor(smoothing_factor)
            lines['D'], = plt.plot(x, spl(x), sns.xkcd_rgb["denim blue"], lw=3, label='Transport')
        else:
            lines['D'], = plt.plot(x, y, sns.xkcd_rgb["denim blue"], lw=3, label='Transport')
        try:
            b = results['Bioirrigation_fx_zt'][0, 0][elem][0, 0][0][start:end]
            if convert_units:
                b = b / (molar_masses[elem] * 10**4 / 365 / 10**6)
            lines['B'], = plt.plot(x, b, sns.xkcd_rgb["medium green"], lw=3, label='Bioirrigation')
            if convert_units:
                x = x / 365
            total['B'] = np.trapz(b, x)
        except:
            pass

        if convert_units:
            lbl_2 = ' $[umol/cm^{2}/y]$'
        else:
            lbl_2 = ' $[mg/m^{2}/y]$'
        leg1 = plt.legend([lines[e] for e in lines.keys()], ["{:.2f} ".format(total[e]) + lbl_2 for e in total.keys()], loc=4)
        ax = plt.gca()
        ax.set_ylabel(lbl)
        ax.ticklabel_format(useOffset=False)
        ax.xaxis.set_major_locator(mdates.MonthLocator(interval=2))
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%b\n%Y'))
        ax.set_xlim([results['days'][0, 0][0][start:end][0] - 366, results['days'][0, 0][0][start:end][-1] - 366])
        ax.grid(linestyle='-', linewidth=0.2)
        ax.legend(loc=1)
        plt.gca().add_artist(leg1)
        legend = plt.legend(frameon=1, loc=1)
        frame = legend.get_frame()
        frame.set_facecolor('white')
        plt.tight_layout()
        plt.show()

    def plot_profile(self, env, elem, convert_units=False):
        results = self.env_getter(env)
        plt.figure(figsize=(6, 4), dpi=192)
        end = int(-365 * (self.years_ago - 1) - 1)
        z = results['z'][0, 0][:, -1]
        mass_per_area = {}
        lines = {}
        for e in elem:
            coef, units = self.unit_converter(convert_units, env, e)
            y = results[e][0, 0][:, -1 + end] * coef
            lines[e], = plt.plot(y, -z, lw=3, label=e[:-2])
            if convert_units and env == 'sediment':
                mass_per_area[e] = np.trapz(y, z / 100)
                lbl = r'$mg / m^2$'
            elif convert_units and env == 'water-column':
                mass_per_area[e] = np.trapz(y, z * 100)
                lbl = r'$umol/cm^{2}$'
            elif not convert_units and env == 'water-column':
                mass_per_area[e] = np.trapz(y, z)
                lbl = r'$mg / m^2$'
            elif not convert_units and env == 'sediment':
                mass_per_area[e] = np.trapz(y, z)
                lbl = r'$umol/cm^{2}$'
        leg1 = plt.legend([lines[e] for e in elem], ["{:.2f} ".format(mass_per_area[e]) + lbl for e in elem], loc=4)
        plt.xlabel(units)
        if env == 'water-column':
            plt.ylabel('Depth, [m]')
        else:
            plt.ylabel('Depth, [cm]')
        ax = plt.gca()
        ax.ticklabel_format(useOffset=False)
        ax.grid(linestyle='-', linewidth=0.2)
        plt.legend(loc=1)
        plt.gca().add_artist(leg1)
        plt.ylim([-results['z'][0, 0][-1], -results['z'][0, 0][0]])
        plt.tight_layout()
        plt.show()

    def plot_rate_profile(self, env, elem):
        results = self.env_getter(env)
        plt.figure(figsize=(6, 4), dpi=192)
        end = int(-365 * (self.years_ago - 1) - 1)
        z = results['z'][0, 0][:, -1]
        rate_per_area = {}
        lines = {}
        for e in elem:
            y = results['dcdt'][0, 0][e[:-2]][0, 0][:, -1 + end]
            lines[e], = plt.plot(y, -z, lw=3, label=e[:-2])
            if env == 'water-column':
                rate_per_area[e] = np.trapz(y, z * 100)
            elif env == 'sediment':
                rate_per_area[e] = np.trapz(y, z)
        lbl = r'$umol/cm^{2} / y$'
        leg1 = plt.legend([lines[e] for e in elem], ["{:.2f} ".format(rate_per_area[e]) + lbl for e in elem], loc=4)
        plt.xlabel('$umol/cm^{2}/y$')
        if env == 'water-column':
            plt.ylabel('Depth, [m]')
        else:
            plt.ylabel('Depth, [cm]')
        ax = plt.gca()
        ax.ticklabel_format(useOffset=False)
        ax.grid(linestyle='-', linewidth=0.2)
        plt.legend(loc=1)
        plt.gca().add_artist(leg1)
        plt.ylim([-results['z'][0, 0][-1], -results['z'][0, 0][0]])
        plt.tight_layout()
        plt.show()

    def contour_plot(self, env, elem, convert_units=False, cmap=ListedColormap(sns.color_palette("Blues", 51))):
        results = self.env_getter(env)
        plt.figure(figsize=(6, 4), dpi=192)
        start = int(-365 * self.years_ago)
        end = int(-365 * (self.years_ago - 1) - 1)
        X, Y = np.meshgrid(results['days'][0, 0][0][start:end] - 366, -results['z'][0, 0][0:end - 1])
        z = 0
        for e in elem:
            coef, units = self.unit_converter(convert_units, env, e)
            z += results[e][0, 0][0:end - 1, start:end] * coef
        CS = plt.contourf(X, Y, z, 51, cmap=cmap, origin='lower')
    #     plt.clabel(CS, inline=1, fontsize=10, colors='w')
        cbar = plt.colorbar(CS)

        plt.ylabel('Depth, [cm]')

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

    def rate_plot(self, env, elem, convert_units=False, cmap=ListedColormap(sns.color_palette("RdBu_r", 101))):
        results = self.env_getter(env)
        plt.figure(figsize=(6, 4), dpi=192)
        start = int(-365 * self.years_ago)
        end = int(-365 * (self.years_ago - 1) - 1)
        X, Y = np.meshgrid(results['days'][0, 0][0][start:end] - 366, -results['z'][0, 0])
        z = 0
        for e in elem:
            z += results['dcdt'][0, 0][e[:-2]][0, 0][:, start:end]
        # CS = plt.contourf(X, Y, z, 51, cmap=cmap, origin='lower')
        lim = np.max(np.abs(z))
        lim = np.linspace(-lim - 1e-16, +lim + 1e-16, 51)
        CS = plt.contourf(X, Y, z, 51, cmap=ListedColormap(sns.color_palette("RdBu_r", 101)), origin='lower', levels=lim, extend='both')
    #     plt.clabel(CS, inline=1, fontsize=10, colors='w')
        cbar = plt.colorbar(CS)

        plt.ylabel('Depth, [cm]')

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
        if elem[0] == 'Tzt':
            cbar.ax.set_ylabel('Temperature, [C]')
        elif elem[0] == 'H_sw_zt':
            cbar.ax.set_ylabel('Light limiting function (group 1), [-]')
        elif elem[0] == 'H_sw_zt_2':
            cbar.ax.set_ylabel('Light limiting function(group 2), [-]')
        else:
            cbar.ax.set_ylabel('Rate ' + lbl + '$[mmol/L/y]$')
        plt.show()

    def bulk_sediment_profile(self, elem, convert_units=False, cmap=ListedColormap(sns.color_palette("RdBu_r", 101))):
        results = self.env_getter('sediment')
        plt.figure(figsize=(6, 4), dpi=192)
        end = int(-365 * (self.years_ago - 1) - 1)
        z = results['z'][0, 0][:, -1]
        mass_per_area = {}
        lines = {}
        fi = results['params'][0, 0]['fi'][0, 0][:, -1]
        for e in elem:
            coef, units = self.unit_converter(convert_units, 'sediment', e)
            if e[:-2] in disolved:
                y = results[e][0, 0][:, -1 + end] * coef * fi
            elif e[:-2] in solid:
                y = results[e][0, 0][:, -1 + end] * coef * (1 - fi)
            lines[e], = plt.plot(y, -z, lw=3, label=e[:-2])
            if convert_units:
                mass_per_area[e] = np.trapz(y, z / 100)
                lbl = r'$mg / m^2$'
            elif not convert_units:
                mass_per_area[e] = np.trapz(y, z)
                lbl = r'$umol/cm^{2}$'
        leg1 = plt.legend([lines[e] for e in elem], ["{:.2f} ".format(mass_per_area[e]) + lbl for e in elem], loc=4)
        plt.xlabel(units)
        plt.ylabel('Depth, [cm]')
        ax = plt.gca()
        ax.ticklabel_format(useOffset=False)
        ax.grid(linestyle='-', linewidth=0.2)
        plt.legend(loc=1)
        plt.gca().add_artist(leg1)
        plt.ylim([-results['z'][0, 0][-1], -results['z'][0, 0][0]])
        plt.tight_layout()
        plt.show()

    def temperature_fit(self):
        results = self.env_getter('water-column', 2)
        Tzt = results['Tzt'][0, 0]
        t = results['days'][0, 0][0] - 366

        df = pd.read_csv('../obs/vanem_obs/temperature.txt', sep=',')
        df.rename(columns={r'Dato': 'date', 'TemperaturC': 'T', 'Depthm': 'z'}, inplace=True)
        df['date'] = pd.to_datetime(df['date'], errors='coerce')
        df = df.dropna()
        unique_depths = df.z.unique()

        fig, axes = plt.subplots(len(unique_depths), 1, sharex='col', figsize=(8, len(unique_depths) * 1.5), dpi=192)

        for i, d in enumerate(np.sort(df.z.unique())):
            inx = np.where(results['z'][0, 0] == d)[0][0]
            axes[i].plot(t, Tzt[inx, :], c=sns.xkcd_rgb["denim blue"], lw=2)
            axes[i].plot(df[df.z == d].date, df[df.z == d]['T'], 'bo', c=sns.xkcd_rgb["pale red"], markersize=4, label=str(d) + ' m')

        axes[3].xaxis.set_major_locator(mdates.MonthLocator(interval=12))
        axes[3].xaxis.set_major_formatter(mdates.DateFormatter('%b\n%Y'))
        axes[1].set_xlim([732313 - 366, 735234 - 366])

        for ax in axes:
            ax.grid(linestyle='-', linewidth=0.2)
            # ax.set_ylim([0, 5])
            ax.set_ylabel(r'$[C]$')
            legend = ax.legend(frameon=1, loc=1)
            frame = legend.get_frame()
            frame.set_facecolor('white')
