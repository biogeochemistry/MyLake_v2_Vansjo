"""
Module with methods required for post processed data analysis
"""
#%%
import datetime
import numpy as np
import metrics
import matplotlib.pyplot as plt
import h5py


def integrated_over_depth_masses(results, start=-365 * 20, end=-1):
    res_wc = {}
    # results = h5py.File('/Volumes/Igor EcoHDD/Scenarios/192ts_Fe_20x_200g_full_scen_base_historical_20y17m_sediment_2015_2100.mat','r')
    # results.close()

    z_wc = np.array(results['MyLake_results']['basin1']['z'])[0]*100

    for elem in ['POP', 'DOP', 'PP', 'Chl', 'C', 'P']:
        y = np.array(results['MyLake_results']['basin1']['concentrations'][elem][start:end, :])/31000
        res_wc[elem] = np.mean(np.trapz(y, z_wc))
    res_wc
    res_sed = {}

    z_sed = np.array(results['Sediment_results']['basin1']['z'])[0]
    y = np.array(results['Sediment_results']['basin1']['concentrations']['PO4'][start:end, :])*0.92
    res_sed['PO4'] = np.mean(np.trapz(y, z_sed))

    for elem in ['POP', 'PO4adsa', 'PO4adsb', 'PO4adsc', 'Fe3PO42', 'Ca3PO42']:
        y = np.array(results['Sediment_results']['basin1']['concentrations'][elem][start:end, :]) * (1 - 0.92)
        res_sed[elem+'_sed'] = np.mean(np.trapz(y, z_sed))

    return res_wc, res_sed



def boundary_P_fluxes(results, start=-365 * 20, end=-1):
    inflow_q = np.array(results['MyLake_results']['basin1']['Inflw'][0, start:end])
    surf_P = np.array(results['MyLake_results']['basin1']['concentrations']['P'][start:end, 0])
    surf_POP = np.array(results['MyLake_results']['basin1']['concentrations']['POP'][start:end, 0])
    surf_DOP = np.array(results['MyLake_results']['basin1']['concentrations']['DOP'][start:end, 0])
    surf_PP = np.array(results['MyLake_results']['basin1']['concentrations']['PP'][start:end, 0])
    surf_Phy = np.array(results['MyLake_results']['basin1']['concentrations']['C'][start:end, 0]) + \
        np.array(results['MyLake_results']['basin1']['concentrations']['Chl'][start:end, 0])
    inflow_POP = np.array(results['MyLake_results']['basin1']['Inflw'][22, start:end])
    inflow_TP = np.array(results['MyLake_results']['basin1']['Inflw'][4, start:end])
    inflow_DOP = np.array(results['MyLake_results']['basin1']['Inflw'][5, start:end])
    area = 2.38e+7

    res = {}
    res['P_outflow'] = np.mean(inflow_q * surf_P) / 31 * 1000 / area * 365 / 1e4
    res['POP_outflow'] = np.mean(inflow_q * surf_POP) / 31 * 1000 / area * 365 / 1e4
    res['PP_outflow'] = np.mean(inflow_q * surf_PP) / 31 * 1000 / area * 365 / 1e4
    res['DOP_outflow'] = np.mean(inflow_q * surf_DOP) / 31 * 1000 / area * 365 / 1e4
    res['Phy_outflow'] = np.mean(inflow_q * surf_Phy) / 31 * 1000 / area * 365 / 1e4
    res['POP_inflow'] = np.mean(inflow_q * inflow_POP) / 31 * 1000 / area * 365 / 1e4
    res['DOP_inflow'] = np.mean(inflow_q * inflow_DOP) / 31 * 1000 / area * 365 / 1e4
    res['TP_inflow'] = np.mean(inflow_q * inflow_TP) / 31 * 1000 / area * 365 / 1e4

    return res


def P_wc_rates(results, start=-365 * 20, end=-1):
    z_wc = np.array(results['MyLake_results']['basin1']['z'])*100
    res =  [\
    np.trapz(np.mean(np.array(results['MyLake_results']['basin1']['rates']['Re'][start:end,:]),axis=0), z_wc[0]), \
    np.trapz(np.mean(np.array(results['MyLake_results']['basin1']['rates']['Rc'][start:end,:]),axis=0), z_wc[0]), \
    -np.trapz(np.mean(1/30973.762*np.multiply(np.array(results['MyLake_results']['basin1']['concentrations']['C'][start:end,:]), np.array(results['MyLake_results']['basin1']['rates']['Growth_bioz'][start:end,:])),axis=0), z_wc[0])+ \
    -np.trapz(np.mean(1/30973.762*np.multiply(np.array(results['MyLake_results']['basin1']['concentrations']['Chl'][start:end,:]), np.array(results['MyLake_results']['basin1']['rates']['Growth_bioz_2'][start:end,:])),axis=0), z_wc[0]), \
    np.trapz(np.mean(1/30973.762*np.multiply(np.array(results['MyLake_results']['basin1']['concentrations']['C'][start:end,:]), np.array(results['MyLake_results']['basin1']['rates']['Loss_bioz'][start:end,:])),axis=0), z_wc[0])+ \
    np.trapz(np.mean(1/30973.762*np.multiply(np.array(results['MyLake_results']['basin1']['concentrations']['Chl'][start:end,:]), np.array(results['MyLake_results']['basin1']['rates']['Loss_bioz_2'][start:end,:])),axis=0), z_wc[0]), \
    -np.trapz(np.mean(np.array(results['MyLake_results']['basin1']['rates']['R31a'][start:end,:]),axis=0), z_wc[0]), \
    -np.trapz(np.mean(np.array(results['MyLake_results']['basin1']['rates']['R32a'][start:end,:]),axis=0), z_wc[0]), \
    -2*np.trapz(np.mean(np.array(results['MyLake_results']['basin1']['rates']['R33a'][start:end,:]),axis=0), z_wc[0]), \
    -2*np.trapz(np.mean(np.array(results['MyLake_results']['basin1']['rates']['R34a'][start:end,:]),axis=0), z_wc[0]), \
    ]


    rate_names = [
        'Re', 'Rc', 'Growth_bioz', 'Loss_bioz', 'R31a', 'R32a', 'R33a', 'R34a'
    ]

    for rn, v in zip(rate_names, res):
        print('{}: {:.3f}'.format(rn, v))

    print('sum: {:.3f}'.format(sum(res)))


def P_rates_sediments(results, start=-365 * 20, end=-1):
    F = (1-0.92)/0.92
    z = np.array(results['Sediment_results']['basin1']['z'])
    R31a = -F*np.trapz(np.mean(np.array(results['Sediment_results']['basin1']['rates']['R31a'][start:end,:]),axis=0), z[0])
    R32a = -F*np.trapz(np.mean(np.array(results['Sediment_results']['basin1']['rates']['R32a'][start:end,:]),axis=0), z[0])
    R33a = -2*np.trapz(np.mean(np.array(results['Sediment_results']['basin1']['rates']['R33a'][start:end,:]),axis=0), z[0])
    R33b = +2*F*np.trapz(np.mean(np.array(results['Sediment_results']['basin1']['rates']['R33b'][start:end,:]),axis=0), z[0])
    R34a = -2*np.trapz(np.mean(np.array(results['Sediment_results']['basin1']['rates']['R34a'][start:end,:]),axis=0), z[0])
    R34b = +2*F*np.trapz(np.mean(np.array(results['Sediment_results']['basin1']['rates']['R34b'][start:end,:]),axis=0), z[0])
    R35a = -F*np.trapz(np.mean(np.array(results['Sediment_results']['basin1']['rates']['R35a'][start:end,:]),axis=0), z[0])
    Ra = +F*np.trapz(np.mean(np.array(results['Sediment_results']['basin1']['rates']['Ra'][start:end,:]),axis=0), z[0])
    Rf = +F*np.trapz(np.mean(np.array(results['Sediment_results']['basin1']['rates']['Rf'][start:end,:]),axis=0), z[0])
    Rc = +np.trapz(np.mean(np.array(results['Sediment_results']['basin1']['rates']['Rc'][start:end,:]),axis=0), z[0])
    R3a_P = +200*4*F*np.trapz(np.mean(np.array(results['Sediment_results']['basin1']['rates']['R3a_P'][start:end,:]),axis=0), z[0])
    R3b_P = +1*4*F*np.trapz(np.mean(np.array(results['Sediment_results']['basin1']['rates']['R3b_P'][start:end,:]),axis=0), z[0])
    R3f_P = +116*4*F*np.trapz(np.mean(np.array(results['Sediment_results']['basin1']['rates']['R3f_P'][start:end,:]),axis=0), z[0])
    R3c_P = +200*4*np.trapz(np.mean(np.array(results['Sediment_results']['basin1']['rates']['R3c_P'][start:end,:]),axis=0), z[0])
    R4a_P = +200*4*F*np.trapz(np.mean(np.array(results['Sediment_results']['basin1']['rates']['R4a_P'][start:end,:]),axis=0), z[0])
    R4b_P = +1*4*F*np.trapz(np.mean(np.array(results['Sediment_results']['basin1']['rates']['R4b_P'][start:end,:]),axis=0), z[0])
    R4f_P = +116*4*F*np.trapz(np.mean(np.array(results['Sediment_results']['basin1']['rates']['R4f_P'][start:end,:]),axis=0), z[0])
    R4c_P = +200*4*np.trapz(np.mean(np.array(results['Sediment_results']['basin1']['rates']['R4c_P'][start:end,:]),axis=0), z[0])
    res = [R31a, R32a, R33a, R33b, R34a, R34b, R35a, Ra, Rf, Rc, R3a_P, R3b_P, R3f_P, R3c_P, R4a_P, R4b_P, R4f_P, R4c_P]

    rate_names = ['R31a', 'R32a', 'R33a', 'R33b', 'R34a', 'R34b', 'R35a', 'Ra', 'Rf', 'Rc', 'R3a_P', 'R3b_P', 'R3f_P', 'R3c_P', 'R4a_P', 'R4b_P', 'R4f_P', 'R4c_P']

    for rn, v in zip(rate_names, res):
        print('{}: {:.3f}'.format(rn, v))

    print('sum: {:.3f}'.format(sum(res)))


def savefig(filename, **kwargs):
    plt.savefig(
        '/Users/imarkelo/Google Drive/GDocuments/Vansjo/latex version/img/{}.pgf'.
        format(filename), **kwargs)
    plt.savefig(
        '/Users/imarkelo/Google Drive/GDocuments/Vansjo/latex version/img/{}.pdf'.
        format(filename), **kwargs)
    plt.savefig(
        '/Users/imarkelo/Google Drive/GDocuments/Vansjo/tables and figures/img/{}.pgf'.
        format(filename), **kwargs)
    plt.savefig(
        '/Users/imarkelo/Google Drive/GDocuments/Vansjo/tables and figures/img/{}.pdf'.
        format(filename), **kwargs)




def convert_timestamp_to_num(timestamp):
    return (timestamp.date() - datetime.date(1, 1, 1)).days + 367


def find_indexes_of_dates(s, o, calibration_end_date='2010-01-01'):
    """
    find indexes of dates simulated vs observed
    Args:
        s: simulated date
        o: observed date
        calibration_end_date (str, optional): calibration/validation date
        border

    Returns:
        idxs_before: indexes of observed values in simulated array before
        calibration_end_date
        idxs_after: indexes of observed values in simulated array after
        calibration_end_date
    """
    idxs_before = np.array([])
    idxs_after = np.array([])
    border_date_num = convert_timestamp_to_num(
        datetime.datetime.strptime(calibration_end_date, '%Y-%m-%d'))
    for md in o:
        idx, = np.where(s == md)
        if idx.size > 0:
            if md < border_date_num:
                idxs_before = np.append(idxs_before, int(idx[0]))
            else:
                idxs_after = np.append(idxs_after, int(idx[0]))
    return idxs_before.astype(int), idxs_after.astype(int)


def elements_at_indexes(sv, idxs):
    """
    takes array and list of indexes and return array elements at these indexes
    Args:
        sv (TYPE): array of elements
        idxs (TYPE): indexes

    Returns:
        values: values at specified indexes
    """
    return sv.take(idxs)


def intersaction_of_sim_with_obs(days_sim, values_sim, days_obs, values_obs):
    """Summary

    Args:
        days_sim (TYPE): Description
        values_sim (TYPE): Description
        days_obs (TYPE): Description
        values_obs (TYPE): Description

    Returns:
        s: Description
        o:
    """
    idxs_sim_before, idxs_sim_after = find_indexes_of_dates(
        days_sim, days_obs, calibration_end_date='2010-01-01')
    idxs_sim = np.append(idxs_sim_before, idxs_sim_after)

    idxs_obs_before, idxs_obs_after = find_indexes_of_dates(
        days_obs, days_sim, calibration_end_date='2010-01-01')
    idxs_obs = np.append(idxs_obs_before, idxs_obs_after)
    return values_sim.take(idxs_sim), values_obs.take(idxs_obs)


def run_metrics(days_sim, values_sim, days_obs, values_obs, calibration_end_date='2010-01-01',
                methods=[metrics.mae, metrics.rmse, metrics.correlation, metrics.rsquared, metrics.pc_bias, metrics.likelihood, metrics.NS]):
    """ run specific metrics for measured and simulated data
    Args:
        days_sim (array): matlabs datenum date
        values_sim (array): simulated values
        days_obs (array): matlabs datenum date
        values_obs (array): observed values
        calibration_end_date (str, optional): date of end of calibration
        methods (list of python methods, optional): methods of metrics
    """
    idxs_sim_before, idxs_sim_after = find_indexes_of_dates(
        days_sim, days_obs, calibration_end_date=calibration_end_date)
    idxs_obs_before, idxs_obs_after = find_indexes_of_dates(
        days_obs, days_sim, calibration_end_date=calibration_end_date)

    print('{0: <35}'.format('Metrics'), end='')
    print('{0: <30}'.format('During calibration'), end='')
    print('{0: <30}'.format('After calibration'))
    for m in methods:
        print('{0: <30}'.format(m.__name__), end='')
        print('\t {0: <30}'.format(
            m(values_sim.take(idxs_sim_before), values_obs.take(idxs_obs_before))), end='')
        print('\t {0: <30}'.format(
            m(values_sim.take(idxs_sim_after), values_obs.take(idxs_obs_after))))


def generate_latex_table(days_sim, values_sim, days_obs, values_obs, calibration_end_date='2010-01-01',
    methods=[metrics.mae, metrics.pc_bias, metrics.rmse, metrics.correlation, metrics.coefficient_of_determination]):
    """ run specific metrics for measured and simulated data
    Args:
        days_sim (array): matlabs datenum date
        values_sim (array): simulated values
        days_obs (array): matlabs datenum date
        values_obs (array): observed values
        calibration_end_date (str, optional): date of end of calibration
        methods (list of python methods, optional): methods of metrics
    """
    idxs_sim_before, idxs_sim_after = find_indexes_of_dates(
        days_sim, days_obs, calibration_end_date=calibration_end_date)
    idxs_obs_before, idxs_obs_after = find_indexes_of_dates(
        days_obs, days_sim, calibration_end_date=calibration_end_date)

    for m in methods:
        #     print('{0: <30}'.format(m.__name__), end='')
        print(' & ', end='')
        print('{0:0.2f}'.format(
            m(values_sim.take(idxs_sim_before), values_obs.take(idxs_obs_before))), end=' / ')
        print('{0:0.2f}'.format(
            m(values_sim.take(idxs_sim_after), values_obs.take(idxs_obs_after))), end='')
    print(' \\\\')




