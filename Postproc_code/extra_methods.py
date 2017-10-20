import numpy as np
import datetime
import metrics


def convert_timestamp_to_num(timestamp):
    return (timestamp.date() - datetime.date(1, 1, 1)).days + 367


def find_indeces_of_dates(s, o, calibration_end_date='2010-01-01'):
    """
    find indexes of dates simulated vs observed
    Args:
        s: simulated date
        o: observed date
        calibration_end_date (str, optional): calibration/validation date border

    Returns:
        idxs_before: indexes of observed values in simulated array before calibration_end_date
        idxs_after: indexes of observed values in simulated array after calibration_end_date
    """
    idxs_before = np.array([])
    idxs_after = np.array([])
    border_date_num = convert_timestamp_to_num(datetime.datetime.strptime(calibration_end_date, '%Y-%m-%d'))
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


def run_metrics(days_sim, values_sim, days_obs, values_obs, calibration_end_date='2010-01-01', methods=[metrics.rmse, metrics.correlation, metrics.pc_bias, metrics.apb, metrics.norm_rmse, metrics.mae, metrics.bias, metrics.NS, metrics.likelihood, metrics.index_agreement, metrics.squared_error, metrics.coefficient_of_determination]):
    """ run specific metrics for measured and simulated data
    Args:
        days_sim (array): matlabs datenum date
        values_sim (array): simulated values
        days_obs (array): matlabs datenum date
        values_obs (array): observed values
        calibration_end_date (str, optional): date of end of calibration
        methods (list of python methods, optional): methods of metrics
    """
    idxs_sim_before, idxs_sim_after = find_indeces_of_dates(days_sim, days_obs, calibration_end_date='2010-01-01')
    idxs_obs_before, idxs_obs_after = find_indeces_of_dates(days_obs, days_sim, calibration_end_date='2010-01-01')

    print('{0: <35}'.format('Metrics'), end='')
    print('{0: <30}'.format('During calibration'), end='')
    print('{0: <30}'.format('After calibration'))
    for m in methods:
        print('{0: <30}'.format(m.__name__), end='')
        print('\t {0: <30}'.format(m(values_sim.take(idxs_sim_before), values_obs.take(idxs_obs_before))), end='')
        print('\t {0: <30}'.format(m(values_sim.take(idxs_sim_after), values_obs.take(idxs_obs_after))))


[metrics.rmse, metrics.correlation, metrics.percentage_deviation, metrics.pc_bias, metrics.apb, metrics.norm_rmse, metrics.mae,
    metrics.bias, metrics.NS, metrics.likelihood, metrics.index_agreement, metrics.squared_error, metrics.coefficient_of_determination]
