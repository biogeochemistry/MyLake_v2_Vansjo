import datetime
from calendar import timegm
from datetime import datetime

import matplotlib.dates as dt
import matplotlib.pyplot as plt
import numpy as np
from pytz import utc


# Convert a unix time u to plot time p, and vice versa
def plottm(u):
    return dt.date2num(datetime.fromtimestamp(u, utc))


def unixtm(p):
    return timegm(dt.num2date(p, utc).utctimetuple())


plottm = np.vectorize(plottm)
unixtm = np.vectorize(unixtm)


def plot_wc(model_results,
            elements,
            depth,
            ax=None,
            dstart='2005-03-07',
            dend='2011-03-07',
            factor=1,
            **kwargs):

    env = 'water'
    results = model_results.env_getter(env)

    inx = np.where(results['z'][0, 0] == depth)[0][0]

    y = 0

    if elements == ['T']:
        y = results['T'][0, 0][inx, :]
        lbl = 'Temperature, C'
    else:
        for e in elements:
            y += results['concentrations'][0, 0][e][0, 0][inx, :] * factor
        lbl = ", ".join(elements) + ' concentration'

    if not ax:
        ax = plt.gca()

    ax.plot(-366 + results['days'][0, 0][0], y, **kwargs)
    # ax.ticklabel_format(useOffset=False)
    ax.grid(linestyle='-', linewidth=0.2)
    plt.tight_layout()
    dend = datetime.datetime.strptime(dend, '%Y-%m-%d')
    dstart = datetime.datetime.strptime(dstart, '%Y-%m-%d')
    plt.xlim(dstart, dend)
    ax.set_xlabel('Date')
    ax.set_ylabel(lbl)
    legend = ax.legend(title='Depth: ' + str(depth) + 'm', frameon=1)
    plt.setp(legend.get_title(), fontsize='xx-small')
    return ax


def get_data_wc(model_results, elements, depth):

    env = 'water'
    results = model_results.env_getter(env)

    inx = np.where(results['z'][0, 0] == depth)[0][0]

    y = 0

    if elements == ['T']:
        y = results['T'][0, 0][inx, :]
    else:
        for e in elements:
            y += results['concentrations'][0, 0][e][0, 0][inx, :]

    return -366 + results['days'][0, 0][0], y
