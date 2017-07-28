import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def read_data(f, h=[1, ]):
    df = pd.read_excel(f, 'WaterChemistry', header=h)
    df['SampleDate'] = pd.to_datetime(df['SampleDate'], format='%d.%m.%Y %H:%M:%S', errors='coerce')
    df = df.convert_objects(convert_numeric=True)
    df = df.replace('NaN', np.NaN)
    return df


def plot_depth(df, y, depth, ax=None):
    x = 'SampleDate'
    if not ax:
        ax = plt.gca()
    ax.plot_date(df[(np.isfinite(df[y])) & (df['Depth1'] == depth)][x], df[(np.isfinite(df[y])) & (df.Depth1 == depth)][y], label=y + ' at ' + str(depth) + ' m')
    return ax
