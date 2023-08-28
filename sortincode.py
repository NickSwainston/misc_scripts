from cProfile import label
from turtle import title
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
import math
from matplotlib.cm import ScalarMappable
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import LinearLocator
import pandas as pd

def condense( infile, sedd):
    p, nd, nf, snr, succ = np.loadtxt(infile, usecols=(1,3,5,7,-1), delimiter=' ', unpack=True)
    trials = np.loadtxt(infile, usecols=(9), delimiter=' ', unpack=True, dtype=str)
    dms = []
    trial_nums = []
    for trial in trials:
        # Remove DM string
        trial = trial.lstrip("DM")
        # dm = float(trial[:6])
        # dms.append(dm)
        trial_num = int(trial[6:])
        trial_nums.append(trial_num)
    # with open('fftsearchreporttolredo.sedd', 'r') as f:

    #     # Loop over each line in the file
    #     for line in f:
    #         # Split the line into fields
    #         fields = line.strip().split()
    #         # Extract the last field
    #         succ = fields[-1]
    #         # Convert the field to a float if needed
    #         # last_value = float(succ)
    #         # Do something with the last value
    #         # print(last_value)
    df = pd.DataFrame(
        {'Period (s)': p, 'Nulling fraction': nf, 'Nulling duration': nd, 'S/N': snr, 'Detected': succ, 'Trial number': trial_nums },
        # dtype={'Period (s)': 'float64', 'Nulling fraction': 'float64', 'Nulling duration': 'float64', 'S/N': 'float64', 'Detected': bool},
    )
    start_row = 0
    end_row = 10
    count = df.loc[start_row:end_row, 'Detected'].sum()
    print(df.loc[start_row:end_row, 'Detected'])
    print(count)
    
    column_names = ['Period (s)', 'Nulling fraction', 'Nulling duration', 'S/N', 'Total detections' ]
    summary_df = pd.DataFrame(columns=column_names)
    for thousand in range(2058):
        print(thousand)
        start_row = thousand * 1000
        end_row = thousand * 1000 + 999
        print(f"Checking row {start_row}: {df.loc[start_row, 'Trial number']} == 1")
        assert df.loc[start_row, 'Trial number'] == 1
        print(f"Checking row {end_row}: {df.loc[end_row, 'Trial number']} == 1000")
        try:
            assert df.loc[end_row, 'Trial number'] == 1000
        except AssertionError:
            for temp_trail in range(start_row, end_row):
                print(temp_trail)
                assert df.loc[temp_trail, 'Trial number'] != df.loc[temp_trail+1, 'Trial number']


        summary_df = summary_df.append(
            {'Period (s)': df.loc[start_row, 'Period (s)'],
             'Nulling fraction': df.loc[start_row, 'Nulling fraction'],
             'Nulling duration': df.loc[start_row, 'Nulling duration'],
             'S/N': df.loc[start_row, 'S/N'],
             'Total detections': df.loc[start_row:end_row, 'Detected'].sum(),
            },
            ignore_index=True,
        )
    print(summary_df)
    summary_df.to_csv(sedd+'.csv', index=False)
    np.save(sedd, summary_df.values)