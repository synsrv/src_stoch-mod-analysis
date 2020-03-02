
import sys, os, itertools, pickle

import numpy as np
import pandas as pd



def get_true_lts(lts, nsp, n_samples, start=0):
    ''' 
    start=0 ---> equal dt method

    start>0 ---> fixed start method, with start_dt=Nsteps/start
    '''

    df = pd.DataFrame(data=lts,
                      columns=['empty', 'c', 't_elim', 't_ins',
                               'elim_during_sim', 'weight', 'pid'])

    df = df.astype({'c': 'int64', 't_elim': 'int64',
                    't_ins': 'int64', 'elim_during_sim': 'int64',
                    'pid': 'int64'})


    if start==0:

        dt = int(nsp['Nsteps']/(n_samples))

        active_at_dt =  np.logical_and(df['t_ins']< dt-1,
                                       df['t_elim'] >= dt-1)
        active_at_2dt = np.logical_and(df['t_ins']< 2*dt-1,
                                       df['t_elim'] >= 2*dt-1)

        df_newins = df[np.logical_and(active_at_2dt,
                                      np.logical_not(active_at_dt))]
        

    else:

        raise NotYetImplementedError



    return np.array(df_newins)
