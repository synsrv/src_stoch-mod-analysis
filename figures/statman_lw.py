
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as pl
from matplotlib import rc

rc('text', usetex=True)
pl.rcParams['text.latex.preamble'] = [
    r'\usepackage{tgheros}',   
    r'\usepackage{sansmath}',  
    r'\sansmath'               
    r'\usepackage{siunitx}',   
    r'\sisetup{detect-all}',   
]  

import argparse, sys, os, itertools, pickle
import numpy as np
import pandas as pd
from scipy import optimize


def powerlaw_func_s(t, gamma, s):
    return (t/s+1)**(-1*gamma)


def statman_lw(dpath, fit=True):

    fig, ax = pl.subplots()
    fig.set_size_inches(5.2,3)

    with open(dpath+'/namespace.p', 'rb') as pfile:
        nsp=pickle.load(pfile)

    with open(dpath+'/lts.p', 'rb') as pfile:
        lts_df=np.array(pickle.load(pfile))


    df = pd.DataFrame(data=lts_df,
                      columns=['empty', 'c', 't_elim', 't_ins',
                               'elim_during_sim', 'weight', 'pid'])

    df = df.astype({'c': 'int64', 't_elim': 'int64',
                    't_ins': 'int64', 'elim_during_sim': 'int64',
                    'pid': 'int64'})


    dt = 4*48 # =4 days!
    n_samples = 13

    active_at_dt = df[(df['t_ins']< dt-1) &
                      (df['t_elim'] >= dt-1)]['pid']
    active_at_2dt = df[(df['t_ins']< 2*dt-1) &
                       (df['t_elim'] >= 2*dt-1)]['pid']

    new_synapse_ids = np.setdiff1d(active_at_2dt, active_at_dt)
    alive_synapse_ids = new_synapse_ids

    dts, subsmp_srvprb = [0.], [1.0]

    for k,j in enumerate(range(3, n_samples+1)):

        alive_until_j =df[(df['t_ins'] < j*dt-1) &
                          (df['t_elim'] >= j*dt-1) &
                          (df['pid'].isin(alive_synapse_ids))]

        alive_synapse_ids = alive_until_j['pid']

        if len(new_synapse_ids) == 0:
            srv_prb=1.0
        else:
            srv_prb = len(alive_synapse_ids)/len(new_synapse_ids)


        dts.append((k+1)*dt)
        subsmp_srvprb.append(srv_prb)


    print(dts, subsmp_srvprb)

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')

    ax.set_xscale('log')
    ax.set_yscale('log')

    ax.set_xlabel('simulation steps')
    ax.set_ylabel('survival probability')


    directory = 'figures/statman_lw/'

    if not os.path.exists(directory):
        os.makedirs(directory)

    ax.legend(frameon=False, loc='lower left',
                  prop={'size': 9})

    fname = dpath[-4:]

    fig.savefig(directory+'/'+fname+'.png', dpi=300,
                bbox_inches='tight')




if __name__ == "__main__":

    data_dirs = sorted(['data/'+pth for pth in next(os.walk("data/"))[1]])

    for dpath in data_dirs:
        statman_lw(dpath)
