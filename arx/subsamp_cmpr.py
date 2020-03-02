
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
from scipy import optimize

from .post_process.lts_subsampling import subsamp_lts

def powerlaw_func_s(t, gamma, s):
    return (t/s+1)**(-1*gamma)


def synsrv_trace_figure(dpath):

    fig, ax = pl.subplots()
    fig.set_size_inches(5.2,3)


    for k in [7, 12, 22]:
        dts, synsrv_prb = subsamp_lts(dpath, k)

        ax.plot(dts, synsrv_prb, 'o', markeredgewidth=1,
                markerfacecolor='None')

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')

    ax.set_xscale('log')
    ax.set_yscale('log')

    ax.set_xlabel('simulation steps')
    ax.set_ylabel('survival probability')

    # ax.set_ylim(10**(-6),1)



    directory = 'figures/subsamp_cmpr/'

    if not os.path.exists(directory):
        os.makedirs(directory)

    ax.legend(frameon=False, loc='lower left',
                  prop={'size': 9})

    fname = dpath[-4:]

    fig.savefig(directory+'/'+fname+'.png', dpi=300,
                bbox_inches='tight')




if __name__ == "__main__":

    data_dirs = sorted(['data/'+pth for pth in next(os.walk("data/"))[1]])

    for dpath in data_dirs[1:]:
        print("excluded dpath 0000!!")
        synsrv_trace_figure(dpath)
