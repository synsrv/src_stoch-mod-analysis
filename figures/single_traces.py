
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


def xt_trace_figure(dpath, fit=True):

    with open(dpath+'/namespace.p', 'rb') as pfile:
        nsp=pickle.load(pfile)

    with open(dpath+'/xt.p', 'rb') as pfile:
        xt=np.array(pickle.load(pfile))


    for k in range(nsp['n_trace_rec']):

        pl.clf()
        

        fig, ax = pl.subplots()
        fig.set_size_inches(5.2,3)

        ax.plot(xt[:,k], lw=0.75)

        if nsp['up_cap'] > 0:
            ax.plot([0, nsp['Nsteps']], [nsp['up_cap']]*2, 'black')

        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.yaxis.set_ticks_position('left')
        ax.xaxis.set_ticks_position('bottom')

        ax.set_xlabel('simulation steps $t$')
        ax.set_ylabel('$X(t)$')

        directory = 'figures/xt_traces/single/' + str(dpath[-4:])

        if not os.path.exists(directory):
            os.makedirs(directory)

        # ax.legend(frameon=False, loc='lower left',
        #               prop={'size': 9})


        fig.savefig(directory+'/'+'%.4d.png' %k, dpi=300,
                    bbox_inches='tight')

        pl.close(fig)


if __name__ == "__main__":

    data_dirs = sorted(['data/'+pth for pth in next(os.walk("data/"))[1]])

    for dpath in data_dirs:
        xt_trace_figure(dpath)
