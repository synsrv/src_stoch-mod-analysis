
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



def kx_trace_figure(dpath):


    with open(dpath+'/namespace.p', 'rb') as pfile:
        nsp=pickle.load(pfile)
    with open(dpath+'/kx.p', 'rb') as pfile:
        kxt_df=np.array(pickle.load(pfile))

    fig, ax = pl.subplots()

    for i in range(len(kxt_df.T)):
        ax.plot(kxt_df[:,i], color='grey')

    directory = 'figures/kxt/'

    if not os.path.exists(directory):
        os.makedirs(directory)

    fname = "%s" % dpath[-4:]

    fig.savefig("./"+directory+fname+'.png', dpi=150,
                bbox_inches='tight')



if __name__ == "__main__":

    data_dirs = sorted(['data/'+pth for pth in next(os.walk("data/"))[1]])

    for dpath in data_dirs:
        kx_trace_figure(dpath)

