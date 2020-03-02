
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

import os, pickle, itertools, powerlaw
import numpy as np
from scipy import optimize



def powerlaw_figures(dpath, fit=False):

    fig, ax = pl.subplots()
    fig.set_size_inches(5.2,3)


    data_dirs = sorted(['data/'+pth for pth in
                        next(os.walk("data/"))[1]])

    with open(dpath+'/namespace.p', 'rb') as pfile:
        nsp=pickle.load(pfile)

    with open(dpath+'/lts.p', 'rb') as pfile:
        lts_df=np.array(pickle.load(pfile))

    # discard synapses present at beginning
    lts_df = lts_df[lts_df[:,1]>0]

    # discard synapses still alive at end of simulation
    # lts_df = lts_df[lts_df[:,4]==1]
    t_split = nsp['Nsteps']/2
    lts_df = lts_df[lts_df[:,3]<t_split]

    lts = lts_df[:,2] - lts_df[:,3]
    lts[lts > t_split] = t_split
    
    fit = powerlaw.Fit(lts, xmin=1, xmax=t_split+1)
    label = '$\gamma = %.4f$, $x_{\mathrm{min}}=%.1f$' %(fit.power_law.alpha,
                                                         fit.power_law.xmin)

    figPDF = powerlaw.plot_pdf(
               lts[lts>fit.power_law.xmin],
               label=label,alpha=1.)

    fit.power_law.plot_pdf(ax=figPDF, linestyle='--')

    lts = lts_df[:,2] - lts_df[:,3]

    fit2 = powerlaw.Fit(lts, xmin=1)
    label = '$\gamma = %.4f$, $x_{\mathrm{min}}=%.1f$' %(fit2.power_law.alpha,
                                                         fit2.power_law.xmin)

    powerlaw.plot_pdf(
               lts[lts>fit2.power_law.xmin],
               label=label,alpha=1.)

    fit2.power_law.plot_pdf(ax=figPDF, linestyle='--')


    pl.legend()

    directory = 'figures/lts_traces_half-split-necessary/'

    if not os.path.exists(directory):
        os.makedirs(directory)
                    
    fname = dpath[-4:]

    fig.savefig(directory+'/'+fname+'.png', dpi=300,
                bbox_inches='tight')



if __name__ == "__main__":

    data_dirs = sorted(['data/'+pth for pth in next(os.walk("data/"))[1]])

    for dpath in data_dirs:
        powerlaw_figures(dpath,fit=True)
