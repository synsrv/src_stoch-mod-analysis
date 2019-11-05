
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
    lts_df = lts_df[lts_df[:,4]==1]

    lts = lts_df[:,2] - lts_df[:,3]


    if fit:

        fit = powerlaw.Fit(lts)
        label = '$\gamma = %.4f$, $x_{\mathrm{min}}=%.1f$' %(fit.power_law.alpha, fit.power_law.xmin)


        figPDF = powerlaw.plot_pdf(
                   lts[lts>fit.power_law.xmin],
                   label=label,alpha=1.)


        # def pwl(x, alph, xmin):
        #     return (alph - 1) * xmin**(alph-1)*x**(-1*alph)

        # figPDF.plot(np.arange(100, 30000, 1),
        #             pwl(np.arange(100, 30000, 1),
        #                 fit.alpha, fit.xmin))

        fit.power_law.plot_pdf(ax=figPDF, linestyle='--')



        # powerlaw.plot_pdf(lts, linear_bins=True,
        #                   ax=figPDF, label=label)




        # manual method
        # bins = np.logspace(np.log10(fit.power_law.xmin),
        #                    np.log10(np.max(lts)),
        #                    15)
        # cts, edgs = np.histogram(
        #                lts[lts>=fit.power_law.xmin],
        #                density=True, bins=bins)

        # centers = (edgs[:-1] + edgs[1:]) / 2

        # figPDF.plot(centers, cts)


    else:
        pass

    pl.legend()

    directory = 'figures/lts_traces/'

    if not os.path.exists(directory):
        os.makedirs(directory)
                    
    fname = dpath[-4:]

    fig.savefig(directory+'/'+fname+'.png', dpi=300,
                bbox_inches='tight')



if __name__ == "__main__":

    data_dirs = sorted(['data/'+pth for pth in next(os.walk("data/"))[1]])

    for dpath in data_dirs:
        powerlaw_figures(dpath,fit=True)
