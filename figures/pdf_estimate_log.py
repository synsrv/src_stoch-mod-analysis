
import argparse, sys, os, itertools, pickle, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as pl
from matplotlib import rc

rc('text', usetex=True)
pl.rcParams['text.latex.preamble'] = [
    r'\usepackage{tgheros}',    # helvetica font
    r'\usepackage{sansmath}',   # math-font matching helvetica
    r'\sansmath'                # actually tell tex to use it!
    r'\usepackage{siunitx}',    # micro symbols
    r'\sisetup{detect-all}',    # force siunitx to use the fonts
]  

import numpy as np
from decimal import Decimal
from scipy.stats import norm, lognorm


def pdf_figure():

    nbins = 125

    pl.close()
    fig = pl.figure()

    s = 150
    ax_lines, ax_cols = 1,1
    axs = {}
    for x,y in itertools.product(range(ax_lines),range(ax_cols)):
        axs['%d,%d'%(x+1,y+1)] = pl.subplot2grid((ax_lines, ax_cols), (x, y))

    fig.set_size_inches(1920/s*ax_lines/4,1080/s*ax_cols/3)

    data_dirs = sorted(['data/'+pth for pth in next(os.walk("data/"))[1]])

    for dpath in data_dirs:

        print('Found ', dpath)

        try:
            with open(dpath+'/namespace.p', 'rb') as pfile:
                nsp=pickle.load(pfile)

            with open(dpath+'/kx.p', 'rb') as pfile:
                kx=np.array(pickle.load(pfile))

                
            log_weights = np.log10(kx[:,0][kx[:,0]>0.000001])

            # axs['1,1'].hist(log_weights, bins=nbins,
            #                 density=True)

            text=r'$N_{\text{steps}} = ' + str(nsp['Nsteps']) +'$'

            floc, fscale = norm.fit(log_weights)
            f_rv = norm(loc=floc, scale=fscale)
            xs = np.linspace(start=np.min(log_weights),
                             stop=np.max(log_weights),
                             num = 1000)
            axs['1,1'].plot(xs, f_rv.pdf(xs), lw=1,
                            linestyle='-', label=text)


            # axs['1,1'].text(1.1, 0.95, text,
            #                 horizontalalignment='left',
            #                 verticalalignment='top',
            #                 linespacing = 1.95,
            #                 fontsize=10,
            #                 bbox={'boxstyle': 'square, pad=0.3',
            #                       'facecolor':'white', 'alpha':1,
            #                       'edgecolor':'none'},
            #                 transform = axs['1,1'].transAxes,
            #                 clip_on=False)

        
        except FileNotFoundError:
            print(dpath[-4:], "reports: Error loading namespace")


 
    
    axs['1,1'].spines['right'].set_visible(False)
    axs['1,1'].spines['top'].set_visible(False)
    axs['1,1'].yaxis.set_ticks_position('left')
    axs['1,1'].xaxis.set_ticks_position('bottom')

    axs['1,1'].set_xlabel(r'$\log_{10}(X)$')
    axs['1,1'].set_ylabel('relative frequency')


    # text = r'$\mu_{b}=' + '%.5E $' % Decimal(nsp['bn_mu']) +\
    #        '\n' + r'$\sigma_{b}=' + '%.2E $' % Decimal(nsp['bn_sig']) +\
    #        '\n ---------------- ' +\
    #        '\n' + r'$N_{\text{process}} = ' + str(nsp['Nprocess']) +'$' +\
    #        '\n' + r'$N_{\text{steps}} = ' + str(nsp['Nsteps']) +'$'
           
    box = axs['1,1'].get_position()
    axs['1,1'].set_position([box.x0, box.y0, box.width * 0.8, box.height])

    axs['1,1'].legend(frameon=False, prop={'size': 7}, loc='center left',
              labelspacing=1.15, borderpad=1.25, bbox_to_anchor=(1, 0.5))


    directory = "figures/weights_log/" 
    if not os.path.exists(directory):
        os.makedirs(directory)
        
    pl.savefig(directory+"/pdf_est_all.png",
               dpi=300, bbox_inches='tight')




    
if __name__ == "__main__":

        
    pdf_figure()
