
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as pl
from matplotlib import rc

matplotlib.rc('text', usetex=True)
pl.rcParams['text.latex.preamble'] = [
    r'\usepackage{tgheros}',    
    r'\usepackage[eulergreek]{sansmath}',   
    r'\sansmath',
    r'\usepackage{siunitx}',    
    r'\sisetup{detect-all}'
]  

import os, pickle
import numpy as np
from scipy import optimize


def powerlaw_func_s(t, gamma, s):
    return (t/s+1)**(-1*gamma)


def npool_effect_figures(fit=False):

    data_dirs = sorted(['data/'+pth for pth in
                        next(os.walk("data/"))[1]])

    df_all = []

    for dpath in data_dirs:

        print('Loading ', dpath)

        try:

            with open(dpath+'/namespace.p', 'rb') as pfile:
                nsp=pickle.load(pfile)

            with open(dpath+'/synsrv_prb_equal.p', 'rb') as pfile:
                synsrv_prbs=np.array(pickle.load(pfile))

            for synsrv in synsrv_prbs:
                concat = {**nsp, **synsrv}
                df_all.append(concat)

        except FileNotFoundError:
            
            print(dpath[-4:], "reports: No namespace or " +\
                              "synsrv_prb data. Skipping.")



    all_k = np.unique([df['k'] for df in df_all])

    for k in all_k:
        
        fig, ax = pl.subplots()
        fig.set_size_inches(5.2,3)

        for df in df_all:

            if df['k']==k:

                label = "$f = "+str(df['Npool']/df['Nprocess'])+"$"
               

                if fit:
                    prm, prm_cov = optimize.curve_fit(powerlaw_func_s,
                                              df['dts'], df['synsrv_prb'], 
                                              p0=[0.5, 0.5])

                    xs = np.arange(df['dts'][0],
                                   df['dts'][-1],
                                   1)
                    
                    bl, = ax.plot(xs, 
                                  powerlaw_func_s(xs,*prm),
                                  linestyle='-', alpha=0.55)

                    label += ', $\gamma = %.4f$ \n $s=%.4f$' %(prm[0], prm[1])

                    ax.plot(df['dts'], df['synsrv_prb'], #  'o',
                            marker='.', linestyle='None',
                            markersize=3,
                            # markeredgewidth=1,
                            # markerfacecolor='None',
                            color=bl.get_color(),
                            label=label)

                    ax_ll = 0.9*df['dts'][1]
                    
                else:
                    ax.plot(df['dts'], df['synsrv_prb'], # 'o',
                            # markeredgewidth=1,
                            # markerfacecolor='None',
                            label=label)


        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.yaxis.set_ticks_position('left')
        ax.xaxis.set_ticks_position('bottom')

        ax.set_xlabel('simulation steps')
        ax.set_ylabel('survival probability')

        # ---------------------------------------------- 
     
        directory = 'figures/npool_effect_equal_linear/'

        if not os.path.exists(directory):
            os.makedirs(directory)


        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

        ax.legend(frameon=False, prop={'size': 7}, loc='center left',
                  labelspacing=1.15, borderpad=1.25, bbox_to_anchor=(1, 0.5))
        

        fname = "npool_effect_k%d" %k

        fig.savefig(directory+'/'+fname+'.png', dpi=300,
                    bbox_inches='tight')

        # ---------------------------------------------- 

        directory = 'figures/npool_effect_equal/'

        if not os.path.exists(directory):
            os.makedirs(directory)

        ax.set_xscale('log')
        ax.set_yscale('log')

        # ax.legend(frameon=False, loc='lower left',
        #           prop={'size': 9})


        fname = "npool_effect_k%d" %k

        fig.savefig(directory+'/'+fname+'.png', dpi=300,
                    bbox_inches='tight')


        # ---------------------------------------------- 
        
        directory = 'figures/npool_effect_equal_trimmed/'

        if not os.path.exists(directory):
            os.makedirs(directory)

        ax.set_xlim(left=ax_ll)

        fname = "npool_effect_k%d_trimmed" %k

        fig.savefig(directory+'/'+fname+'.png', dpi=300,
                    bbox_inches='tight')


        


if __name__ == "__main__":

    npool_effect_figures(fit=True)
