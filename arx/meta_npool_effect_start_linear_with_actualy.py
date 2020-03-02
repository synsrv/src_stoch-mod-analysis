
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

    print("excluded dpath 0000!!")

    for dpath in data_dirs[1:]:

        print('Loading ', dpath)

        try:

            with open(dpath+'/namespace.p', 'rb') as pfile:
                nsp=pickle.load(pfile)

            with open(dpath+'/synsrv_prb_10.p', 'rb') as pfile:
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
        fig.set_size_inches(4.8,3)

        df_saved = None

        for df in df_all:

            if df['k']==k:

                df_saved = df

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

                    label += ',\n $\gamma = %.4f$' %(prm[0])

                    ax.plot(df['dts'], df['synsrv_prb'], '.', #  'o',
                            # markeredgewidth=1,
                            # markerfacecolor='None',
                            color=bl.get_color(),
                            label=label)

                    # ax.set_xlim(left=0.95*df['dts'][1])

                else:
                    ax.plot(df['dts'], df['synsrv_prb'], # 'o',
                            # markeredgewidth=1,
                            # markerfacecolor='None',
                            label=label)


        bin_w = 1
        with open(df_saved['dpath']+'/lts.p', 'rb') as pfile:
            lts_df=np.array(pickle.load(pfile))

        # discard synapses present at beginning
        lts_df = lts_df[lts_df[:,1]>0]

        # only take synapses grown in first half of simulation
        t_split = df_saved['Nsteps']/2
        lts_df = lts_df[lts_df[:,3]<t_split]

        lts = lts_df[:,2] - lts_df[:,3]

        assert np.min(lts) > 0

        lts[lts>t_split]=t_split

        bins = np.arange(1,t_split+bin_w,bin_w)

        counts, edges = np.histogram(lts,
                                     bins=bins,
                                     density=False)

        srv = 1. - np.cumsum(counts)/float(np.sum(counts))

        centers = (edges[:-1] + edges[1:])/2.
        ax.plot(centers, srv, '.')


        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.yaxis.set_ticks_position('left')
        ax.xaxis.set_ticks_position('bottom')

        # ax.set_xscale('log')
        # ax.set_yscale('log')

        fig.suptitle('$k= %d $' %k)
        fig.tight_layout(rect=[0., 0., 1, 0.95])

        ax.set_xlabel('simulation steps')
        ax.set_ylabel('survival probability')

        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

        ax.legend(frameon=False, prop={'size': 8}, loc='center left',
                  labelspacing=1.15, borderpad=1.25,
                  bbox_to_anchor=(1, 0.5))

        # ax.legend(frameon=False, loc='lower left',
        #               prop={'size': 9})


        directory = 'figures/meta_npool_effect_start_linear/'

        if not os.path.exists(directory):
            os.makedirs(directory)

        fname = "npool_effect_k%d_with-actual" %k

        fig.savefig(directory+'/'+fname+'.png', dpi=300,
                    bbox_inches='tight')




if __name__ == "__main__":

    npool_effect_figures(fit=True)
