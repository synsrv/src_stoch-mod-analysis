
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


def powerlaw_func_s(t, gamma, s):
    return (t/s+1)**(-1*gamma)


def powerlaw_figures(fit=False, manual=False):

    data_dirs = sorted(['data/'+pth for pth in
                        next(os.walk("data/"))[1]])

    df_all = []

    for dpath in data_dirs:

        print('Loading ', dpath)

        try:

            with open(dpath+'/namespace.p', 'rb') as pfile:
                nsp=pickle.load(pfile)

            with open(dpath+'/true_lts_equal.p', 'rb') as pfile:
                true_lts_df=np.array(pickle.load(pfile))

            for true_lts in true_lts_df:
                concat = {**nsp, **true_lts}
                df_all.append(concat)

        except FileNotFoundError:
            
            print(dpath[-4:], "reports: No namespace or " +\
                              "synsrv_prb data. Skipping.")



    all_Npool = np.unique([df['Npool'] for df in df_all])
    all_k = np.unique([df['k'] for df in df_all])

    for k in all_k[all_k >= 10]:
        
        # fig, ax = pl.subplots()
        # fig.set_size_inches(5.2,3)

        for df in df_all:

            if df['pl_alpha']==1.5 and df['k']==k: 
            # if df['Npool']==npool and df['k'] in [4, 5, 7, 8, 10]:
            # if df['Npool']==npool and df['k'] in [100, 1000]:
            # if df['Npool']==npool and df['k'] in [1000]:
            

                label = "$k = "+str(df['k'])+"$"

                dt = df['dts'][1]-df['dts'][0]

                # print(dt)
                # print(df['synsrv_prb'])

                lts_dat = np.array(df['df_newins'][:,2] - df['df_newins'][:,3])
                # lts_dat[lts_dat > (df['Nsteps']-dt)] = df['Nsteps']-dt
                # print(lts_dat)

                # lts_dat=np.trim_zeros(lts_dat)

                # with open(dpath+'/lts.p', 'rb') as pfile:
                #     lts_df=np.array(pickle.load(pfile))

                # # discard synapses present at beginning
                # lts_df = lts_df[lts_df[:,1]>0]
    
                # # only take synapses grown in first half of simulation
                # t_split = nsp['Nsteps']/2
                # lts_df = lts_df[lts_df[:,3]<t_split]
                # lts_dat = lts_df[:,2] - lts_df[:,3]



                if fit:
                    # prm, prm_cov = optimize.curve_fit(powerlaw_func_s,
                    #                           df['dts'], df['synsrv_prb'], 
                    #                           p0=[0.5, 0.5])

                    # xs = np.arange(df['dts'][0],
                    #                df['dts'][-1],
                    #                1)


                    
                    # bl, = ax.plot(xs, 
                    #               powerlaw_func_s(xs,*prm),
                    #               linestyle='-', alpha=0.55)

                    # label += ', $\gamma = %.4f$' %(prm[0])

                    fit = powerlaw.Fit(lts_dat, xmin=df['dts'][1])
                    label += ', $\gamma = %.4f$, $x_{\mathrm{min}}=%.1f$' %(fit.power_law.alpha, fit.power_law.xmin)
                    print(df['k'],df['Npool'])


                    figPDF = powerlaw.plot_pdf(
                               lts_dat[lts_dat>fit.power_law.xmin],
                               label=label,alpha=0.2)


                    # def pwl(x, alph, xmin):
                    #     return (alph - 1) * xmin**(alph-1)*x**(-1*alph)
                    
                    # figPDF.plot(np.arange(100, 30000, 1),
                    #             pwl(np.arange(100, 30000, 1),
                    #                 fit.alpha, fit.xmin))
                    
                    fit.power_law.plot_pdf(ax=figPDF, linestyle='--')

                                 
                   
                    # powerlaw.plot_pdf(lts_dat, linear_bins=True,
                    #                   ax=figPDF, label=label)

                    
                    if manual:
                        bins = np.logspace(np.log10(fit.power_law.xmin),
                                           np.log10(np.max(lts_dat)),
                                           15)
                        cts, edgs = np.histogram(
                                       lts_dat[lts_dat>=fit.power_law.xmin],
                                       density=True, bins=bins)

                        centers = (edgs[:-1] + edgs[1:]) / 2

                        figPDF.plot(centers, cts)


                    # ax_ll = 0.9*df['dts'][1]
                    
                else:
                    pass
                    # print(df['synsrv_prb'])
                    # ax.plot(df['dts'], df['synsrv_prb'], '.', # 'o',
                    #         # markeredgewidth=1,
                    #         # markerfacecolor='None',
                    #         label=label)
                    # ax_ll = 0.9*df['dts'][1]


        directory = 'figures/true_lts/'

        if not os.path.exists(directory):
            os.makedirs(directory)
                    
        pl.legend()
        pl.savefig(directory+'k%d.png' %k,
                            bbox_inches='tight')

        pl.clf()


        # ax.spines['right'].set_visible(False)
        # ax.spines['top'].set_visible(False)
        # ax.yaxis.set_ticks_position('left')
        # ax.xaxis.set_ticks_position('bottom')

        # ax.set_xlabel('simulation steps')
        # ax.set_ylabel('survival probability')

        # # ---------------------------------------------- 

        # directory = 'figures/lifetimes_powerlaw_linear/' 

        # if not os.path.exists(directory):
        #     os.makedirs(directory)

        # box = ax.get_position()
        # ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

        # ax.legend(frameon=False, prop={'size': 7}, loc='center left',
        #           labelspacing=1.15, borderpad=1.25, bbox_to_anchor=(1, 0.5))


        # fname = "lifetimes_powerlaw_Npool%d" %(npool)

        # fig.savefig(directory+'/'+fname+'.png', dpi=300,
        #             bbox_inches='tight')


        # # ---------------------------------------------- 

        # directory = 'figures/lifetimes_powerlaw/' 

        # if not os.path.exists(directory):
        #     os.makedirs(directory)
        
        # ax.set_xscale('log')
        # ax.set_yscale('log')

        # fname = "lifetimes_powerlaw_Npool%d" %(npool)

        # fig.savefig(directory+'/'+fname+'.png', dpi=300,
        #             bbox_inches='tight')

        # # ---------------------------------------------- 

        # directory = 'figures/lifetimes_powerlaw_trimmed/' 

        # if not os.path.exists(directory):
        #     os.makedirs(directory)
        
        # ax.set_xlim(left=ax_ll)

        # fname = "lifetimes_powerlaw_Npool%d" %(npool)

        # fig.savefig(directory+'/'+fname+'.png', dpi=300,
        #             bbox_inches='tight')



        # pl.close(fig)



if __name__ == "__main__":

   powerlaw_figures(fit=True)
