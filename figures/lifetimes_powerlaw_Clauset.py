
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

import os, pickle, itertools
import numpy as np
from scipy import optimize

import matlab.engine
eng = matlab.engine.start_matlab()

eng.addpath('/home/hoffmann/lab/stoch-mod/analysis-dev/m/',nargout=0)



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

            with open(dpath+'/synsrv_prb_equal.p', 'rb') as pfile:
                synsrv_prbs=np.array(pickle.load(pfile))

            for synsrv in synsrv_prbs:
                concat = {**nsp, **synsrv}
                df_all.append(concat)

        except FileNotFoundError:
            
            print(dpath[-4:], "reports: No namespace or " +\
                              "synsrv_prb data. Skipping.")



    all_Npool = np.unique([df['Npool'] for df in df_all])
    all_k = np.unique([df['k'] for df in df_all])

    for k in all_k[all_k > 10]:
        

        for df in df_all:

            if df['pl_alpha']==1.5 and df['k']==k: 

                label = "$k = "+str(df['k'])+"$"

                dt = df['dts'][1]-df['dts'][0]

                lts_dat = np.array(df['lts'])


                if fit:


                    alpha, bmin, L = eng.bplfit(matlab.double(list(xxx),
                                                matlab.double([1, 10, 100, 1000]),
                                                nargout=3)


                    # def pwl(x, alph, xmin):
                    #     return (alph - 1) * xmin**(alph-1)*x**(-1*alph)
                    

                    
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


        directory = 'figures/lifetimes_powerlaw/'

        # if not os.path.exists(directory):
        #     os.makedirs(directory)
                    
        # pl.legend()
        # pl.savefig('figures/lifetimes_powerlaw/k%d.png' %k,
        #                     bbox_inches='tight')

        # pl.clf()



if __name__ == "__main__":

   powerlaw_figures(fit=True)
