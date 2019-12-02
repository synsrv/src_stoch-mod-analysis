
import os, pickle
import numpy as np

from .lts_subsampling import subsamp_lts, log_subsamp_lts
from .true_lts_subsmp_pop import get_true_lts


def subsamp_equal_dt(nsp, lts, dpath):

    df, df_nw = [], []
    

    # for k in [4, 5, 7, 8, 10, 12, 15, 17, 20, 22, 25, 27, 30, 32, 100,
    #           250, 500, 1000]:
    for k in [5, 10, 25, 50, 100, 250, 2500]:
        
        dts, synsrv_prb, lts_subsmp = subsamp_lts(lts, nsp, k, start=0)

        log_dts, log_synsrv_prb, log_lts_subsmp = log_subsamp_lts(lts,
                                                                  nsp, k,
                                                                  start=0)

        df_entry = {'k' : k, 'dts' : dts,
                    'synsrv_prb' : synsrv_prb,
                    'lts' : lts_subsmp,
                    'log_dts': log_dts,
                    'log_synsrv_prb': log_synsrv_prb,
                    'log_lts_subsmp': log_lts_subsmp,
                    'start' : 0,
                    'dpath': dpath}

        df.append(df_entry)

        
        df_newins = get_true_lts(lts, nsp, k, start=0)
        
        df_entry = {'k' : k, 'dts' : dts,
                    'start' : 0,
                    'df_newins' : df_newins}

        df_nw.append(df_entry)

        
    fpath = 'data/' + dpath[-4:] + '/synsrv_prb_equal.p'
    with open(fpath, 'wb') as pfile:
        pickle.dump(df, pfile)

    fpath = 'data/' + dpath[-4:] + '/true_lts_equal.p'
    with open(fpath, 'wb') as pfile:
        pickle.dump(df_nw, pfile)




if __name__ == "__main__":

    data_dirs = sorted(['data/'+pth for pth in next(os.walk("data/"))[1]])
    
    for dpath in data_dirs:

        print("Loding lts data...")
  
        with open(dpath+'/namespace.p', 'rb') as pfile:
            nsp=pickle.load(pfile)

        with open(dpath+'/lts.p', 'rb') as pfile:
            lts=np.array(pickle.load(pfile))

        print("Finished loading.")

       
        subsamp_equal_dt(nsp, lts, dpath)
