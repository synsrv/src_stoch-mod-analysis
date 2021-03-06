
import os, pickle
import numpy as np

from .lts_subsampling import subsamp_lts


def subsamp_fixed_start_dt(nsp, lts, dpath):

    df = []

    # for k in [4, 5, 7, 8, 10, 12, 15, 17, 20, 22, 25, 27, 30, 32]:
    for start in [3, 4, 5, 10, 25, 50, 100, 250]:

        for k in [4, 8, 16, 32, 64, 128, 256, 1024]:
                       
            dts, synsrv_prb = subsamp_lts(lts, nsp, k, start=start)

            df_entry = {'k' : k, 'dts' : dts,
                        'synsrv_prb' : synsrv_prb,
                        'start' : start,
                        'dpath': dpath}

            df.append(df_entry)


        
    fpath = 'data/' + dpath[-4:] + '/synsrv_prb_start.p'

    with open(fpath, 'wb') as pfile:
        pickle.dump(df, pfile)




if __name__ == "__main__":

    data_dirs = sorted(['data/'+pth for pth in next(os.walk("data/"))[1]])
    
    for dpath in data_dirs:
        
        print("Loding lts data...")

        with open(dpath+'/namespace.p', 'rb') as pfile:
            nsp=pickle.load(pfile)

        with open(dpath+'/lts.p', 'rb') as pfile:
            lts=np.array(pickle.load(pfile))

        print("Finished loading.")
        
        subsamp_fixed_start_dt(nsp, lts, dpath)
