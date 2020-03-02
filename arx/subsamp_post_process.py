
import os, pickle
import numpy as np

from .post_process.lts_subsampling import subsamp_lts


def subsamp_and_save(dpath):

    df = []

    for k in [5, 7, 10, 12, 15, 17, 20, 22, 25, 27, 30, 32]:
    # for k in [5,10]:

        dts, synsrv_prb = subsamp_lts(dpath, k)

        df_entry = {'k' : k, 'dts' : dts,
                    'synsrv_prb' : synsrv_prb,
                    'dpath': dpath}

        df.append(df_entry)

        print(df_entry)

        
    fpath = 'data/' + dpath[-4:] + '/synsrv_prb.p'

    with open(fpath, 'wb') as pfile:
        pickle.dump(df, pfile)




if __name__ == "__main__":

    data_dirs = sorted(['data/'+pth for pth in next(os.walk("data/"))[1]])

    for dpath in data_dirs[1:]:
        print("excluded dpath 0000!!")
        subsamp_and_save(dpath)
