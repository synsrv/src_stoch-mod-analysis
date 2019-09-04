
import sys, os, itertools, pickle

import numpy as np
import pandas as pd



def subsamp_lts(lts, nsp, n_samples, start=0):
    ''' 
    start=0 ---> equal dt method

    start>0 ---> fixed start method, with start_dt=Nsteps/start
    '''

    df = pd.DataFrame(data=lts,
                      columns=['empty', 'c', 't_elim', 't_ins',
                               'elim_during_sim', 'weight', 'pid'])

    df = df.astype({'c': 'int64', 't_elim': 'int64',
                    't_ins': 'int64', 'elim_during_sim': 'int64',
                    'pid': 'int64'})

    if start==0:

        dt = int(nsp['Nsteps']/(n_samples))

        active_at_dt = df[(df['t_ins']< dt-1) &
                          (df['t_elim'] >= dt-1)]['pid']
        active_at_2dt = df[(df['t_ins']< 2*dt-1) &
                           (df['t_elim'] >= 2*dt-1)]['pid']

        new_synapse_ids = np.setdiff1d(active_at_2dt, active_at_dt)
        alive_synapse_ids = new_synapse_ids

        dts, subsmp_srvprb = [0.], [1.0]

        for j in range(3, n_samples+1):

            alive_until_j =df[(df['t_ins'] < j*dt-1) &
                              (df['t_elim'] >= j*dt-1) &
                              (df['pid'].isin(alive_synapse_ids))]

            alive_synapse_ids = alive_until_j['pid']

            if len(new_synapse_ids) == 0:
                srv_prb=0.
            else:
                srv_prb = len(alive_synapse_ids)/len(new_synapse_ids)

                
            dts.append(j*dt-1)
            subsmp_srvprb.append(srv_prb)

        return dts, subsmp_srvprb
    

    else:

        dt_init = int(nsp['Nsteps']/start)
        dt = int((nsp['Nsteps']-(2*dt_init))/(n_samples))

        active_at_1dt_init = df[(df['t_ins']< dt_init) &
                                (df['t_elim'] >= dt_init)]['pid']
        
        active_at_2dt_init = df[(df['t_ins']< 2*dt_init) &
                                (df['t_elim'] >= 2*dt_init)]['pid']

        new_synapse_ids = np.setdiff1d(active_at_2dt_init,
                                       active_at_1dt_init)
        alive_synapse_ids = new_synapse_ids

        dts, subsmp_srvprb = [0.], [1.0]

        for j in range(1, n_samples+1):

            # missing the previous ones!

            alive_until_j =df[(df['t_ins'] < (2*dt_init)+j*dt-1)
                              & (df['t_elim'] >= (2*dt_init)+j*dt-1)
                              & (df['pid'].isin(alive_synapse_ids))]

            alive_synapse_ids = alive_until_j['pid']

            if len(new_synapse_ids) == 0:
                srv_prb=1.0
            else:
                srv_prb = len(alive_synapse_ids)/len(new_synapse_ids)


            subsmp_srvprb.append(srv_prb)
            dts.append(j*dt-1)


        return dts, subsmp_srvprb
    

    



