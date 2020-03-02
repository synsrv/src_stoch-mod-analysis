
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
        n_alive = len(alive_synapse_ids)

        dts, subsmp_srvprb, subsmp_lts = [1.], [1.0], []

        for j in range(3, n_samples+1):

            alive_until_j = df[(df['t_ins'] < j*dt-1) &
                               (df['t_elim'] >= j*dt-1) &
                               (df['pid'].isin(alive_synapse_ids))]

            subsmp_lts.extend((n_alive - len(alive_until_j)) * [(j-2)*dt])
            

            alive_synapse_ids = alive_until_j['pid']
            n_alive = len(alive_synapse_ids)

            if len(new_synapse_ids) == 0:
                srv_prb=1.0
            else:
                srv_prb = len(alive_synapse_ids)/len(new_synapse_ids)

                
            dts.append((j-2)*dt+1)
            subsmp_srvprb.append(srv_prb)

        return dts, subsmp_srvprb, subsmp_lts, []
    

    else:

        # start is the value of the first two
        # initial sampling steps

        dt = start

        active_at_dt = df[(df['t_ins']< dt-1) &
                          (df['t_elim'] >= dt-1)]['pid']
        active_at_2dt = df[(df['t_ins']< 2*dt-1) &
                           (df['t_elim'] >= 2*dt-1)]['pid']

        new_synapse_ids = np.setdiff1d(active_at_2dt, active_at_dt)
        alive_synapse_ids = new_synapse_ids
        n_alive = len(alive_synapse_ids)

        # also store true lts
        active_syn = df[(df['t_ins'] >= dt-1) &
                        (df['t_ins'] < 2*dt-1) &
                        (df['t_elim'] >= 2*dt-1) &
                        (df['pid'].isin(alive_synapse_ids))]

        np.testing.assert_array_equal(np.sort(active_syn['pid']),
                                      np.sort(new_synapse_ids))
        
        true_lts = active_syn['t_elim']-active_syn['t_ins']
        

        dts, subsmp_srvprb, subsmp_lts = [1.], [1.0], []

        # now use dt as derived from n_samples
        dt = int((nsp['Nsteps']-2*start)/(n_samples-2))

        for j in range(1, n_samples-1):

            alive_until_j = df[(df['t_ins'] < 2*start+j*dt-1) &
                               (df['t_elim'] >= 2*start+j*dt-1) &
                               (df['pid'].isin(alive_synapse_ids))]

            subsmp_lts.extend((n_alive - len(alive_until_j)) * [j*dt])
            

            alive_synapse_ids = alive_until_j['pid']
            n_alive = len(alive_synapse_ids)

            if len(new_synapse_ids) == 0:
                srv_prb=1.0
            else:
                srv_prb = len(alive_synapse_ids)/len(new_synapse_ids)

                
            dts.append(j*dt+1)
            subsmp_srvprb.append(srv_prb)

        return dts, subsmp_srvprb, subsmp_lts, true_lts
    

    



def log_subsamp_lts(lts, nsp, n_samples, start=0):
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
        n_alive = len(alive_synapse_ids)

        subsmp_srvprb, subsmp_lts = [1.0], []

        dts = np.logspace(start=np.log10(1),
                          stop=np.log10(nsp['Nsteps']+1-2*dt),
                          num=n_samples-1)
        
        for j in range(len(dts)-1):

            alive_until_j = df[(df['t_ins']-2*dt < dts[j]) &
                               (df['t_elim']-2*dt >= dts[j]) &
                               (df['pid'].isin(alive_synapse_ids))]

            subsmp_lts.extend((n_alive - len(alive_until_j)) * [np.round(dts[j+1], 3)])
            

            alive_synapse_ids = alive_until_j['pid']
            n_alive = len(alive_synapse_ids)

            if len(new_synapse_ids) == 0:
                srv_prb=1.0
            else:
                srv_prb = len(alive_synapse_ids)/len(new_synapse_ids)
       
            subsmp_srvprb.append(srv_prb)

        return np.round(dts,3), subsmp_srvprb, subsmp_lts
    

    else:

        raise NotYetImplementedError

    

    



