import numpy as np
import random as rn
import math
import time
import sys
import copy
from itertools import combinations
import os
from tqdm import tqdm

def FunAnneal(n, m, alignment_r, alignment_c):
    m=m[alignment_r,:]
    m=m[:,alignment_c]
    E=np.sum(((n-m)**2)/(n+m))
    #E=1-np.corrcoef(n.ravel(),m.ravel())[0,1]
    m=m[:,np.argsort(alignment_c)]
    m=m[np.argsort(alignment_r),:]
    return E

def Change_alignment(alignment):
    idx, idy = np.random.choice(alignment, size=2, replace=False)
    alignment[idx], alignment[idy] = alignment[idy], alignment[idx]

def update(start, step, T, E, acceptance, improvement, IterMax):

    elapsed = time.time() - start
    if step == 0:
        print(' Temperature        Energy    Accept   Improve     Elapsed   Remaining')
        sys.stdout.write('\r%12.2f  %12.2f                      %s            ' % \
        (T, E, time_string(elapsed)))
        sys.stdout.flush()
    else:
        remain = (IterMax - step) * (elapsed / step)
        sys.stdout.write('\r%12.2f  %12.2f  %7.2f%%  %7.2f%%  %s  %s' % \
            (T, E, 100.0 * acceptance, 100.0 * improvement, time_string(elapsed), time_string(remain))),
        sys.stdout.flush()

def time_string(seconds):
    """Returns time in seconds as a string formatted HHHH:MM:SS."""
    s = int(round(seconds))  # round to nearest second
    h, s = divmod(s, 3600)   # get hours and remainder
    m, s = divmod(s, 60)     # split remainder into minutes and seconds
    return '%4i:%02i:%02i' % (h, m, s)


def align(n, m, Tmax=10, Tmin=0.01, IterMax=100000, updates=1000):
    n = n.astype(np.float)
    m = m.astype(np.float)
    step=0
    start = time.time()
    alignment_r = np.arange(0,np.shape(n)[0])
    alignment_c = np.arange(0,np.shape(n)[1])

    if Tmin<=0.0:
        raise Exception('Exponential cooling requires a minimum temperature greater than zero.')

    Tfactor = -math.log(float(Tmax) / Tmin)
    T = Tmax
    E=FunAnneal(n, m, alignment_r, alignment_c)
    prev_alignment_r = copy.copy(alignment_r)
    prev_alignment_c = copy.copy(alignment_c)
    prevEnergy = E
    best_alignment_r = copy.copy(alignment_r)
    best_alignment_c = copy.copy(alignment_c)
    bestEnergy = E
    trials, accepts, improves = 0, 0, 0

    if updates>0:
        updateWavelength = float(IterMax) / updates
        update(start, step, T, E, None, None, IterMax)

    # Attempt moves to new states
    while step < IterMax:
        step += 1
        T = Tmax * math.exp(Tfactor * step / IterMax)

        if bool(rn.getrandbits(1)):
            Change_alignment(alignment_r)
        else:
            Change_alignment(alignment_c)

        E = FunAnneal(n, m, alignment_r, alignment_c)
        dE = E - prevEnergy
        trials += 1

        if dE > 0.0 and math.exp(-dE / T) < rn.random():
            # Restore previous state
            alignment_r = copy.copy(prev_alignment_r)
            alignment_c = copy.copy(prev_alignment_c)

            E = prevEnergy
        else:
            # Accept new state and compare to best state
            accepts += 1
            if dE < 0.0:
                improves += 1
            prev_alignment_r = copy.copy(alignment_r)
            prev_alignment_c = copy.copy(alignment_c)
            prevEnergy = E
            if E < bestEnergy:
                best_alignment_r = copy.copy(alignment_r)
                best_alignment_c = copy.copy(alignment_c)
                bestEnergy = E
        if updates > 1:
            if step // updateWavelength > (step - 1) // updateWavelength:
                update(start,step, T, E, float(accepts) / trials, float(improves) / trials, IterMax)
                trials, accepts, improves = 0, 0, 0

    if updates>1:
        print('\n')

    # Return best state and energy
    return bestEnergy, best_alignment_r, best_alignment_c

def run_alignments(reference=1, rows=10, columns=10, nsamples=500, perspective_column=0):
    folderin = "../Data/output/raw/proba_mat/"
    folderout = "../Data/output/processed/alignments/"

    nmat=[x.strip().split() for x in open(folderin+str(reference)+"proba10gp-10gp.dat", "r").readlines()]
    nmat=np.asarray(nmat)[:,perspective_column].reshape((rows,columns)).astype(np.float)

    results=np.zeros((1+rows+columns,nsamples)).astype(int)
    results[1:,reference]=np.arange(0, rows+columns)
    results[0,reference]=0

    for iii in tqdm(range(0, nsamples)):
        if iii==reference:
            continue

        mmat=[x.strip().split() for x in open(folderin+"/"+str(iii)+"proba10gp-10gp.dat", "r").readlines()]
        mmat=np.asarray(mmat)[:,perspective_column].reshape((rows,columns)).astype(np.float)

        bestEnergy, best_alignment_r, best_alignment_c = align(nmat, mmat, updates=0)

        results[0,iii]=bestEnergy
        results[1:,iii]=np.concatenate((best_alignment_r, best_alignment_c), axis=0)

    np.save(folderout+str(reference)+"proba10gp-10gp.npy", results)
