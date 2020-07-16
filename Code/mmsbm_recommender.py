#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Mixed membership stochastic block model (Godoy-Lorite et al. 2016)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#-----------------------------------------------------------------------
# Following modifications has been done to the original code
# 1) Script has been broken down into functions, can now be uploaded as module
# 2) Order of vectors not the same as in the original code
# 3) Names of variables have been changed
# 4) theta.dat; eta.dat; proba-K-L.dat are saved for each sampling
#-----------------------------------------------------------------------

#importing the different modules
import sys
import numpy as np
from math import *
import copy
import random
import csv

def read_files(training, test, zeros_as_null):

    user_dict = {}
    place_dict = {}
    visits = {}
    linksr = []

    file = open(training,"r")
    for line in file:
        about = line.strip().split("\t")
        if int(about[2])!=0 or not zeros_as_null:
            try:
                x = user_dict[about[0]][0]
                user_dict[about[0]][1] += 1
            except KeyError:
                x = len(user_dict)
                user_dict[about[0]] = [x, 1]

            try:
                y = place_dict[about[1]][0]
                place_dict[about[1]][1] += 1
            except KeyError:
                y = len(place_dict)
                place_dict[about[1]] = [y, 1]

            try:
                v = visits[about[2]]
            except KeyError:
                v = len(visits)
                visits[about[2]] = v

            linksr.append([x, y, v])
    file.close()

    file = open(test,"r")
    test=[]
    for line in file:
        about = line.strip().split("\t")
        if int(about[2])!=0 or not zeros_as_null:
            test.append([about[0], about[1], about[2]])
    file.close()

    # TODO: run tests for data

    return test, user_dict, place_dict, linksr, visits

def sampling(c, ofolder, linksr, nsampling, iterations, K, L, R,
             user_dict, n_users, users_denom, place_dict, n_places, places_denom, rat, verbose):

    if verbose:
        sys.stderr.write(" ".join(['sampling',str(c),'\n']))

    #theta = vector containing the different groups to which each user belongs to
    theta = np.random.rand(n_users,K) / users_denom[:,np.newaxis]

    #eta = vector containing the different groups to which each place belongs to
    eta = np.random.rand(n_places,L) / places_denom[:,np.newaxis]

    # 3d matrix containing random probabilities of ratings across user-group and place-group combos
    # NOTE: I have changed the structure of this and related variables!!!
    pr = np.random.rand(R, K, L)

    # normalize the probabilities across ratings
    # should divide by: sum of all ratings corresponding to a group-group combo
    pr = pr / pr.sum(axis=0)

    # create empty containers for the calculations that are made during each iteration
    ntheta = np.zeros((n_users, K))
    neta = np.zeros((n_places, L))
    npr = np.zeros((R, K, L))

    #########################################################################################
    for g in range(iterations):
        if verbose:
            sys.stderr.write(" ".join(['iteration',str(g),'\n']))

        # update the parameters using each observed 'rating'
        for n, m, ra in linksr:

            # calculate the sum of all mixtures for rating ra by
            # multiplying rating probabilities rowwise by the user group
            # membership probabilities and columnwise by the place group
            # membership probabilities

            D = (pr[ra].T * theta[n]).T * eta[m]

            # normalize these values
            a = D / D.sum()

            # update the new (n) parameter estimates
            npr[ra] = npr[ra] + a
            ntheta[n] = ntheta[n] + a.sum(axis=1)
            neta[m] = neta[m] + a.sum(axis=0)

        # normalize the users' membership probabilities across groups
        ntheta = ntheta / users_denom[:,np.newaxis]

        # normalize the places' membership probabilities across groups
        neta = neta / places_denom[:,np.newaxis]

        # normalize the probabilities across ratings
        npr = npr / npr.sum(axis=0)

        # create copies of previous values and zero'd estimates as placeholders
        theta = copy.deepcopy(ntheta)
        eta = copy.deepcopy(neta)
        pr = copy.deepcopy(npr)

        # restart arrays
        ntheta = ntheta*0
        neta = neta*0
        npr = npr*0

    # calculate the likelihood given the probabilities
    Like = 0.
    for n, m, ra in linksr:
        D = 0.
        for l in range(L):
            for k in range(K):
                D = D+theta[n][k]*eta[m][l]*pr[ra][k][l]
        for l in range(L):
            for k in range(K):
                Like = Like+(theta[n][k]*eta[m][l]*pr[ra][k][l])*log(D)/D

    if verbose:
        print Like

    fout = open(ofolder+'likelihood/'+str(c)+'_Likelihood'+str(K)+'gp'+str(L)+'gp.dat',"w")
    fout.write('%s\n' % Like)
    fout.close()

    #theta : gp membership vector of visitors
    theta_vector = open(ofolder+'visitor_gp/'+str(c)+'_visitor_gp'+str(K)+'.dat',"w")
    for i in user_dict.keys():
        row_num = user_dict[i][0]
        for k in range(K):
            theta_vector.write('%s ' % i)
            theta_vector.write('%s ' % theta[row_num][k])
            theta_vector.write('\n')
    theta_vector.close()

    #eta : gp membership vectors of places
    eta_vector = open(ofolder+'place_gp/'+str(c)+'_place_gp'+str(L) +'.dat',"w")
    for j in place_dict.keys():
        row_num = place_dict[j][0]
        for l in range(L):
            eta_vector.write('%s ' % j)
            eta_vector.write('%s ' % eta[row_num][l])
            eta_vector.write('\n')
    eta_vector.close()

    #pr : probability matrix of observed links between gp of users & gp of places
    proba_visit = open(ofolder+'proba_mat/'+str(c)+'proba'+str(K)+'gp-'+str(L)+'gp.dat',"w")
    for k in range(K):
        for l in range(L):
            for r in range(R):
                proba_visit.write('%s ' % pr[r][k][l])
            proba_visit.write('\n')
    proba_visit.close()

    #save the probability distribution for c of the links in rat
    ##i.e. save the probability dist when the rating is 0 or 1.

    """nl = 0

    #This is the test from older versions of the code, but here is not being used
    for about in test:
        for r in range(R):
            pra = 0.
            for k in range(K):
                for l in range(L):
                    pra = pra + theta[user_dict[about[0]]][k]*eta[place_dict[about[1]]][l]*pr[r][k][l]
            rat[nl][r] = rat[nl][r] + pra/nsampling
        nl=nl+1"""


def run_sampling(training,
             test,
             ofolder="../Data/output/raw/",
             K=10,
             L=10,
             nsampling=1,
             iterations=200,
             zeros_as_null=False,
             verbose=True):

    if not zeros_as_null:
        sys.stderr.write("\nCareful! Zeros in your data will represent a type of interaction and will be used in the sampling.\n\n")

    test, user_dict, place_dict, linksr, visits = read_files(training, test, zeros_as_null)

    n_users = len(user_dict)
    n_places = len(place_dict)
    R = len(visits)

    #TODO: figure if I need to clean this or get random sample
    rat = np.zeros((len(test),R))

    users_denom = np.asarray(user_dict.values())
    users_denom = users_denom[users_denom[:,0].argsort(),1]
    places_denom = np.asarray(place_dict.values())
    places_denom = places_denom[places_denom[:,0].argsort(),1]

    for c in range(nsampling):
        sampling(c, ofolder, linksr, nsampling, iterations, K, L, R,
                 user_dict, n_users, users_denom, place_dict, n_places, places_denom, rat, verbose)
    return "done!"
