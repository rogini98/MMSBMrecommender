#----------------------------------------------------
# run alignment function on the 500 independent runs
# to find the permuted labels
# then permute the labels in raw mmsbm outputs
#----------------------------------------------------

import os
import numpy as np
import sys
import copy
from itertools import combinations
from tqdm import tqdm
import pandas as pd
from collections import OrderedDict
#import alignment_cutre as bbm

# run alignment function on the 500 independent runs to find the permuted labels
#results = bbm.run_alignments(reference = 149, rows = 6, columns = 8, perspective_column = 1)

def permute_proba_mat(
	reference = 41,
    	folderin = "../../../results/2.process-selected-mmsbm/results-region/rc-10by8/raw/prob-gp-gp/", 
	folderout = "../../../results/2.process-selected-mmsbm/results-region/rc-10by8/processed/reordered/proba-gp-gp/", 
	rows = 10, 
	columns = 8, 
	nsamples = 500):
    
    # read in the best permutations from alignments
    results = np.load('../../../results/2.process-selected-mmsbm/results-region/rc-10by8/processed/alignments/41proba10gp-8gp.dat.npy')
    
    for iii in tqdm(range(0, nsamples)):
        
        # if sampling no == reference, do not permute labels
        if iii==reference:
            continue

        # read in the sample and reshape it in matrix format
        mmat = [x.strip().split() for x in open(folderin+str(iii)+"proba10gp-8gp.dat", "r").readlines()]
        ## select the 2nd column and format it
        mmat = np.asarray(mmat)[:,0].reshape((rows,columns)).astype(np.float)
        print mmat
        #print mmat.shape

        # extract permuted labels of sample from results
        permuted_labels = results[1:, iii]
        #print permuted_labels
        print permuted_labels[:rows] # permuted labels of rows
        print permuted_labels[rows:] # permuted labels of columns

        # now reorder sample 
        reordered_mat = mmat[permuted_labels[:rows],:] # reorder rows of 
        reordered_mat = mmat[:, permuted_labels[rows:]] # reorder columns of proba mat

        #print reordered_mat
        #reordered_mat = pd.DataFrame(reordered_mat)
        #reordered_mat.to_csv(folderout+str(iii)+"reordered_proba10gp-8gp.csv", index = True, header = True, sep = "\t")

        np.savetxt(folderout+str(iii)+"reordered_proba10gp-8gp.txt", reordered_mat)

def permute_user_gp(
    reference = 41,
	folderin = "../../../results/2.process-selected-mmsbm/results-region/rc-10by8/raw/user-gp/", 
	folderout = "../../../results/2.process-selected-mmsbm/results-region/rc-10by8/processed/reordered/user-gp/", 
	rows = 213484, 
	columns = 10, 
	nsamples = 500):
    
    # read in the best permutations from alignments
    results = np.load('../../../results/2.process-selected-mmsbm/results-region/rc-10by8/processed/alignments/41proba10gp-8gp.dat.npy')
    #print results
    
    for iii in tqdm(range(0, nsamples)):
    	#print iii
        if iii==reference:
            continue
        
        # read in the sample and reshape it in matrix format
        mmat = [x.strip().split() for x in open(folderin+str(iii)+"_usergp10.dat", "r").readlines()]

        # get the user id
        user_id = np.asarray(mmat)[:,0]
        user_id = list(OrderedDict.fromkeys(user_id))
        
        ## select the 2nd column [contains proba of gp membership] 
        mmat = np.asarray(mmat)[:,1].reshape((rows,columns)).astype(np.float)
        #print mmat.shape
        #print mmat[1,:]

        # extract permuted labels of sample from results
        permuted_labels = results[1:, iii]
        #print permuted_labels
        #print permuted_labels[:columns] # permuted labels of columns
        
        # now reorder sample 
        reordered_mat = mmat[:, permuted_labels[:columns]] # reorder columns of proba mat
        #print reordered_mat[1,:]
	np.savetxt(folderout+str(iii)+"reordered_usergp10.txt", reordered_mat)
        # as data frame to add the rownames       
        #reordered_mat = pd.DataFrame(reordered_mat, index= user_id)

        # save all the reordered matrices  
        #reordered_mat.to_csv(folderout+"/"+str(iii)+"reordered_usergp10.csv", sep= "\t", index = True)

def permute_place_gp(
    reference = 41,
	folderin = "../../../results/2.process-selected-mmsbm/results-region/rc-10by8/raw/place-gp/", 
	folderout = "../../../results/2.process-selected-mmsbm/results-region/rc-10by8/processed/reordered/place-gp/", 
	rows = 17, 
	columns = 8, 
	nsamples = 500):
    
    # read in the best permutations from alignments
    results = np.load('../../../results/2.process-selected-mmsbm/results-region/rc-10by8/processed/alignments/41proba10gp-8gp.dat.npy')

    for iii in tqdm(range(0, nsamples)):
        #print iii
        
        if iii==reference:
            continue

        # read in the sample and reshape it in matrix format
        mmat = [x.strip().split() for x in open(folderin+str(iii)+"_placegp8.dat", "r").readlines()]
        
        # get the place id
        place_id = np.asarray(mmat)[:,0]
        place_id = list(OrderedDict.fromkeys(place_id))
        print len(place_id) # yay ca marche poulette!
        print place_id 
        

        ## select the 1st column and format it
        mmat = np.asarray(mmat)[:,1].reshape((rows,columns)).astype(np.float)
        print mmat.shape

        # extract permuted labels of sample from results
        permuted_labels = results[1:, iii]
        print permuted_labels
        print permuted_labels[-columns:] # permuted labels of columns
        #print permuted_labels[columns:]
        # now reorder sample 
        reordered_mat = mmat[:, permuted_labels[-columns:]] # reorder columns of proba mat
        #reordered_mat = pd.DataFrame(reordered_mat)
        #print reordered_mat
        #reordered_mat.shape

        # as data frame to add the rownames       
        #reordered_mat2 = pd.DataFrame(reordered_mat.values, index= place_id)
        #print reordered_mat
        #print reordered_mat2
        # save as csv 
        #np.save(folderout+str(iii)+"reordered_placegp8", reordered_mat)
        #reordered_mat2.to_csv(folderout+"/"+str(iii)+"rc_reordered_placegp8.csv", sep= "\t", index = True)
        # w/0 place name
        #reordered_mat.to_csv(folderout+"/"+str(iii)+"reordered_placegp8.csv", index = False, header = True, sep= "\t", )
        np.savetxt(folderout+str(iii)+"reordered_placegp8.txt", reordered_mat)
