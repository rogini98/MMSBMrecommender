# MMSBMrecommender for visitation networks 
To characterise the visitor--place network, we used the Mixed Membership Stochastic Block Model inference approach of Godoy-Lorite et al. (2016). 
The original code for mixed-membership stochastic block model algorithm was taken from Godoy-Lorite et al. 2016 (article https://www.pnas.org/content/113/50/14207).

This model uses an Expectation Maximisation algorithm to infer the model parameters visitor groups membership, place groups membership and the probability of interaction between visitor and place groups which maximise the likelihood of the observed travelling patterns (Godoy-Lorite et al. 2016). 
As iterating the Expectation Maximisation algorithm with different initial conditions can lead to different solutions, we performed 500 independent runs ('sampling' in the code). We suggest distributing the sampling process, given that they are independent processes which speed up the computation.

Additionally, as we are particularly interested in characterising the group memberships of visitors and places from the observed visitation network, in our current analysis we fit the MMSBM to identify group membership of visitors and places respectively (using the *mmsbm_recommender.py* script). We then used the Simulated Annealing algorithm (using the *SA_algorithm.py* script) to find the best permuation of the order of the visitor and place group labels (refer to the sections below for further information) and reorganise our model fits (using *permute_labels.py* on output from the *SA_algorithm.py* script) obtained from the MMSBM.

We have also provided an example of the analysis and output generated from the test data set provided.

## 1. Identify group membership of visitors and places by fitting MMSBM
*mmsbm_recommender.py*

Changes made to original code:
Note that we made the following changes to original code **mmsbm_recommender.py**: 
- changed the names of the variables
- split the original script into functions

The code is written in python. To run the code:
python -c "import mmsbm_recommender as mm; mm.run_sampling(nsampling = 1, iterations = 200, training='./Data/raw/training_dataset.dat', test='./Data/raw/test_dataset.dat', ofolder = './Data/output/', K=10, L=10)"; where K is the number of group-memberships for visitors and L is the number of group memberships for places, and both "nsampling" and "iterations" are optional parameters such that default values are 1 for nsampling and 200 for iterations.

Note that, the parameters (theta/eta/probability matrix) are saved in separate files for each sampling. 
The following output files 'visitor_gpK.dat', 'place_gpL.dat', 'proba-Kgp-Lgp.dat' are generated for each of the 500 runs of the model.

## 2. Find the best matching of group labels using Simulated Annealing algorithm
*SA_algorithm.py*

Though the MMSBM allows us to identify the optimal number of visitor groups (*K*) and place groups (*L*), the model does not distinguish between group labels---which is crucial in our case. For instance distinguishing whether a particular visitor group 2 was consistently observed to travel to places of group 1 in the 500 independent runs of the model is important---especially to characterise the behaviour of the different visitor groups. However, as we had no prior knowledge of "true" labels of the visitor and place groups, we considered the solution having the maximum likelihood as our "true" labels of the visitor groups and place groups. Using the aforementioned assumed "true" probability matrix, we used a Simulated Annealing algorithm---which is an optimisation function---to find the best matching of the labels for each of the probability matrix obtained from the 500 independent fits of our model. 

The code is written in python. To run the code:
python -c "import SA_algorithm as sa; sa.run_alignments(reference = 1, rows = 10, columns = 10)"; where *reference* refers to the probability matrix assumed to have the "true" labels; rows and columns refer to the number of rows and columns of the probability matrix.

The following output file is produced "reference_proba_K_L.npy".

## 3. Reordering the group labels of posteriors
*permute_labels.py*

Using the best permutation of the group labels obtained from the previous step, we then reorder the labels of the raw outputs obtained from the MMSBM -- i.e. the visitor group and place groups across the 500 independent model fits. 

The code is written in python. To run the code:
python -c "import permute_labels as perm; perm.permute_proba_mat(); perm.permute_place_gp(); perm.permute_user_gp()"

The following output file is produced "X_reordered_proba_K_L.txt", "X_visitor_gp_K.txt", "X_place_gp_L.txt" where X is the sampling number.
