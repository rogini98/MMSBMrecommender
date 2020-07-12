# MMSBMrecommender for visitation networks 

## Identify group membership of visitors and places
To characterise the visitor--place network, we used the Mixed Membership Stochastic Block Model inference approach of Godoy-Lorite et al. (2016). 
The original code for mixed-membership stochastic block model algorithm was taken from Godoy-Lorite et al. 2016 (article https://www.pnas.org/content/113/50/14207).
The code is written in Pypy (a Python interpreter)---you can use python just changing 'import _numpypy as np' by 'import numpy as np'. To run the code:
pypy mmsbm_recommender.py training_dataset test_dataset K L (sampling iterations) ; where K is the number of group-memberships for visitors and L is the number of group memberships for places, and "(sampling iterations)" are optional parameters such that default values are 1 for sampling and 200 for iterations.

This model uses an Expectation Maximisation algorithm to infer the model parameters visitor groups membership, place groups membership and the probability of interaction between visitor and place groups which maximise the likelihood of the observed travelling patterns. 
As iterating the Expectation Maximisation algorithm with different initial conditions can lead to different solutions, we performed 500 independent runs ('sampling' in the code). We suggest distributing the sampling process, given that they are independent processes which speed up the computation.

Note that to characterise the visitor groups (K) and place groups (L) from the visitation network we made the following modifications to the original code: 
- length of dictionnaries are used instead of taking the max number of the "visitors"/"places" columns as per original code
- parameters (theta/eta/probability matrix) are saved in separate files for each sampling 

The following output files 'theta.dat', 'eta.dat', 'proba-K-L.dat' are saved for each sampling.

## Find the best matching of group labels
Though the MMSBM allows us to identify the optimal number of visitor groups ($K$) and place groups ($L$), the model does not distinguish between group labels---which is crucial in our case. For instance distinguishing whether a particular visitor group 2 was consistently observed to travel to places of group 1 in the 500 independent runs of the model---especially to characterise the behaviour of the different visitor groups. However, as we had no prior knowledge of "true" labels of the visitor and place groups, we considered the solution having the maximum likelihood as our "true" labels of the visitor groups and place groups. Using our assumed "true" probability matrix, we used a Simulated Annealing algorithm---which is an optimisation function---to find the best matching of the labels for each of the probability matrix obtained from the 500 independent fits of our model. 

## Reordering the probability matrices based on the best matching of group labels





