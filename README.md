# MMSBMrecommender for visitation networks 
Original code for mixed-membership stochastic block model algorithm was taken from Godoy-Lorite et al. 2016 (article https://www.pnas.org/content/113/50/14207). 
To characterise the visitor groups (K) and place groups (L) from the visitation network we made the following modifications to the original code: 

The code is written in Pypy (a Python interpreter)---you can use python just changing 'import _numpypy as np' by 'import numpy as np'. To run the code:
pypy mmsbm_recommender.py training_dataset test_dataset K L (sampling iterations) ; where K is the number of group-memberships for visitors and L is the number of group memberships for places, and "(sampling iterations)" are optional parameters such that default values are 1 for sampling and 200 for iterations.

The algorithm gives slightly different solutions depending on the initialization of the parameters. As none of these solutions is significantly better than the other we perform different random initializations ('sampling' in the code), and the final probability distribution over the ratings is the average over all of them.
We suggest distributing the sampling process, given that they are independent processes which speed up the computation. For the iterations, fixing the iterations save computational time, but to fix the value ensure that likelihood is in a plateau. You can compute the likelihood from time to time to stop the iterations when the likelihood is not growing, but it also takes computational time.

The output file 'predictions.dat' includes:
userID itemID real_rating prediction_rating ratings_probability_distribution
