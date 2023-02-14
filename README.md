# Scaling-of-whole-brain-resting-state-dynamics

- PRG_function.m:

This function implements the Phenomenological Renormalization Group
(PRG) method, introduced by Meshulam et al. (2018, 2019). It also can
implement the extended version of PRG based on connectivity.

Within this method, the collective activity is iteratively coarse-grained 
by grouping maximally correlated variables (or maximally coupled variables 
in the connectivity-based case). At each coarse-graining step 
k=0,1,…,kmax, clusters of size K=2^k are built, resulting in a system of 
N/K coarse-grained variables and successively ignoring degrees of
freedom. This code computes several observables of the coarse-grained 
variables as a function of K.   


Refs:
Meshulam, L., Gauthier, J.L., Brody, C.D., Tank, D.W. & Bialek, W. 
Coarse graining, fixed points, and scaling in a large population of neurons. 
Phys. Rev. Lett. 123, 178103 (2019).

Meshulam, L., Gauthier, J.L., Brody, C.D., Tank, D.W. & Bialek, W. 
Coarse-graining and hints of scaling in a population of 1000+ neurons. 
arXiv, 1812.11904 (2018).

- metropolis_spin_model.m

Simulates the spin model using the Metropolis algorithm 
and applies the PRG method to the model's activity.

- Connectivity_matrices.mat

Contains connectivity and distances matrices
