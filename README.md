## Simulation Versions

**Mean-density FBM via Binned Gaussian**: `mdfbm_binned_gaussians.f90` - Uses a binned shard mean-density distribution to calculate the gradient, but each walker contributes a little Gaussian over multiple bins in the shared distribution rather than to onlya single bin. 

**Mean-density FBM via Analytic (Unbinned) Gaussian:** `wfBM-prl-genbin.f90` - For a given time and walker in the set, a loop over all walkers evalutating `sum_gaus_derivs(x0=x(it-1), traj(i), it-1)` which computes analytically the derivative at point `x0` will be performed. The gradient will be the sum of derivative contributions from all walkers. In this way, the force a walker feels will be proportional to this analytically computed gradient term effected by all walkers in the set. 

**Mean-density FBM via Simple Binning:** `soft_fbm_parallel1-mill.f90` - Walkers contribute a single entry to the shared-history distribution.


## Mill Cluster Reference

**Compilation settings**
```
module load openmpi/4.1.6/gcc/12.2.0
mpif90 *.f90 *.f -o fbm -mcmodel=large -cpp -fallow-argument-mismatch -ffree-line-length-512
```

**Reference commands**

`sinfo` - lists summary statistics of available nodes 

`squeue` - shows specific jobs on each node - can add `-p vojtat-lab` or `-u [username]` flags to narrow search results 

