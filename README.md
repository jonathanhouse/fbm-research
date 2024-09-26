`wfBM-prl-genbin.f90` - outer-loop is time-increments and inner-loop is walkers in a set. For every walker in a set, take the `sum_gaus_derivs(x0=x(it-1), traj(i), it-1)` which computes analytically the derivates at point x0 for a given trajectory. This function is placed inside of a for-loop over all walkers in a set, and incremented to a gradient variable. In that way, the loop will produce a gradient term that is the sum of all the gradients for every walker in the set evaluated at the point the current walker is sitting at. 

`wfbm-bin-to-unbin.f90` - this code attempts to study the behavior of switching from a binned gradient calculation (like ASYM5050) to an unbinned version (like Gaussians, etc.). 

`soft_fbm_parallel1-mill.f90` - this code whould calculate for binned versions of the gradient. I want to verify that this code reverts back to that 4/3 exponent that was observed in the past. 

`gaus-unbinned-wfbm.f90` - work in progress code that attempts to still use a binned version of the gradient, but also fills in bins for each walker proportional to the weights of some Gaussian 


On Mill settings: 

```
module load openmpi/4.1.6/gcc/12.2.0
mpif90 *.f90 *.f -o fbm -mcmodel=large -cpp -fallow-argument-mismatch -ffree-line-length-512
```

`sinfo` - lists summary statistics of available nodes 
`squeue` - shows specific jobs on each node - can add `-p vojtat-lab` or `-u username` to narrow search results 

