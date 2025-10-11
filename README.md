# aggregatedbounds
Replication code for "Debiased Machine Learning of Aggregated Intersection Bounds and Other Causal Parameters"

This repository includes replication files for simulations and empirical application sections.
In both cases there will be separate folders as of now, one for non-orthogonal approach and the other for orthogonal.

In the orthogonal approach - there will be two use cases: once for binary outcomes and the other for discrete non-binary outcomes.
The necessary packages are included in the beginning of the code. 
Empirical code should take less than 15 min to run. 

To run the discrete non-bonary regime, make sure to use this in the main file:

ortho_leebounds(
  leedata = myleedata,
  s.hat   = s.hat.sub,
  y.hat   = y.hat.sub,     # unused in discrete path
  flag_binary   = FALSE,
  flag_discrete = TRUE,
  flag_helps    = FALSE,
  ortho_d       = TRUE     # selection-orthogonal denominator
)


As well as uncomment the following versions instead of the binary ones:

res <- basic_lee_bound_discrete(leedata)
...
sharp_res_basic <- basic_lee_bound_discrete(myleedata)


