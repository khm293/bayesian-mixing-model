# bayesian-mixing-model
stan model for endmember mixing analysis

Contents

1. mixing_model_generic.R
This is a generic code that can be modified for different numbers of tracers/endmembers (always satisfying n+1 endmembers and n tracers). User must specifiy priors and likelihood if they don't want the defaults. Inputs include endmember data distribution (right now default is normal (mean plus sd) and vector of stream concentrations for which the mixing model is performed. 

2. 3comp_2tracer_student.R
This is the code used in Markovich et al. (2019) which has two tracers and three endmembers, includes a varying effect for time (weekly samples), and uses a student-t distribution for the Si tracer to capture the thick tails. 
