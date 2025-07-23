Sample code for the paper *Improving Sampling By Modifying The Effective Diffusion*, T. Lelièvre, R. Santet, G. Stoltz

The folder includes:
- System and collective variable definitions in `input*.jl` ;
- Thermodynamic integration in `thermodynamic_integration.jl` ;
- MALA, adaptive MALA, RMHMC and RMGHMC algorithms in their associated Julia scripts, relying on the numerical schemes defined in `schemes.jl`. Note that RMHMC and RMHMC include numerical reversibility checks to remove asymptotic bias because of the use of implicit schemes ;
- Rejection probabilities computations in `compute_rejection_probabilities.jl` ;
- Transition times computations in `*_transition_times.jl`. Some hyperthreading on the time step is included in the `*_transition_times_alpha_threading.jl` scripts.