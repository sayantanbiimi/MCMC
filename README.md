# Metropolis–Hastings Sampler demo 
Deployed at: https://sayantanbiimi.github.io/MCMC/Metropolis-Hastings/

An interactive Shiny app for demonstrating both general and random-walk
Metropolis–Hastings.

The app uses only `shiny` and base R, making it suitable for Shinylive and
GitHub Pages. It includes four targets:

- standard normal;
- a separated bimodal normal mixture;
- a highly correlated bivariate normal;
- a banana-shaped bivariate target.
- a gamma target with an asymmetric log-normal random walk;
- a beta target with a Beta(2,2) independence proposal.

The last two examples retain the complete Hastings correction and show why the
proposal-density ratio cannot generally be discarded.

Diagnostics include trace plots, acceptance rate, autocorrelation functions,
effective sample size, retained-sample plots, and running summaries.

## Classwork

1. For the standard normal target, compare proposal SDs 0.1, 1, and 6.
2. Explain why the chain with SD 0.1 accepts often but has low ESS.
3. For the bimodal mixture, start at -4 and try proposal SDs 1 and 5.
4. Use the correlated target to show that acceptance rate alone misses poor geometry.
5. Use the banana target to motivate adaptive, gradient-based, or geometry-aware MCMC.

# Gibbs Sampling demo
Deployed at: https://sayantanbiimi.github.io/MCMC/gibbs-sampling/

This focused Shiny app contains two examples:

1. A bivariate normal target with editable means, standard deviations and correlation.
2. A conjugate normal hierarchical model with editable data-generating parameters,
   group structure, and prior hyperparameters.

For the hierarchical model, users may choose between:

- an independent Normal prior for the grand mean with an inverse-gamma prior for the between-group variance;
- a coupled Normal–inverse-gamma prior.

The app provides trace plots, ACF, effective sample sizes, group-level recovery,
hierarchical shrinkage, and prior-versus-posterior comparisons. It uses only
`shiny` and base R and is suitable for Shinylive deployment.

```r
install.packages("shiny")
shiny::runApp("gibbs-sampling-lab")
```

