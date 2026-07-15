# Metropolis–Hastings Laboratory

An interactive Shiny app for demonstrating both general and random-walk
Metropolis–Hastings.

## Run locally

```r
install.packages("shiny")
shiny::runApp("metropolis-hastings-lab")
```

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

## Suggested classroom sequence

1. For the standard normal target, compare proposal SDs 0.1, 1, and 6.
2. Explain why the chain with SD 0.1 accepts often but has low ESS.
3. For the bimodal mixture, start at -4 and try proposal SDs 1 and 5.
4. Use the correlated target to show that acceptance rate alone misses poor geometry.
5. Use the banana target to motivate adaptive, gradient-based, or geometry-aware MCMC.
