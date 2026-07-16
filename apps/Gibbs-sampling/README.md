# Gibbs Sampling Laboratory

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
