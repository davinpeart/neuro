# neuro

Helper functions for behavioural neuroscience analyses in R — plotting themes
and palettes, posterior summaries for `brms` fits, correlation-matrix tools, and
a few growth-curve helpers.

This is a working toolkit rather than a general-purpose package. It exists so
that collaborators can install it instead of receiving script files by email.
Parts of it duplicate better-maintained packages, and it is not on CRAN.
Where something here is a thinner version of an established tool, this README
says so.

## Installation

```r
# install.packages("remotes")
remotes::install_github("davinpeart/neuro")
```

`brms`, `rstan`, `bridgesampling`, `patchwork`, `Bessel` and `matrixStats` are
suggested rather than required — install them only if you use the functions that
need them. Everything else is pulled in automatically.

## What is here

Twenty-four exported functions. Each one runs on a clean `library(neuro)` with
nothing else attached.

### Themes and palettes

| Function | Does |
| --- | --- |
| `theme_nature()` | ggplot2 theme, Nature-journal proportions |
| `theme_cell()` | ggplot2 theme, Cell-journal proportions |
| `nature_palette(colour, n)` | 15 hues x 6 shades. `colour = NULL` returns the full list |
| `springer(N, gradn)` | Ordered categorical palettes of length 1–8 |
| `colours5(colour)` | Five saturated hues |
| `seaborn(n)` | Muted four-colour palette |

### Posterior summaries and predictive checks

For `brms` fits.

| Function | Does |
| --- | --- |
| `ribbon_epreds(fit, x, y, ...)` | Posterior expected values by group, as a ribbon plot. Returns the plot, the draws, and the summary |
| `savage.dickey(fit, ...)` | Savage–Dickey density ratio as a Bayes factor, with a prior/posterior plot. **See the note below** |
| `pp_stat_dens(yrep, y)` | Two-statistic posterior predictive check as a 2-D density |
| `enframe_descriptives(y, stat1, stat2)` | Two summary statistics of a vector or draws matrix, as a data frame |
| `enframe_prop_integer(y)` | Proportion of each integer value in a vector or draws array |

`ribbon_epreds` overlaps `tidybayes::add_epred_draws()` with
`ggdist::stat_lineribbon()`, and `pp_stat_dens` overlaps
`bayesplot::ppc_stat_2d()`. Both of those are more general and better tested.
`savage.dickey` overlaps `bayestestR::bayesfactor_parameters()`; what it adds is
an expression interface, so you can test a contrast such as `b_a - b_b` directly.

### Growth and decay curves

| Function | Does |
| --- | --- |
| `nlnb_formula(dv, t, curve, t0)` | Builds a non-linear `brms` formula. `curve` is `"exp"`, `"double exp"` or `"logistic"` |
| `asymptotic_exponential()` | Exponential approach to an asymptote |
| `double_exponential()` | Difference of two exponentials |
| `generalized_logistic()` | Richards curve |
| `logit(p)` | Log odds |

The formulae from `nlnb_formula` are parameterised so every non-linear parameter
is estimated in log space, keeping the mean non-negative. That makes them usable
as-is with Poisson, negative binomial and gamma responses — which is the point
of them, and the reason they are not simply `y ~ SSasymp(...)`.

### Correlations

| Function | Does |
| --- | --- |
| `cor_test_mat(x, method)` | Matrix of correlation p-values |
| `plot_cor_mat(rmat, pmat, ...)` | Correlation matrix as a tile plot with coefficients and significance marks |
| `fisher_r_to_z(r1, r2, n1, n2, return)` | Fisher z test for the difference between two correlations |
| `plot_volcano(cor1, cor2, n1, n2, ...)` | Volcano plot of the differences between two correlation matrices |
| `k_means_scree_plot(matrix, max.k)` | Within-cluster sum of squares against k. Returns the plot and the fitted clusterings |

`plot_volcano` is the one with no close equivalent elsewhere — volcano plots of
*between-group correlation differences* are common in c-Fos network papers and
not packaged anywhere I know of. For the rest, `corrplot`, `ggcorrplot` and
`corrr::rplot()` are more capable than `plot_cor_mat`, and
`factoextra::fviz_nbclust()` is more capable than `k_means_scree_plot`.

### Principal components

| Function | Does |
| --- | --- |
| `plot_importance(pc, ncomp, nret)` | Scree plot of individual and cumulative variance from a `prcomp` object |

`factoextra::fviz_eig()` does the same job with more options.

### Skellam

| Function | Does |
| --- | --- |
| `rskellam(n, mu, sigma)` | Draws from the Skellam distribution, parameterised by mean and excess variance |
| `rskellamzi(n, mu, sigma, zeta)` | Draws from the zero-inflated Skellam |

Both are parameterised by the mean `mu` and `sigma`, the amount by which the
variance exceeds `|mu|` — rather than by the two underlying Poisson rates. The
Skellam requires variance > `|mean|`, and this parameterisation satisfies that
for every value the parameters can take.

## Sharp edges

Things that will bite, in rough order of how likely you are to hit them.

> **`savage.dickey()` needs a proper prior, not just `sample_prior = "yes"`.**
> The Savage–Dickey ratio is the prior density at the null divided by the
> posterior density at the null, so a flat prior gives it nothing to divide by.
> brms's default prior on `class = "b"` is flat, and the failure surfaces as
> `Error: Not enough unique values`, which does not point at the cause. Fit with
> something like `prior = prior(normal(0, 5), class = "b")`.

**`plot_cor_mat()` has `rsize` and `psize` the wrong way round.** The
documentation says `rsize` sizes the correlation coefficients and `psize` sizes
the significance marks. The code does the opposite. The arguments have been left
as they are so that existing calls keep producing the same figures — if your
numbers come out the wrong size, swap the two.

**`k_means_scree_plot(matrix, max.k)` needs `max.k` strictly less than
`nrow(matrix)`.** `kmeans()` cannot fit as many centres as there are rows, and
the resulting error mentions `nrow(x)` without saying which argument to change.

**The `brms` families for the Skellam are not exported.** `skellam()` and
`zero_inflated_skellam()` exist in `R/skellam_functions.R` but are deliberately
unreachable. Their log-likelihood evaluates the Bessel term through two
different routines depending on whether `y == 0`, which makes the returned
values inconsistent between zero and non-zero observations and distorts `loo()`
comparisons; the zero-inflated version additionally collapses all posterior
draws to a single value at `y == 0`. Both problems are being fixed properly in a
separate package. `rskellam()` and `rskellamzi()` are exported and fine.

**`R/stan_functions.R` is unexported and does not work.** `neuro_model()` and
its negative-binomial code generator predate the rest of the package and fail on
their own default arguments. They are kept only because they are the reference
for a rewrite.

## Requirements

R >= 4.1. Tested on R 4.5.3 with ggplot2 4.0.3 and brms 2.23.0.

## Licence

MIT — see `LICENSE`.
