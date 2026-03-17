# CausalSpline-Stata

**Nonlinear Causal Dose-Response Estimation via Splines**

Stata 14.1 port of the [CausalSpline R package](https://github.com/causalfragility-lab/CausalSpline).

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

---

## Overview

**causalspline** estimates the causal dose-response function E[Y(t)] for a
continuous treatment T under unconfoundedness using restricted cubic splines.

Most causal inference tools assume a linear treatment effect. Real policy and
health problems often involve thresholds, diminishing returns, and nonmonotone
relationships. causalspline recovers these nonlinear structures without
parametric assumptions.

Three identification strategies are implemented:

- **IPW**: Inverse Probability Weighting via the Generalised Propensity Score
- **gcomp**: G-computation (outcome regression)
- **dr**: Doubly Robust (combines both)

A key feature is **geometric fragility diagnostics** -- tools for evaluating
the structural stability of the estimated dose-response curve using first and
second derivatives. These complement traditional sensitivity analyses by
identifying treatment regions where causal relationships are flat, rapidly
changing, or structurally ambiguous.

---

## Installation

### From SSC

```stata
* Coming soon to SSC
* ssc install causalspline
```

### Manual installation

Download the ZIP from this repository and copy all `.ado` and `.sthlp` files
to your personal ado folder:

```stata
sysdir          /* shows your PERSONAL path */
```

Copy all files there, then verify:

```stata
which causalspline
help causalspline
```

---

## Quick Example

```stata
* Simulate threshold dose-response data
cs_simulate 500, dgp(threshold) seed(1) clear

* Fit causal dose-response curve via IPW
causalspline, outcome(Y) treatment(T) confounders(X1 X2 X3) ///
    method(ipw) dfexposure(5) evalgrid(100) bootreps(200) verbose

* Plot dose-response curve
cs_plot

* Overlap diagnostics
cs_overlap

* Derivatives
cs_gradient

* Fragility diagnostics
cs_fragility, type(curvature_ratio)
cs_fragility, type(inverse_slope)

* Regional fragility
cs_region, a(2) b(4)
cs_region, a(4) b(8)
```

---

## Commands

| Command | Description |
|---|---|
| `causalspline` | Estimate the causal dose-response curve |
| `cs_fragility` | Geometric fragility curve with dual-panel plot |
| `cs_gradient` | First and second derivatives of estimated curve |
| `cs_region` | Regional fragility integral over treatment interval |
| `cs_overlap` | Overlap / positivity diagnostics (ESS, histogram) |
| `cs_plot` | Dose-response curve plot with confidence band |
| `cs_simulate` | Simulate nonlinear dose-response datasets |

---

## Supported DGPs (for simulation)

```stata
cs_simulate 500, dgp(threshold)    /* flat then linear rise */
cs_simulate 500, dgp(diminishing)  /* concave / diminishing returns */
cs_simulate 500, dgp(nonmonotone)  /* inverted-U / hump shape */
cs_simulate 500, dgp(linear)       /* linear baseline */
cs_simulate 500, dgp(sinusoidal)   /* oscillatory / complex */
```

---

## Supported Methods

| Method | `method()` | Consistent if ... |
|---|---|---|
| Inverse Probability Weighting | `ipw` | GPS model correct |
| G-computation | `gcomp` | Outcome model correct |
| Doubly Robust | `dr` | At least one model correct |

---

## Requirements

- Stata 14.1 or later
- All files are pure ASCII (no Unicode characters)
- No user-written dependencies required

---

## Companion R Package

This Stata package is a port of the
[CausalSpline R package](https://github.com/causalfragility-lab/CausalSpline).
Both packages implement identical methodology and produce comparable results.

---

## References

- Hirano, K. & Imbens, G.W. (2004). The propensity score with continuous
  treatments. In *Applied Bayesian Modeling and Causal Inference* (pp. 73-84).
  Wiley.
- Imbens, G.W. (2000). The role of the propensity score in estimating
  dose-response functions. *Biometrika*, 87(3), 706-710.
- Robins, J.M., Hernan, M.A. & Brumback, B. (2000). Marginal structural
  models and causal inference in epidemiology. *Epidemiology*, 11(5), 550-560.

---

## Citation

```stata
* In Stata:
which causalspline
```

Or cite as:

> Hait, S. (2026). *causalspline: Nonlinear Causal Dose-Response Estimation
> via Splines*. Stata package.
> https://github.com/causalfragility-lab/CausalSpline-Stata

---

## License

GPL (>= 3)

Maintainer: Subir Hait <haitsubi@msu.edu>
Michigan State University
