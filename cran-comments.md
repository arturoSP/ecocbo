## Resubmission

This is a resubmission of `ecocbo`, updating it to version 0.13.0.

## Summary of Changes

This update includes one major change, several bug fixes, and internal improvements. The key changes are:

1.  **Updated Feature**: `sim_cbo()` is redesigned as a new function that evaluates for the sampling designs that best approximate to $(1-\alpha)$ and then finds the cost for each one of the selected sampling efforts. Lastly, the function marks the best options. `underwood_cbo()`, has been added to keep the possibility of performing Underwood's optimization.
2.  **Bug Fixes**: Addressed a `NOTE` present in previous checks concerning undefined global variables.
3.  **Documentation**: Updated the package vignette to demonstrate the new functionality.

The `NEWS.md` file provides a more detailed list of changes.

## Addressing a Previous `NOTE`

The previous version of the package generated a `NOTE` during `R CMD check` regarding "no visible binding for global variables". This issue has been fully resolved in this version by:

* Explicitly importing functions from other packages (e.g., `stats::rchisq`) using the `NAMESPACE` file via `roxygen2` tags (`@importFrom`).
* Declaring variables used in non-standard evaluation with `utils::globalVariables()`.

## Test Environments

This package has been checked locally using `devtools::check(args = "--as-cran")` on:

* Platform: x86_64-pc-linux-gnu
* R version: 4.5.1 (2025-06-13)

── R CMD check results ─────────────────────────────────────────── ecocbo 0.13.0 ────
Duration: 3m 32.6s

0 errors ✔ | 0 warnings ✔ | 0 notes ✔

It has also been successfully checked on:

* `win-builder` (R-devel)
* `github actions` (covering various distributions, including Linux and macOS)

No new issues were reported in these checks. There are no reverse dependencies for this package.

Thank you for your time and for maintaining CRAN.
