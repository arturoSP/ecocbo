# ecocbo 0.13.1

# ecocbo 0.13.0

This is a maintenance and feature-update release, preparing the package for resubmission to CRAN.

## NEW FEATURES

- Updated the function `sim_cbo()` by changing the structure to an empirical optimization based on precision and then cost. The function evaluates for the sampling designs that best approximate to $(1-\alpha)$ and then finds the cost for each one of the selected sampling efforts. Lastly, the function marks the best options. The function is documented and includes examples.
- Added the function `underwood_cbo()` to keep the possibility of using Underwood's optimization that was used in previous versions of the package.  The function is documented and includes examples.

## BUG FIXES

- Resolved `NOTE`s related to "no visible binding for global variable" and "no visible global function definition". This was addressed by explicitly importing functions with `@importFrom` and declaring global variables with `utils::globalVariables()`, making the package more robust and compliant with CRAN policies.
- Corrected potential miscalculations in internal simulation functions by ensuring all variable scopes are handled correctly.

## IMPROVEMENTS

- The main vignette has been updated to include a detailed example and workflow for the new `permanova_twoway()` function.
- Internal code has been refactored to reduce dependencies and improve clarity, for instance, by favoring base R functions where appropriate.

# ecocbo 0.12.0

- 'ecocbo' now can work with either single-factor and nested-symmetric experiments. 

- A new function, `prep_data()`, is added. This allows the user to select which model to use and then prepares the data for using with the rest of the functions. 

- `plot_power()` was updated. Readability of the power curve is improved by differentiating the optimal, or user selected, experimental design. 

# ecocbo 0.11.0

# ecocbo 0.10.2

This is a resubmission, in this version we have:

- corrected the presentation for function names in the description texts by removing single quotes and adding () after each name.

- added ISBN to the reference for Underwood (1997) in the DESCRIPTION, as well as in the manual files and vignette.

- changed \dontrun{} to \donttest{} for examples involving 'sim_beta()' as it takes more than 5 seconds to run them as they are. Changing the parameter values to make the function run faster would not be instructional to the final user, as it would not demonstrate the function's actual runtime and functionality. The example code will not be run automatically, but it can still be run manually if desired.

# ecocbo 0.10.1
- We added '\\donttest{}' to the 'sim_beta()' example in the DESCRIPTION file because it takes more than 5 seconds to run. Changing the parameter values to make the function run faster would not be instructional to the final user, as it would not demonstrate the function's actual runtime and functionality. The example code will not be run automatically, but it can still be run manually if desired.

# ecocbo 0.9.1

* Added a `NEWS.md` file to track changes to the package.
