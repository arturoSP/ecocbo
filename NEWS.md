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
