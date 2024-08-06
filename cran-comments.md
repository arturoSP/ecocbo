## Update to ecocbo 0.12.0

- 'ecocbo' now can work with either single-factor and nested-symmetric experiments. 

- A new function, `prep_data()`, is added. This allows the user to select which model to use and then prepares the data for using with the rest of the functions. 

- `plot_power()` was updated. Readability of the power curve is improved by differentiating the optimal, or user selected, experimental design. 

### R CMD check results 

Duration: 4m 3.5s

0 errors ✔ | 0 warnings ✔ | 0 notes ✔

## Resubmission

This is a resubmission, in this version we have:

- corrected the presentation for function names in the description texts by removing single quotes and adding () after each name.

- added ISBN to the reference for Underwood (1997) in the DESCRIPTION, as well as in the manual files and vignette.

- changed \dontrun{} to \donttest{} for examples involving 'sim_beta()' as it takes more than 5 seconds to run them as they are. Changing the parameter values to make the function run faster would not be instructional to the final user, as it would not demonstrate the function's actual runtime and functionality. The example code will not be run automatically, but it can still be run manually if desired.

### R CMD check results

0 errors ✔ | 0 warnings ✔ | 1 notes ✖

- This is a new release

  > checking CRAN incoming feasibility ... [164s] NOTE <br/>
  > - Maintainer: 'Arturo Sanchez-Porras <sp.arturo@gmail.com>'
  
  > Possibly mis-spelled words in DESCRIPTION: <br/>
  > -  Permutational (14:70)
  
  - New submission
  - The term "PERMANOVA" was coined by [Anderson (2017)](https://doi.org/10.1002/9781118445112.stat07841) to refer to Permutational Multivariate Analysis of Variance. This is part of the statistical theory that underlies the calculations performed in the function sim_beta(). As you can see, the word "permutational" is used correctly in the title of the paper.


## Test environments

- Local:
  - Linux (PopOS), R 4.1.2(x86_64-pc-linux-gnu)
  - Windows 10, R 4.3.0(x86_64-w64-mingw32/x64 (64-bit))
- win-builder:
  - Windows Server 2022 x64 (build 20348)
- R-hub builder (https://builder.r-hub.io):
  - Fedora Linux, R-devel, clang, gfortran
  - Windows Server 2022, R-devel, 64 bit
  - Ubuntu Linux 20.04.1 LTS, R-release, GCC
- macOS builder (https://mac-r-project.org/macbuilder/submit.html):
  - aarch64-apple-darwin20 (64-bit)
