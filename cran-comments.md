## Resubmission

This is a resubmission, in this version I have:
- added \donttest{} to the 'sim_beta' example because it takes more than 5 seconds to run. Changing the parameter values to make the function run faster would not be instructional to the final user, as it would not demonstrate the function's actual runtime. The example code will not be run automatically, but it can still be run manually if desired.

## R CMD check results

0 errors ✔ | 0 warnings ✔ | 1 notes ✖

- This is a new release

  > checking CRAN incoming feasibility ... [164s] NOTE <br/>
  > - Maintainer: 'Arturo Sanchez-Porras <sp.arturo@gmail.com>'
  
  > Possibly mis-spelled words in DESCRIPTION: <br/>
  > -  Permutational (14:49)
  
  - New submission
  - Part of this work follows the ideas developed by Anderson (2017) in her paper "Permutational Multivariate Analysis of Variance (PERMANOVA)", therefore "Permutational" is not a mis-spelled word.

## Test environments

- Local:
  - Linux (PopOS), R 4.1.2(x86_64-pc-linux-gnu)
  - Windows 10, R 4.3.0(x86_64-w64-mingw32/x64 (64-bit))
- win-builder:
  - Windows Server 2022 x64 (build 20348)
- R-hub builder (https://builder.r-hub.io):
  - Fedora Linux, R-devel, clang, gfortran
    - Gave the same notes that were reported locally by the Linux session.
- macOS builder (https://mac-r-project.org/macbuilder/submit.html):
  - aarch64-apple-darwin20 (64-bit)
    - Status: OK


