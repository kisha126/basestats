
# basestats

<figure>
<img src="https://img.shields.io/badge/License-MIT-yellow.svg"
alt="alt text" />
<figcaption aria-hidden="true">alt text</figcaption>
</figure>

A collection of hypothesis test functions implemented in R, designed to
be used directly as modules rather than a traditional package. This
project focuses on providing clear implementations and demonstrating
modular code organization with `box`.

## Overview

`basestats` provides recreations and potential extensions of common
statistical hypothesis tests. Instead of bundling these into an R
package, it leverages the `box` package for a modular approach. This
allows users to directly use specific functions without installing a
full package, promoting code transparency and easier modification.

The functions are currently organized, e.g. into parametric test
functions (`par_funcs`) and non-parametric test functions (`npar_funcs`)
categories.

## Dependencies

- **`box`**: You can’t use the function without `box` package installed
  in your system. This is the main requirement to use the modular
  functions in `basestats`.

  This is how you install the package:

  - Installing this package from CRAN:

    ``` r
    install.packages('box')
    ```

  - From R-universe:

    ``` r
    install.packages('box', repos = 'https://klmr.r-universe.dev')
    ```

  - Or a development version:

    ``` r
    devtools::install_github("klmr/box")
    ```

- **`tabhelpers`**: My newly created package, used for displaying
  formatted statistical outputs from a data frame or a matrix. Take note
  that this is still currently developing, so install this from GitHub:

  ``` r
  devtools::install_github("kisha126/tabhelpers")
  ```

- **Core `tidyverse` packages (e.g., `dplyr`, `tidyr`)**: While the full
  `tidyverse` meta-package is convenient
  (`install.packages("tidyverse")`), the examples and potentially the
  `tabhelpers` package might rely specifically on packages like `dplyr`
  for data manipulation (e.g., using the `|>` pipe). Ensure relevant
  `tidyverse` components are installed.

- R version must be \>4.2.0

# Getting Started

1.  Get the code: Clone or download this repository to your local
    machine.

    ``` bash
    git clone https://github.com/kisha126/basestats.git
    cd basestats # Optional
    ```

2.  **Use in R:** Navigate your R session’s working directory to outside
    the `basestats` folder (or adjust the path in `box::use`
    accordingly). You can then load modules using `box::use()`. I wanna
    recommend you to not use `setwd` to set up a project with this
    repository.

## Project Structure

Here is the project structure:

    .
    ├── basestats
    │   ├── __init__.R
    │   ├── npar_funcs
    │   │   ├── __init__.R
    │   │   └── wilcoxtest.R
    │   ├── par_funcs
    │   │   ├── __init__.R
    │   │   └── ttest.R
    │   └── utils.R
    ├── README.md
    └── README.rmd

*(Note: You’ll need appropriate **init**.R files to export functions and
sub-modules correctly for box::use)*

## Example Usage

Since `basestats` is not a package, but rather a structural directory,
load the directory as a module with `bs` as an alias of it, and don’t
forget the `./` as the prefix. Ensure your R working directory is set
correctly relative to the `basestats` directory. If `basestats` is a
subdirectory of your current working directory, use `bs = ./basestats`.

``` r
box::use(bs = ./basestats)
data("iris")

iris |> 
    bs$par$t_test(Sepal.Length, Sepal.Width, Petal.Length)
#> 
#> ─────────────────────────────── One-Sample t-test ──────────────────────────────
#> 
#> 
#> 1. Sepal.Length (μ = 0) 
#> 
#> 
#>                            ─────────────────────────
#>                                   Test Result       
#>                            ─────────────────────────
#>                              t_statistic   p-value  
#>                            ─────────────────────────
#>                                86.425      <0.001   
#>                            ─────────────────────────
#> 
#> 
#> 
#>                                Other Information                
#>                ─────────────────────────────────────────────────
#>                  Alternative Hypothesis   :   different from 0
#>                  True Population          :                  0
#>                  Estimate                 :              5.843
#>                ─────────────────────────────────────────────────
#> 2. Sepal.Width (μ = 0) 
#> 
#> 
#>                            ─────────────────────────
#>                                   Test Result       
#>                            ─────────────────────────
#>                              t_statistic   p-value  
#>                            ─────────────────────────
#>                                85.908      <0.001   
#>                            ─────────────────────────
#> 
#> 
#> 
#>                                Other Information                
#>                ─────────────────────────────────────────────────
#>                  Alternative Hypothesis   :   different from 0
#>                  True Population          :                  0
#>                  Estimate                 :              3.057
#>                ─────────────────────────────────────────────────
#> 3. Petal.Length (μ = 0) 
#> 
#> 
#>                            ─────────────────────────
#>                                   Test Result       
#>                            ─────────────────────────
#>                              t_statistic   p-value  
#>                            ─────────────────────────
#>                                26.073      <0.001   
#>                            ─────────────────────────
#> 
#> 
#> 
#>                                Other Information                
#>                ─────────────────────────────────────────────────
#>                  Alternative Hypothesis   :   different from 0
#>                  True Population          :                  0
#>                  Estimate                 :              3.758
#>                ─────────────────────────────────────────────────

box::unload(bs)
```

## Contributing

\[Optional: Add guidelines here if you welcome contributions, e.g.,
“Pull requests are welcome. For major changes, please open an issue
first to discuss what you would like to change.”\]

## License

This project is licensed under the MIT License. See the
[LICENSE](https://opensource.org/license/MIT) file or link for details.
