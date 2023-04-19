# Welcome to variogam

The variogam folder contains all code necessary to recreate the simulation study found in Duskey and Sullivan 2023.  We will update with the full citation when the article has been accepted for publication.  Briefly, the code can simulate data with a mean trend and autocorrelated residuals, and then run a Bayesian regression kriging model in Stan through CmdStanR.  The regression is a penalized spline, commonly associated with GAMs, and the variogram is a simple exponential model.

## Requirements

The Bayesian models contained in variogam run in Stan via the R package CmdStanR.  See the following website for instructions on how to install CmdStan, CmdStanR, and Stan:

https://mc-stan.org/cmdstanr/

Please follow all instructions therein to ensure that CmdStanR has been properly installed and configured.

## Usage

Follow these steps prior to first time usage to ensure the code runs properly on your machine:

1. Download the "variogam-main" zipped folder from Github in its entirety
2. Unzip and place "variogam-main" in your preferred directory
3. Navigate to the Code folder and open the file "depend.R" in R or RStudio
4. Confirm that all listed packages are installed on your machine; uncomment and run lines corresponding to packages you have not yet installed, then re-comment
5. Change the "mypath" variable to the file path which contains the variogam folder
6. Save and close "depend.R"

Upon all subsequent uses, we recommend running the .R files in the following order:

1. depend.R -- to set mypath; alternatively, change setwd() within other R files
2. allstan.R -- to simulate data and run all models on autocorrelated data
3. nosim.R -- to simulate data and run all models on non-autocorrelated data
4. gamkrige.R -- to run frequentist GAM + variogram on all simulated data
5. looic.R -- to calculate log likelihood and compute LOOIC for all models
6. proplm.R -- to run linear models on median variogram parameter estimates as a function of maximum median estimates
7. allplot.R -- to recreate all simulation figures in Duskey and Sullivan 2021

All files not appearing here are sourced in one or more of the files listed above, and typically contain functions necessary to run the code in each script.

NOTE: There are several nearly empty output folders contained in variogam, containing only a textfile called ".keep".  These are meant to contain and organize your own output.  Feel free to delete the ".keep" files after downloading.  If you would like our original model output, it is available upon request.

## Contact us

If you have questions or concerns regarding this code, or would like help in re-formatting it for your own use, please do not hesitate to contant the corresponding author at:

elizabeth.duskey@slu.se

## License

Copyright 2022 Elizabeth Duskey

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
