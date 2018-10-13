# README

This file is here to document the replication package for the choice
model estimated in Johanna Mollerstrom, Bjørn-Atle Reme, Erik
Ø. Sørensen, Luck, choice and responsibility — An experimental study
of fairness views, Journal of Public Economics, Volume 131, 2015,
Pages 33-40, https://doi.org/10.1016/j.jpubeco.2015.08.010.

data/DataErik.dta and data/situations.dta: Data

preprocess.do (and .log): Stata dataprocessing script which
   takes the files created above and creates data/DataErik_long.dta and
   data/DataErik_players.dta.
   
luck-model.Rmd: R file with definition of data structures, log
   likelihood functions, standard error calculations and such. No
   actual estimation is done by this file.  luck-model.{R,html} are
   products of this file. The html-file has explanations of how the
   likelihood function is calculated as well.
   An explicit design goal for the code is to be value explicit typing
   and design of data structures and straightforward calculation above
   computational speed. But an understanding of the R S4 object system
   might be necessary to fully understand the code here. 

luck-main-table.Rmd: R file with estimation code for Table 6 in the paper.
   luck-main-table.html is the file with calculated output.

insurance_no_insurance: R file with estimation code for Table 8 in the paper.
   insurance_no_insurance.html is the file with calculated output.

pooling_with_extension.Rmd: R file with estimation code for Table 9 in the
   paper. pooling_with_extension.html is the file with calculated output.

predictions.do: Stata file that do predictions to generate Table 7. Generates
   predictions.log


The html files are readable with standard browsers. The three files
with estimation code relies on the luck-model.R file for its
functions. It might be wise to study one of the generated .html files with
calculated results from the .Rmd files to see what installed packages
are needed. The library file luck-model.Rmd also relies on setting R to
work with 8 processor cores. For a different computer, this should
be changed to reflect the computer it is running on.








 
