[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/fowler-lab/cryptic-ecoffs/HEAD?labpath=reproduce-tables-and-figures.ipynb)

# Reproducing the Tables and Figures in the CRyPTIC ECOFF/ECV paper

This repository contains a [jupyter notebook](https://jupyter-notebook-beginner-guide.readthedocs.io/en/latest/index.html) that reads in the raw data from the [CRyPTIC](http://www.crypticproject.org/) project and reproduces the vast majority of the Tables and Figures in the below peer-reviewed paper.

> Epidemiological cutoff values for a 96-well broth microdilution plate for high-throughput research antibiotic susceptibility testing of M. tuberculosis
> 
> The CRyPTIC Consortium
> 
> European Respitatory Journal [doi:10.1183/13993003.00239-2022](https://doi.org/10.1183/13993003.00239-2022)

## Running the notebook

### Using MyBinder

The easiest way by far is to click the `launch binder` button at the top of this document. This will open a MyBinder session which will install all the required software on a remote container (this can take a few minutes) and furnish you the notebook in browser. You can then run each cell in turn, or add new cells and write your own Python/Pandas code to interrogate the data and analysis. The graphs are, however, not written to the `graphs/` folder when using MyBinder. For that you need to install it locally.

### Installing locally

This requires a little more knowledge, but if you have any Python experience should be familiar. In brief, in a terminal

```
$ git clone git@github.com:fowler-lab/cryptic-ecoffs.git
$ cd cryptic-ecoffs/
$ pip install .
$ jupyter-lab reproduce-tables-and-figures.ipynb
```
The last command should open a window in your default browser and load the notebook. Again, like with MyBinder you can then step through and run all the cells. All the graphs are written as PDF files to `graphs/` and are identical to the ones in the manuscript.

