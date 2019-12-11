# MonteCarlo-DiffusionGrowth
 Molecular Thin Film Growth by Diffusion using Monte Carlo

## Model

An empty slate is populated with dot molecules that can move up / down or (model b) sideways randomly for Tdifussion steps. More diffusion will create larger clusters as expected. 

## Simulation Images

Here are two images generated for a 200x200 cluster size with diffusion 10 and 1000, respectively:

![200x200 with high diffusion](images/MolDiff-Side-200x200-D1000.png)

![200x200 with low diffusion](images/MolDiff-Side-200x200-D10.png)

## Measuring the height-height correlation length

We have discussed the merits of [height-height correlation function to determine grain size](https://iopscience.iop.org/article/10.1088/1742-6596/417/1/012069/pdf), so we are applying this method to these images.

![Correlation Length Fit](images/MolDiff-Side-400x400-D2000-CorrLen.png)

In order to limit the correlation length fit standard deviation to about 3%, at least 5 million iterations for the HHCF are needed for a 400x400 grid, see `test.hhcf.R`. The correlation length dependence with iterations for 3 different diffusion lengths are shown here:

![Correlation length versus number of iterations](images/test.hhcf-400x400-HHCF1000.png)

## Grain Size Dependence on Diffusion Steps

