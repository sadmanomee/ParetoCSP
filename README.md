# ParetoCSP
Github repository for our manuscript - **"Crystal structure prediction using neural network potential and age-fitness Pareto genetic algorithm"**

Authors: Sadman Sadeed Omee, Lai Wei, and Jianjun Hu.


[Machine Learning and Evolution Laboratory,](http://mleg.cse.sc.edu)<br />
Department of Computer Science and Engineering, <br />
University of South Carolina,<br/>
SC, 29201, United States.

**ParetoCSP** is a crystal structure prediction (CSP) algorithm that given a chemical composition of a crystal, utilizes a combination of a multi-objective genetic algorithm (GA) and a neural network inter-atomic potential (IAP) model to search for the energetically optimal structure in the structure space. Specifically, we have enhanced the NSGA-III to use the genotypic age of the population as a separate optimization criterion so that it is treated as an independent dimension in the multi-objective Pareto front, and the M3GNet universal IAP to guide the age fitness Pareto optimized GA search function towards identifying the structure with the lowest final energy per atom, which is then relaxed using M3GNet IAP.

# Table of Contents
* [Necessary Installations](#installation)
* [How to run](#usage)
* [Contributors](#contributors)
* [Acknowledgement](#acknowledgement)

<a name="installation"></a>
## Necessary Installations
Please install the following packages if not already installed. We show how to install them using **pip** only, but you can also use **conda** for the installation purpose. Also you can a virtual environment using conda or pip for this purpose (recommended).

Use the following commands to install the necessary packages:
```bash
git clone https://github.com/sadmanomee/ParetoCSP.git
cd ParetoCSP
pip install -r requirements.txt
```

<a name="usage"></a>
## How to run
The default configuration are mentioned in the ```config.py``` file. So simply running the following command will run the algorithm will all default configuration (including compostion and algorithm name):
```bash
python predict_structure.py
```
To specify the crystal composition and other arguments, use the following command (recommended):
```bash
python predict_structure.py --comp=Sr1Ti1O3 --alg=afpo --pop=100 --max_step=5000
```
This will run ParetoCSP for the composition SrTiO<sub>3</sub> (you need to specify the atom number after each compound symbol) for the PartoCSP algorithm with a population size of 100 for the NSGA-III and the algorithm will run for a total 5000 generations. The output crystal structure (CIF file) path will be printed at the end of the algorithm run. Other arguments that can be passed with the command are mentioned in ```predict_structure.py``` file.

<a name="contributors"></a>
## Contributors

1. Sadman Sadeed Omee (<https://www.sadmanomee.com/>)
2. Dr. Jianjun Hu (<http://www.cse.sc.edu/~jianjunh/>)

## Acknowledgement

Our code is based on the [GN-OA](http://www.comates.group/links?software=gn_oa) algorithm's repository, which has a well-developed CSP code. We also used the [Pymoo](https://github.com/anyoptimization/pymoo) code's repository for implementing the AFPO-enhanced NSGA-III.
