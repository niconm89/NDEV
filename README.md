# NDEV v01: Nucleotide Diversity Estimation and Visualization
This program calculate the Genetic Diversity of sequences present in a nucleotide aligment.
The file must be in FASTA format. If you need to convert any of your alignments to FASTA format, these [python scripts](https://github.com/niconm89/scripts_bioinformatics) may help you.

------
## Requieriments
NDEV is written in python3 and requiere the following python modules/libraries:
* os
* argparse
* math
* random
* pandas
* numpy
* matplotlib

------
## Install
To install the required modules, you can use pip3 and install all in one command:

`pip3 install os argparse math random pandas numpy matplotlib`

However, I recommend you to create a [conda environment](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html).

```bash
#create the environment named NDEV
conda create -n NDEV
#list all environments I have
conda env list
#activate the environment we just created
conda activate NDEV
#install modules
conda install -c jmcmurray os
conda install -c conda-forge argparse matplotlib
conda install -c anaconda numpy pandas
#the above steps can be done in one simple command, but I like to split steps as more as possible
```

## Usage & Examples
* Program help

`python NDEV.py -h`

* Program version

`python NDEV.py -v`

* Program 

`python NDEV.py -v`



## Outputs


## Contact

Nicolás Nahuel Moreyra, PhD
* :e-mail: niconm89@gmail.com
