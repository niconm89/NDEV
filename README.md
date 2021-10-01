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
* matplotlib (in the next version ggplot2 will be used to generate plots)

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

------
## Usage & Examples
* Program help

`python NDEV.py -h`

* Program version

`python NDEV.py -v`

* Loading an alignment to estimate PI values and generate a line plot.

`python GenDivEst.py --mode full --infile data_example/ingroups.fasta --prefix test1_full --wsize 500 --woverlap 100 --title 'This tis the title' --title_fontsize 32 --labels cluster_buzzatii --xlabel 'Windows in alignment' --xlimit 0 39 --ylabel 'Nucleotide Diversity (PI)' --ylimit 0 0.02 --axis_fontsize 28 28 --ticks_size 22 22 --legend 'upper right' --legend_fontsize 26 --colors red --linestyles solid --linewidth 2.5 --dpi 300`

* Loading an alignment to estimate PI values and generate a line plot and a boxplot of PI values between all pairwise comparisons.

`python GenDivEst.py --mode full --infile data_example/ingroups.fasta --prefix test_ingroups --wsize 500 --woverlap 100 --title 'This tis the title' --title_fontsize 32 --labels cluster_buzzatii --xlabel 'Windows in alignment' --xlimit 0 39 --ylabel 'Nucleotide Diversity (PI)' --ylimit 0 0.02 --axis_fontsize 28 28 --ticks_size 18 18 --legend 'upper right' --legend_fontsize 26 --colors red --linestyles solid --linewidth 2.5 --dpi 300 --pairwise_plot --rotation 45 0`

* Loading an alignment to estimate PI values without generating any plot

`python GenDivEst.py --mode full --infile data_example/ingroups.fasta --prefix test2_full-noplots --wsize 500 --woverlap 100 --noplot`

* Using plot mode with two files with PI values

`python GenDivEst.py --mode plot --pifile data_example/ingroups.pi_valores.txt data_example/outgroups.pi_valores.txt --prefix test3_plot --title 'This tis the title' --title_fontsize 32 --labels ingroups outgroups --xlabel 'Windows in alignment' --xlimit 0 39 --ylabel 'Nucleotide Diversity (PI)' --ylimit 0 0.02 --axis_fontsize 28 28 --ticks_size 22 22 --legend 'upper right' --legend_fontsize 26 --colors red green --linestyles solid dashed --linewidth 4.0 2.5 --dpi 300`

------
## Outputs

Acá van los ejemplos de salidas

------
## Contact

Nicolás Nahuel Moreyra, PhD
* :e-mail: niconm89@gmail.com
