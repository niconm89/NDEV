#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Nucleotide Diversity Estimation and Visualization (NDEV)

@author: Nicolás Nahuel Moreyra (niconm89@gmail.com)
"""
#%% Imports
import os
import argparse
from math import factorial
from random import choices as rchoices
import pandas as pd
import numpy as np
import matplotlib.colors as mcolors
#import matplotlib.pyplot as plt
from plotnine import *
#from math import comb #only for python3.8+

#%% Functions
def load_PIvalues_from_file(PI_file):
    'Receive a path to numpy file and return an numpy array'
    return np.loadtxt(PI_file).copy()
#
def is_fasta(first_line_file): #checked
    'Receive the first line of a file to check if starst with the > symbol.'
    if first_line_file[0] != '>':
        print("Input file does not have fasta format.\nLeaving...")
        exit(1)
    return True
#
def read_aln_fasta(input_file):
    'Receive a path to FASTA file and return an dictionary of sequences'
    dict_seqs = {}
    with open(input_file,'rt') as infile:
        FASTA = infile.readlines()
        is_fasta(FASTA[0]) #check if the first line of the file starts with '>', otherwise it will exit
        for line in FASTA:
            line = line.rstrip() #let's discard the newline at the end (if any)
            if line[0] == '>': #identify the id of each sequence
                name = line[1:] #avoid the '>'
                dict_seqs[name] = ''
            else:
                dict_seqs[name] = dict_seqs[name] + line
    return dict_seqs
#
def read_aln_to_array(input_file):
    'Receive a path to FASTA file and return two numpy arrays with ids and sequences'
    alignment_dict = read_aln_fasta(input_file)
    ids = np.array([ idseq for idseq in alignment_dict.keys() ])
    sequences = np.array([ list(seq) for seq in alignment_dict.values() ])
    return ids.copy(),sequences.copy() # para enviar las secuencias sin los IDs
#
def calc_pi(list_of_mutations, list_of_valid_sites):
    'Receive two numpy arrays with the mutations and the valid sites (no gaps or N) per position of a pairwise comparison, and return the proportion of different sites.'
    nvalid = list_of_valid_sites.sum()
    if nvalid == 0: #no sites to compare in this window
        return None
    return list_of_mutations.sum() / nvalid
#
def compare_two_sequence(seq1,seq2):
    'Receive two numpy arrays with one sequence in each one and return two numpy arrays: one with 1 in positions where the sequences are different, and other with 1 in positions where both sequences have a valid base (A, C, G or T).'
    seqlen = len(seq1)
    list_of_differences = np.zeros(seqlen)
    list_of_valid_sites = np.zeros(seqlen)
    for i in range(0,seqlen):
        base1 = seq1[i].upper()
        base2 = seq2[i].upper()
        if base1 in 'ACGT' and base2 in 'ACGT': #gaps, N and ambiguous bases are excluded
            list_of_valid_sites[i] = 1
            if base1 != base2:
                list_of_differences[i] = 1
    return list_of_differences, list_of_valid_sites
#
def number_of_combinations(n, m):
    'Calculate the number of possible combinations or n elements grouping by m elements.'
    return factorial(n) // (factorial(m) * factorial(n - m))
#
def analyze_window(window_aln):
    'Receive a window of an alignment and return the pi value (mean proportion of different sites across all pairwise comparisons).'
    number_seqs = len(window_aln)
    pi_values_pairs = []
    for i in range(number_seqs):
        for j in range(i+1,number_seqs):
            mutations, valid_sites = compare_two_sequence(window_aln[i],window_aln[j])
            pi_pair = calc_pi(mutations, valid_sites)
            if pi_pair is not None: #pairs without valid sites are not considered
                pi_values_pairs.append(pi_pair)
    if len(pi_values_pairs) == 0: #no valid sites in the whole window
        return np.nan
    pi = sum(pi_values_pairs) / len(pi_values_pairs)
    return pi
#
def number_of_windows(length_aln, windows_size, windows_overlap):
    'Receive the length of an aligment, the window size and the overlap length between windows, and return the number of windows for the alignment to generate.'
    return sum([ 1 for i in range(0, length_aln - windows_overlap, windows_size - windows_overlap) ])
#
def generate_windows(alignment, wsize, woverlap):
    'Receive an alignment and windows size and overlap, and return a list of numpy arrays that correspond to windows of the aligments.'
    length_aln = len(alignment[0])
    return [ alignment[:,i:i+wsize] for i in range(0, length_aln - woverlap, wsize - woverlap) ]
#
def calculate_pi_per_window(alignment, wsize, woverlap):
    'Receive an alignment and windows size and overlap, and return a numpy array with pi values per windows generated.'
    lenseq = len(alignment[0]) #sequence length
    numwin = number_of_windows(lenseq, wsize, woverlap) #number of windows
    pi_values = np.zeros(numwin) #array of pi values (per window)
    for n,window_aln in enumerate(generate_windows(alignment,wsize,woverlap)): #for each window, we analyze it (calc pi)
        pi_values[n] = analyze_window(window_aln)
    return pi_values
#
def check_file_exists(full_path_to_file):
    'Check if the file exist. Otherwise, the program will exit.'
    if not os.path.exists(full_path_to_file):
        print("Error: File does not exists! " + full_path_to_file)
        exit(1)
    return True
#
def generate_pairwise_alingments(ids_seqs, alignment, windows_size, windows_overlap):
    'Receive an alignment in a numpy array (with windows size and overlap) and return a numpay array (matrix) with PI values per each pairwise comparison in the alignment.'
    nseqs = len(alignment)
    num_comparisons = number_of_combinations(nseqs, 2)
    numwin = number_of_windows(len(alignment[0]), windows_size, windows_overlap)
    pairwise_PIvalues = np.zeros(num_comparisons*numwin)

    ids = []
    count = 0
    for i in range(nseqs):
        for j in range(i+1,nseqs):
            pairwise_aln = np.array([alignment[i],alignment[j]])
            pairwise_PIvalues[count:count+numwin] = calculate_pi_per_window(pairwise_aln, windows_size, windows_overlap)
            ids.append(ids_seqs[i]+'_'+ids_seqs[j])
            count += numwin

    labels = []
    windows = []
    for ID in ids:
        for w in range(0,numwin):
            labels.append(ID)
            windows.append(w+1)
    
    return windows, labels, pairwise_PIvalues
#
def check_prefix(prefix):
    'Receive a prefix and extract (if present) the full or relative path and check if exists.'
    out_prefix = prefix.split('/')[-1]
    full_path = '/'.join(prefix.split('/')[:-1])
    if not full_path: #if no path (/) was passed
        full_path = os.getcwd() #current working directory
    if not os.path.isdir(full_path): #if path does not exists
        raise RuntimeError("Output folder does not exists. Check the PATH!. Aborting...")
    return full_path,out_prefix
#
def generate_default_colors(nfiles):
    'Receive a number and return this amount of different colors.'
    CSS_colors = mcolors.CSS4_COLORS.copy()
    if 'white' in CSS_colors:
        del CSS_colors['white']
    return rchoices(list(CSS_colors.values()),k=nfiles)
#
def plot_pi_pairwise(data, full_path, out_prefix, plot_params):
    'Receive a dataframe with pairwise PI values and a dictionary with plot parameters, and plot (save) a figure with one box per comparison.'
    if plot_params['violin']:
        distribution = geom_violin(alpha=plot_params['alpha'])
    else:
        distribution = geom_boxplot(alpha=plot_params['alpha'])
    boxplot = ggplot(data, aes(x='Comparison', y='PI_value', fill='Comparison')) + \
    distribution + \
    ggtitle(plot_params['title']) + \
    theme(panel_grid_major=element_blank(),panel_grid_minor=element_blank()) + \
    theme(text=element_text(family="Arial", size=plot_params['labels_fontsize'])) + \
    theme(plot_title=element_text(size=plot_params['title_fontsize'], family="Arial", face="bold")) + \
    theme(axis_text_x = element_text(color="black", angle=plot_params['rotation'][0], hjust=.5, vjust=1, size=plot_params['axis_fontsize'][0])) + \
    theme(axis_text_y = element_text(color="black", angle=plot_params['rotation'][1], hjust=1, vjust=0.5, size=plot_params['axis_fontsize'][1])) + \
    theme(axis_line_x = element_line(size=plot_params['axis_linesize'], colour=plot_params['axis_linecolor'])) + \
    theme(axis_line_y = element_line(size=plot_params['axis_linesize'], colour=plot_params['axis_linecolor'])) + \
    theme(legend_position='none') + \
    xlab(plot_params['xlabel']) + ylab(plot_params['ylabel'])
    if plot_params['jittler']: #adding points over the boxes
        boxplot = boxplot + geom_jitter(width=0.2, size=1)
    if plot_params['boxes_colors']: #colors passed by the user
        boxplot = boxplot + scale_fill_manual(values=plot_params['boxes_colors'])
    if plot_params['ylimit']:
        boxplot = boxplot + scale_y_continuous(limits=plot_params['ylimit'])
    #saving the figure
    out_figure = os.path.join(full_path, out_prefix + '.pairwise_boxplot.png') #creating full output path
    boxplot.save(filename=out_figure, dpi=plot_params['dpi'], height=plot_params['figsize'][0], width=plot_params['figsize'][1])
#
def plot_pi_values(data, full_path, out_prefix, plot_params):
    'Receive a list of path to files with pi values and a dictionary with plot parameters, and plot (save) a figure with all curves (one per file).'
    #curves
    lineplot = ggplot(data, aes(x='Num_Win', y='PI_value', group='Comparison', color='Comparison', fill='Comparison')) + \
    geom_line(size = plot_params['linewidth'], linetype=plot_params['linestyle']) + \
    ggtitle(plot_params['title']) + \
    theme(panel_grid_major=element_blank(),panel_grid_minor=element_blank()) + \
    theme(text=element_text(family="Arial", size=plot_params['labels_fontsize'])) + \
    theme(plot_title=element_text(size=plot_params['title_fontsize'], family="Arial", face="bold")) + \
    theme(axis_text_x = element_text(color="black", angle=plot_params['rotation'][0], hjust=.5, vjust=1, size=plot_params['axis_fontsize'][1])) + \
    theme(axis_text_y = element_text(color="black", angle=plot_params['rotation'][1], hjust=1, vjust=0.5, size=plot_params['axis_fontsize'][1])) + \
    theme(axis_line_x = element_line(size=plot_params['axis_linesize'], colour=plot_params['axis_linecolor'])) + \
    theme(axis_line_y = element_line(size=plot_params['axis_linesize'], colour=plot_params['axis_linecolor'])) + \
    theme(legend_direction="vertical",legend_box_spacing=0.2) + \
    xlab(plot_params['xlabel']) + ylab(plot_params['ylabel'])
    #saving the figure
    out_figure = os.path.join(full_path, out_prefix + '.pivalues.png') #creating full output path
    lineplot.save(filename=out_figure, dpi=plot_params['dpi'], height=plot_params['figsize'][0], width=plot_params['figsize'][1])
#
def manage_plots(plot_mode, pi_files, full_path, out_prefix, plot_params):
    'Receive a list of files and a dictionary with plot parameters and admin the reformat variables to make the plot.'

    if not plot_params['legend']:
        plot_params['legend'] = ''

    columns = ['Num_Win','PI_value','Comparison']
    if plot_mode == 'boxplot':
        data = pd.read_csv(pi_files, sep='\t')
        for col in data.columns:
            if col not in columns:
                raise RuntimeError('Inconsistent column names. Aborting...')
        if plot_params['pairs']: #keep only the comparisons selected by the user
            data = data[ data['Comparison'].isin(plot_params['pairs']) ].copy()
            if len(data) == 0:
                raise RuntimeError('None of the selected comparisons were found in the pairwise file. Aborting...')
        if plot_params['boxes_colors']:
            ncomp = len(data['Comparison'].unique())
            if len(plot_params['boxes_colors']) != ncomp:
                print('Warning: the number of colors is different from the number of comparisons. Default colors will be used.')
                plot_params['boxes_colors'] = None
        #plotting
        plot_pi_pairwise(data, full_path, out_prefix, plot_params)

    else:#lineplot
        dfs = []
        nfiles = len(pi_files)
        for i,pi_file in enumerate(pi_files):
            df = pd.read_csv(pi_file,sep='\t')
            nrows = len(df['Num_Win'])
            cols = df.columns
            for col in cols:
                if col not in columns:
                    raise RuntimeError('Inconsistent column names. Aborting...')
            dfs.append(df)
        data = pd.concat(dfs, ignore_index = True)
        #plotting
        plot_pi_values(data, full_path, out_prefix, plot_params)
#
def manage_calc(args, full_path, out_prefix):
    aln_ids,aln_seqs = read_aln_to_array(args.infile) #read the aln file
    pi_values_in_windows = calculate_pi_per_window(aln_seqs, args.wsize, args.woverlap) #get pi values per window (array of values)
    pi_values_out_fname = out_prefix + '.nuc-diversity.txt'
    full_pi_values_out_fname = os.path.join(full_path,pi_values_out_fname)
    windows = [ i for i in range(1,len(pi_values_in_windows)+1) ]
    df_pi = pd.DataFrame(windows, columns = ['Num_Win'])
    df_pi['PI_value'] = pd.Series(pi_values_in_windows)
    df_pi['Comparison'] = pd.Series([args.pilabel]*len(windows))
    df_pi.to_csv(full_pi_values_out_fname, sep='\t', index=False)
    #pairwise_PIvalues is a list of nucleotide diversity values among pairwise sequence comparisons
    windows, labels, pairwise_PIvalues = generate_pairwise_alingments(aln_ids, aln_seqs, args.wsize, args.woverlap)
    pairwise_out_fname = out_prefix + '.pairwise.tsv'
    full_pairwise_out_fname = os.path.join(full_path,pairwise_out_fname)
    df_pairwise = pd.DataFrame(windows, columns = ['Num_Win'])
    df_pairwise['PI_value'] = pd.Series(pairwise_PIvalues)
    df_pairwise['Comparison'] = pd.Series(labels)
    df_pairwise.to_csv(full_pairwise_out_fname, sep='\t', index=False)        
#
def check_usage(args):
    'This function check the correctness of arguments passed to the program. If any problem is found, the program will exit.'
    if args.mode == 'calculate': #ask if mode = calculate
        if args.pifiles: #no pifile allowed in calculate mode.
            print("It is not allowed to supply a file with PI values in calculate mode. File will be omitted.")
        #--prefix argument is mandatory and in calculate mode --infile must be passed
        required_together = ('infile', 'pilabel')
        if not all( [ getattr(args,x) for x in required_together]):
            raise RuntimeError("Cannot supply --mode `calculate` without --infile and --pilabel")

        variable_to_plot = ['pifiles','title','labels','xlabel','xlimit','ylabel','ylimit','legend','legend_fontsize']
        if any( [ getattr(args, argument) for argument in variable_to_plot ] ): #if any argument is present in args, it will return True
            print("Plot parameters are ignored in `calculate` mode.")
    else: #ask for plot mode (line or boxplot)
        if args.mode == 'lineplot':
            if not args.pifiles:
                raise RuntimeError("No PI file was passed. Aborting...")
        if args.mode == 'boxplot':
            if not args.pairwise:
                raise RuntimeError("No pairwise PI file was passed (--pairwise). Aborting...")
#
def usage():
    parser = argparse.ArgumentParser(
    prog='NDEV.py',
    description='''Nucleotide Diversity Estimation and Visualization.\nEstimation of PI values (Nei, 1987) along an alignment using overlapping windows.''', epilog="""NDEV was developed using easy to install libraries/modules to facilitate its use. This is the first version of this program so if you have any suggestion or you find any bug/issue, please send me an email at niconm89@gmail.com. Thanks!""")
    parser.add_argument('-v', '--version', action='version', version='%(prog)s 1.0')
    parser.add_argument('--mode', type=str, required=True, choices=['calculate', 'lineplot', 'boxplot'], help='Working mode. In `full mode` the program will expect a FASTA format alignment file to estimate the nucleotide diversity in windows. In `plot mode`, the program will expect at least one file of PI values, following the format of one per line (as the output in `calculate mode`). [default: calculate]')
    parser.add_argument('--prefix', metavar='<str>', type=str, required=True, help='Prefix with full path (if necessary) for the outputs. By default, when full path is not present, the output files will be placed in the current working directory (from where the script is run).')
    #calculate mode
    parser.add_argument('--infile', metavar='<str>', type=str, required=False, help='Path to the FASTA alignment file to analyze. All sequences must have the same length (number of characters). In case you need to convert your alignment from PHYLIP or CLUSTAL format to FASTA, you can try some of the python scripts available at https://github.com/niconm89/scripts_bioinformatics.')
    parser.add_argument('--pilabel', metavar='<str>', type=str, required=False, help='Prefix with full path (if necessary) for the outputs. By default, when full path is not present, the output files will be placed in the current working directory (from where the script is run).')
    parser.add_argument('--wsize', metavar='<int>', type=int, required=False, default = 500, help='Windows size to split the alignment into parts and to calculate a PI value for each one. [default = 500]')
    parser.add_argument('--woverlap', metavar='<int>', type=int, required=False, default = 100, help='Windows overlap size. We recommend to try different sizes according to the windows size value. [default = 100] .')
    #shared plot parameters
    parser.add_argument('--title', metavar='<str>', type=str, required=False, default = '', help='Title of the final figure.')
    parser.add_argument('--title_fontsize', metavar='size', type=float, required=False, default = 28, help='Font size of Title. In addition to numeric values, the following strings are allowed: xx-small, x-small, small, medium, large, x-large, xx-large, larger, smaller and None. [ default = x-large ]')
    parser.add_argument('--xlabel', metavar='<str>', type=str, required=False, default = '', help='Label of the x axis.')
    parser.add_argument('--xlimit', metavar=('x0', 'x1'), type=float, nargs=2, required=False, help='X limits values for x0 >= 0 and x1 > x0.')
    parser.add_argument('--ylabel', metavar='<str>', type=str, required=False, default = '', help='Label of the y axis.')
    parser.add_argument('--ylimit', metavar=('y0', 'y1'), type=float, nargs=2, required=False, help='Y limits values for y0 >= 0 and y1 > x0.')
    parser.add_argument('--labels_fontsize', metavar='size', type=float, required=False, default = 23, help='Font size of x and y labels. In addition to numeric values, the following strings are allowed: xx-small, x-small, small, medium, large, x-large, xx-large, larger, smaller and None. [ default = x-large ]')
    parser.add_argument('--linestyle', metavar='style', type=str, required=False, default = 'solid', help='Line styles to use for the curves. By default all curves will have the same style. The styles must be in the same order as files. Available options are: solid,dashed,dashdot,dotted,None (quoted symbols also seems to work but we do not guarantee the functionality of the program). [default = solid]')
    parser.add_argument('--linewidth', metavar='width', type=float, required=False, default = 1.0, help='Width to use for all the curves (they will have the same). [default = 1.0]')
    parser.add_argument('--dpi', metavar='<int>', type=int, required=False, default = 72, help='Dot Per Inch (DPI) value for PNG format figure. [default = 72]')
    parser.add_argument('--alpha', metavar='<float>', type=float, required=False, default = 0.7, help='Dot Per Inch (DPI) value for PNG format figure. [default = 72]')
    parser.add_argument('--figsize', metavar=('height','width'), type=float, nargs='*', required=False, default = [15,20], help='Width to use for all the curves (they will have the same). [default = 1.0]')
    parser.add_argument('--axis_fontsize', metavar=('xsize','ysize'), type=float, nargs=2, required=False, default = [15,15], help='Font size of x and y axis labels. Allowed values are equal to option --title_fontsize.')
    parser.add_argument('--axis_linesize', metavar='size', type=float, required=False, default = 1.5, help='Integer value for degree of rotation of xtick labels in the boxplot.')
    parser.add_argument('--axis_linecolor', metavar='color', type=str, required=False, default = 'black', help='Integer value for degree of rotation of xtick labels in the boxplot.')
    parser.add_argument('--rotation', metavar=('x ticks rot', 'y ticks rot'), type=int, required=False, nargs=2, default = [0,0], help='Integer value for degree of rotation of xtick labels in the boxplot.')
    parser.add_argument('--theme_color', metavar='color', type=str, required=False, default = 'lightgrey', help='Color of the theme.')
    parser.add_argument('--grid', action='store_true', required=False, default=False, help='Draw grid on theme.')
    #boxplot
    parser.add_argument('--pairwise', metavar='pairwise PI values file', type=str, required=False, help='Mandatory in `boxplot mode`. Path to the file with pairwise PI values. File must be in tsv format and have the same name columns as the pairwise output of the program.')
    parser.add_argument('--select_comp', metavar=('file1', 'file2'), type=str, nargs='*', required=False, help='Select the pairwise comparisons (labels in `Comparison` column) to plot.')
    parser.add_argument('--violin', action='store_true', required=False, default=False, help='Generate a violin plot instead of a boxplot.')
    parser.add_argument('--jittler', action='store_true', required=False, default=False, help='Adding points to boxplot.')
    parser.add_argument('--boxes_colors', metavar=('color1', 'color2'), type=str, nargs='*', required=False, help='Color for each box (pairwise comparison). The number of colors introduced must be the same as the pairwise comparison to plot (it considers --select_comp).')
    #lineplot
    parser.add_argument('--pifiles', metavar=('file1', 'file2'), type=str, nargs='*', required=False, help='Mandatory in `plot mode`. Path to the file with PI values. File must be in txt format and have one PI value per line. In case you want to compare PI values across different alignments, it is possible to load more than one file. Files must have same number of PI values (lines), otherwise the program will failed.')
    parser.add_argument('--labels', metavar=('label1', 'label2'), type=str, action="extend", nargs='*', required=False, help='Labels for curve/s, they must have the same order as file names introduced. If plotting two or more curves, the number of labels must be the same as curves plotted.')
    parser.add_argument('--legend', metavar='<str>', type=str, required=False, choices=['upper left', 'upper right', 'lower left', 'lower right','upper center', 'lower center', 'center left', 'center right', ''], help='Position in the figure where the legend will be placed. Options (pass between quotes): upper left, upper right, lower left, lower right,upper center, lower center, center left, center right, None.')
    parser.add_argument('--legend_fontsize', metavar='size', type=str, required=False, help='Legend font size. Allowed values are equal to option --title_fontsize. [ default = x-large ]')
    parser.add_argument('--lines_colors', metavar=('color1','color2'), type=str, nargs='*', required=False, default = 'red' ,help='[string] Colors to use for the curves. The program will assign color in the same order as files introduced. By default each curve will have a randomly assigned color. In case you need more information, read the documentation of the matplotlib.colors module.')
    args=parser.parse_args()

    check_usage(args)
    return args
#
#%% Main program
def main():
    args = usage()
    full_path,out_prefix = check_prefix(args.prefix)

    pi_files = []
    if args.mode == 'calculate': #calculate pi values from alignment
        if check_file_exists(args.infile): #if alignment file exists
            manage_calc(args, full_path, out_prefix)
    else: #line or boxplot mode
        plot_params = { 'title':args.title, 'title_fontsize':args.title_fontsize, 'labels':args.labels, 
                       'xlabel':args.xlabel, 'ylabel':args.ylabel, 'axis_fontsize':args.axis_fontsize, 
                       'legend':args.legend, 'legend_fontsize':args.legend_fontsize, 
                       'lines_colors':args.lines_colors, 'linestyle':args.linestyle, 'linewidth':args.linewidth, 
                       'dpi':args.dpi, 'rotation':args.rotation, 'pairs':args.select_comp,  'figsize':args.figsize, 
                       'axis_linecolor':args.axis_linecolor, 'axis_linesize':args.axis_linesize, 
                       'theme_color':args.theme_color, 'labels_fontsize':args.labels_fontsize,
                       'violin':args.violin, 'jittler':args.jittler, 'boxes_colors':args.boxes_colors,
                       'ylimit':args.ylimit, 'alpha':args.alpha }
        if args.mode == 'boxplot':
            pairwise_file = args.pairwise
            check_file_exists(pairwise_file) #check if file passed exists
            manage_plots(args.mode, pairwise_file, full_path, out_prefix, plot_params)
        else: #'lineplot'
            pi_files = args.pifiles
            for pf in pi_files: #check if files passed exists
                check_file_exists(pf)
            #call to plot
            manage_plots(args.mode, pi_files, full_path, out_prefix, plot_params)
# End of main program

# Name of the program
if __name__ == '__main__':
	main()
	
#%% Testing
'''
Testing here!
'''
