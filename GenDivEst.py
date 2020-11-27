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
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
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
def calc_pi(list_of_mutations, number_of_comparisons):
    'Receive a numpy array with aligment mutations per positions in the aligment and return the pi value'
    pi = list_of_mutations.sum() / number_of_comparisons / len(list_of_mutations)
    return pi
#
def compare_two_sequence(seq1,seq2):
    'Receive two numpy arrays with one sequence in each one and return a numpy array with 0 or 1 in alignment positions where the sequences are equal or not, respectively.'
    seqlen = len(seq1)
    list_of_differences = np.zeros(seqlen)
    for i in range(0,seqlen):
        if seq1[i] != seq2[i]:
            list_of_differences[i] = 1
    return list_of_differences
#
def number_of_combinations(n, m):
    'Calculate the number of possible combinations or n elements grouping by m elements.'
    return factorial(n) // (factorial(m) * factorial(n - m))
#
def analyze_window(window_aln):
    'Receive a window of an alignment and return the pi value.'
    ncomparisons = number_of_combinations(len(window_aln),2)
    #ncomparisons = comb(len(window_aln),2) #only for python3.8+
    number_seqs = len(window_aln)
    wsize = len(window_aln[0])
    mutations_window = np.zeros(wsize)
    for i in range(number_seqs):
        for j in range(i+1,number_seqs):
            mutations_window += compare_two_sequence(window_aln[i],window_aln[j])
    pi = calc_pi(mutations_window,ncomparisons)
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
def save_pi_values(output_file,list_of_pi):
    'Receive the name of an output file and a list with pi values per windows, and save the pi values in a numpy txt format.'
    np.savetxt(output_file, list_of_pi)
#
def generate_pairwise_alingments(alignment, windows_size, windows_overlap):
    'Receive an alignment in a numpy array (with windows size and overlap) and return a numpay array (matrix) with PI values per each pairwise comparison in the alignment.'
    nseqs = len(alignment)
    num_comparisons = number_of_combinations(nseqs, 2)
    numwin = number_of_windows(len(alignment[0]), windows_size, windows_overlap)
    pairwise_PIvalues = np.zeros((num_comparisons,numwin))

    count = 0
    for i in range(nseqs):
        for j in range(i+1,nseqs):
            pairwise_aln = np.array([alignment[i],alignment[j]])
            pairwise_PIvalues[count] = calculate_pi_per_window(pairwise_aln, windows_size, windows_overlap)
            count += 1
    
    return pairwise_PIvalues
#
def plot_pi_pairwise(aln_ids, pi_file, full_path, out_prefix, plot_params, selected_columns=False):
    ''''
    Receive a numpy array with PI values (pi_file) and return a bloxplot figure for all comparisons.
    '''
    data = load_PIvalues_from_file(pi_file) #load the file
    ncomp = number_of_combinations(len(aln_ids),2) #number of comparisons (boxplots)
    
    #setting the figure
    plt.figure()
    width1 = 20
    height1 = 10
    width_height_1 = (width1, height1)
    colors = generate_default_colors(ncomp) #generate colors    
    plt.figure(figsize=width_height_1, facecolor='w')    

    #making the boxplot
    indices = [i+1 for i in range(38)]
    columns = [ 'Comp'+str((i+1)) for i in range(ncomp) ]    
    df = pd.DataFrame(data, columns = columns, index = indices)
    
    if not selected_columns:
        selected_columns = columns
    boxplot = df.boxplot(column=selected_columns)
    if not plot_params['rotation']:
        plot_params['rotation'] = 0
    if plot_params['axis_fontsize']:
        boxplot = df.boxplot(grid=True, rot=plot_params['rotation'], fontsize=plot_params['axis_fontsize'])
    else:
        boxplot = df.boxplot(grid=True, rot=plot_params['rotation'], fontsize=15)
    
    if plot_params['ticks_size']:
        plt.xticks(fontsize=plot_params['ticks_size'][0], rotation=90)
        plt.yticks(fontsize=plot_params['ticks_size'][1], rotation=90)
    else:
        plt.xticks(fontsize=14, rotation=90)
        plt.yticks(fontsize=14, rotation=90)

    '''
    boxplot 
    color=dict(boxes='r', whiskers='r', medians='r', caps='r'),
             boxprops=dict(linestyle='-', linewidth=1.5),
             flierprops=dict(linestyle='-', linewidth=1.5),
             medianprops=dict(linestyle='-', linewidth=1.5),
             whiskerprops=dict(linestyle='-', linewidth=1.5),
             capprops=dict(linestyle='-', linewidth=1.5),
             showfliers=False, grid=True, rot=0
    '''
    #adding colors
    if 'xlabel_fontsize' in plot_params:
        boxplot.set_xlabel('Pairwise comparisons', fontsize = plot_params['xlabel_fontsize'][0])
        boxplot.set_ylabel('Nucleotide Diversity', fontsize = plot_params['xlabel_fontsize'][1])
    else:
        boxplot.set_xlabel('Pairwise comparisons', fontsize = 25)
        boxplot.set_ylabel('Nucleotide Diversity', fontsize = 25)

    out_figure = os.path.join(full_path, out_prefix+'.pairwise_boxplot.png') #creating full output path
    plt.savefig(out_figure, bbox_inches='tight', dpi = plot_params['dpi']) #save figure
    
    return     
#
def plot_pi_values(list_of_files, full_path, out_prefix, plot_params):
    'Receive a list of path to files with pi values and a dictionary with plot parameters, and plot (save) a figure with all curves (one per file).'
    #create a list of numpy arrays with PI values to plot
    pi_curves = [ load_PIvalues_from_file(pi_file) for pi_file in list_of_files ]
    
    nwin_max = max([ len(pi_curve) for pi_curve in pi_curves ])
    x_axis = np.linspace(1, nwin_max, num=nwin_max)
    
    width1 = 20
    height1 = 10
    width_height_1 = (width1, height1)
    plt.figure(figsize=width_height_1, facecolor='w')
    
    labels = plot_params['labels']
    colors = plot_params['colors']
    linestyles = plot_params['linestyles']
    linewidths = plot_params['linewidth']
    
    for i,pi in enumerate(pi_curves):
        plt.plot(x_axis, pi, label=labels[i], color=colors[i], linewidth=linewidths[i], linestyle=linestyles[i])
    
    plt.xlim(0, nwin_max+1) # x axis range
    if plot_params['ticks_size']:
        plt.xticks(np.linspace(0, nwin_max, 10), fontsize = plot_params['ticks_size'][0]) # Ponemos marcas (ticks) en el eje x
    else:
        plt.xticks(np.linspace(0, nwin_max, 10)) # Ponemos marcas (ticks) en el eje x
    pi_max =  max([max(windows) for windows in pi_curves ]) #len(list_of_pi)
    limit_factor = pi_max * 0.05 #upper bound of y-axis
    plt.ylim(0, pi_max+limit_factor) # y axis range
    if plot_params['ticks_size']:
        plt.yticks(np.linspace(0, pi_max+limit_factor, 10), fontsize = plot_params['ticks_size'][1]) # Ponemos marcas (ticks) en el eje y
    else:
        plt.yticks(np.linspace(0, pi_max+limit_factor, 10)) # Ponemos marcas (ticks) en el eje y

    #setting titles in figure
    if plot_params['title']:
        if plot_params['title_fontsize']:
            plt.title(plot_params['title'], fontsize=plot_params['title_fontsize']) # It is used to give title to the graph. Title is passed as the parameter to this funtion.
        else:
            plt.title(plot_params['title'], fontsize='x-large') # It is used to give title to the graph. Title is passed as the parameter to this funtion.
    if plot_params['xlabel']:
        if plot_params['axis_fontsize']:
            plt.xlabel(plot_params['xlabel'], fontsize=plot_params['axis_fontsize'][0])
        else:
            plt.xlabel(plot_params['xlabel'], fontsize='x-large') # It sets the label name at X-axis. Name of X-axis is passed as argument to this funtion.
    if plot_params['ylabel']:
        if plot_params['axis_fontsize']:
            plt.ylabel(plot_params['ylabel'], fontsize=plot_params['axis_fontsize'][1])
        else:
            plt.ylabel(plot_params['ylabel'], fontsize='x-large') # It sets the label name at Y-axis. Name of Y-axis is passed as argument to this funtion.
    if plot_params['legend']:#show legend
        if plot_params['legend_fontsize']:
            plt.legend(loc=plot_params['legend'], fontsize=plot_params['legend_fontsize']) 
        else:
            plt.legend(loc=plot_params['legend'], fontsize='x-large')

    #saving the figure
    out_figure = os.path.join(full_path, out_prefix) #creating full output path
    plt.savefig(out_figure, bbox_inches='tight', dpi = plot_params['dpi']) #save figure
#
def generate_default_colors(nfiles):
    'Receive a number and return this amount of different colors.'
    CSS_colors = mcolors.CSS4_COLORS.copy()
    if 'white' in CSS_colors:
        del CSS_colors['white']
    return rchoices(list(CSS_colors.values()),k=nfiles)
#
def manage_plots(pi_files, full_path, out_prefix, plot_params):
    'Receive a list of files and a dictionary with plot parameters and admin the reformat variables to make the plot.'
    number_of_files = len(pi_files)
    ##reformat labels
    #check that if the number of labels coincide with number of files, otherwise we abort!
    if plot_params['labels']:
        if number_of_files != len(plot_params['labels']):
            raise RuntimeError("Different numbers of files and labels to plot. Aborting...")
    else: #generate empty labels
        plot_params['labels'] = [''] * number_of_files
    ##reformat linestyles
    #check that if the number of linestyles coincide with number of files, otherwise all curves will have the same style
    if isinstance(plot_params['linestyles'],list):
        linestyles = plot_params['linestyles']
    else:
        linestyles = [ plot_params['linestyles'] ]
        plot_params['linestyles'] = linestyles

    if number_of_files != len(plot_params['linestyles']):
        if len(plot_params['linestyles']) == 1:
            plot_params['linestyles'] = plot_params['linestyles'] * number_of_files
        else:
            print("All the curves will have the same linestyle (first supplied).")
            plot_params['linestyles'] = plot_params['linestyles'][0] * number_of_files
    ##reformat colors
    if not plot_params['colors']:
        print("No color/s passed. We will generate for you " + str(number_of_files) + " random colors.")
        plot_params['colors'] = generate_default_colors(number_of_files)
    else:
        if number_of_files != len(plot_params['colors']):
            if number_of_files == 1:
                print("Number of colors passed is 1 and the number of files passed is " + str(number_of_files) + ". All the curves will have the same color (" + "".join(plot_params['colors']) + ")")
                plot_params['colors'] = plot_params['colors'] * number_of_files
            else:
                print("Number of colors passed is " + str(len(plot_params['colors'])) + " and the number of files passed is: " + str(number_of_files) + ". We will generate for you " + str(number_of_files) + " random colors.")
                plot_params['colors'] = generate_default_colors(number_of_files)
    ##reformat linewidth
    if isinstance(plot_params['linewidth'],list):
        linewidths = plot_params['linewidth']
    else:
        linewidths = [ plot_params['linewidth'] ]
        plot_params['linewidth'] = linewidths
    if number_of_files != len(plot_params['linewidth']):
        if len(plot_params['linewidth']) == 1:
            plot_params['linewidth'] = plot_params['linewidth'] * number_of_files
        else:
            print("All the curves will have the same linewidth (first supplied).")
            plot_params['linewidth'] = plot_params['linewidth'][0] * number_of_files
    #reformat font sizes

    plot_pi_values(pi_files, full_path, out_prefix, plot_params)
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
def check_usage(args):
    'This function check the correctness of arguments passed to the program. If any problem is found, the program will exit.'
    if args.mode == 'full': #ask if mode = full
        if args.pifile: #no pifile allowed in full mode.
            print("It is not allowed to supply a file with PI values in full mode. The file will be omitted.")
        #--prefix argument is mandator and in full mode --pifile must be passed
        required_together = ('infile','prefix')
        if not all( [ getattr(args,x) for x in required_together]):
            raise RuntimeError("Cannot supply --mode `full` without --infile and --prefix")
        #arguments only necessary for making the figure. Any of them will be omitted if --noplot argument is passed to the program.
        used_for_plot = ['title','labels','xlabel','ylabel','legend','colors','title_fontsize','axis_fontsize','legend_fontsize']
        if args.noplot: #if args.noplot = True, no plot will be created
            if any( [ getattr(args, argument) for argument in used_for_plot ] ): #if any argument is present in args, it will return True
                print("--noplot supplied. No figure will be create. One or some of the following parameters was also supplied and this is not compatible, they will be omitted:")
                print(used_for_plot)
    else: #ask for plot mode
        if args.noplot:#plot mode and --noplot not compatible
            raise RuntimeError("Cannot supply --mode plot with --noplot")

        used_for_pi_calc = ['infile','numwin']
        for param in used_for_pi_calc:
            if getattr(args,param):
                print(param + ' is not allowed in --mode plot and will be omitted.')
#
def usage():
    parser = argparse.ArgumentParser(
    prog='NDEV.py',
    description='''Nucleotide Diversity Estimation and Visualization.\nEstimation of PI values (Nei, 1987) along an alignment using overlapping windows.''', epilog="""NDEV was developed using easy to install libraries/modules to facilitate its use. This is the first version of this program so if you have any suggestion or you find any bug/issue, please send me an email at niconm89@gmail.com. Thanks!""")
    parser.add_argument('-v', '--version', action='version', version='%(prog)s 1.0')
    parser.add_argument('--mode', type=str, required=True, choices=['full', 'plot'], default='full', help='Working mode. In `full mode` the program will expect a FASTA format alignment file to estimate the nucleotide diversity in windows. In `plot mode`, the program will expect at least one file of PI values, following the format of one per line (as the output in `full mode`). [default: full]')
    #full mode
    parser.add_argument('--infile', metavar='<str>', type=str, required=False, help=' Path to the FASTA alignment file to analyze. All sequences must have the same length (number of characters). In case you need to convert your alignment from PHYLIP or CLUSTAL format to FASTA, you can try some of the python scripts available at https://github.com/niconm89/scripts_bioinformatics.')
    parser.add_argument('--prefix', metavar='<str>', type=str, required=True, help='Prefix with full path (if necessary) for the outputs. By default, when full path is not present, the output files will be placed in the current working directory (from where the script is run).')
    parser.add_argument('--wsize', metavar='<int>', type=int, required=False, default = 500, help='Windows size to split the alignment into parts and to calculate a PI value for each one. [default = 500]')
    parser.add_argument('--woverlap', metavar='<int>', type=int, required=False, default = 100, help='Windows overlap size. We recommend to try different sizes according to the windows size value. [default = 100] .')
    parser.add_argument('--pairwise_plot', action='store_true', required=False, default=False, help='Generate a boxplot of nucleotide diversity values for pairwise comparisons (sequences in the alignment)')
    parser.add_argument('--rotation', action='store_true', required=False, default=False, help='Integer value for degree of rotation of xtick labels in the boxplot.')
    parser.add_argument('--noplot', action='store_true', required=False, default=False, help='Disable plotting in `full mode`.')
    #plot mode & options for plotting
    parser.add_argument('--pifile', metavar=('file1', 'file2'), type=str, nargs='*', required=False, help='Mandatory in `plot mode`. Path to the file with PI values. File must be in txt format and have one PI value per line. In case you want to compare PI values across different alignments, it is possible to load more than one file. Files must have same number of PI values (lines), otherwise the program will failed.')
    parser.add_argument('--title', metavar='<str>', type=str, required=False, help='Title of the final figure.')
    parser.add_argument('--title_fontsize', metavar='size', type=str, required=False, help='Font size of Title. In addition to numeric values, the following strings are allowed: xx-small, x-small, small, medium, large, x-large, xx-large, larger, smaller and None. [ default = x-large ]')
    parser.add_argument('--labels', metavar=('label1', 'label2'), type=str, action="extend", nargs='*', required=False, help='Labels for curve/s, they must have the same order as file names introduced. If plotting two or more curves, the number of labels must be the same as curves plotted.')
    parser.add_argument('--labels_fontsize', metavar='<str>', type=str, required=False, help='Font size of Title. In addition to numeric values, the following strings are allowed: xx-small, x-small, small, medium, large, x-large, xx-large, larger, smaller and None. [ default = x-large ]')
    parser.add_argument('--xlabel', metavar='<str>', type=str, required=False, help='Label of the x axis.')
    parser.add_argument('--xlimit', metavar=('x0', 'x1'), type=float, nargs=2, required=False, help='X limits values for x0 >= 0 and x1 > x0.')
    parser.add_argument('--ylabel', metavar='<str>', type=str, required=False, help='Label of the y axis.')
    parser.add_argument('--ylimit', metavar=('y0', 'y1'), type=float, nargs=2, required=False, help='Y limits values for y0 >= 0 and y1 > x0.')
    parser.add_argument('--axis_fontsize', metavar=('xsize','ysize'), type=str, action="extend", nargs=2, required=False, help='Font size of x and y axis labels. Allowed values are equal to option --title_fontsize.')
    parser.add_argument('--ticks_size', metavar=('xsize','ysize'), type=str, action="extend", nargs=2, required=False, help='Font size of x and y axis ticks labels. Allowed values are equal to option --title_fontsize.')
    parser.add_argument('--legend', metavar='<str>', type=str, required=False, choices=['upper left', 'upper right', 'lower left', 'lower right','upper center', 'lower center', 'center left', 'center right', ''], help='Position in the figure where the legend will be placed. Options (pass between quotes): upper left, upper right, lower left, lower right,upper center, lower center, center left, center right, None.')
    parser.add_argument('--legend_fontsize', metavar='size', type=str, required=False, help='Legend font size. Allowed values are equal to option --title_fontsize. [ default = x-large ]')
    parser.add_argument('--colors', metavar=('color1','color2'), type=str, nargs='*', required=False, help='[string] Colors to use for the curves. The program will assign color in the same order as files introduced. By default each curve will have a randomly assigned color. In case you need more information, read the documentation of the matplotlib.colors module.')
    parser.add_argument('--linestyles', metavar=('style1','style2'), type=str, nargs='*', required=False, default = 'solid', help='Line styles to use for the curves. By default all curves will have the same style. The styles must be in the same order as files. Available options are: solid,dashed,dashdot,dotted,None (quoted symbols also seems to work but we do not guarantee the functionality of the program). [default = solid]')
    parser.add_argument('--linewidth', metavar=('width1','width2'), type=float, nargs='*', required=False, default = 1.0, help='Width to use for all the curves (they will have the same). [default = 1.0]')
    parser.add_argument('--dpi', metavar='<int>', type=int, required=False, default = 72, help='Dot Per Inch (DPI) value for PNG format figure. [default = 72]')
    args=parser.parse_args()

    check_usage(args)
    return args
#
#%% Main program
def main():
    args = usage()
    full_path,out_prefix = check_prefix(args.prefix)

    pi_files = []
    if args.mode == 'full': #calculate pi values from alignment
        if check_file_exists(args.infile): #if alignment file exists
            aln_ids,aln_seqs = read_aln_to_array(args.infile) #read the aln file
            pi_values_in_windows = calculate_pi_per_window(aln_seqs, args.wsize, args.woverlap) #get pi values per window (array of values)
            pi_values_out_fname = out_prefix + '.pi-values.txt'
            full_pi_values_out_fname = os.path.join(full_path,pi_values_out_fname)
            save_pi_values(full_pi_values_out_fname,pi_values_in_windows)
            pi_files.append(full_pi_values_out_fname)
            #pairwise_PIvalues is a list of nucleotide diversity values among pairwise sequence comparisons
            pairwise_PIvalues = generate_pairwise_alingments(aln_seqs, args.wsize, args.woverlap)
            pairwise_out_fname = out_prefix + '.pairwise.pi-values.txt'
            full_pairwise_out_fname = os.path.join(full_path,pairwise_out_fname)
            save_pi_values(full_pairwise_out_fname,pairwise_PIvalues.transpose())
            #manage_plot
            if not args.noplot:
                plot_params = { 'title':args.title, 'title_fontsize':args.title_fontsize, 'labels':args.labels, 'labels_fontsize':args.labels_fontsize, 'xlabel':args.xlabel, 'ylabel':args.ylabel, 'axis_fontsize':args.axis_fontsize, 'ticks_size':args.ticks_size, 'legend':args.legend, 'legend_fontsize':args.legend_fontsize, 'colors':args.colors, 'linestyles':args.linestyles, 'linewidth':args.linewidth, 'dpi':args.dpi, 'rotation':args.rotation }
                manage_plots(pi_files, full_path, out_prefix, plot_params)
                if args.pairwise_plot:
                    plot_pi_pairwise(aln_ids, full_pairwise_out_fname, full_path, out_prefix, plot_params)
    else: #plot mode
        pi_files = args.pifile
        for pf in pi_files:
            check_file_exists(pf)
        plot_params = { 'title':args.title, 'title_fontsize':args.title_fontsize, 'labels':args.labels, 'labels_fontsize':args.labels_fontsize, 'xlabel':args.xlabel, 'ylabel':args.ylabel, 'axis_fontsize':args.axis_fontsize, 'ticks_size':args.ticks_size, 'legend':args.legend, 'legend_fontsize':args.legend_fontsize, 'colors':args.colors, 'linestyles':args.linestyles, 'linewidth':args.linewidth, 'dpi':args.dpi, 'rotation':args.rotation }
        manage_plots(pi_files, full_path, out_prefix, plot_params)
#        if args.pairwise_plot:
#           plot_pi_pairwise(aln_ids, full_pairwise_out_fname, full_path, out_prefix, plot_params)
# End of main program

# Name of the program
if __name__ == '__main__':
	main()
	
#%% Testing
'''
Testing here!
'''
