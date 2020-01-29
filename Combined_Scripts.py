import matplotlib.pyplot as plt
import random
import csv
import math
import statistics
import pandas as pd
import matplotlib
import numpy as np
from statistics import mean
import sklearn.metrics
from math import sqrt
import rpy2.robjects as ro
from rpy2.robjects import r, pandas2ri
from rpy2.robjects.packages import importr
import re
from scipy.optimize import minimize_scalar
import argparse
import os
import glob
from matplotlib.backends.backend_pdf import FigureCanvasPdf, PdfPages
from matplotlib.figure import Figure
from matplotlib.pyplot import cm
import randomcolor
#This program is split into three scripts. 
#Script_1 aligns the template with your data, and creates a 
#table that is readable by script 2. Script 2 pumps out the paramters needed for curve fitting.
#Script three then uses those parameters and plots your data. 

with open('IC50_and_Scatter_Index.csv', mode='a') as ic50_SI:
    fnames = ['drug', 'IC50', 'Scatter Index', 'b', 'c', 'd', 'e']
    writer = csv.DictWriter(ic50_SI, fieldnames=fnames)    
    writer.writeheader()

class Script_1: 
    def __init__(self):
        
        ''' Constructor for this class. '''

    def format_script(self, template_name, raw_data_file): 
        '''
        Args: template name: Name of the plate template, 
              raw_data_file: data from plate reader
        
        Returns: None 


        #The template should be modified by the user to fit the actual template of their 
        #experiment. For example, it should be correctly labelled which wells have DSMO
        #and which have what drug. Use the given text inside the template as a guide. 
        #####                                        #####
        ##### If you have empty wells, label it NA so that the script can remove it from the dataframe when it gets processed#####
        #####


        '''
       
        with open(template_name) as df:
            lines = df.readlines()
        drugs = lines[1:17] #these are the lines for each section of the plate template.  
        #There are three templates in the plate template file. The drug template, replicate template,
        # and dosage template. We are
        #going to make these three different templates three different columns. We then will put 
        #those three columns into one dataframe, and then add 1 more column, containing the luminosity of each well--the luminosity is the raw data, which will come from the raw_data_file.
        replicates = lines[19:35]
        dosage = lines[37:53]

        drug_list = []
        for i in drugs: 
            lsp = i.strip().split(",")
            drug_list.append(lsp)

        replicates_list = []
        for i in replicates:
            lsp = i.strip().split(",")
            replicates_list.append(lsp)

        dosage_list = []
        for i in dosage:
            lsp = i.strip().split(",")
            dosage_list.append(lsp)


        drugdf = pd.DataFrame(drug_list)
        drugdf.columns = range(1,25)
        repdf = pd.DataFrame(replicates_list)
        repdf.columns = range(1,25)
        dosedf= pd.DataFrame(dosage_list)
        dosedf.columns = range(1,25)


        ##stacked stacks all of the columns together##
        drugstacked = drugdf.stack().reset_index()
        drugstacked.columns = ['row','col','drug']
        repstacked = repdf.stack().reset_index()
        repstacked.columns = ['row','col','replicate']
        dosestacked = dosedf.stack().reset_index()
        dosestacked.columns = ['row','col','dosage']

        merge_1 = pd.merge(drugstacked,repstacked)
        merge_2 = pd.merge(merge_1, dosestacked)
        #now, you should have three columns, one for drug_#, one for replicate #, and one for dose.

        ###Merging template with data ###
        with open(raw_data_file) as df2:
            lines2 = df2.readlines()

        lum = lines2[1:17]
        lum_list = []
        for i in lum:
            lsp2 = i.strip().split(",")
            lum_list.append(lsp2)
        lumdf = pd.DataFrame(lum_list)
        lumdf.drop(lumdf.columns[0], axis=1, inplace=True) #unneseccary columns
        lumdf.drop(lumdf.columns[24], axis=1, inplace=True) #unnessecary columns 
        lum_stack =lumdf.stack().reset_index()

        final_merge = pd.merge(merge_2, lum_stack, left_on=['row','col'], right_on=['level_0','level_1'])
        final_merge.drop(final_merge.columns[5], axis=1, inplace=True)
        final_merge.drop(final_merge.columns[5], axis=1, inplace=True)
        final_merge.columns = ['row','col','drug', 'replicate', 'dosage', 'response']


        ###Editing for final export###
        final_merge = final_merge[~final_merge['replicate'].isin(['NA'])] #Dropping all rows that contain NA
        final_merge.drop(final_merge.columns[0], axis=1, inplace=True)
        final_merge.drop(final_merge.columns[0], axis=1, inplace=True)
        final_merge.drop(final_merge.columns[1], axis=1, inplace=True)
        final_merge = final_merge[['dosage','response','drug']]

        final_merge['response']=final_merge.response.astype('int64')
        final_merge['response'] = final_merge['response'].div(1000).round(4)
        #the template contains 1,2,3, as a way to mark the drugs. We need to change these labels to drug_#
        final_merge['drug'] = final_merge['drug'].str.replace('1','drug_1') 
        final_merge['drug'] = final_merge['drug'].str.replace('2','drug_2')
        final_merge['drug'] = final_merge['drug'].str.replace('3','drug_3')
        final_merge['drug'] = final_merge['drug'].str.replace('4','drug_4')
        final_merge['drug'] = final_merge['drug'].str.replace('5','drug_5')
        final_merge['drug'] = final_merge['drug'].str.replace('6','drug_6')
        final_merge['drug'] = final_merge['drug'].str.replace('7','drug_7')
        final_merge['drug'] = final_merge['drug'].str.replace('8','drug_8')




        final_merge.to_csv(r'reformatteddata_for_scriptII.csv', header=True, index = False)


#This script calculates the parameters by sourcing in an R script. 
class Script_2:
    
    def __init__(self):
    
        ''' Constructor for this class. '''
    def find_parameters(self, drug_name): 
        '''
        Args: drug_name: this is the name of the drug and should be in the format 'drug_#'.
        See the template for an example. The drug_name will be read in from a drug_list, which is
        given in the following function.

        Returns: None 
        '''
        r_source = ro.r['source']
        fitting_func = r_source('Fitting_single_drug.R')[0] #this is the R code. 
        dataframe = pd.read_csv("reformatteddata_for_scriptII.csv")
        df1 = dataframe[dataframe['drug'].str.contains(drug_name)] 
        #print(df1)
        f = open("parameters.txt", "a")
        #convert dataframe to R (psudocode)
        pandas2ri.activate()
        r_df = pandas2ri.py2ri(df1)


        #run drc's drm function from BK code 
        result = fitting_func(r_df)

        with open("output" + "_" + str(drug_name) + ".txt", "a") as f:
            print(result, file=f) #output textfile. 
    
    #we are going to append all the drug names (e.g., drug_#) to a list.
    #this list will be read by the find_parameters to determine which drugs are to be plotted/get the parameters from.  
    def drug_list_def(self):
        '''
        Args: None

        Returns: a drug list that will be read by find_parameters. 
        '''
        drug_list = []
        with open('reformatteddata_for_scriptII.csv') as df: #we take in our modified dataframe from Script 1 and basically 
            #get all of the different 'Drug_#'s and append them to this list. This should give us a list of all of our drugs
            #that were tested. 
                    split_list = []
                    
                    lines = df.read()
                    lines2 = re.split(', |, |\*|\n',lines)
                    for element in lines2:
                        element2 = element.split(',')
                        split_list.append(element2)
                    for element in split_list:
                        if len(element) == 3 and 'drug' not in element:
                            drug_list.append(element[2])
                    drug_list = list(dict.fromkeys(drug_list))
        return drug_list



class Script_3:
    def __init__(self):
      ''' Constructor for this class. '''

    #the text files contained the paramters from Script 2 vary; sometimes 
    #the desired parameters are at @ index 0,1,2,3 and other times, they are at
    #index 1,2,3,4. This is based on the length of the particular line where the 
    #paramters are located. The parameters will always be on line 8; it just depends
    #where on line 8.
    '''
    Args: parameters_test_file: your 'output_drug_#.txt' from Script 2. drug_x: whatever drug is being read in from the drug_list

    Returns: None 
    '''
    def format_parameters(self, parameters_test_file, drug_x):
        desired_lines = []
        with open(parameters_test_file) as df:
            lines = df.readlines()
            lines2 = lines[8]
            lines3 = lines2.split()
        if len(lines3) == 4:
            desired_lines.append(float(lines3[0]))
            desired_lines.append(float(lines3[1]))
            desired_lines.append(float(lines3[2]))
            desired_lines.append(float(lines3[3]))

        elif len(lines3) == 5:
            desired_lines.append(float(lines3[1]))
            desired_lines.append(float(lines3[2]))
            desired_lines.append(float(lines3[3]))
            desired_lines.append(float(lines3[4]))

        drug_dict_para[drug_x] = desired_lines
        return drug_dict_para

    def find_mean(self, formatted_datafile, drug_x): #here, the first parameter is reformatteddata_for_scriptII.csv
        '''
        Args: formatted_datafile: your datafile that was created in Script 1
              drug_x: whatever drug is being read in from the drug_list
          
        Returns: drug_dict_x, drug_dict_y, drug_dict_std. Each one is a dictionary; the
         first two are dictionaries that carry means of the x and y data, the third one caries the standard deviation for each drug @ a 
         specific concentration. 
        '''
        x=[] #all of the x  data
        y=[] #all of the y data
        std_graph = []
        mean_new_x = [] #the mean of each three replicates. Since this is the concentration, and the concentration is constant 
        #throughout the replicates, this shouldn't be any different than the original x values.
        mean_new_y = [] #the mean of the luminosities amongst the three replicates for each drug @ each specific concentration. 
        with open(formatted_datafile, 'r') as csvfile:
            data = pd.read_csv(csvfile)
            data = data.sort_values(['drug', 'dosage'], ascending=[True, True])
            data_list = data.values.tolist()
            for row in data_list[1:]:
                if drug_x in row: 
                    x.append(float(row[0]))
                    y.append(float(row[1]))
        n = 3 #we are clustering every three rows because those should be the three
        #replicates for each drug @ each specific concentration. 
        for i in range(0, len(x), n):
            mean_new_x.append(mean(x[i:i+n]))
        for i in range(0, len(y), n):
            mean_new_y.append(mean(y[i:i+n]))
        for i in range(0,len(y),n):
            std_graph.append(statistics.stdev(y[i:i+n]))

        xdata_drug1 = np.array(mean_new_x)
        ydata_drug1 = np.array(mean_new_y)
        std_data = np.array(std_graph)

        drug_dict_x[drug_x] = (xdata_drug1)
        drug_dict_y[drug_x] = (ydata_drug1)
        drug_dict_std[drug_x] = (std_data)

        return drug_dict_x, drug_dict_y, drug_dict_std



    def sigmoid_new(self, x, b, c, d, e): #sigmoid function 
        '''
        args: x (x data), b (parameter of sigmoid function), c (parameter of sigmoid function), d (parameter of sigmoid function), 
        e (parameter of sigmoid function)

        returns: equation of sigmoid function
        '''
        y = c + ((d-c) / (1 + np.exp(b*(np.log(x)-np.log(e)))))
        return y

    def objective_new(self, x, b, c, d, e):
        '''
        args: same as sigmoid_new
        returns: data used for IC50 calcuation
        '''
        return ((((d-c)/2)+c) - self.sigmoid_new(x, b, c, d, e)) ** 2

    def plotDRC_new(self, xdata,ydata, std, b, c, d, e, err = None, IC50lines=True,
                xlabel="Concentration (micromol/microliter)",ylabel="Luminosity",title="Drug Name",shade=False,
            shadeColor = 'white', labelpos = "lower left",logscale=True, ax = None, drug = '', color = None, write = False, label = False):
        
        
        popt = (b, c, d, e)
        res = minimize_scalar(self.objective_new, args=tuple(popt))
        x = np.logspace(np.min(np.log(xdata)), np.max(np.log(xdata)), 10000)
        y = self.sigmoid_new(x, *popt)
        y_pred_rmse = self.sigmoid_new(xdata, *popt)
        ##If you want them on the same graph, then change ax = ax in the function call.##
        if ax is None:
            fig, ax = plt.subplots(1,1)
        ax.plot(np.ndarray.flatten(xdata), np.ndarray.flatten(ydata), 'o', label='data points for' + ' ' + drug, color = color) #plotting of data
        ax.errorbar(xdata, ydata, yerr=std, fmt='o', capsize=3, color = color)
        ax.plot(np.ndarray.flatten(x),np.ndarray.flatten(y), label= None, color = color) #plotting of curve
        ax.legend(bbox_to_anchor=(1,0), loc="lower right", prop={'size': 6})
        ax.set_title('Drug Response')
        ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(True); ax.spines["bottom"].set_linewidth(1.5)
        ax.spines['left'].set_visible(True); ax.spines["left"].set_linewidth(1.5)
        # if IC50lines:
        #     p = np.polyfit(np.log(x), y, 1)
        #     d = np.polyval(p,np.log(x))
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        
        ax.xaxis.set_tick_params(width=1.5,size=10)
        ax.yaxis.set_tick_params(width=1.5)
        ax.set_xticks([0,0.0001,0.001,0.01,0.1, 1, 10,100,1000])
        if logscale: #setting the x axis to a log scale. 
            ax.set_xscale('log')
            ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
            ax.xaxis.set_tick_params(which='minor', width=1.5,size=5)
        #ax.set_xlim(xlim[0]+0.0000001,xlim[1])
        for tick in ax.get_xticklabels():
            tick.set_rotation(45)
        ax.set_xticklabels(["0","0.0001","0.001","0.01","0.1", "1", "10","100","1000"])
        rmse = sklearn.metrics.mean_squared_error(ydata, y_pred_rmse, sample_weight=None, multioutput='uniform_average', squared=True)
        si = (np.mean(ydata)/rmse)
        if label is True: #depending on the nature of the graph, we will want a label that tells the IC50 or scatter index
            ax.text(0.4, 0.1, "IC50 for" + " " + drug + " =" + str(round(res.x,4)), horizontalalignment='center', verticalalignment='center', transform = ax.transAxes)
            ax.text(0.4, 0.06, "SI for" + " " + drug + " =" + str(round(si,6)), horizontalalignment='center', verticalalignment='center', transform = ax.transAxes)
        if label is False:
            pass
        if write is True:
             #write basically writes down the IC50/Scatter index into a seperate CSV.
            with open('IC50_and_Scatter_Index.csv', mode='a') as ic50_SI:
                ic50_writer = csv.writer(ic50_SI, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
                ic50_writer.writerow([drug, str(round(res.x,4)), str(round(si,6)), str(b), str(c), str(d), str(e)])
                print(d)

        
        if write is False:
            pass
        if ax is None:
            return fig 
        else:
            return ax

if __name__ == '__main__':
    #definition for one of the parsers#
    
    # set up parser
    parser = argparse.ArgumentParser("This is Ari's helper string")
        # add arguments/options
    parser.add_argument("--input_file", 
            #required=True,
            type=str,
            help="This is the raw data file.")
    parser.add_argument("--path", 
            #required=True,
            type=str,
            help="This is the raw data folder.")
    parser.add_argument('--input_file_combination',
            type=str,
            help="This is the raw data file.")
    parser.add_argument('--combination_drugs', 
            type=str, nargs='+', action='append', 
            help="Type which drugs you want in a single graph in the format drug_a drug_b --combination_drugs  drug_c drug_d // WHITESPACE BETWEEN DRUG NAMES NESSECARY //. Use '--combination_drugs' to seperate graphs.")
    parser.add_argument("--plate_map", 
            default="./plate_map.csv", 
            type=str,
            help="This is the plate map template file.")
    parser.add_argument('-e', '--extension', default='', help='File extension to filter by.')
        # get provided user options
    args = parser.parse_args()
    folder = os.listdir(args.path)
    ###for graphs with specific drug combinations ####
    if args.combination_drugs and args.input_file_combination:
        Script1 = Script_1()
        Script1.format_script('plate_map_template.csv', args.input_file_combination)  
        for graph_list in args.combination_drugs:
            for element in graph_list:
                Script2 = Script_2()   
                Script2.find_parameters(element)

        Script3 = Script_3()
        drug_dict_para = {}
        drug_dict_x = {}
        drug_dict_y = {}
        drug_dict_std = {}
            #fig, ax = plt.subplots(1,1)
        err = np.array([0.1,0.1,0.1,0.1,0.1,0.1])
        
        pp = matplotlib.backends.backend_pdf.PdfPages('graphs_combinations.pdf')
        for graph_list in args.combination_drugs:
            color_single_graph = ['#2BBA7B', '#7BB3F1', '#7A51CC', '#FFC300', '#E8732B', '#2BDFE8', '#990000', '#000000', '#FFCCFF', '#E0BBE4'
                '#957DAD', '#D291BC', '#FEC8D8', '#FFDFD3', '#5918D0']
            fig, ax = plt.subplots(1,1)
            for element in graph_list:
                color = random.choice(color_single_graph)
                Script3.format_parameters('output_%s.txt' % element, element)
                Script3.find_mean('reformatteddata_for_scriptII.csv', element)
                ax = Script3.plotDRC_new(drug_dict_x[element], drug_dict_y[element], 
                drug_dict_std[element], drug_dict_para[element][0], drug_dict_para[element][1], 
                drug_dict_para[element][2], drug_dict_para[element][3], shade=True, labelpos=[0.0001,4], 
                ax = ax, drug = element, color = color, write = False)
                color_single_graph.remove(color)
            pp.savefig()      
        pp.close() 
        for element in sorted(Script2.drug_list_def()): 
            os.remove('output_%s.txt' % element)    
        os.remove('parameters.txt')
        os.remove('reformatteddata_for_scriptII.csv')


    ###if you only have one input file of raw data###  
    if args.input_file:
        Script1 = Script_1()
        Script1.format_script('plate_map_template.csv', args.input_file)
        Script2 = Script_2()    
        for element in Script2.drug_list_def():
            Script2.find_parameters(element)

        Script3 = Script_3()
        drug_dict_para = {}
        drug_dict_x = {}
        drug_dict_y = {}
        drug_dict_std = {}
            #fig, ax = plt.subplots(1,1)
        err = np.array([0.1,0.1,0.1,0.1,0.1,0.1])
        color_single_graph = ['#2BBA7B', '#7BB3F1', '#7A51CC', '#FFC300', '#E8732B', '#2BDFE8', '#990000', '#000000', '#FFCCFF', '#E0BBE4'
                '#957DAD', '#D291BC', '#FEC8D8', '#FFDFD3', '#5918D0']
        color_multi_graph = ['#2BBA7B', '#7BB3F1', '#7A51CC', '#FFC300', '#E8732B', '#2BDFE8', '#990000', '#000000']
        #with matplotlib.backends.backend_pdf.PdfPages('multipage.pdf') as pages:
        pp = matplotlib.backends.backend_pdf.PdfPages('graphs_input_file.pdf')
        for element in sorted(Script2.drug_list_def()): 
            color = random.choice(color_single_graph)
            Script3.format_parameters('output_%s.txt' % element, element)
            Script3.find_mean('reformatteddata_for_scriptII.csv', element)
            print(drug_dict_para)
            fig = Script3.plotDRC_new(drug_dict_x[element], drug_dict_y[element], 
                drug_dict_std[element], drug_dict_para[element][0], drug_dict_para[element][1], 
                drug_dict_para[element][2], drug_dict_para[element][3], shade=True, labelpos=[0.0001,4], 
                ax = None, drug = element, color = color, write = True, label = True)
            color_single_graph.remove(color)
            pp.savefig()
            
        #pp.close()
                #canvas = FigureCanvasPdf(fig)
                #canvas.print_figure(pages)
                #matplotlib.pyplot.close()
            
        ###multiple figures on one plot###
        fig, ax = plt.subplots(1,1)
        for element in sorted(Script2.drug_list_def()): 
            color = random.choice(color_multi_graph)
            Script3.format_parameters('output_%s.txt' % element, element)
            Script3.find_mean('reformatteddata_for_scriptII.csv', element)
            ax = Script3.plotDRC_new(drug_dict_x[element], drug_dict_y[element], 
                drug_dict_std[element], drug_dict_para[element][0], drug_dict_para[element][1], 
                drug_dict_para[element][2], drug_dict_para[element][3], shade=True, labelpos=[0.0001,4], 
                ax = ax, drug = element, color = color, write = False)
            color_multi_graph.remove(color)
        pp.savefig()
        pp.close()
        #plt.savefig('multigraphs.pdf')
                #canvas = FigureCanvasPdf(fig)
                #canvas.print_figure(pages)
                #matplotlib.pyplot.close()
        for element in sorted(Script2.drug_list_def()): 
            os.remove('output_%s.txt' % element)
        os.remove('parameters.txt')
        os.remove('reformatteddata_for_scriptII.csv')

    ###if you have a path/folder containing multiple CSVs of raw data####    
    if args.path:
        for name in folder:
            if name.startswith('.'):
                continue
            name = os.path.join(args.path, name)
            Script1 = Script_1()
            Script1.format_script('plate_map_template.csv', name)
            Script2 = Script_2()    
            for element in Script2.drug_list_def():
                Script2.find_parameters(element)

            Script3 = Script_3()
            drug_dict_para = {}
            drug_dict_x = {}
            drug_dict_y = {}
            drug_dict_std = {}
            #fig, ax = plt.subplots(1,1)
            err = np.array([0.1,0.1,0.1,0.1,0.1,0.1])
 
        #with matplotlib.backends.backend_pdf.PdfPages('multipage.pdf') as pages:
            pp = matplotlib.backends.backend_pdf.PdfPages('graphs_%s.pdf' % os.path.basename(name).split('.')[0])
            with open('IC50_and_Scatter_Index.csv', mode='a') as ic50_SI:
                ic50_writer = csv.writer(ic50_SI, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
                ic50_writer.writerow([name])
            color_single_graph = ['#2BBA7B', '#7BB3F1', '#7A51CC', '#FFC300', '#E8732B', '#2BDFE8', '#990000', '#000000', '#FFCCFF', '#E0BBE4'
                '#957DAD', '#D291BC', '#FEC8D8', '#FFDFD3', '#5918D0']
            color_multi_graph = ['#2BBA7B', '#7BB3F1', '#7A51CC', '#FFC300', '#E8732B', '#2BDFE8', '#990000', '#000000']
            for element in sorted(Script2.drug_list_def()): 
                color = random.choice(color_single_graph)
                Script3.format_parameters('output_%s.txt' % element, element)
                Script3.find_mean('reformatteddata_for_scriptII.csv', element)
                fig = Script3.plotDRC_new(drug_dict_x[element], drug_dict_y[element], 
                    drug_dict_std[element], drug_dict_para[element][0], drug_dict_para[element][1], 
                    drug_dict_para[element][2], drug_dict_para[element][3], shade=True, labelpos=[0.0001,4], 
                    ax = None, drug = element, color = random.choice(colors, 1, replace=False), write = True, label = True)
                color_single_graph.remove(color)
                pp.savefig()
            #pp.close()
                    #canvas = FigureCanvasPdf(fig)
                    #canvas.print_figure(pages)
                    #matplotlib.pyplot.close()
                
            ###multiple figures on one plot###
            fig, ax = plt.subplots(1,1)
            for element in sorted(Script2.drug_list_def()): 
                color = random.choice(color_multi_graph)
                Script3.format_parameters('output_%s.txt' % element, element)
                Script3.find_mean('reformatteddata_for_scriptII.csv', element)
                ax = Script3.plotDRC_new(drug_dict_x[element], drug_dict_y[element], 
                    drug_dict_std[element], drug_dict_para[element][0], drug_dict_para[element][1], 
                    drug_dict_para[element][2], drug_dict_para[element][3], shade=True, labelpos=[0.0001,4], 
                    ax = ax, drug = element, color = x, write= False)
                color_multi_graph.remove(color)
            pp.savefig()
            pp.close()
            #plt.savefig('multigraph_%s.pdf' % os.path.basename(name).split('.')[0])
                    #canvas = FigureCanvasPdf(fig)
                    #canvas.print_figure(pages)
                    #matplotlib.pyplot.close()
            Script2 = Script_2()  
            for element in sorted(Script2.drug_list_def()): 
                os.remove('output_%s.txt' % element)   
            os.remove('reformatteddata_for_scriptII.csv')
            os.remove('parameters.txt')
        #removing the .txt files to clear up unnessecary, inbetween files. 
        

           
