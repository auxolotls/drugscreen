This is a script that can give you graphs of drug dosage versus luminosity, along with the corresponding IC50s. You can plot multiple drugs on a single graph, or just a single drug on a graph, depending on your arg parse command.

Here’s what all the different files and folders do:
[SCROLL DOWN FOR INSTRUCTIONS ON USING ARGPARSE TO RUN THE SCRIPT] 

Sample_data: This is the sample raw data to be fed into the script. This data should be obtained from a plate reader that can print out luminosity readings.

Fitting_single_drug.R: This is the R script that will fit your data to a sigmoid curve. It gives out the four parameters needed to fit your curve. Your data may not fit using this four parameter method—if this is the case, the scatter index will be sufficiently high or low for you to tell; you will also be able to tell because your curve will not be properly fit. 

Combined_scripts.Py: Inside are the python scripts needed to graph your graphs. 

Plate_map_templates: This is a .csv file containing the plate map templates. ‘Combined_scripts.py’ will combine the data in this file and the data in ‘Sample_data.csv’ into one dataframe called ‘reformatteddata_for_scriptII.csv’. This csv file can be read by the R script because of its specific formatting. 

The first “plate” in the template is the drug labelling, in the sense that 1 = drug_1, 2 = drug_2…8 = drug_8.

The second “plate” are the replicates for each drug.

The third plate are the concentrations of each drug in each well.

#######
NOTE: You will have to edit this template if your well setup looks different from the one given in this template. 

Some things to note: Either label your drugs as simply the ‘#’ or ‘drug_#’. (Without the quotes). The script cannot read them otherwise. 

If you have empty wells, label them ‘NA’ in the template (without the quotes). The script will later remove those wells. 
#######

Sample_Folder_for_ArgParsing: You have the choice of either reading through a folder that contains multiple files of raw data, or just one data file. This folder is a sample path that you can use to demo how to plot multiple files of data. 

#######################
How to use argparse to run your python script. 

First, make sure you have all of these packages installed: 
Matplotlib
Statistics
Pandas
Numpy 
Scipy.Optimize
SKlearn.metrics
Rpy2
Randomcolor 
rDRC

TO INPUT A SINGLE FILE, PUT INTO THE TERMINAL:


python Combined_Scripts.py --input_file <raw data file> 

TO INPUT A PATH, PUT INTO THE TERMINAL:


python Combined_Scripts.py --input_path <folder path> 

TO INPUT A PATH, PUT INTO THE TERMINAL:


python Combined_Scripts.py --combination_drugs drug_1 drug_2 drug_3  --input_file <raw data file>

^^^ Make sure there is space between the drug names; it is not separated by commas.

If you have multiple graphs with multiple drug combinations, simply go: 

python Combined_Scripts.py --combination_drugs drug_1 drug_2 drug_3 --combination_drugs drug_6 drug_5 —input_file_combination <raw data file>

OUTPUTS
You will get a PDF/PDFs of your graph(s). 
You will also get a textfile containing your IC50s and Scatter Index value for each drug.
Some files may be created, and then deleted after the script is done running. That is part of the script; it is used to prevent clutter buildup. 


TODO:
Test!
