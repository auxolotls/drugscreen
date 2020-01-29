#!/usr/bin/env bash
#make sure you’re in the right directory containing the Script_Externship_2020#

##to input a single raw data file, and output graphs of one drug per graph##

python combined_scripts.py --input_file sample_data.csv

##to input a single raw data file, and output combinations of drugs on graphs##

python combined_scripts.py --combination_drugs drug_1 drug_2 drug_3 --combination_drugs drug_6 drug_5 --input_file_combination sample_data.csv

##to input a path, and output graphs of one drug per graph##

Python combined_scripts.py --path Sample_Folder_for_ArgParsing



