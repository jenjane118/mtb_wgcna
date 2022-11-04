#!/usr/bin python3
# coding=utf-8

#iterate_david.py

# wrapper script to get gene lists and submit to David web service
from termClusteringReport import run_termClusteringReport
from chartReport import get_chartReport

taxon_id = 83332


with open('module_list.txt', 'r') as f:
    mod_list = f.read().split('\n')
for module in mod_list:
    not_mods = [module !="", module !="grey"]
    if all(not_mods):
        mod_file = "gene_strings/" + module + ".txt"
        with open(mod_file, 'r') as genefile:
            genes = genefile.read()
            #run_termClusteringReport(module, genes)
            get_chartReport(module, genes)

print("DONE!")
