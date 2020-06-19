# TRNstability

This repository contains various codes for simulating the gene expression model and for exploring the stability of transcriptional regulatory networks.

### System Requirements/Installation:
These codes were ran with MATLAB R2019a. 

### Instructions for use:
In each of the folders, there is a file name ending with 'forcluster'. These are the main scripts that will reproduce and save the data for how the stability of the network varies with parameters such as N. (Within each of these files, there is a section 'Define parameters' in which one can specify the values of the different parameters of the model.) In each of the folders, the file 'Script_AnalyzeData.m' can then be used to analyze the saved data and output the corresponding figures. 

The folder 'fullyRandomNetworks' also contains a third Script file which can be used to run a single instance of a quenched network for different values of maximum fold-change (maxfc). Running this will produce figures for how the concentrations change over time for these different values of maxfc.

### Demo:
As an example, within the folder 'randomDAGs', running 'Script_eigvaluescalingwithN_withfc_forcluster.m' will produce a data output file similar to the data file that is currently in that folder. 

With the current data, running 'Script_AnalyzeData.m' in that folder will produce figures showing how the maximum eigenvalue of the Jacobian matrix vary with N and maxfc.`
