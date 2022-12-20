%% check individual people

clear all; close all; clc
%% 1. input data
cd("C:\Users\u0117545\Documents\GitHub\ULIFT_BC")
addpath("C:\Users\u0117545\OneDrive - KU Leuven\2.Dataprocessing\Matlab\addons")

Timepoint   = 'T1';
movement    = "F";
path.root   = 'C:\Users\u0117545\KU Leuven\An De Groef - DATA';
path.out    = fullfile(path.root,'Output','Database_MovQual.mat');
fs          = 60;
plot_or_not = 0;

Affected_table = readtable(fullfile(path.root,"Aangedane zijde.xlsx"));

subj_name   = ['BC_00' num2str(subj)];
affected = Affected_table(strcmp(Affected_table.ppID, subj_name), "involved");

load 