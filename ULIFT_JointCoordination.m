%% ULIFT joint coordination V1

clear all; clc; close all
%% 1. load data
cd("C:\Users\u0117545\Documents\GitHub\ULIFT_BC")
addpath("C:\Users\u0117545\OneDrive - KU Leuven\2.Dataprocessing\Matlab\addons")
addpath("C:\Users\u0117545\Documents\GitHub\ULIFT_BC")
addpath('C:\Users\u0117545\Documents\GitHub\biomechZoo-master\Toolbox\inter_joint_coordination')
path.root   = 'C:\Users\u0117545\KU Leuven\An De Groef - DATA';
path.in     = fullfile(path.root, 'Output', 'Database_ULIFT.mat');
path.out    = fullfile(path.root,'Output','Database_SPM.mat');

if exist(path.in,'file')
    load(path.in)
end

affected_table = readtable(fullfile(path.root, "Aangedane zijde.xlsx"));
Phase = 'phase1';
%% data setup

proximal    = 'Scapula_abbuction';
distal      = 'Glenohumeraal_flexion';

subj = 6;
if subj < 10
    subj_name   = ['BC_00' num2str(subj)];
elseif subj < 100
    subj_name   = ['BC_0' num2str(subj)];
else
    subj_name   = ['BC_', num2str(subj)];
end
%%
% find the affected side
involved = string(affected_table(find(strcmp(affected_table.ppID, subj_name)),:).involved);
% setup row number for exporting the data
d = strfind(subj_name,'_');
rownr = str2double(subj_name(d+1:end));

% pre-op
proximal_angle(:,1) = mean(Data.(subj_name).T0.(Phase).(involved).(proximal), 2);
distal_angle(:,1) = mean(Data.(subj_name).T0.(Phase).(involved).(distal),2);

%post-op
proximal_angle(:,2) = mean(Data.(subj_name).T1.(Phase).(involved).(proximal), 2);
distal_angle(:,2) = mean(Data.(subj_name).T1.(Phase).(involved).(distal),2);
%% Hilbert Transform per joint angle
% Pre-op
[PA_data_prox(1,:)] = Hilbert_PA(proximal_angle(:,1));  
[PA_data_dist(1,:)] = Hilbert_PA(distal_angle(:,1)); 

% Post-op 
[PA_data_prox(2,:)] = Hilbert_PA(proximal_angle(:,2));  
[PA_data_dist(2,:)] = Hilbert_PA(distal_angle(:,2)); 

%% Continuous Relative Phase 
CRP_data(1,:) = CRP(PA_data_dist(1,:),PA_data_prox(1,:));
CRP_data(2,:) = CRP(PA_data_dist(2,:),PA_data_prox(2,:));

%% 
%[MARP,DP] = CRP_stats(CRP_data);

%% figures  
%figure;

subplot(2,3,1); 
plot(proximal_angle)
xlim([0, 100])
title('(a) Proximal angle' )

subplot(2,3,2); 
plot(PA_data_prox', '.')
ylim([-180, 180])
xlim([0, 100])
title('(c) Proximal Phase Angles')

subplot(2,3,4); 
plot(distal_angle)
xlim([0, 100])
title('(b) Distal Joint Angle')

subplot(2,3,5); 
plot(PA_data_dist', '.')
ylim([-180, 180])
xlim([0, 100])
title('(d) Distal Phase Angle')

subplot(2,3,[3,6]); 
plot(CRP_data')
xlim([0, 100])
title('(e) Proximal-Distal CRP Angle')

