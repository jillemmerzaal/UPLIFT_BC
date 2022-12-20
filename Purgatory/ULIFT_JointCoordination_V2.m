%% ULIFT joint coordination V2

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

for subj = [1:3, 5:10]
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
    colnr = str2double(subj_name(d+1:end));

    % pre-op
    proximal_angle.Pre_op(:,colnr) = mean(Data.(subj_name).T0.(Phase).(involved).(proximal), 2);
    distal_angle.Pre_op(:,colnr) = mean(Data.(subj_name).T0.(Phase).(involved).(distal),2);

    %post-op
    proximal_angle.Post_op(:,colnr) = mean(Data.(subj_name).T1.(Phase).(involved).(proximal), 2);
    distal_angle.Post_op(:,colnr) = mean(Data.(subj_name).T1.(Phase).(involved).(distal),2);

end% end of subjects

%% 
figure; 
temp = mean(proximal_angle.Pre_op(:, [1:3 5:10]), 2)./mean(distal_angle.Pre_op(:, [1:3 5:10]), 2)
plot(temp)
hold on 
plot(mean(proximal_angle.Post_op(:, [1:3 5:10]), 2), mean(distal_angle.Post_op(:, [1:3 5:10]), 2), 'r')
legend('Pre-op', 'Post-op')
xlabel(proximal)
ylabel(distal)

