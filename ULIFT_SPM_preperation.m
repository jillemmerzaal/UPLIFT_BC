%% ULIFT SMP preperation
% this code will prepare the data to be used in the SPM analysis in
% python. It will:
    % 1) read the data from the preprocessing step
    % 2) structure the data per group (Healthy control vs breast cancer/
    % pain vs no pain)
    % 3) check the data per participant. i.e. the phases of the ULIFT task
    % are autmatically extracated based on change points, peak detection,
    % and zerro crossing, but this could be faulty for some trials. Thus
    % the correlation coefficient of each trail is compared to the average
    % for that participant. When the R is lower than 0.7 that particular
    % trial will be flagged and removed. 
    % 3) data will be saved and stored in a seperate struct to be used in
    % python. 
clear all; close all; clc
cd("C:\Users\u0117545\Documents\GitHub\ULIFT_BC")


path.root   = 'C:\Users\u0117545\KU Leuven\An De Groef - DATA';
load(fullfile(path.root,'Output','Database_ULIFT.mat'));



timepoint = 'T0';

subjects = fieldnames(Data);
nsubjects = size(subjects,1);

IK_angles = fieldnames(Data.(subjects{1}).ULIFT.(timepoint).IK.L.Phase1.normalised);
nangles = size(IK_angles,1);

for subj = 1:nsubjects
%% check left
    temp = mean(Data.(subjects{subj}).ULIFT.(timepoint).IK.L.Phase1.normalised.Glenohumeraal_flexion, 2);

    nReps = size(Data.(subjects{subj}).ULIFT.(timepoint).IK.L.Phase1.normalised.Glenohumeraal_flexion,2);
    for reps = 1:nReps
        R = corrcoef(Data.(subjects{subj}).ULIFT.(timepoint).IK.L.Phase1.normalised.Glenohumeraal_flexion(:,reps), temp);
        BC.IK.L.general(reps).corrcoef = R(1,2)';
        if BC.IK.L.general(reps).corrcoef < 0.8
            for ang = 1:nangles
                Data.(subjects{subj}).ULIFT.(timepoint).IK.L.Phase1.normalised.(IK_angles{ang})(:,reps) = nan;
            end
        end
    end
    clear temp
%% check right 
temp = mean(Data.(subjects{subj}).ULIFT.(timepoint).IK.R.Phase1.normalised.Glenohumeraal_flexion, 2);

    nReps = size(Data.(subjects{subj}).ULIFT.(timepoint).IK.R.Phase1.normalised.Glenohumeraal_flexion,2);
    for reps = 1:nReps
        R = corrcoef(Data.(subjects{subj}).ULIFT.(timepoint).IK.R.Phase1.normalised.Glenohumeraal_flexion(:,reps), temp);
        BC.IK.R.general(reps).corrcoef = R(1,2)';
        if BC.IK.R.general(reps).corrcoef < 0.7
            for ang = 1:nangles
                Data.(subjects{subj}).ULIFT.(timepoint).IK.R.Phase1.normalised.(IK_angles{ang})(:,reps) = nan;
            end
        end
    end
clear temp

%% create averages 
    for ang = 1:nangles
        BC.IK.L.(IK_angles{ang})(:,subj) = nanmean(Data.(subjects{subj}).ULIFT.(timepoint).IK.L.Phase1.normalised.(IK_angles{ang}), 2);
        BC.IK.R.(IK_angles{ang})(:,subj) = nanmean(Data.(subjects{subj}).ULIFT.(timepoint).IK.R.Phase1.normalised.(IK_angles{ang}), 2);
    end
end

% %% Combine left and right of GC subjects
% for ang = 1:nangles
%     BC.IK.B.(IK_angles{ang})   = cat(2,BC.IK.R.(IK_angles{ang}),BC.IK.L.(IK_angles{ang}));
% end

%% Kinematics
fig1 = figure;
set(fig1,'Position',[1 31.4000 1536 758.4])

color_BC = [0 1 0];
color_GC = [0 0 0];
x = [1:101];

p = numSubplots(nangles);
for ang = 1:nangles
    subplot(p(1),p(2),ang)
    %PREOP
    %---------------------------------
    Av_bc = nanmean(BC.IK.R.(IK_angles{ang}),2);
    SD_bc = nanstd(BC.IK.R.(IK_angles{ang}),0,2);
    
    L_BC = shadedErrorBar(x,Av_bc,SD_bc,{'Color',color_BC},[0.5]);


    set(L_BC.mainLine,'LineWidth',1.5)
    set(gca,'XLim',[1 100])

    title(IK_angles{ang},'FontSize',13);
    xlabel('% of ULIFT cycle','FontSize',13);
end




