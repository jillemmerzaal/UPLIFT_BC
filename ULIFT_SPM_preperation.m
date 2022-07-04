%% ULIFT SPM preperation
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
    %% LEFT arm
    %----------
    % phase 1
    temp = mean(Data.(subjects{subj}).ULIFT.(timepoint).IK.L.Phase1.normalised.Glenohumeraal_flexion, 2);

    nReps = size(Data.(subjects{subj}).ULIFT.(timepoint).IK.L.Phase1.normalised.Glenohumeraal_flexion,2);
    for reps = 1:nReps
        R = corrcoef(Data.(subjects{subj}).ULIFT.(timepoint).IK.L.Phase1.normalised.Glenohumeraal_flexion(:,reps), temp);
        BC.IK.L.general(reps).corrcoef_ph1 = R(1,2)';
        rounded(reps) = round(R(1,2)*100)/100;
        if rounded(reps) < 0.8
            for ang = 1:nangles
                Data.(subjects{subj}).ULIFT.(timepoint).IK.L.Phase1.normalised.(IK_angles{ang})(:,reps) = nan;
            end
        end
    end
    clear temp rounded

    % phase 4
    temp = mean(Data.(subjects{subj}).ULIFT.(timepoint).IK.L.Phase4.normalised.Glenohumeraal_flexion, 2);

    nReps = size(Data.(subjects{subj}).ULIFT.(timepoint).IK.L.Phase4.normalised.Glenohumeraal_flexion,2);
    for reps = 1:nReps
        R = corrcoef(Data.(subjects{subj}).ULIFT.(timepoint).IK.L.Phase4.normalised.Glenohumeraal_flexion(:,reps), temp);
        BC.IK.L.general(reps).corrcoef_ph4 = R(1,2)';
        rounded(reps) = round(R(1,2)*100)/100;
        if rounded(reps) < 0.8
            for ang = 1:nangles
                Data.(subjects{subj}).ULIFT.(timepoint).IK.L.Phase4.normalised.(IK_angles{ang})(:,reps) = nan;
            end
        end
    end
    clear temp rounded


    %% RIGHT
    %-------
    % phase 1
    temp = mean(Data.(subjects{subj}).ULIFT.(timepoint).IK.R.Phase1.normalised.Glenohumeraal_flexion, 2);

    nReps = size(Data.(subjects{subj}).ULIFT.(timepoint).IK.R.Phase1.normalised.Glenohumeraal_flexion,2);
    for reps = 1:nReps
        R = corrcoef(Data.(subjects{subj}).ULIFT.(timepoint).IK.R.Phase1.normalised.Glenohumeraal_flexion(:,reps), temp);
        BC.IK.R.general(reps).corrcoef_ph1 = R(1,2)';
        rounded(reps) = round(R(1,2)*100)/100;
        if rounded(reps) < 0.8
            for ang = 1:nangles
                Data.(subjects{subj}).ULIFT.(timepoint).IK.R.Phase1.normalised.(IK_angles{ang})(:,reps) = nan;
            end
        end
    end
    clear temp rounded

    % phase 4
    temp = mean(Data.(subjects{subj}).ULIFT.(timepoint).IK.R.Phase4.normalised.Glenohumeraal_flexion, 2);

    nReps = size(Data.(subjects{subj}).ULIFT.(timepoint).IK.R.Phase4.normalised.Glenohumeraal_flexion,2);
    for reps = 1:nReps
        R = corrcoef(Data.(subjects{subj}).ULIFT.(timepoint).IK.R.Phase4.normalised.Glenohumeraal_flexion(:,reps), temp);
        BC.IK.R.general(reps).corrcoef_ph4 = R(1,2)';
        rounded(reps) = round(R(1,2)*100)/100;
        if rounded(reps) < 0.8
            for ang = 1:nangles
                Data.(subjects{subj}).ULIFT.(timepoint).IK.R.Phase4.normalised.(IK_angles{ang})(:,reps) = nan;
            end
        end
    end
    clear temp rounded

    %% create averages
    for ang = 1:nangles
        %LEFT
        BC.IK.L_phase1.(IK_angles{ang})(:,subj) = nanmean(Data.(subjects{subj}).ULIFT.(timepoint).IK.L.Phase1.normalised.(IK_angles{ang}), 2);
        BC.IK.L_phase4.(IK_angles{ang})(:,subj) = nanmean(Data.(subjects{subj}).ULIFT.(timepoint).IK.L.Phase4.normalised.(IK_angles{ang}), 2);

        %RIGHT
        BC.IK.R_phase1.(IK_angles{ang})(:,subj) = nanmean(Data.(subjects{subj}).ULIFT.(timepoint).IK.R.Phase1.normalised.(IK_angles{ang}), 2);
        BC.IK.R_phase4.(IK_angles{ang})(:,subj) = nanmean(Data.(subjects{subj}).ULIFT.(timepoint).IK.R.Phase4.normalised.(IK_angles{ang}), 2);


    end
end

% %% Combine left and right of GC subjects
% for ang = 1:nangles
%     BC.IK.B.(IK_angles{ang})   = cat(2,BC.IK.R.(IK_angles{ang}),BC.IK.L.(IK_angles{ang}));
% end

%% Kinematics

jointsOfInterst ={'Scapula_abbuction'      
    'Scapula_rotation'       
    'Scapula_flexion'        
    'Glenohumeraal_abbuction'
    'Glenohumeraal_rotation' 
    'Glenohumeraal_flexion'         
    'Elbow_flexion'          
    'Trunk_lateroflexion'    
    'Trunk_rotation'         
    'Trunk_flexion'      };    

nangles = size(jointsOfInterst,1);


% phase 1
fig1 = figure;
set(fig1,'Position',[1 31.4000 1536 758.4])

color_BC_l = [0 1 0];
color_BC_r = [0 0 0];
x = [1:101];

p = numSubplots(nangles);

for ang = 1:nangles
    subplot(p(1),p(2),ang)
    %PREOP
    %---------------------------------
    Av_bc_l = nanmean(BC.IK.L_phase1.(jointsOfInterst{ang}),2);
    Sd_bc_l = nanstd(BC.IK.L_phase1.(jointsOfInterst{ang}),0,2);

    Av_bc_r = nanmean(BC.IK.R_phase1.(jointsOfInterst{ang}),2);
    Sd_bc_r = nanstd(BC.IK.R_phase1.(jointsOfInterst{ang}),0,2);

    L_BC = shadedErrorBar(x,Av_bc_l,Sd_bc_l,{'Color',color_BC_l},[0.5]);
    hold on
    R_BC = shadedErrorBar(x,Av_bc_r, Sd_bc_r, {'Color', color_BC_r}, [0.5]);

    set(L_BC.mainLine,'LineWidth',1.5)
    set(R_BC.mainLine, 'LineWidth', 1.5)
    set(gca,'XLim',[1 100])

    title(jointsOfInterst{ang},'FontSize',13);
    xlabel('% of ULIFT cycle','FontSize',13);
end


% phase 4
fig2 = figure;
set(fig2,'Position',[1 31.4000 1536 758.4])

color_BC_l = [0 1 0];
color_BC_r = [0 0 0];
x = [1:101];

p = numSubplots(nangles);

for ang = 1:nangles
    subplot(p(1),p(2),ang)
    %PREOP
    %---------------------------------
    Av_bc_l = nanmean(BC.IK.L_phase4.(jointsOfInterst{ang}),2);
    Sd_bc_l = nanstd(BC.IK.L_phase4.(jointsOfInterst{ang}),0,2);

    Av_bc_r = nanmean(BC.IK.R_phase4.(jointsOfInterst{ang}),2);
    Sd_bc_r = nanstd(BC.IK.R_phase4.(jointsOfInterst{ang}),0,2);

    L_BC = shadedErrorBar(x,Av_bc_l,Sd_bc_l,{'Color',color_BC_l},[0.5]);
    hold on
    R_BC = shadedErrorBar(x,Av_bc_r, Sd_bc_r, {'Color', color_BC_r}, [0.5]);

    set(L_BC.mainLine,'LineWidth',1.5)
    set(R_BC.mainLine, 'LineWidth', 1.5)
    set(gca,'XLim',[1 100])

    title(jointsOfInterst{ang},'FontSize',13);
    xlabel('% of ULIFT cycle','FontSize',13);
end


