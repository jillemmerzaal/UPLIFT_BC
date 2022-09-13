%% plot joint angle data

% clear all; close all; clc
%% 1. load data
cd("C:\Users\u0117545\Documents\GitHub\ULIFT_BC")
addpath("C:\Users\u0117545\OneDrive - KU Leuven\2.Dataprocessing\Matlab\addons")
addpath("C:\Users\u0117545\Documents\GitHub\ULIFT_BC")
path.root   = 'C:\Users\u0117545\KU Leuven\An De Groef - DATA';
path.out    = fullfile(path.root,'Output','Database_ULIFT.mat');

if exist(path.out,'file')
    load(path.out)
end

Timepoint = 'T1';


%% setup data for figure
IK_angles = fieldnames(Data.BC_001.(Timepoint).phase1.L);

jointsOfInterst ={'Scapula_abbuction'      
    'Scapula_rotation'       
    'Scapula_flexion'        
    'Glenohumeraal_abbuction'
    'Glenohumeraal_rotation' 
    'Glenohumeraal_flexion'         
    'Elbow_flexion'          
    'Trunk_lateroflexion'    
    'Trunk_rotation'         
    'Trunk_flexion'};    

nangles = size(jointsOfInterst,1);

% phase 1 Left
Subj = fieldnames(Data);
nSubj = size(Subj,1);

fig1 = figure;
set(fig1,'Position',[1 31.4000 1536 758.4])

color_BC_L = [0 1 0];
color_BC_R = [0 0 0];
x = [1:101];
p = numSubplots(nSubj);

for s = 1:nSubj
    subplot(p(1), p(2), s)
    plot(Data.(Subj{s}).(Timepoint).phase1.L.(jointsOfInterst{6}));
    hold on
    plot(mean(Data.(Subj{s}).(Timepoint).phase1.L.(jointsOfInterst{6}),2 ), 'LineWidth',2)

    set(gca,'XLim',[1 100])

    title(Subj{s},'FontSize',13);
    xlabel('% of ULIFT cycle','FontSize',13);
end

sgtitle('ULIFT phase 1 Left arm') 



%% averages per subject
try
for s = 1:nSubj
    for ang = 1:nangles
        BC.L.phase1.(Timepoint).(jointsOfInterst{ang})(:,s) = mean(Data.(Subj{s}).(Timepoint).phase1.L.(jointsOfInterst{ang}),2);
        if strcmp(jointsOfInterst{ang}, 'Trunk_lateroflexion') || strcmp(jointsOfInterst{ang}, 'Trunk_rotation')
            BC.R.phase1.(Timepoint).(jointsOfInterst{ang})(:,s) = mean(-1 .* Data.(Subj{s}).(Timepoint).phase1.R.(jointsOfInterst{ang}),2);
        else
            BC.R.phase1.(Timepoint).(jointsOfInterst{ang})(:,s) = mean(Data.(Subj{s}).(Timepoint).phase1.R.(jointsOfInterst{ang}),2);
        end
    end
end
catch ME
end

% plot(mean(Data.(Subj{s}).(Timepoint).phase1.L.(jointsOfInterst{6}),2 ))

% phase 1
fig2 = figure;
set(fig1,'Position',[1 31.4000 1536 758.4])

color_BC_L = [0 1 0];
color_BC_R = [0 0 0];
x = [1:101];

p = numSubplots(nangles);

for ang = 1:nangles
    subplot(p(1),p(2),ang)
    %PRE-OP
    %---------------------------------
    Av_bc_L = nanmean(BC.L.phase1.(Timepoint).(jointsOfInterst{ang}),2);
    Sd_bc_L = nanstd(BC.R.phase1.(Timepoint).(jointsOfInterst{ang}),0,2);

    Av_bc_R = nanmean(BC.R.phase1.(Timepoint).(jointsOfInterst{ang}),2);
    Sd_bc_R = nanstd(BC.R.phase1.(Timepoint).(jointsOfInterst{ang}),0,2);

    L_BC = shadedErrorBar(x, Av_bc_L,Sd_bc_L,{'Color',color_BC_L},[0.5]);
    hold on
    R_BC = shadedErrorBar(x, Av_bc_R, Sd_bc_R, {'Color', color_BC_R}, [0.5]);

    set(L_BC.mainLine,'LineWidth',1.5)
    set(R_BC.mainLine, 'LineWidth', 1.5)
    set(gca,'XLim',[1 100])

    title(jointsOfInterst{ang},'FontSize',13);
    xlabel('% of ULIFT cycle','FontSize',13);
    
end
sgtitle('ULIFT phase 1 Left (green)/Right(black) comparison') 

