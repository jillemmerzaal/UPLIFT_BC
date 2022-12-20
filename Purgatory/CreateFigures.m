%% plot joint angle data

clear all; clc; close all
%% 1. load data
cd("C:\Users\u0117545\Documents\GitHub\ULIFT_BC")
addpath("C:\Users\u0117545\OneDrive - KU Leuven\2.Dataprocessing\Matlab\addons")
addpath("C:\Users\u0117545\Documents\GitHub\ULIFT_BC")
path.root   = 'C:\Users\u0117545\KU Leuven\An De Groef - DATA';
path.in     = fullfile(path.root, 'Output', 'Database_ULIFT.mat');
path.out    = fullfile(path.root,'Output','Database_SPM.mat');
Timpoints = {'T0', 'T1'};

if exist(path.in,'file')
    load(path.in)
end

affected_table = readtable(fullfile(path.root, "Aangedane zijde.xlsx"));


Phase       = 'phase1';
%% setup data for figure
% IK_angles = fieldnames(Data.BC_001.((Timpoints{t})).(Phase).L);

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
ntime = size(Timpoints,2);


for t = 1%:ntime
    %% phase 1
    Subj = fieldnames(Data);
    nSubj = size(Subj,1);

    fig1 = figure;
    set(fig1,'Position',[1 31.4000 1536 758.4])

    color_BC_aff = [0 1 0];
    color_BC_unaff = [0 0 0];
    x = [1:101];
    p = numSubplots(nSubj);

    try
        for s = [1:13, 15:18, 20:nSubj] % remove subj 14 and 20
            subplot(p(1), p(2), s)
            plot(Data.(Subj{s}).(Timpoints{t}).(Phase).L.(jointsOfInterst{6}));
            hold on
            plot(mean(Data.(Subj{s}).(Timpoints{t}).(Phase).L.(jointsOfInterst{6}),2 ), 'LineWidth',2)

            plot(Data.(Subj{s}).((Timpoints{t})).(Phase).R.(jointsOfInterst{6}), '--');
            plot(mean(Data.(Subj{s}).(Timpoints{t}).(Phase).R.(jointsOfInterst{6}),2 ), '--', 'LineWidth',2)

            set(gca,'XLim',[1 100])

            title(Subj{s},'FontSize',13);
            xlabel('% of ULIFT cycle','FontSize',13);
        end

        sgtitle(['ULIFT ' Phase, ' Left (solid) and Right(dashed) at ' (Timpoints{t})])
    catch ME
    end


    % averages per subject

    %try
        for s = [1:13, 15:18, 20:nSubj]
            involved = affected_table(find(strcmp(affected_table.ppID, Subj{s})),:).involved;
            for ang = 1:nangles
                %% if LEFT is affected
                if strcmp(involved, 'L')
                    BC.aff.(Phase).((Timpoints{t})).(jointsOfInterst{ang})(:,s) = mean(Data.(Subj{s}).((Timpoints{t})).(Phase).L.(jointsOfInterst{ang}),2);   
                    if strcmp(Phase, 'phase1')
                        if strcmp(jointsOfInterst{ang}, 'Trunk_lateroflexion') || strcmp(jointsOfInterst{ang}, 'Trunk_rotation')
                            BC.unaff.(Phase).((Timpoints{t})).(jointsOfInterst{ang})(:,s) = mean(-1 .* Data.(Subj{s}).((Timpoints{t})).(Phase).R.(jointsOfInterst{ang}),2);
                        else
                            BC.unaff.(Phase).((Timpoints{t})).(jointsOfInterst{ang})(:,s) = mean(Data.(Subj{s}).((Timpoints{t})).(Phase).R.(jointsOfInterst{ang}),2);
                        end
                    elseif strcmp(Phase, 'phase4')
                        BC.unaff.(Phase).((Timpoints{t})).(jointsOfInterst{ang})(:,s) = mean(Data.(Subj{s}).((Timpoints{t})).(Phase).R.(jointsOfInterst{ang}),2);

                    end
                    %% if RIGHT is affected
                elseif strcmp(involved, 'R')
                    BC.unaff.(Phase).((Timpoints{t})).(jointsOfInterst{ang})(:,s) = mean(Data.(Subj{s}).((Timpoints{t})).(Phase).L.(jointsOfInterst{ang}),2);

                    if strcmp(Phase, 'phase1')
                        if strcmp(jointsOfInterst{ang}, 'Trunk_lateroflexion') || strcmp(jointsOfInterst{ang}, 'Trunk_rotation')
                            BC.aff.(Phase).((Timpoints{t})).(jointsOfInterst{ang})(:,s) = mean(-1 .* Data.(Subj{s}).((Timpoints{t})).(Phase).R.(jointsOfInterst{ang}),2);
                        else
                            BC.aff.(Phase).((Timpoints{t})).(jointsOfInterst{ang})(:,s) = mean(Data.(Subj{s}).((Timpoints{t})).(Phase).R.(jointsOfInterst{ang}),2);
                        end
                    elseif strcmp(Phase, 'phase4')
                        BC.aff.(Phase).((Timpoints{t})).(jointsOfInterst{ang})(:,s) = mean(Data.(Subj{s}).((Timpoints{t})).(Phase).R.(jointsOfInterst{ang}),2);
                    end

                end
                


                %% save to table 
                 ppID(s,:)    = string(Subj{s});
                 time    = string(Timpoints{t});
                switch Phase      
                    case 'phase1'
                         phase1.(Timpoints{t}).(jointsOfInterst{ang})(s,:) = table(ppID, time, BC.aff.(Phase).((Timpoints{t})).(jointsOfInterst{ang})(:,s)')
%                         phase1.(Timpoints{t}).(jointsOfInterst{ang})(s,:) = table(ppID, time, BC.aff.(Phase).((Timpoints{t})).(jointsOfInterst{ang})(:,s)');
                    case 'phase4'
                        phase4.(Timpoints{t}).(jointsOfInterst{ang})(s,:) = table(ppID, time);

%                         phase4.(Timpoints{t}).(jointsOfInterst{ang})(s,:) = table(ppID, time, BC.aff.(Phase).((Timpoints{t})).(jointsOfInterst{ang})(:,s)');
                end
            end
        end
%     catch ME
%     end

    fig2 = figure;
    set(fig1,'Position',[1 31.4000 1536 758.4])

    color_BC_aff = [0 1 0];
    color_BC_unaff = [0 0 0];
    x = [1:101];

    p = numSubplots(nangles);

    for ang = 1:nangles
        subplot(p(1),p(2),ang)
        %PRE-OP
        %---------------------------------
        Av_bc_aff = nanmean(BC.aff.(Phase).((Timpoints{t})).(jointsOfInterst{ang}),2);
        Sd_bc_aff = nanstd(BC.aff.(Phase).((Timpoints{t})).(jointsOfInterst{ang}),0,2);

        Av_bc_unaff = nanmean(BC.unaff.(Phase).((Timpoints{t})).(jointsOfInterst{ang}),2);
        Sd_bc_unaff = nanstd(BC.unaff.(Phase).((Timpoints{t})).(jointsOfInterst{ang}),0,2);

        aff_BC = shadedErrorBar(x, Av_bc_aff,Sd_bc_aff,{'Color',color_BC_aff},[0.5]);
        hold on
        unaff_BC = shadedErrorBar(x, Av_bc_unaff, Sd_bc_unaff, {'Color', color_BC_unaff}, [0.5]);

        set(aff_BC.mainLine,'LineWidth',1.5)
        set(unaff_BC.mainLine, 'LineWidth', 1.5)
        set(gca,'XLim',[1 100])

        title(jointsOfInterst{ang},'FontSize',13);
        xlabel('% of ULIFT cycle','FontSize',13);

    end
    sgtitle(['ULIFT ' Phase ' Affected (green)/Unaffected(black) comparison ' (Timpoints{t})])
end
%% save to table

switch Phase 
    case 'phase1'
        for ang = 1:nangles
            temp.t0 = sortrows(phase1.T0.(jointsOfInterst{ang}),'ppID'); % sort the rows of t0 based on subject id
            temp.t1 = sortrows(phase1.T1.(jointsOfInterst{ang}),'ppID'); % sort the rows of t1 based on subject id
            temp.df = [temp.t0; temp.t1]; % concatinate the two tables into 1
            temp.df.Properties.VariableNames = {'ppID', 'Time', jointsOfInterst{ang}}; % set table header names
            writetable(temp.df, 'Phase1.xlsx', Sheet=jointsOfInterst{ang}) % write the table to an excel file
        end

    case 'phase4'
        for ang = 1:nangles
            temp.t0 = sortrows(phase4.T0.(jointsOfInterst{ang}),'ppID'); % sort the rows of t0 based on subject id
            temp.t1 = sortrows(phase4.T1.(jointsOfInterst{ang}),'ppID'); % sort the rows of t1 based on subject id
            temp.df = [temp.t0; temp.t1]; % concatinate the two tables into 1
            temp.df.Properties.VariableNames = {'ppID', 'Time', jointsOfInterst{ang}}; % set table header names
            writetable(temp.df, 'Phase4.xlsx', Sheet=jointsOfInterst{ang}) % write the table to an excel file
        end
end


%% purgatory 
% Phase1_Gleno_flex = BC.aff.phase1.T1.Glenohumeraal_flexion'


% Phase1_Gleno_flex = table(T, Phase1_Gleno_flex)

% writetable(Phase1_Gleno_flex, 'Phase1.xlsx', Sheet='Gleno_flex')


% for ang = 1:nangles
% 
%     Time1 = repmat((Timpoints{1}), size(BC.aff.phase1.(Timpoints{1}).(jointsOfInterst{ang})',1),1);
%     Time2 = repmat((Timpoints{2}), size(BC.aff.phase1.(Timpoints{2}).(jointsOfInterst{ang})',1),1);
%     T = [Time1;Time2];
%     temp_df = [BC.aff.phase1.(Timpoints{1}).(jointsOfInterst{ang})'; BC.aff.phase1.(Timpoints{2}).(jointsOfInterst{ang})'];
%     temp = table(T, temp_df);
% 
%     temp.Properties.VariableNames = {'Time', jointsOfInterst{ang}};
% 
%     writetable(temp, 'Phase1.xlsx', Sheet=jointsOfInterst{ang})
% 
%     clear temp Time1 Time2 T temp_df
% end