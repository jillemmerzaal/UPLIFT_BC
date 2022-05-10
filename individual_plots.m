%individual plots
clear all; close all; clc

subj_id = 'BCT_004';
Timepoint   = 'T0';
movement    = "ULIFT";
path.root   = 'C:\Users\u0117545\KU Leuven\An De Groef - 5. Data';
path.out    = fullfile(path.root,'Output','Database_ULIFT.mat');
if exist(path.out,'file')
    load(path.out)
end



df.L.Phase1 = Data.(subj_id).ULIFT.(Timepoint).IK.L.Phase1.normalised;
df.L.Phase4 = Data.(subj_id).ULIFT.(Timepoint).IK.L.Phase4.normalised;

df.R.Phase1 = Data.(subj_id).ULIFT.(Timepoint).IK.R.Phase1.normalised;
df.R.Phase4 = Data.(subj_id).ULIFT.(Timepoint).IK.R.Phase4.normalised;

JointAngles = fieldnames(df.R.Phase1);

%% check left phase 1
    temp = mean(df.L.Phase1.Glenohumeraal_flexion, 2);

    nReps = size(df.L.Phase1.Glenohumeraal_flexion,2);
    for reps = 1:nReps
        R = corrcoef(df.L.Phase1.Glenohumeraal_flexion(:,reps), temp);
        df.L.Phase1.general(reps).corrcoef = R(1,2)';
        if df.L.Phase1.general(reps).corrcoef < 0.7
            for ang = 1:length(JointAngles)
                df.L.Phase1.(JointAngles{ang})(:,reps) = nan;
            end
        end
    end
    clear temp

%% check right phase 1
    temp = mean(df.R.Phase1.Glenohumeraal_flexion, 2);

    nReps = size(df.R.Phase1.Glenohumeraal_flexion,2);
    for reps = 1:nReps
        R = corrcoef(df.R.Phase1.Glenohumeraal_flexion(:,reps), temp);
        df.R.Phase1.general(reps).corrcoef = R(1,2)';
        if df.R.Phase1.general(reps).corrcoef < 0.7
            for ang = 1:length(JointAngles)
                df.R.Phase1.(JointAngles{ang})(:,reps) = nan;
            end
        end
    end
    clear temp

    %% check left phase 4
    temp = mean(df.L.Phase4.Glenohumeraal_flexion, 2);

    nReps = size(df.L.Phase4.Glenohumeraal_flexion,2);
    for reps = 1:nReps
        R = corrcoef(df.L.Phase4.Glenohumeraal_flexion(:,reps), temp);
        df.L.Phase4.general(reps).corrcoef = R(1,2)';
        if df.L.Phase4.general(reps).corrcoef < 0.7
            for ang = 1:length(JointAngles)
                df.L.Phase4.(JointAngles{ang})(:,reps) = nan;
            end
        end
    end
    clear temp

    %% check right phase 4
    temp = mean(df.R.Phase4.Glenohumeraal_flexion, 2);

    nReps = size(df.R.Phase4.Glenohumeraal_flexion,2);
    for reps = 1:nReps
        R = corrcoef(df.R.Phase4.Glenohumeraal_flexion(:,reps), temp);
        df.R.Phase4.general(reps).corrcoef = R(1,2)';
        if df.R.Phase4.general(reps).corrcoef < 0.7
            for ang = 1:length(JointAngles)
                df.R.Phase4.(JointAngles{ang})(:,reps) = nan;
            end
        end
    end
    clear temp



figure
p = tiledlayout('flow', 'TileSpacing','compact');
title(p,['Phase 1 Left ', (Timepoint)])
xlabel(p,'% movement')
ylabel(p,['Joint angle', char(176)])

for jnt = 1:length(JointAngles)
    nexttile
    plot(df.L.Phase1.(JointAngles{jnt}))
    xlim([0, 100])

    [PlotTitle] = regexprep(JointAngles{jnt}, '_', ' ');
    title(PlotTitle)
end

figure
p = tiledlayout('flow', 'TileSpacing','compact');
title(p,['Phase 4 Left ', (Timepoint)])
xlabel(p,'% movement')
ylabel(p,['Joint angle', char(176)])

for jnt = 1:length(JointAngles)
    nexttile
    plot(df.L.Phase4.(JointAngles{jnt}))
    xlim([0, 100])

    [PlotTitle] = regexprep(JointAngles{jnt}, '_', ' ');
    title(PlotTitle)
end

figure
p = tiledlayout('flow', 'TileSpacing','compact');
title(p,['Phase 1 Right ', (Timepoint)])
xlabel(p,'% movement')
ylabel(p,['Joint angle', char(176)])

for jnt = 1:length(JointAngles)
    nexttile
    plot(df.R.Phase1.(JointAngles{jnt}))
    xlim([0, 100])

    [PlotTitle] = regexprep(JointAngles{jnt}, '_', ' ');
    title(PlotTitle)
end

figure
p = tiledlayout('flow', 'TileSpacing','compact');
title(p,['Phase 4 Right ', (Timepoint)])
xlabel(p,'% movement')
ylabel(p,['Joint angle', char(176)])

for jnt = 1:length(JointAngles)
    nexttile
    plot(df.R.Phase4.(JointAngles{jnt}))
    xlim([0, 100])

    [PlotTitle] = regexprep(JointAngles{jnt}, '_', ' ');
    title(PlotTitle)
end


%% Averages
figure
p = tiledlayout('flow', 'TileSpacing','compact');
title(p,['Phase 1 Average left and right ', (Timepoint)])
xlabel(p,'% movement')
ylabel(p,['Joint angle', char(176)])

for jnt = 1:length(JointAngles)
    df_avg.L.Phase1.(JointAngles{jnt}) = mean(Data.(subj_id).ULIFT.(Timepoint).IK.L.Phase1.normalised.(JointAngles{jnt}), 2);
    df_avg.R.Phase1.(JointAngles{jnt}) = mean(Data.(subj_id).ULIFT.(Timepoint).IK.R.Phase1.normalised.(JointAngles{jnt}), 2);

    nexttile
    plot(df_avg.L.Phase1.(JointAngles{jnt}), 'DisplayName', "Avg Left")
    hold on
    plot(df_avg.R.Phase1.(JointAngles{jnt}), 'DisplayName', "Avg Right")
    xlim([0, 100])

    [PlotTitle] = regexprep(JointAngles{jnt}, '_', ' ');
    title(PlotTitle)
end



figure
p = tiledlayout('flow', 'TileSpacing','compact');
title(p,['Phase 4 Average left and right ', (Timepoint)])
xlabel(p,'% movement')
ylabel(p,['Joint angle', char(176)])

for jnt = 1:length(JointAngles)
    df_avg.L.Phase4.(JointAngles{jnt}) = mean(Data.(subj_id).ULIFT.(Timepoint).IK.L.Phase4.normalised.(JointAngles{jnt}), 2);
    df_avg.R.Phase4.(JointAngles{jnt}) = mean(Data.(subj_id).ULIFT.(Timepoint).IK.R.Phase4.normalised.(JointAngles{jnt}), 2);

    nexttile
    plot(df_avg.L.Phase4.(JointAngles{jnt}), 'DisplayName', "Avg Left")
    hold on
    plot(df_avg.R.Phase4.(JointAngles{jnt}), 'DisplayName', "Avg Right")
    xlim([0, 100])

    [PlotTitle] = regexprep(JointAngles{jnt}, '_', ' ');
    title(PlotTitle)
end




