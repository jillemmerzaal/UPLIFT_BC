%% Workflow_CNN
% this code determines the wear and non wear sections in the data ussing the
% Covolutional neural network developed by Syed et al., (2022).
% IMPORTANT NOTE: if you have first run the extract_gt3x3df script. Restart
% matlab to read in the new environment.

% funtions needed:
% 1. raw_non_wear_functions.py
% 2. numpy arrays extracted from Extract_gt3x2df.m

clear; clc; close all;

cd('C:\Users\u0117545\Documents\GitHub\ULIFT_BC\DailyLife\')
path.root   = 'C:\Users\u0117545\Documents\GitHub\ULIFT_BC\DailyLife\';
path.data   = 'C:\Users\u0117545\KU Leuven\An De Groef - DATA';
timepoint   = 'T1';
plot_or_not = 1;
fs          = 30;


%% set python environment
pe = pyenv("Version", "C:\GBW_MyPrograms\Anaconda3\envs\MATLAB_PYTHON\python.exe");

pathToFunc = fileparts(which('raw_non_wear_functions.py'));
if count(py.sys.path,pathToFunc) == 0
    insert(py.sys.path,int64(0),pathToFunc);
end

%%
for subj = [1 3]%6:10
    if subj < 10
        subj_name = ['BC_00', num2str(subj)];
    elseif subj < 100
        subj_name = ['BC_0', num2str(subj)];
    else
        subj_name = ['BC_', num2str(subj)];
    end


    disp(' ')
    disp(['Processing: ' subj_name])
    disp(['     at timepoint: ', timepoint])

    path.subj   = fullfile(path.data, subj_name, 'Accelerometrie', timepoint);
    check_subj  = exist(path.subj, "dir");
    content_check = isfile(fullfile(path.subj, 'rawdata_hip.parq'));

    if check_subj == 7 && content_check == 1
        clearvars -except path timepoint fs plot_or_not pe pathToFunc subj_name subj content nfiles file_counter file_name

        disp('     Reading in the hip data......')

        rawdata.hip_t = parquetread(fullfile(path.subj, 'rawdata_hip.parq'));
        rawdata.hip = table2array(rawdata.hip_t);
        hip_data = mat2np(rawdata.hip);

        %% obtain raw acceleration data non-wear-periods from the hip sensor
        disp('     Segment wear vs non-wear periods......')
        sample_rate = py.int(fs);
        pyOut.wear_vector_hees = py.raw_non_wear_functions.hees_2013_calculate_non_wear_time(hip_data, sample_rate);
        wear_vector_hees = np2mat(pyOut.wear_vector_hees);

        figure;
        plot(rawdata.hip);
        hold on
        plot(wear_vector_hees, 'LineWidth',1.5)

        start_indices = find([0; diff(wear_vector_hees) == 1]);
        end_indices = find([0; diff(wear_vector_hees) == -1]);

        if ~isempty(start_indices) && ~isempty(end_indices)
            if end_indices(1) - start_indices(1) < 0
                start_indices = [1; start_indices];
            end

            if end_indices(end) - start_indices(end) < 0
                end_indices = [end_indices; size(rawdata.hip,1)];
            end

            % second threshold
            min_non_wear_time_window = 135;

            second_threshold = (fs * 60) * min_non_wear_time_window;

            mask = find(end_indices - start_indices <= second_threshold);
            start_indices(mask) = [];
            end_indices(mask) = [];

            xline(start_indices)

            %% calculate wear time
            disp('     Calulate wear time......')

            wear_blocks(1,:) = end_indices - start_indices; % number of samples in a block
            wear_blocks(2,:) = wear_blocks(1,:) / fs; % number of seconds in a block
            wear_blocks(3,:) = wear_blocks(2,:) / 60; % number of minutes in a block
            wear_blocks(4,:) = wear_blocks(3,:) / 60; % number of hours in a block

            disp('     Done')
            disp('     __________________________')

            %if size(wear_blocks,2) >= 5 && sum(wear_blocks(4,:) >= 12) >= 5
            %% finish rest of the analysis!
            j = 0;
            for idx = 1:size(wear_blocks,2)
                if wear_blocks(4,idx) >= 12 && wear_blocks(4,idx) <= 24
                    j=j+1;
                    blocks(j,:) = [start_indices(idx), end_indices(idx)];
                end
            end

            if size(blocks,1) >= 5
                disp('     Load left arm data:......')
                rawdata.left_t = parquetread(fullfile(path.subj, 'rawdata_left.parq'));
                rawdata.left = table2array(rawdata.left_t);
                %left_data = mat2np(rawdata.left);


                disp('     Load right arm data:......')
                rawdata.right_t = parquetread(fullfile(path.subj, 'rawdata_right.parq'));
                rawdata.right = table2array(rawdata.right_t);
                %right_data = mat2np(rawdata.right);

                if ~isempty(rawdata.right) && ~isempty(rawdata.left)
                    disp(['     Wear time and data check completed......', newline, '     Conclusion: Data will be analysed'])
                    Conclusion = 'Analysed';
                    basic_table = table(string(subj_name), string(timepoint), wear_blocks(4,:), string(Conclusion));
                    basic_table.Properties.VariableNames = {'ppID', 'Timepoint', 'wear time', 'Conclusion'};
                    writetable(basic_table, 'Weartime.xlsx', Sheet=timepoint, WriteMode='append')
                    disp('     __________________________')

                    if size(rawdata.left,1) >= (size(rawdata.hip,1) / 100) * 85 && size(rawdata.right,1) >= (size(rawdata.hip,1) / 100) * 85
                        disp('     Shape check completed......')
                        for b = 1:size(blocks,1)
                            disp(['     Analyzing wear block ', num2str(b), ' of ', num2str(size(blocks,1))])
                            data.left{b,1} = rawdata.left(blocks(b,1):blocks(b,2), :);
                            data.left{b,2} = mean(rawdata.left(blocks(b,1):blocks(b,2), :));
                            data.right{b,1} = rawdata.right(blocks(b,1):blocks(b,2), :);
                            data.right{b,2} = mean(rawdata.right(blocks(b,1):blocks(b,2), :));
                            data.hip{b,1} = rawdata.hip(blocks(b,1):blocks(b,2), :);

                            %% correct sensor tilt and axis definition to match Lum et al.,

                            disp('     Redefine axis defenition......')
                            % Left
                            acc.x = data.left{b,1}(:,2);
                            acc.y = data.left{b,1}(:,1);
                            acc.z = data.left{b,1}(:,3) * -1;

                            data.L{b,1} = [acc.x, acc.y,  acc.z];

                            clear acc

                            % Right
                            acc.x = data.right{b,1}(:,2) * -1;
                            acc.y = data.right{b,1}(:,1) * -1;
                            acc.z = data.right{b,1}(:,3) * -1;

                            data.R{b,1} = [acc.x, acc.y,  acc.z];

                            clear acc

                            %% resample data from 30Hz (our sensors) to 50Hz (needed for the model)
                            disp('     Resample data from 30 to 50Hz......')
                            fs = 30;
                            fs_new = 50;

                            x = 1:length(data.L{b,1}); % old time axis of the data
                            xq = 1:fs/fs_new:length(data.L{b,1}); % new time axis for the data

                            data.L_sp{b,1} = interp1(x, data.L{b,1}, xq, 'spline');
                            data.R_sp{b,1} = interp1(x, data.R{b,1}, xq, 'spline');

                            %% prepare for model
                            disp('     Calulate features from data block......')

                            load('model.mat')
                            n=200;
                            feature_l = featurecalc1(data.L_sp{b,1},n);
                            feature_r = featurecalc1(data.R_sp{b,1},n);

                            yyfit_l = trainedModel.predictFcn(feature_l);
                            yyfit_r = trainedModel.predictFcn(feature_r);

                            %% calculate minutes functionally active
                            disp(['     Calculate the total functional activity for block ', num2str(b)])
                            data.wear_time(b,1) = ((size(data.L_sp{b,1},1) / 50) / 3600);
                            data.pred_L(b,1) = (size(yyfit_l(yyfit_l == 1), 1) * 4) / 60; %minuten
                            data.pred_L(b,2) = data.pred_L(b,1)/60; % hours

                            data.pred_R(b,1) = (size(yyfit_r(yyfit_r == 1), 1) * 4) / 60; % minutes
                            data.pred_R(b,2) = data.pred_R(b,1) / 60; % hours
                        end
                        %% Output table
                        disp('      Write data to output table.........')
                        no_days = size(blocks,1);
                        avg_weartime = mean(data.wear_time);
                        avg_act_l = mean(data.pred_L(:,2));
                        avg_perc_l = avg_act_l / avg_weartime;
                        avg_act_r = mean(data.pred_R(:,2));
                        avg_perc_r = avg_act_r / avg_weartime;

                        tbl = table(string(subj_name), string(timepoint),no_days, avg_weartime, avg_act_l, avg_perc_l, avg_act_r, avg_perc_r);
                        tbl.Properties.VariableNames = {'subject id', 'Timepoint', 'wear days', 'wear time', 'mean use L', 'percentage use L', 'mean use R', 'percentage use R'};
                        writetable(tbl, 'output.xlsx', 'WriteMode','append','AutoFitWidth',true, Sheet=timepoint)
                        %% figures
                        if plot_or_not
                            maximum_signal = max(max(abs(rawdata.hip)));
                            ymax = maximum_signal + ((maximum_signal/100) * 10);
                            ymin = ymax * -1;

                            figure;
                            hold on

                            x_values = [blocks(:,1) blocks(:,2) blocks(:,2) blocks(:,1)];
                            y_values = [ymin ymin ymax ymax];

                            fill(x_values, y_values, [0.78, 0.81, 0.84], 'EdgeColor','none');
                            plot(rawdata.hip);
                            hold off

                        end
                    else
                        disp(['     Shape check completed......', newline, '     Revision of Conclusion: Data will not be analysed'])

                        disp('      Write data to output table.........')
                        no_days = size(blocks,1);
                        avg_weartime = "shape check failed, check de data in actigraph software for errors";
                        avg_act_l = nan;
                        avg_perc_l = nan;
                        avg_act_r = nan;
                        avg_perc_r = nan;

                        tbl = table(string(subj_name), string(timepoint),no_days, avg_weartime, avg_act_l, avg_perc_l, avg_act_r, avg_perc_r);
                        tbl.Properties.VariableNames = {'subject id', 'Timepoint', 'wear days', 'wear time', 'mean use L', 'percentage use L', 'mean use R', 'percentage use R'};
                        writetable(tbl, 'output.xlsx', 'WriteMode','append','AutoFitWidth',true, Sheet=timepoint)
                    end% shape check completed

                else
                    disp(['     Wear time and data check completed......', newline, '     Conclusion: Data will not be analysed'])
                    Conclusion = 'Corrupt files';
                    basic_table = table(string(subj_name), string(timepoint), wear_blocks(4,:), string(Conclusion));
                    basic_table.Properties.VariableNames = {'ppID', 'Timepoint', 'wear time', 'Conclusion'};
                    writetable(basic_table, 'Weartime.xlsx', Sheet=timepoint, WriteMode='append')

                end %left and right data is not empty
            else
                disp(['     Wear time check completed......', newline, '     Conclusion: Data will not be analysed'])
                Conclusion = 'Not Analysed';
                basic_table = table(string(subj_name), string(timepoint), wear_blocks(4,:), string(Conclusion));
                basic_table.Properties.VariableNames = {'ppID', 'Timepoint', 'wear time', 'Conclusion'};
                writetable(basic_table, 'Weartime.xlsx', Sheet=timepoint, WriteMode='append')

                disp(['Processing: ' subj_name newline ' at timepoint: ' timepoint ' COMPLETED'])
                disp('====================================================')
            end% wear blocks are sufuicuently long
        else
            disp(['     Wear time check completed......', newline, '     Conclusion: Data will not be analysed'])
            Conclusion = "Not Analysed";

            basic_table = table(string(subj_name), string(timepoint), "no wear blocks", Conclusion);
            basic_table.Properties.VariableNames = {'ppID', 'Timepoint', 'wear blocks',  'Conclusion'};
            writetable(basic_table, 'Weartime.xlsx', Sheet=timepoint, WriteMode='append')

            disp(['Processing: ' subj_name newline ' at timepoint: ' timepoint ' COMPLETED'])
            disp('====================================================')
        end % start/end indices are not empty
    end % subject path exists
end % subject loop
