%% Extract_gt3x
% Extract the information from the .gt3x files. and saves the numpy arrays
% to file in the subject folder.
 
% functions needed: 
% 1. gt3x2df.py this extract the raw acceleration data from the nativ
% actigraph files
% 2. gt3x_functions.py this extracts the meta information from the sensors
% to determine the accelerometer site.

clear; clc; close all;

cd('C:\Users\u0117545\Documents\GitHub\ULIFT_BC\DailyLife\')
path.root 		= 'C:\Users\u0117545\Documents\GitHub\ULIFT_BC\DailyLife\';
path.data       = 'C:\Users\u0117545\KU Leuven\An De Groef - DATA';
timepoint       = 'T1';
plot_or_not     = 1;


%% set python environment
pe = pyenv("Version", "C:\GBW_MyPrograms\Anaconda3\python.exe");

pathToFunc = fileparts(which('gt3x2df.py'));
if count(py.sys.path,pathToFunc) == 0
    insert(py.sys.path,int64(0),pathToFunc);
end

%%
for subj = 1:3%6:10
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
    content_check = dir(path.subj);

    if check_subj == 7 && size(content_check,1) > 3
        clearvars -except path timepoint plot_or_not pe pathToFunc subj_name subj content nfiles file_counter file_name

        content = dir(path.subj);
        nfiles = size(content,1);
        file_counter = 0;
        for file = 1:nfiles
            if contains(content(file).name, '.gt3x')
                file_counter = file_counter + 1;
                file_name{file_counter} = content(file).name;
            end
        end


        %% obtain raw acceleration data non-wear-periods from the hip sensor
        %try
        disp('     Extracting meta data .gt3x files......')
        for f = 1:size(file_name,2)
            temp = py.gt3x_functions.unzip_gt3x_file(fullfile(path.subj,file_name{f}));
            info = py.gt3x_functions.extract_info(temp{2});
            acceleration_scale = py.float(info{'Acceleration_Scale'});
            sample_rate = py.int(info{'Sample_Rate'});

            rawdata.acceleration_scale = double(acceleration_scale);
            rawdata.sample_rate = double(sample_rate);

            if info{'Limb'} == 'Wrist'
                if info{'Side'} == 'Left'
                    left_data = file_name{f};
                    disp(['         ', file_name{f}, newline, '         Side: Left Wrist'])
                    log.left = py.gt3x_functions.extract_log(temp{1}, acceleration_scale, sample_rate);
                elseif info{'Side'} == 'Right'
                    right_data = file_name{f};
                    disp(['         ', file_name{f}, newline, '         Side: Right Wrist'])
                    log.right = py.gt3x_functions.extract_log(temp{1}, acceleration_scale, sample_rate);
                end
            elseif info{'Limb'} == 'Waist'
                hip_data = file_name{f};
                disp(['         ', file_name{f}, newline, '         Side: Hip'])
                log.hip = py.gt3x_functions.extract_log(temp{1}, acceleration_scale, sample_rate);
            end
            clear temp info
        end

        disp('     Done')
        disp('     __________________________')
        %% obtain raw acceleration data non-wear-periods from the hip sensor
        disp(['     Reading in the hip data: ', hip_data '......'])
        path.subj = [path.subj, '\\'];
        pyOut.hip = py.gt3x2df.gt3x2df(path.subj, hip_data);
        rawdata.hip = double(pyOut.hip{1});
        T_hip = table(rawdata.hip(:,1), rawdata.hip(:,2), rawdata.hip(:,3));
        T_hip.Properties.VariableNames = {'X', 'Y', 'Z'};
        parquetwrite([path.subj, 'rawdata_hip.parq'],T_hip)


        disp(['     Reading in the left data: ', left_data '......'])
        pyOut.left = py.gt3x2df.gt3x2df(path.subj, left_data);
        rawdata.left = double(pyOut.left{1});
        T_left = table(rawdata.left(:,1), rawdata.left(:,2), rawdata.left(:,3));
        T_left.Properties.VariableNames = {'X', 'Y', 'Z'};
        parquetwrite([path.subj, 'rawdata_left.parq'], T_left)


        disp(['     Reading in the right data: ', right_data '......'])
        pyOut.right = py.gt3x2df.gt3x2df(path.subj, right_data);
        rawdata.right = double(pyOut.right{1});
        T_right = table(rawdata.right(:,1), rawdata.right(:,2), rawdata.right(:,3));
        T_right.Properties.VariableNames = {'X', 'Y', 'Z'};
        parquetwrite([path.subj, 'rawdata_right.parq'],T_right)


        disp('     Done')

        %save([path.subj, 'rawdata.mat'], 'rawdata', '-v7.3')
        


%         if ~isempty(double(pyOut.hip{1,1}))
%             %% from np.array to double
%            
%             raw_data = pyOut.hip{1};
% 
%             disp('     Done')
%             disp('     __________________________')
            
%             return
%             %% Wear periods
%             disp('     Segment wear vs non-wear periods......')
%             start_indices = find([0; diff(hip.data(:,4)) == 1]);
%             end_indices = find([0; diff(hip.data(:,4)) == -1]);
% 
%             if ~isempty(start_indices) && ~isempty(end_indices)
%                 if end_indices(1) - start_indices(1) < 0
%                     start_indices = [1; start_indices];
%                 end
% 
%                 if end_indices(end) - start_indices(end) < 0
%                     end_indices = [end_indices; size(hip.data(:,4),1)];
%                 end
% 
%                 hz = 30;
%                 min_non_wear_time_window = 135;
% 
%                 second_threshold = (hz * 60) * min_non_wear_time_window;
% 
%                 mask = find(end_indices - start_indices <= second_threshold);
%                 start_indices(mask) = [];
%                 end_indices(mask) = [];
% 
%                 %% calculate wear time
%                 disp('     Calulate wear time......')
%                 wear_blocks(1,:) = end_indices - start_indices; % number of samples in a block
%                 wear_blocks(2,:) = wear_blocks(1,:) / hz; % number of seconds in a block
%                 wear_blocks(3,:) = wear_blocks(2,:) / 60; % number of minutes in a block
%                 wear_blocks(4,:) = wear_blocks(3,:) / 60; % number of hours in a block
% 
%                 disp('     Done')
%                 disp('     __________________________')
%                 if size(wear_blocks,2) >= 5 && sum(wear_blocks(4,:) >= 12) >= 5
%                     %% finish rest of the analysis!
% 
%                     blocks_of_interest = find(wear_blocks(4,:) >= 12);
%                     for idx = 1:size(blocks_of_interest,2)
%                         blocks(idx,:) = [start_indices(idx), end_indices(idx)];
%                     end
% 
%                     disp(['     Reading in the left arm data: ', left_data '......'])
%                     pyOut.left = py.gt3x2df.gt3x2df(path.subj, left_data);
%                     disp(['     Reading in the right arm data: ', right_data '......'])
%                     pyOut.right = py.gt3x2df.gt3x2df(path.subj, right_data);
% 
%                     if ~isempty(double(pyOut.left{1,1})) && ~isempty(double(pyOut.right{1,1}))
%                         disp(['     Wear time and data check completed......', newline, '     Conclusion: Data will be analysed'])
%                         Conclusion = 'Analysed';
%                         basic_table = table(string(subj_name), string(timepoint), wear_blocks(4,:), string(Conclusion));
%                         basic_table.Properties.VariableNames = {'ppID', 'Timepoint', 'wear time', 'Conclusion'};
%                         writetable(basic_table, 'Weartime.xlsx', Sheet=timepoint, WriteMode='append')
%                         disp('     __________________________')
% 
%                         left.data(:,1) = double(py.array.array('d', py.numpy.nditer(pyOut.left{2})))';
%                         left.data(:,2) = double(py.array.array('d', py.numpy.nditer(pyOut.left{3})))';
%                         left.data(:,3) = double(py.array.array('d', py.numpy.nditer(pyOut.left{4})))';
% 
% 
%                         right.data(:,1) = double(py.array.array('d', py.numpy.nditer(pyOut.right{2})))';
%                         right.data(:,2) = double(py.array.array('d', py.numpy.nditer(pyOut.right{3})))';
%                         right.data(:,3) = double(py.array.array('d', py.numpy.nditer(pyOut.right{4})))';
% 
% 
% 
%                         if size(left.data,1) >= (size(hip.data,1) / 100) * 85 && size(right.data,1) >= (size(hip.data,1) / 100) * 85
%                             disp('     Shape check completed......')
%                             for b = 1:size(blocks,1)
%                                 disp(['     Analyzing wear block ', num2str(b)])
%                                 data.left{b,1} = left.data(blocks(b,1):blocks(b,2), :);
%                                 data.left{b,2} = mean(left.data(blocks(b,1):blocks(b,2), :));
%                                 data.right{b,1} = right.data(blocks(b,1):blocks(b,2), :);
%                                 data.right{b,2} = mean(right.data(blocks(b,1):blocks(b,2), :));
%                                 data.hip{b,1} = hip.data(blocks(b,1):blocks(b,2), :);
% 
%                                 %% correct sensor tilt and axis definition to match Lum et al.,
% 
%                                 disp('     Redefine axis defenition......')
%                                 % Left
%                                 acc.x = data.left{b,1}(:,2);
%                                 acc.y = data.left{b,1}(:,1);
%                                 acc.z = data.left{b,1}(:,3) * -1;
% 
%                                 data.L{b,1} = [acc.x, acc.y,  acc.z];
% 
%                                 clear acc
% 
%                                 % Right
%                                 acc.x = data.right{b,1}(:,2) * -1;
%                                 acc.y = data.right{b,1}(:,1) * -1;
%                                 acc.z = data.right{b,1}(:,3) * -1;
% 
%                                 data.R{b,1} = [acc.x, acc.y,  acc.z];
% 
%                                 clear acc
% 
% 
%                                 %% resample data from 30Hz (our sensors) to 50Hz (needed for the model)
%                                 disp('     Resample data from 30 to 50Hz......')
%                                 fs = 30;
%                                 fs_new = 50;
% 
%                                 x = 1:length(data.L{b,1}); % old time axis of the data
%                                 xq = 1:fs/fs_new:length(data.L{b,1}); % new time axis for the data
% 
%                                 data.L_sp{b,1} = interp1(x, data.L{b,1}, xq, 'spline');
%                                 data.R_sp{b,1} = interp1(x, data.R{b,1}, xq, 'spline');
% 
%                                 %% prepare for model
%                                 disp('     Calulate features from data block......')
% 
%                                 load('model.mat')
%                                 n=200;
%                                 feature_l = featurecalc1(data.L_sp{b,1},n);
%                                 feature_r = featurecalc1(data.R_sp{b,1},n);
% 
%                                 yyfit_l = trainedModel.predictFcn(feature_l);
%                                 yyfit_r = trainedModel.predictFcn(feature_r);
% 
%                                 %% calculate minutes functionally active
%                                 disp(['     Calculate the total functional activity for block ', num2str(b)])
%                                 data.wear_time(b,1) = ((size(data.L_sp{b,1},1) / 50) / 3600);
%                                 data.pred_L(b,1) = (size(yyfit_l(yyfit_l == 1), 1) * 4) / 60; %minuten
%                                 data.pred_L(b,2) = data.pred_L(b,1)/60; % hours
% 
%                                 data.pred_R(b,1) = (size(yyfit_r(yyfit_r == 1), 1) * 4) / 60; % minutes
%                                 data.pred_R(b,2) = data.pred_R(b,1) / 60; % hours
% 
%                             end
% 
%                             %% Output table
%                             disp('      Write data to output table.........')
%                             no_days = size(blocks,1);
%                             avg_weartime = mean(data.wear_time);
%                             avg_act_l = mean(data.pred_L(:,2));
%                             avg_perc_l = avg_act_l / avg_weartime;
%                             avg_act_r = mean(data.pred_R(:,2));
%                             avg_perc_r = avg_act_r / avg_weartime;
% 
%                             tbl = table(string(subj_name), string(timepoint),no_days, avg_weartime, avg_act_l, avg_perc_l, avg_act_r, avg_perc_r);
%                             tbl.Properties.VariableNames = {'subject id', 'Timepoint', 'wear days', 'wear time', 'mean use L', 'percentage use L', 'mean use R', 'percentage use R'};
%                             writetable(tbl, 'output.xlsx', 'WriteMode','append','AutoFitWidth',true, Sheet=timepoint)
% 
% 
%                             %% figures
%                             if plot_or_not
%                                 maximum_signal = max(max(abs(hip.data(:,1:3))));
%                                 ymax = maximum_signal + ((maximum_signal/100) * 10);
%                                 ymin = ymax * -1;
% 
%                                 figure;
%                                 hold on
% 
%                                 x_values = [blocks(:,1) blocks(:,2) blocks(:,2) blocks(:,1)];
%                                 y_values = [ymin ymin ymax ymax];
% 
%                                 fill(x_values, y_values, [0.78, 0.81, 0.84], 'EdgeColor','none');
%                                 plot(hip.data(:,1:3));
%                                 hold off
%                                
%                             end
%                         else
%                             disp(['     Shape check completed......', newline, '     Revision of Conclusion: Data will not be analysed'])
% 
% 
%                             disp('      Write data to output table.........')
%                             no_days = size(blocks,1);
%                             avg_weartime = "shape check failed, check de data in actigraph software for errors";
%                             avg_act_l = nan;
%                             avg_perc_l = nan;
%                             avg_act_r = nan;
%                             avg_perc_r = nan;
% 
%                             tbl = table(string(subj_name), string(timepoint),no_days, avg_weartime, avg_act_l, avg_perc_l, avg_act_r, avg_perc_r);
%                             tbl.Properties.VariableNames = {'subject id', 'Timepoint', 'wear days', 'wear time', 'mean use L', 'percentage use L', 'mean use R', 'percentage use R'};
%                             writetable(tbl, 'output.xlsx', 'WriteMode','append','AutoFitWidth',true, Sheet=timepoint)
% 
% 
%                         end %Shape check
% 
%                     else
%                         disp(['     Wear time and data check completed......', newline, '     Conclusion: Data will not be analysed'])
%                         Conclusion = 'Corrupt files';
%                         basic_table = table(string(subj_name), string(timepoint), wear_blocks(4,:), string(Conclusion));
%                         basic_table.Properties.VariableNames = {'ppID', 'Timepoint', 'wear time', 'Conclusion'};
%                         writetable(basic_table, 'Weartime.xlsx', Sheet=timepoint, WriteMode='append')
% 
%                     end % python data is not empty
%                 else
%                     disp(['     Wear time check completed......', newline, '     Conclusion: Data will not be analysed'])
%                     Conclusion = 'Not Analysed';
%                     basic_table = table(string(subj_name), string(timepoint), wear_blocks(4,:), string(Conclusion));
%                     basic_table.Properties.VariableNames = {'ppID', 'Timepoint', 'wear time', 'Conclusion'};
%                     writetable(basic_table, 'Weartime.xlsx', Sheet=timepoint, WriteMode='append')
% 
% 
% 
%                     disp(['Processing: ' subj_name newline ' at timepoint: ' timepoint ' COMPLETED'])
%                     disp('====================================================')
%                     %catch
% 
% 
%                 end % if the wear time are long enough
%             else
%                 disp(['     Wear time check completed......', newline, '     Conclusion: Data will not be analysed'])
%                 Conclusion = "Not Analysed";
% 
% 
%                 basic_table = table(string(subj_name), string(timepoint), "no wear blocks", Conclusion);
%                 basic_table.Properties.VariableNames = {'ppID', 'Timepoint', 'wear blocks',  'Conclusion'};
%                 writetable(basic_table, 'Weartime.xlsx', Sheet=timepoint, WriteMode='append')
% 
% 
% 
%                 disp(['Processing: ' subj_name newline ' at timepoint: ' timepoint ' COMPLETED'])
%                 disp('====================================================')
% 
%             end % wear time could be detected
%         end % if hip data is not empty
    end % subject path exists and folder is not empty
end% loop through subjects