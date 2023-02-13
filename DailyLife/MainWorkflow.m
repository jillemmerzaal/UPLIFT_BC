%% Functional activities in daily life.

%%
%
% * UPLIFT (and EduCan) project
% * at least 5 days of ware time
% * at least 12 hours per day.
% * minutes active as an average per day
% * sample


clear; clc; close all;

cd('C:\Users\u0117545\Documents\GitHub\ULIFT_BC\DailyLife\')
path.root 		= 'C:\Users\u0117545\Documents\GitHub\ULIFT_BC\DailyLife\';
path.data       = 'C:\Users\u0117545\KU Leuven\An De Groef - DATA';
timepoint       = 'T0';


%% set python environment
pe = pyenv("Version", "C:\GBW_MyPrograms\Anaconda3\python.exe");

pathToFunc = fileparts(which('gt3x2df.py'));

if count(py.sys.path,pathToFunc) == 0
    insert(py.sys.path,int64(0),pathToFunc);
end



%% 
for subj = 2%1:5
    if subj < 10
        subj_name = ['BC_00', num2str(subj)];
    elseif subj < 100
        subj_name = ['BC_0', num2str(subj)];
    else
        subj_name = ['BC_', num2str(subj)];
    end


    disp(' ')
    disp(['Processing ' subj_name])
    disp(['     at timepoint: ', timepoint])

    path.subj   = fullfile(path.data, subj_name, 'Accelerometrie', timepoint);
    check_subj  = exist(path.subj, "dir");
    content_check = dir(path.subj);

    if check_subj == 7 && size(content_check,1) > 3
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
        %TODO: automate reading in the hip data

        hip_data = file_name{2};

        disp(['     Reading in the hip data: ', hip_data])
        path.subj = [path.subj, '\\'];
        pyOut = py.gt3x2df.gt3x2df(path.subj, hip_data);

        disp('Done')

        %% from np.array to double
       disp('   from numpy to double')
        x               = double(py.array.array('d', py.numpy.nditer(pyOut{2})))';
        y               = double(py.array.array('d', py.numpy.nditer(pyOut{3})))';
        z               = double(py.array.array('d', py.numpy.nditer(pyOut{4})))';
        non_wear_vector = double(py.array.array('d', py.numpy.nditer(pyOut{5})))';

        disp('  Done')

        %% Wear periods
        disp('  Segment wear vs non-wear periods')

       
        start_indices = find([0; diff(non_wear_vector) == 1]);
        end_indices = find([0; diff(non_wear_vector) == -1]);

        if end_indices(1) - start_indices(1) < 0
            start_indices = [1; start_indices];
        end

        if end_indices(end) - start_indices(end) < 0
            end_indices = [end_indices; size(non_wear_vector,1)];
        end

        hz = 30;
        min_non_wear_time_window = 60;

        second_threshold = (hz * 60) * min_non_wear_time_window;

        mask = find(end_indices - start_indices <= second_threshold);
        start_indices(mask) = [];
        end_indices(mask) = [];

        figure; plot([x, y , z]); hold on
        plot(diff(non_wear_vector), 'LineWidth', 2)
        plot(non_wear_vector, 'LineStyle', '--')
        xline(start_indices, 'Color', "#7E2F8E", LineWidth=2)
        xline(end_indices, 'Color',"#A2142F", LineWidth=2)

        disp('  Done')
        %% calculate wear time
        wear_blocks(1,:) = end_indices - start_indices; % number of samples in a block
        wear_blocks(2,:) = wear_blocks(1,:) / hz; % number of seconds in a block
        wear_blocks(3,:) = wear_blocks(2,:) / 60; % number of minutes in a block
        wear_blocks(4,:) = wear_blocks(3,:) / 60; % number of hours in a block

        if size(wear_blocks,2) >= 5 && sum(wear_blocks(4,:) >= 12) >= 5
            Conclusion = 'Analysed';
            basic_table = table(string(subj_name), string(timepoint), wear_blocks(4,:), string(Conclusion));
            basic_table.Properties.VariableNames = {'ppID', 'Timepoint', 'wear time', 'Conclusion'};
            writetable(basic_table, 'Weartime.xlsx', Sheet=timepoint, WriteMode='append')

            %% finish rest of the analysis! 
            % import left and right arm data and analyse functionaly active
            % the start and end indices from the hip data should be used to
            % analyse the different days. 

        else
            Conclusion = 'Not Analysed';
            basic_table = table(string(subj_name), string(timepoint), wear_blocks(4,:), string(Conclusion));
            basic_table.Properties.VariableNames = {'ppID', 'Timepoint', 'wear time', 'Conclusion'};
            writetable(basic_table, 'Weartime.xlsx', Sheet=timepoint, WriteMode='append')
        end % wear time is long enough
        clear wear_blocks
    end % subject path exists and folder is not empty
end% loop through subjects