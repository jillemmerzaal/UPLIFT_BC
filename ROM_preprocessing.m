%% maximal range of motion task
%

%% 1. input data
cd("C:\Users\u0117545\Documents\GitHub\ULIFT_BC")
addpath("C:\Users\u0117545\OneDrive - KU Leuven\2.Dataprocessing\Matlab\addons")


Timepoint   = 'T0';
movement    = 'Abductie'; %Abductie Anteflexie Exorotatie
path.root   = 'C:\Users\u0117545\KU Leuven\An De Groef - UPLIFT-BC\INVESTIGATOR SITE FILE\5. Data';
path.out    = fullfile(path.root,'Output','Database_ULIFT.mat');

plot_or_not = 1;


%% 2. load data
for subj = 2%1:3
    if subj < 10
        subj_name   = ['BCT_00' num2str(subj)];
    elseif subj < 100
        subj_name   = ['BCT_0' num2str(subj)];
    else
        subj_name   = ['BCT_', num2str(subj)];
    end

    disp(' ')
    disp(['Processing ' subj_name ': ' Timepoint '.....'])

    path.subj   = fullfile(path.root, subj_name, 'Xsens', Timepoint);
    check_subj  = exist(path.subj);

    if check_subj == 7
        %initialize counters
        counterR = 0;
        counterL = 0;

        content = dir(path.subj);
        nfiles = size(content,1);

        % Start loop through ULIFT files per subject
        for file = 4%1:nfiles


            if contains(content(file).name, movement) && contains(content(file).name, '.mvnx')
                number  = str2num(content(file).name(13:end-5));
                file_ik = fullfile(path.subj, content(file).name);

                [~,name, ~] = fileparts(content(file).name);
                [fileName] = regexprep(name, '-', '_');

                if contains(content(file).name, 'R')
                    arm = 'right';
                else
                    arm = 'left';
                end

                %% 2.1 Load xsens data
                % Change the filename here to the name of the file you would like to import
                disp(['    ' content(file).name ': read xsens file'])
                tree = load_mvnx(file_ik);

                % Read some basic data from the file
                mvnxVersion = tree.metaData.mvnx_version; % version of the MVN Studio used during recording

                if (isfield(tree.metaData, 'comment'))
                    fileComments = tree.metaData.comment; % comments written when saving the file
                end

                % Read some basic properties of the subject;

                frameRate = tree.metaData.subject_frameRate;
                suitLabel = tree.metaData.subject_label;
                originalFilename = tree.metaData.subject_originalFilename;
                recDate = tree.metaData.subject_recDate;
                segmentCount = tree.metaData.subject_segmentCount;

                % Retrieve sensor label
                % creates a struct with sensor data
                if isfield(tree,'sensorData') && isstruct(tree.sensorData)
                    sensorData = tree.sensorData;
                end

                % Retrieve segment labels
                % creates a struct with segment definitions
                if isfield(tree,'segmentData') && isstruct(tree.segmentData)
                    segmentData = tree.segmentData;
                end

                % Retrieve joint labels
                % creates a struct with segment definitions
                if isfield(tree,'jointData') && isstruct(tree.jointData)
                    jointData = tree.jointData;
                end

                if strcmp(arm, 'left')
                    jointno     = 12;
                    segmentno   = 12;
                else
                    jointno     = 8;
                    segmentno   = 8;
                end


                %% display results
                figure;
                stackedplot(jointData(jointno).jointAngleXZY)

                figure;
                stackedplot(jointData(jointno).jointAngle)

            end % end if movement && .mvnx
        end % loop of number of files
    end %check if subject exists
end %end subjects
