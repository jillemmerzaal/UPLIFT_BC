%% ULIFT_preprocessing_SingleRep

% this code processess the ULIFT data, but only a single, middel repetition from
% phase 1 and a single middel repetition from phase 4

clear all; close all; clc
%% 1. input data
cd("C:\Users\u0117545\Documents\GitHub\ULIFT_BC")
addpath("C:\Users\u0117545\OneDrive - KU Leuven\2.Dataprocessing\Matlab\addons")


Timepoint   = 'T0';
movement    = "ULIFT";
path.root   = 'C:\Users\u0117545\KU Leuven\An De Groef - DATA';
path.out    = fullfile(path.root,'Output');

figure;
tiledlayout('flow')

plot_or_not = 1;
 %initialize counters
 counter=0;
%% 2. load data
for subj = 3
    if subj < 10
        subj_name   = ['BC_00' num2str(subj)];
    elseif subj < 100
        subj_name   = ['BC_0' num2str(subj)];
    else
        subj_name   = ['BC_', num2str(subj)];
    end

    disp(' ')
    disp(['Processing ' subj_name ': ' Timepoint '.....'])

    path.subj   = fullfile(path.root, subj_name, 'Xsens', Timepoint, 'Reproces');
    check_subj  = exist(path.subj);
    if check_subj == 7
       

        content = dir(path.subj);
        nfiles = size(content,1);
        % Start loop through ULIFT files per subject
        for file = 1:nfiles
            if contains(content(file).name, movement) && contains(content(file).name, '.mvnx')
                number  = str2num(content(file).name(13:end-5));
                file_ik = fullfile(path.subj, content(file).name);

                [~,name, ~] = fileparts(content(file).name);
                [fileName] = regexprep(name, '-', '_');

                d = strfind(name,'_');
                if size(d,2) == 1
                    arm = content(file).name(d+1);

                elseif size(d,2) == 2
                    temp = content(file).name(d(2)+1);
                    arm = [temp, '_SSS']; % SSS = self-selected speed
                    clear temp
                end

                disp(['     ' 'Analysing: ' fileName '.....'])
                disp(['   ' 'Arm of interst: ' arm '.....'])

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

                if contains(arm, 'L')
                    jointno     = 14;
                    segmentno   = 14;
                else
                    jointno     = 10;
                    segmentno   = 10;
                end


                %% set counters
                
                counter         = counter +1 ;
                
                %% 2.2 D efine start & end points of each repetition based on velocity of lower arm
                %---------------------------------------------------------------------------------
                disp(['    ' content(file).name ': define change points'])

                % change points
                % Find the change points of the position data and the minima and maxima in the
                % velocity data of the lower arm.
                %
                % Hypothesis is that on average the position of the high and the low sections
                % can be seperated using the position data. And a better idea of the precice start
                % end end points of the different target peaks
                %
                % Using both the changepoints and the arm velocity data a more robust segmenttation
                % of the ULIFT task can be achieved
                %----------------------------------

                posdata = segmentData(jointno).position(:,3);
                [changeIndices,segmentMean] = ischange(posdata,"MaxNumChanges",2);
                x = find(changeIndices);


                thresh_ph1 = mean(segmentMean(1:x(1))) + mean(segmentMean(1:x(1)))*0.05;
                [maxIndices_ph1, ~] = peakfinder(posdata(1:x(1)), [], thresh_ph1, 1, []);

                thresh_ph4 = mean(segmentMean(x(2):end)) + mean(segmentMean(x(2):end))*0.025;
                [maxIndices_ph4, ~] = peakfinder(posdata(x(2):end), [], thresh_ph4, 1, []);
                maxIndices_ph4 = maxIndices_ph4 + x(2);

                % display the results of the change points
                nexttile
                plot(posdata,"Color",[77 190 238]/255,"DisplayName","Input data")
                hold on

                %               Plot change points
                x = repelem(find(changeIndices),3);
                y = repmat([ylim(gca) missing]',nnz(changeIndices),1);
                plot(x,y,"Color",[51 160 44]/255,"LineWidth",1,"DisplayName","Change points")
                title("Number of change points: " + nnz(changeIndices))
                % Plot segments between change points
                plot(segmentMean,"Color",[64 64 64]/255,"DisplayName","Segment mean")

                % Plot local maxima phase 1
                plot(maxIndices_ph1,posdata(maxIndices_ph1),"^","Color",[217 83 25]/255,...
                    "MarkerFaceColor",[217 83 25]/255,"DisplayName","Local maxima")
                % Plot local maxima phase 4
                plot(maxIndices_ph4,posdata(maxIndices_ph4),"^","Color",[217 83 25]/255,...
                    "MarkerFaceColor",[217 83 25]/255,"DisplayName","Local maxima")


                ppID{counter, 1} = subj_name;
                Phase1(counter,1) = size(maxIndices_ph1,1);
                filename{counter,1} = fileName;
                Phase4(counter,1) = size(maxIndices_ph4,1);
                
                if Phase1(counter,1) == 3 && Phase4(counter,1) == 3
                    run(counter,1) = 1;
                else
                    run(counter,1 ) = 0;
                end



               
            end % if contains relevant names
        end% loop through number of files
        T = table(ppID, filename, Phase1, Phase4, run);

        writetable(T, fullfile(path.out, [subj_name, '.xlsx']))
    end% if subject exists
end %number of subjects

