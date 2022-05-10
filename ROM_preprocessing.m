%% maximal range of motion task
%

%% 1. input data
cd("C:\Users\u0117545\Documents\GitHub\ULIFT_BC")
addpath("C:\Users\u0117545\OneDrive - KU Leuven\2.Dataprocessing\Matlab\addons")


Timepoint   = 'T0';
movement    =  'AB'; %Abductie'; %Abductie Anteflexie Exorotatie
path.root   = 'C:\Users\u0117545\KU Leuven\An De Groef - 5. Data';
path.out    = fullfile(path.root,'Output','Database_ULIFT.mat');

plot_or_not = 1;


%% 2. load data
for subj = 5%1:3
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

        % Start loop through ROM files per subject
        for file = 1:nfiles


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
                mvnxVersion = tree.metaData.mvnx_version; %version of the MVN Studio used during recording

                if (isfield(tree.metaData, 'comment'))
                    fileComments = tree.metaData.comment; %comments written when saving the file
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
                    segmentno   = 14;
                else
                    jointno     = 8;
                    segmentno   = 10;
                end


                %% display results
                figure;
                stackedplot(jointData(jointno).jointAngleXZY)
                title("Shoulder XZY joint angle")
                xlabel('Sample no')

                figure;
                stackedplot(segmentData(segmentno).velocity)
                title('Velocity lowerarm segment')

                figure;
                stackedplot(segmentData(segmentno).position)
                title('Position lowerarm segment')
                %% 2.2 Define start & end points of each repetition based on change in variance of shoulder angles
                %-------------------------------------------------------------------------------------------------
                disp(['    ' content(file).name ': define start and end points'])

                b = mean(jointData(jointno).jointAngleXZY(1:50, 1)) + mean(jointData(jointno).jointAngleXZY(1:50, 1))* 0.05;
                
               start = find(jointData(jointno).jointAngleXZY(:,1) > b, 1, "first");
               eind = find(jointData(jointno).jointAngleXZY(:,1) > b, 1, "last");

               figure; 
               tiledlayout('flow')
               nexttile
                plot(jointData(jointno).jointAngleXZY(:,1))
                xline(start)
                xline(eind)
                yline(b)
                title('Repetitions on baseline')
                ylabel("Shoulder XZY joint angle")
                xlabel('Sample no')

              
                % find change points (Method, mean)
                df = jointData(jointno).jointAngleXZY(:,1);
                [changeIndices,segmentMean] = ischange(df,"MaxNumChanges",6);

                x = find(changeIndices);
                nexttile; plot(jointData(jointno).jointAngleXZY(:,1))
                xline(x)
                title('Repetitions on change points (Method, mean)')
                ylabel("Shoulder XZY joint angle")
                xlabel('Sample no')

                % find change points (Method: linear)
                [changeIndices,segmentMean] = ischange(df, "linear", 'MaxNumChanges',6);

                x = find(changeIndices);
                nexttile; plot(jointData(jointno).jointAngleXZY(:,1))
                xline(x)
                title('Repetitions on change points (Method: linear)')
                ylabel("Shoulder XZY joint angle")
                xlabel('Sample no')


                % find change points (Method: variance)
                [changeIndices,segmentMean] = ischange(df, 'variance', 'MaxNumChanges',6);

                x = find(changeIndices);
                
                
                nexttile; plot(jointData(jointno).jointAngleXZY(:,1))
                xline(x)
                xline(start, 'g')
                title('Repetitions on change points (Method: variance)')
                ylabel("Shoulder XZY joint angle")
                xlabel('Sample no')

                if start-x(1) > 10 
                    x(1) = start;
                end

                %% 2.3 Extract the relevant kinematics
                %-------------------------------------
                disp(['    ' content(file).name ': extract relevant Kinematics'])
                if strcmp(arm, 'right')
                    counterR    = counterR + 1;
                    counter     = counterR;
                    jointNo     = 7:10; % right upper extremity
                else
                    counterL    = counterL + 1;
                    counter     = counterL;
                    jointNo     = 11:14; % left upper extremity
                end

                

                % Save the timecurves to struct, per participant, per
                % repetition.
                % Struct holds:
                    % full timecurve
                    % start and end points of repetitions
                    % timecurves for each repetition
                    % time normalised timecurves for each repetition

                % set General information per participant, per trial. 
                Data_out.(movement).(Timepoint).General.CutIndices.(fileName) = x;

                % initialise joint names
                jointNames = ['Scapula', "Glenohumeraal", "Elbow", "Wrist"];

                
                nf = 101;
                for jnt = 1:4
                    for rep = 1:3
                        % temporary time curves with the phase data
                        %------------------------------------------
                        s2 = rep*2; 
                        s1 = s2-1;
                        
                        if strcmp(jointNames{jnt}, 'Glenohumeraal')
                            temp.X = jointData(jointNo(jnt)).jointAngleXZY(x(s1):x(s2) ,1);
                            temp.Y = jointData(jointNo(jnt)).jointAngleXZY(x(s1):x(s2) ,2);
                            temp.Z = jointData(jointNo(jnt)).jointAngleXZY(x(s1):x(s2) ,3);
                        else
                            temp.X = jointData(jointNo(jnt)).jointAngle(x(s1):x(s2) ,1);
                            temp.Y = jointData(jointNo(jnt)).jointAngle(x(s1):x(s2), 2);
                            temp.Z = jointData(jointNo(jnt)).jointAngle(x(s1):x(s2), 3);
                        end

                        % initiation of joint names
                        if strcmp(jointNames{jnt}, 'Scapula')
                            IK_X = [jointNames{jnt}, '_lateralRotation'];
                            IK_Y = [jointNames{jnt}, '_Protraction'];
                            IK_Z = [jointNames{jnt}, '_AnterieureTilt'];

                        else
                            IK_X = [jointNames{jnt}, '_abbuction'];
                            IK_Y = [jointNames{jnt}, '_rotation'];
                            IK_Z = [jointNames{jnt}, '_flexion'];
                        end

                        % full ROM timeseries
                        %--------------------

                        % Timedata of the repetitions
                        %----------------------------
                        repetition = ['rep', num2str(rep)];
                        Data_out.(movement).(Timepoint).IK.(arm).(repetition).(fileName).(IK_X) = temp.X;
                        Data_out.(movement).(Timepoint).IK.(arm).(repetition).(fileName).(IK_Y) = temp.Y;
                        Data_out.(movement).(Timepoint).IK.(arm).(repetition).(fileName).(IK_Z) = temp.Z;

                        % Time normalised repetitions
                        %----------------------------
                        Data_out.(movement).(Timepoint).IK.(arm).normalised.(IK_X)(:,rep) = interp1([1:size(temp.X,1)], ...
                            temp.X', [1:(size(temp.X,1))/nf:size(temp.X,1)], 'spline');
                        Data_out.(movement).(Timepoint).IK.(arm).normalised.(IK_Y)(:,rep) = interp1([1:size(temp.Y,1)], ...
                            temp.Y', [1:(size(temp.Y,1))/nf:size(temp.Y,1)], 'spline');
                        Data_out.(movement).(Timepoint).IK.(arm).normalised.(IK_Z)(:,rep) = interp1([1:size(temp.Z,1)], ...
                            temp.Z', [1:(size(temp.Z,1))/nf:size(temp.Z,1)], 'spline');

                        clear temp
                    end
                end





            end % end if movement && .mvnx
        end % loop of number of files
    end %check if subject exists
     %% Save data
    %--------------

    if exist('Data_out','var')
        if exist(path.out,'file')
            load(path.out)
        end


        [Data.(subj_name).(movement).(Timepoint)] = Data_out.(movement).(Timepoint);
        save(path.out,'Data')
        clear Data Data_out
    end

    disp(['*********Finished ' subj_name '**********'])
    disp(' ')
end %end subjects
