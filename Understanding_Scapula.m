%% Scapula angles check
% This script has been created to help understand the joint angles of the
% scapula and genohumaral joint in a controlled environment and with
% controlled momvement.

clear all; clc;
%% 1. input data
cd("C:\Users\u0117545\Documents\GitHub\ULIFT_BC")
addpath("C:\Users\u0117545\OneDrive - KU Leuven\2.Dataprocessing\Matlab\addons")

movement    =  'ABD';
path.root   = 'C:\Users\u0117545\KU Leuven\An De Groef - DATA\XsensTest';
%path.out    = fullfile(path.root,'Output','Database_ULIFT.mat');

plot_or_not = 1;

T = readtable("C:\Users\u0117545\KU Leuven\An De Groef - DATA\XsensTest\coordinatesystem.xlsx");


path.subj   = fullfile(path.root);



content = dir(path.subj);
nfiles = size(content,1);

% Start loop through isolated ROM files per subject
for file = 1:nfiles
    if contains(content(file).name, '.mvnx')
        number  = str2num(content(file).name(13:end-5));
        file_ik = fullfile(path.subj, content(file).name);

        [~,name, ~] = fileparts(content(file).name);
        [fileName] = regexprep(name, '-', '_');

        div = strfind(name, '-');
        fileroot = name(1:div-1);

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

        %% 2.2 Extract the relevant kinematics
        %-------------------------------------
        disp(['    ' content(file).name ': extract relevant Kinematics'])
        if strcmp(arm, 'right')
            jointNo     = 7:10; % right upper extremity
        else
            jointNo     = 11:14; % left upper extremity
        end

        % initialise joint names
        jointNames = ['Scapula', "Glenohumeraal", "Elbow", "Wrist"];
        nf = 101;

        coord = cell2mat(T.CoordinateSystem(find(strcmp(T.Naam, fileroot))));



        for jnt = 1:4

            % temporary time curves
            %----------------------

            if strcmp(coord, 'XYZ')

                if strcmp(jointNames{jnt}, 'Glenohumeraal')
                    temp.X = jointData(jointNo(jnt)).jointAngle(: ,1);
                    temp.Y = jointData(jointNo(jnt)).jointAngle(: ,2);
                    temp.Z = jointData(jointNo(jnt)).jointAngle(: ,3);

                elseif strcmp(jointNames{jnt}, "Scapula")
                    temp.X = jointData(jointNo(jnt)).jointAngle(: ,1);
                    temp.Y = jointData(jointNo(jnt)).jointAngle(: ,2);
                    temp.Z = jointData(jointNo(jnt)).jointAngle(: ,3);

                else
                    temp.X = jointData(jointNo(jnt)).jointAngle(: ,1);
                    temp.Y = jointData(jointNo(jnt)).jointAngle(:, 2);
                    temp.Z = jointData(jointNo(jnt)).jointAngle(:, 3);
                end
            elseif strcmp(coord, 'XZY')
                if strcmp(jointNames{jnt}, 'Glenohumeraal')
                    temp.X = jointData(jointNo(jnt)).jointAngleXZY(: ,1);
                    temp.Y = jointData(jointNo(jnt)).jointAngleXZY(: ,2);
                    temp.Z = jointData(jointNo(jnt)).jointAngleXZY(: ,3);

                elseif strcmp(jointNames{jnt}, "Scapula")
                    temp.X = jointData(jointNo(jnt)).jointAngle(: ,1);
                    temp.Y = jointData(jointNo(jnt)).jointAngle(: ,2);
                    temp.Z = jointData(jointNo(jnt)).jointAngle(: ,3);

                else
                    temp.X = jointData(jointNo(jnt)).jointAngle(: ,1);
                    temp.Y = jointData(jointNo(jnt)).jointAngle(:, 2);
                    temp.Z = jointData(jointNo(jnt)).jointAngle(:, 3);
                end

            end
            % initiation of joint names

            if strcmp(jointNames{jnt}, "Elbow") || strcmp(jointNames{jnt}, "Wrist")
                IK_X = [jointNames{jnt}, '_FrontalPlaneRadial'];
                IK_Y = [jointNames{jnt}, '_TransversalPlanePronation'];
                IK_Z = [jointNames{jnt}, '_SagitalPlaneFlexion'];
            else
                IK_X = [jointNames{jnt}, '_FrontalPlaneAbduction'];
                IK_Y = [jointNames{jnt}, '_TransversalPlaneInternal'];
                IK_Z = [jointNames{jnt}, '_SagitalPlaneFlexion'];
            end



            % Time normalised repetitions
            %----------------------------
            df.(fileName).(IK_X) = interp1([1:size(temp.X,1)], ...
                temp.X', [1:(size(temp.X,1))/nf:size(temp.X,1)], 'spline');
            df.(fileName).(IK_Y) = interp1([1:size(temp.Y,1)], ...
                temp.Y', [1:(size(temp.Y,1))/nf:size(temp.Y,1)], 'spline');
            df.(fileName).(IK_Z) = interp1([1:size(temp.Z,1)], ...
                temp.Z', [1:(size(temp.Z,1))/nf:size(temp.Z,1)], 'spline');

            clear temp

        end


        %% plot the individual normalised timecurves
        JointAngles = fieldnames(df.(fileName));

        figure('units','normalized','outerposition',[0 0 1 1])
        p = tiledlayout(4,3);
        [figureTitle] = regexprep(name, '_', ' ');
        title(p, [figureTitle, ' -- ', coord])
        xlabel(p,'% movement')
        ylabel(p,['Joint angle', char(176)])

        axislabel = ["AD(-)/AB(+)", "Ext(-)/Int(+)", "Ex(-)/Fl(+)",...
            "AD(-)/AB(+)", "Ext(-)/Int(+)", "Ex(-)/Fl(+)",...
            "Ulnar(-)/Ratial(+)", "Sup(-)/Pro(+)", "Ex(-)/Fl(+)",...
            "Ulnar(-)/Ratial(+)", "Sup(-)/Pro(+)", "Ex(-)/Fl(+)"];

        for jnt = 1:length(JointAngles)
            nexttile
            plot(df.(fileName).(JointAngles{jnt}))
            xlim([0, 100])
            ylabel(axislabel{jnt})

            [PlotTitle] = regexprep(JointAngles{jnt}, '_', ' ');
            title(PlotTitle)
        end

        filepath = fullfile('C:\Users\u0117545\KU Leuven\An De Groef - DATA\XsensTest\Figures', [fileName, '.tif']);

        saveas(p,filepath);


close all
    end% contains filename and .mvnx
end% number of fils
