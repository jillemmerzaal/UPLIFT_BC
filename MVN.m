function [sensorData, segmentData, jointData]= MVN(file_ik)

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


