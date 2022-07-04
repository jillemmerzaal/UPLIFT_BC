function [sensorData, segmentData, jointData, tree]= MVN(file_ik)

tree = load_mvnx(file_ik);

% Read some basic data from the file
tree.mvnxVersion = tree.metaData.mvnx_version; %version of the MVN Studio used during recording

if (isfield(tree.metaData, 'comment'))
    tree.fileComments = tree.metaData.comment; %comments written when saving the file
end

% Read some basic properties of the subject;

tree.frameRate = tree.metaData.subject_frameRate;
tree.suitLabel = tree.metaData.subject_label;
tree.originalFilename = tree.metaData.subject_originalFilename;
tree.recDate = tree.metaData.subject_recDate;
tree.segmentCount = tree.metaData.subject_segmentCount;

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


