%% call python test

clear; clc

%% set python environment
pe = pyenv("Version", "C:\GBW_MyPrograms\Anaconda3\python.exe");

pathToFunc = fileparts(which('gt3x2df.py'));

if count(py.sys.path,pathToFunc) == 0
    insert(py.sys.path,int64(0),pathToFunc);
end

%% set data file
path.root = "C:\Users\u0117545\KU Leuven\An De Groef - DATA\BC_002\Accelerometrie\T0\";
file = "BC_002_H_T0.gt3x";



%% run python code
pyOut = py.gt3x2df.gt3x2df(path.root, file);  

%% from np.array to double
x               = double(py.array.array('d', py.numpy.nditer(pyOut{2})))';
y               = double(py.array.array('d', py.numpy.nditer(pyOut{3})))';
z               = double(py.array.array('d', py.numpy.nditer(pyOut{4})))'; 
non_wear_vector = double(py.array.array('d', py.numpy.nditer(pyOut{5})))';


%% find uninterupted wear periods

start_indices = find([0; diff(non_wear_vector) == 1]);
end_indices = find([0; diff(non_wear_vector) == -1]);

% if length(start_indices) ~= length(end_indices)
%     if (start_indices(1) - end_indices(1)) > 0
%         start_indices(1) = [];
%     elseif (start_indices(end) - end_indices(end)) > 0
%         end_indices(end+1) = length(non_wear_vector);
% 
%     end
% end


if end_indices(1) - start_indices(1) < 0
    end_indices(1) = [];
end

if end_indices(end) - start_indices(end) < 0 
    start_indices(end) = [];
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

%% calculate wear time
wear_blocks(:,1) = end_indices - start_indices; % number of samples in a block
wear_blocks(:,2) = wear_blocks(:,1) / hz; % number of seconds in a block
wear_blocks(:,3) = wear_blocks(:,2) / 60; % number of minutes in a block
wear_blocks(:,4) = wear_blocks(:,3) / 60; % number of hours in a block

% wear time should be at least:
%%
% 
%  5 days
%  12 hour per day
% 

if size(wear_blocks,1) >= 5 && sum(wear_blocks(:,4) > 12) >= 5



end


