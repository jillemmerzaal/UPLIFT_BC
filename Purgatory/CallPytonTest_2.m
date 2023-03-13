%% call python test

clear; clc

%% set python environment
pe = pyenv("Version", "C:\GBW_MyPrograms\Anaconda3\python.exe");

pathToFunc = fileparts(which('gt3x_functions.py'));

if count(py.sys.path,pathToFunc) == 0
    insert(py.sys.path,int64(0),pathToFunc);
end

%% set data file
path.root = "C:\Users\u0117545\KU Leuven\An De Groef - DATA\BC_001\Accelerometrie\T0\";
content = dir(path.root)


%%
pyOut = py.gt3x_functions.unzip_gt3x_file(fullfile(path.root, file), save_location = None, delete_source_file = False);

%%
pyOut = py.gt3x_functions.extract_info(pyOut{2})