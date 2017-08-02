function RunProjectList(varargin)
% Load the project list
load('projectlist.mat')

% If no arguments are given, then continue from last stopping point in list
if nargin == 0
    next_dir_number = allprojects{1};
% Start running the project list at a specific point
else
    next_dir_number = varargin{1}+1;
end

% Run the project (unless at end of project list)
if next_dir_number < numel(allprojects)+1
    nextproject = allprojects{next_dir_number};
    cd(nextproject)
    % system(['nohup sh ' nextproject '/launch_all.sh &'])
    system('nohup ./launch_all.sh &');
    allprojects{1} = allprojects{1} + 1;
    cd ..
    save('projectlist.mat', 'allprojects')
end