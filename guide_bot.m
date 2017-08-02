function [ output_args ] = guide_bot(name,demands,requirements,defaults,input,options)
% Looping over the inputs and running initialize
clc

filename_guide_bot = mfilename('fullpath');
char_name = 'guide_bot';
path_variable = filename_guide_bot(1:end-length(char_name));

addpath([path_variable '/guide_bot_source'])

options.guide_bot_path = path_variable;

initialize(name)

for i=1:length(input)
mcstas_bot(input{i},demands,requirements,defaults,name,options)
end 

disp('guide_bot function finished. Please look through the above command window output for warnings.')
%mfilename('fullpath')
end

