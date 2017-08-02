function initialize(name)
% initialize the records and sh scripts

[status,message,messageid]=mkdir(name);


if ispc
fid = fopen(['./' name '/windows_support.sh'],'w');
fprintf(fid,'#!/bin/bash\n# on cluster, excecute "chmod +x windows_support.sh\n# then "./windows_support.sh"\n');
fclose(fid);
end

l{1}='#!/bin/bash';
l{end+1}='# Made by guide_bot';
l{end+1}='# Will compile all the .instr files ';
l{end+1}='';

CompileStr='';
for i=1:length(l)
  CompileStr=[CompileStr l{i} '\n'];
end
clear l;

fid = fopen(['./' name '/compile_all.sh'],'w');
fprintf(fid,CompileStr);
fclose(fid);
if ~ispc
unix(['chmod 744 ./' name '/compile_all.sh']);
else
fid = fopen(['./' name '/windows_support.sh'],'a');
fprintf(fid,['chmod 744 ./compile_all.sh\n']);
fclose(fid);    
end

fid = fopen(['./' name '/compile_all_py.sh' ],'w');
fprintf(fid,CompileStr);
fclose(fid);
if ~ispc
unix(['chmod 744 ./' name '/compile_all_py.sh']);
else
fid = fopen(['./' name '/windows_support.sh'],'a');
fprintf(fid,['chmod 744 ./compile_all_py.sh\n']);
fclose(fid);    
end

% Lunch script for the cluster

l{1}='#!/bin/bash';
l{end+1}='# Made by guide_bot';
l{end+1}='# Will launch all the batch files (works on DMSC cluster)';
l{end+1}='';

LaunchStr='';
for i=1:length(l)
  LaunchStr=[LaunchStr l{i} '\n'];
end
clear l;

fid = fopen(['./' name '/launch_all.sh' ],'w');
fprintf(fid,LaunchStr);
fclose(fid);
if ~ispc
unix(['chmod 744 ./' name '/launch_all.sh']);
else
fid = fopen(['./' name '/windows_support.sh'],'a');
fprintf(fid,['chmod 744 ./launch_all.sh\n']);
fclose(fid);    
end


% Suggested reruns

l{1}='#!/bin/bash';
l{end+1}='# Made by guide_bot';
l{end+1}='# Will relaunch the runs which failed for some reason';
l{end+1}='# Commands will be added to this script if something fails';
l{end+1}='# When this runs, it will make a new suggested_reruns-fail for further fails';
l{end+1}='';
l{end+1}='Num=$(ls | grep suggested_reruns-fails | wc -l)';
l{end+1}='let Num+=1';
l{end+1}='NumStr=$(echo $Num)';
l{end+1}='#echo "testfile$NumStr.txt"';
l{end+1}='head -13 suggested_reruns-fails.sh > "suggested_reruns-fails$NumStr.sh"';
l{end+1}='chmod 744 "suggested_reruns-fails$NumStr.sh"';
l{end+1}='';

rerunStr='';
for i=1:length(l)
  rerunStr=[rerunStr l{i} '\n'];
end
clear l;

fid = fopen(['./' name '/suggested_reruns-fails.sh'],'w');
fprintf(fid,rerunStr);
fclose(fid);
if ~ispc
unix(['chmod 744 ./' name '/suggested_reruns-fails.sh']);
else
fid = fopen(['./' name '/windows_support.sh'],'a');
fprintf(fid,['chmod 744 ./suggested_reruns-fails.sh\n']);
fclose(fid);    
end

end



