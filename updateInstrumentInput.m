% Script that allows you to update the input of instrument files with the 
% values from parameter files from guide_bot output. This does for example 
% allow you to use mcgui with guide_bot output. Just insters the relevant
% filenames and run the script.
% Jonas Okkels Birk 2017

% insert your own file paths+names here:
parfileName='mypath/myfile.par';
instrFileName='mypath/myfile.par.instr'; 

%The rest should not need any input from you, just run the script.

% Read instrument file into cell array
fid = fopen(instrFileName,'r');
instrLines = 1;
tline = fgetl(fid);
instrFile{instrLines} = tline;
while ischar(tline)
    instrLines = instrLines+1;
    tline = fgetl(fid);
    instrFile{instrLines} = tline;
end
fclose(fid); %finished reading instrument file

% read Par file into cell array
fid = fopen(parfileName,'r');
parLines = 1;
tline = fgetl(fid);
parFile{parLines} = tline;
while ischar(tline)
    parLines = parLines+1;
    tline = fgetl(fid);
    parFile{parLines} = tline;
end
fclose(fid); %finished reading par file
missingLines=0;
fprintf('Looking through the Instrument file. No equivalents \nfor these par file lines was found:\n\n');
for plNr=1:parLines-1 % loop over par cell entries
     if length(parFile{plNr})>0
        CurrentText = strsplit(parFile{plNr},'=');
            
        flag=1;
        for lNr=1:instrLines % check if corresponding line in instrument cell exist
            if length(instrFile{1,lNr})>length(CurrentText{1})+3
                fprintf('%d    %d')
                if strcmp(instrFile{1,lNr}(1:length(CurrentText{1})),CurrentText{1})
                        if strcmp(CurrentText{1},'WaveMin') % skip if waveMin as this has supoptimal value in par file
                            instrFile{1,lNr}=instrFile{1,lNr};
                            flag=0;
                        elseif strcmp(CurrentText{1},'WaveMax') % skip if waveMax as this has supoptimal value in par file
                            flag=0;
                        else % update line in instrument cell array with values from par file
                            flag=0;
                        end
                    break
                end
            end   
        end
        if flag==1 % print lines that was not found in instrument cell
            fprintf([CurrentText{1} '=' CurrentText{2} '\n']);
            missingLines=missingLines+1;
        end
    end
end
if missingLines==0
    fprintf('All lines was found, the script is happy!\n')
end

% Write amended instruemnt cell array into file
fprintf('\nWriting to updated instrument file...')
fid = fopen(instrFileName, 'w');
for i = 1:numel(instrFile)
   if instrFile{i+1} == -1
       fprintf(fid,'%s', instrFile{i});
       break
   else
       fprintf(fid,'%s\n', instrFile{i});
   end
end
fclose(fid);
fprintf('aaand done!\n')
