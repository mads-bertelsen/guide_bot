function brilliance_ref(McStasStr,source_for_bt_refference,t,input,filename,Project_name,requirements,demands,endinput,enddeclare,NumSnaps,scan,options_general, modules, modulelist, defaults, globalinfo)

%%%%%%%%%%%%%%%%%%%%% INPUTS MODIFIED BY LELAND %%%%%%%%%%%%%%%%%%%%
% Added modules, modulelist, defaults and globalinfo to the inputs.
% These are needed to calculate wavelength dependent snapwidths when a
% monochromator is present. See my other modifications below for more
% details.

% Original function definition:
% function brilliance_ref(McStasStr,source_for_bt_refference,t,input,filename,Project_name,requirements,demands,endinput,enddeclare,NumSnaps,scan,options_general)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for scanj=1:scan.length(2)
for scani=1:scan.length(1)
    %if scan.length(scanj)==1;
    if scan.length(1)==1 && scan.length(2)==1;
      scanname='';    
    else
      if scan.dimension == 2;
        scanname=['_' num2str(scani) scan.names{1} '_' num2str(scanj) scan.names{2}];
      else
        scanname=['_' num2str(scani) scan.names{1}];
      end
    end
    
    % what if a 2d scan is one demand and one requirement?
    demandnames=fieldnames(demands);
    for i=1:length(demandnames)
        if sum(ismember(demandnames{i},scan.names))>0.5
            if strcmp(demandnames{i},scan.names{1})
                demandscan.(demandnames{i})=demands.(demandnames{i})(scani);
            else
                demandscan.(demandnames{i})=demands.(demandnames{i})(scanj);
            end
        else
            demandscan.(demandnames{i})=demands.(demandnames{i});
        end
    end
    reqnames=fieldnames(requirements);
    for i=1:length(reqnames)
        if sum(ismember(reqnames{i},scan.names))>0.5
            if strcmp(reqnames{i},scan.names{1})
                reqscan.(reqnames{i})=requirements.(reqnames{i})(scani);
            else
                reqscan.(reqnames{i})=requirements.(reqnames{i})(scanj);
            end
        else
            reqscan.(reqnames{i})=requirements.(reqnames{i});
        end
    end
    % I could just change reqscan / demand scan with the locked parameter
    % in order to support locked scans here.
    % locked as in, scan sample_size_x and y together for only squares
    if scan.locked_mode == 1
       for i = 1:length(demandnames)
            members = ismember(scan.locked.names,demandnames{i});
            if sum(members) > 0.5
                num_locked = find(members);
                if strcmp(scan.locked.names_parent{num_locked},scan.names{1})
                  demandscan.(demandnames{i}) = scan.locked.values{num_locked}(scani);
                else
                  demandscan.(demandnames{i}) = scan.locked.values{num_locked}(scanj);
                end
            end
       end
       for i = 1:length(reqnames)
            members = ismember(scan.locked.names,reqnames{i});
            if sum(members) > 0.5
                num_locked = find(members);
                if strcmp(scan.locked.names_parent{num_locked},scan.names{1})
                  reqscan.(reqnames{i}) = scan.locked.values{num_locked}(scani);
                else
                  reqscan.(reqnames{i}) = scan.locked.values{num_locked}(scanj);
                end
            end
       end
    end


run=0;
if (exist(['./' Project_name '/brilliance_refference'])>0.5 && exist(['./' Project_name '/brilliance_refference/input_used' scanname '.txt'])>0.5)
% check if this have allready been written for this exact run

wrong=0;
% demands
demandnames=fieldnames(demands);
reqnames=fieldnames(requirements);
fid=fopen(['./' Project_name '/brilliance_refference/input_used' scanname '.txt']);
testd=textscan(fid,'%s = %f',1,'HeaderLines',3);
if ~(strcmp(testd{1},demandnames{1}) && testd{2}==demandscan.(demandnames{1})); wrong=1; end;


for i=1:length(fieldnames(demands))-1;
testd=textscan(fid,'%s = %f',1);
if ~(strcmp(testd{1},demandnames{i+1}) && testd{2}==demandscan.(demandnames{i+1})); wrong=1; end;
end;

% requirements
testr=textscan(fid,'%s = %f',1,'HeaderLines',2);
if ~(strcmp(testr{1},'source'))
if ~(strcmp(testr{1},reqnames{1}) && testr{2}==reqscan.(reqnames{1})); wrong=1; end;
end
for i=1:length(fieldnames(requirements))-1 
testr=textscan(fid,'%s = %f',1);
if ~(strcmp(testr{1},'source'))
if ~(strcmp(testr{1},reqnames{i+1}) && testr{2}==reqscan.(reqnames{i+1})); wrong=1; end;
end
end
%save('debug.mat')
if wrong==1; run=0; end;

else
    run=1;
end

if run==1;
disp('Overwriting brilliance files')
% Make the directory and copy files if needed
[status,message,messageid]=mkdir(['./' Project_name '/brilliance_refference']);
[status,message,messageid]=copyfile([options_general.guide_bot_path '/components/*comp'],['./' Project_name '/brilliance_refference']);
[status,message,messageid]=copyfile([options_general.guide_bot_path '/components/*c'],['./' Project_name '/brilliance_refference']);

% write a file describing what input this file describes:

l{1}=filename;
l{end+1}=input;

l{end+1}='demands:';
% demands
names=fieldnames(demandscan);
for i=1:length(names)
   l{end+1}=[names{i} ' = ' num2str(demandscan.(names{i})) ]; 
end
clear names;

l{end+1}='requirements:';
% requirements
names=fieldnames(requirements);
for i=1:length(names)
   l{end+1}=[names{i} ' = ' num2str(reqscan.(names{i})) ]; 
end

input_used='';
for i=1:length(l)
    input_used=[input_used l{i} '\n'];
end

fid = fopen(['./' Project_name '/brilliance_refference/input_used' scanname '.txt' ],'w');
fprintf(fid,input_used);
fclose(fid);
clear l;

% write the McStas file
% A line in the source is changed
%source_for_bt_refference{6}='dist = guide_start, focus_xw = 6*sizeX, focus_yh = 6*sizeY,';

l{1}='COMPONENT Origin = Progress_bar()';
l{2}=' AT (0,0,0) ABSOLUTE';
l{3}='';
l{4}='COMPONENT source = Source_div(';
l{5}='yheight = mod_y, xwidth = mod_x,';
l{6}=['focus_aw = divreq_x*6.1, focus_ah = divreq_y*6.1,'];
l{7}='lambda0 = Lambda0, dlambda = dLambda, flux = 100, gauss = 0)';
l{8}='AT (0, 0, 0) RELATIVE Origin';
l{9}='ROTATED (0,0,0) RELATIVE Origin';
l{10}='';
l{11}='COMPONENT StartOfGuide = Arm()';
l{12}='AT (0,0,guide_start) RELATIVE source';

source_for_bt_refference=l;

brilliance.trace='';
l=source_for_bt_refference;
tracestring='';
for i=1:length(l)
    tracestring=[tracestring l{i} '\n'];
end
brilliance.trace=[tracestring '\n\n' brilliance.trace '\n\n'];
clear l;

% using input t which is the detector bank
tracestring_analyze='';
for i=1:length(t)
    tracestring_analyze=[tracestring_analyze t{i} '\n'];
end
brilliance.trace=['TRACE\n' brilliance.trace tracestring_analyze '\n\n'];
clear l;

brilliance.declarestring='\nDECLARE\n%%{\n';
for i=1:enddeclare 
    brilliance.declarestring = [brilliance.declarestring 'double ' McStasStr.declare{i} ';\n' ];
end
for i=1:length(McStasStr.declareint) 
    brilliance.declarestring = [brilliance.declarestring 'int ' McStasStr.declareint{i} ';\n' ];
end
brilliance.declarestring = [brilliance.declarestring 'double endx1;\n'];
brilliance.declarestring = [brilliance.declarestring 'double endy1;\n'];
brilliance.declarestring = [brilliance.declarestring '%%}\n'];

% the things done to avoid accesing inputvalue(i) where i is to large are
% stupid. Instead, check if inputvalue(imax) exists, if not, create it with
% a zero to indicate that it does not have an input value.
brilliance.inputstring=['DEFINE INSTRUMENT brilliance_ref(\n'];
for i=1:endinput
    if (i<=numel(McStasStr.inputvalue))
        if (McStasStr.inputvalue(i)~=0)
        brilliance.inputstring = [brilliance.inputstring McStasStr.input{i} '=' num2str(McStasStr.inputvalue(i)) ',\n' ];
        else
        brilliance.inputstring = [brilliance.inputstring McStasStr.input{i} ',\n' ];
        end
    else
    brilliance.inputstring = [brilliance.inputstring McStasStr.input{i} ',\n' ];
    end
end
brilliance.inputstring = [brilliance.inputstring 'guide_start = 0.001) \n']; 


l{1}='dLambda = 0.5*(WaveMax - WaveMin);';
l{end+1}='Lambda0 = dLambda+WaveMin;';
l{end+1}='';
l{end+1}='var_divreq_x = divreq_x;';
l{end+1}='var_divreq_y = divreq_y;';
l{end+1}='';
l{end+1}='endx1 = mod_x;';
l{end+1}='endy1 = mod_y;';
l{end+1}='';
l{end+1}='u=1e-6;';
l{end+1}='';

% these are now input
%l{end+1}=['mod_x=' num2str(requirements.moderator_size_x) ';'];
%l{end+1}=['mod_y=' num2str(requirements.moderator_size_y) ';'];

initstring='';
for i=1:length(l)
    initstring=[initstring l{i} '\n'];
end
clear l;
brilliance.initialize=['\nINITIALIZE\n%%{\n' initstring '\n%%}\n'];
        
        
fid = fopen(['./' Project_name '/brilliance_refference/brilliance.instr'], 'w');
write=[brilliance.inputstring brilliance.declarestring brilliance.initialize brilliance.trace '\nEND'];
fprintf(fid,write);
fclose(fid);


% Write ifit file

clear l;

%l{1}=['instrument_name=''brilliance'';'];
l{1}=['filenamebrill=''brilliance'';'];
%l{end+1}=['inputstring=''' input ''';'];
%l{end+1}='';
if scan.mode==0
  l{end+1}=['scanname='''';'];
else
  l{end+1}=['scanname=''' scanname ''';'];
end
l{end+1}='bpath=pwd;';
l{end+1}='i=0;done=''1''';
l{end+1}='%% creates a new filename';
l{end+1}='while(str2num(done)==1); i=1+i; [a,done]=unix([''test -d \'' bpath ''/'' filenamebrill scanname num2str(i) '' && echo "1" || echo "0"'']); end;';
l{end+1}='filenamebrill=[filenamebrill scanname num2str(i)];';
l{end+1}='';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% LELAND MODIFICATION %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Expand to pc case
l{end+1} = 'if ispc==1';
l{end+1} = '    select=1;';
l{end+1} = 'else';
l{end+1}='    %% checks if the script is on the ESS cluster';
l{end+1}=['    [a,check]=unix(''test -f ' options_general.cluster_path 'clusterid && echo "2" || echo "1"'')'];
l{end+1}='    select=str2num(check);';
l{end+1} = 'end';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% END OF LELAND MODIFICATION %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Old code
%l{end+1}='%% checks if the script is on the ESS cluster';
%l{end+1}=['[a,check]=unix(''test -f ' options_general.cluster_path 'clusterid && echo "2" || echo "1"'')'];
%l{end+1}='select=str2num(check);';

l{end+1}='';
l{end+1}='%% primary parameters';

iFitStr='';
for i=1:length(l)
    iFitStr=[iFitStr l{i} '\n'];
end
clear l

% Making the optimize array the same size as the input array.
optimize=zeros(1,length(McStasStr.input));
for i=1:length(McStasStr.optimize)
    optimize(i)=McStasStr.optimize(i);
end
inputvalue=zeros(1,length(McStasStr.input));
for i=1:length(McStasStr.inputvalue)
    inputvalue(i)=McStasStr.inputvalue(i);
end

% Writing the input parameters to the ifit string
for i=1:endinput
    if scan.mode==0
        if strcmp(McStasStr.input{i},'sample_dist')
            iFitStr=[iFitStr 'p.' McStasStr.input{i} '=''' num2str(0.001) ''';\n']; %? Do i need a value for everything?
        elseif (strcmp(McStasStr.input{i},'mod_x') || strcmp(McStasStr.input{i},'mod_y'))
            if strcmp(McStasStr.input{i},'mod_x')
                compare_to = 'sizeX';
            else
                compare_to = 'sizeY'; 
            end
            
            for test = 1:length(McStasStr.input)
               if strcmp(McStasStr.input{test},compare_to)
                  compare_index = test; 
               end
            end
            
            if inputvalue(i) < 4.5*inputvalue(compare_index)
                iFitStr=[iFitStr 'p.' McStasStr.input{i} '=''' num2str(4.5*inputvalue(compare_index)) ''';\n'];  
            else
               % old code    
               iFitStr=[iFitStr 'p.' McStasStr.input{i} '=''' num2str(inputvalue(i)) ''';\n'];
            end
        else
            iFitStr=[iFitStr 'p.' McStasStr.input{i} '=''' num2str(inputvalue(i)) ''';\n']; %? Do i need a value for everything?
        end
    else
%        match=0;
%         for stest=1:length(scan.names)
%            if strcmp(scan.namesmcstas{stest},McStasStr.input{i})
%                 match=stest;
%            end
%         end
%         if match==0
%             if strcmp(McStasStr.input{i},'sample_dist')
%                 iFitStr=[iFitStr 'p.' McStasStr.input{i} '=''' num2str(0.1) ''';\n']; %? Do i need a value for everything?
%             else
%                 iFitStr=[iFitStr 'p.' McStasStr.input{i} '=''' num2str(inputvalue(i)) ''';\n']; %? Do i need a value for everything?
%             end
%         else
%             %save('debug.mat')
%             if (strcmp(scan.names{stest},'Hsize') || strcmp(scan.names{stest},'Vsize'))
%               iFitStr=[iFitStr 'p.' McStasStr.input{i} '=''' num2str(demands.(scan.names{stest})(scani)/100) ''';\n'];    
%             elseif (strcmp(scan.names{stest},'moderator_size_x') || strcmp(scan.names{stest},'moderator_size_y'))
%               iFitStr=[iFitStr 'p.' McStasStr.input{i} '=''' num2str(requirements.(scan.names{stest})(scani)) ''';\n'];  
%             else
%               iFitStr=[iFitStr 'p.' McStasStr.input{i} '=''' num2str(demands.(scan.names{stest})(scani)) ''';\n'];
%             end
%         end
        % Above, 1 version
%         if strcmp(scan.namesmcstas{1},McStasStr.input{i})
%                % Add to the ifit string in the three cases
%                 identifier = 1;
%                 if (strcmp(scan.names{identifier},'Hsize') || strcmp(scan.names{identifier},'Vsize'))
%                   iFitStr=[iFitStr 'p.' McStasStr.input{i} '=''' num2str(demands.(scan.names{identifier})(scani)/100) ''';\n'];    
%                 elseif (strcmp(scan.names{identifier},'moderator_size_x') || strcmp(scan.names{identifier},'moderator_size_y'))
%                   iFitStr=[iFitStr 'p.' McStasStr.input{i} '=''' num2str(requirements.(scan.names{identifier})(scani)) ''';\n'];  
%                 else
%                   iFitStr=[iFitStr 'p.' McStasStr.input{i} '=''' num2str(demands.(scan.names{identifier})(scani)) ''';\n'];
%                 end
%         elseif scan.dimension == 2 && strcmp(scan.namesmcstas{2},McStasStr.input{i})
%                % There is another layer to check for, if yes, add
%                % If this fails to, either because there is none or
%                % dimensionality is to low, write normally
%                identifier = 2;
%                 if (strcmp(scan.names{identifier},'Hsize') || strcmp(scan.names{identifier},'Vsize'))
%                   iFitStr=[iFitStr 'p.' McStasStr.input{i} '=''' num2str(demands.(scan.names{identifier})(scanj)/100) ''';\n'];    
%                 elseif (strcmp(scan.names{identifier},'moderator_size_x') || strcmp(scan.names{identifier},'moderator_size_y'))
%                   iFitStr=[iFitStr 'p.' McStasStr.input{i} '=''' num2str(requirements.(scan.names{identifier})(scanj)) ''';\n'];  
%                 else
%                   iFitStr=[iFitStr 'p.' McStasStr.input{i} '=''' num2str(demands.(scan.names{identifier})(scanj)) ''';\n'];
%                 end
%         elseif strcmp(McStasStr.input{i},'sample_dist')
%                 iFitStr=[iFitStr 'p.' McStasStr.input{i} '=''' num2str(0.001) ''';\n']; %? Do i need a value for everything?
%         else
%                % write normally
%                iFitStr=[iFitStr 'p.' McStasStr.input{i} '=''' num2str(inputvalue(i)) ''';\n'];
%         end
        % Above, 2 version, below current version
        
        clear members_scan;
            members_scan = ismember(scan.namesmcstas,McStasStr.input{i});
            
            locked_parameter = 0;
            if scan.locked_mode == 1
                members_child  = ismember(scan.locked.namesmcstas,McStasStr.input{i});
                if sum(members_child) > 0.5
                    child_index = find(members_child);
                    locked_parameter = 1;
                end
            end
            
            if sum(members_scan)>0.5 || locked_parameter == 1;
                if locked_parameter == 1
                    % if locked, it can not also be scanned.
                    variable_name = scan.locked.names{child_index};
                else
                    variable_name = scan.names{find(members_scan,1)};
                end    
                
                % Needs attention.
                if (strcmp(variable_name,'Hsize') || strcmp(variable_name,'Vsize'))
                  iFitStr=[iFitStr 'p.' McStasStr.input{i} '=''' num2str(demandscan.(variable_name)/100) ''';\n'];    
                elseif (strcmp(variable_name,'moderator_size_x') || strcmp(variable_name,'moderator_size_y'))
                  % Need to check for moderator size, if below 4.5xSize, increase to 4.5xSize.
                  if strcmp(variable_name,'moderator_size_x')
                     compare_to = 'Hsize';
                  else
                     compare_to = 'Vsize'; 
                  end
                  
                  if reqscan.(variable_name) < 4.5*demandscan.(compare_to)/100
                    iFitStr=[iFitStr 'p.' McStasStr.input{i} '=''' num2str(4.5*demandscan.(compare_to)/100) ''';\n'];  
                  else
                    % old code    
                    iFitStr=[iFitStr 'p.' McStasStr.input{i} '=''' num2str(reqscan.(variable_name)) ''';\n'];  
                  end
                else
                  iFitStr=[iFitStr 'p.' McStasStr.input{i} '=''' num2str(demandscan.(variable_name)) ''';\n'];
                end
            elseif strcmp(McStasStr.input{i},'sample_dist')
                 iFitStr=[iFitStr 'p.' McStasStr.input{i} '=''' num2str(0.001) ''';\n']; %? Do i need a value for everything?
            elseif (strcmp(McStasStr.input{i},'mod_x') || strcmp(McStasStr.input{i},'mod_y'))
                if strcmp(McStasStr.input{i},'mod_x')
                    compare_to = 'sizeX';
                else
                    compare_to = 'sizeY'; 
                end

                for test = 1:length(McStasStr.input)
                   if strcmp(McStasStr.input{test},compare_to)
                      compare_index = test; 
                   end
                end

                if inputvalue(i) < 4.5*inputvalue(compare_index)
                    iFitStr=[iFitStr 'p.' McStasStr.input{i} '=''' num2str(4.5*inputvalue(compare_index)) ''';\n'];  
                else
                   % old code    
                   iFitStr=[iFitStr 'p.' McStasStr.input{i} '=''' num2str(inputvalue(i)) ''';\n'];
                end
            else
               % write normally
               iFitStr=[iFitStr 'p.' McStasStr.input{i} '=''' num2str(inputvalue(i)) ''';\n'];
               % Need to check for moderator size, if below 4.5xSize, increase to 4.5xSize.
            end
    end
end
iFitStr=[iFitStr 'p.guide_start = ''0.001'';\n']; %? Do i need a value for everything?



% calculate resonable centers of wavelength snapshots.
% Always 5 snapshots

deltaL=demandscan.WaveLmax-demandscan.WaveLmin;
List=0:(NumSnaps-1);
centerscalc=demandscan.WaveLmin+deltaL/(NumSnaps-1)*List;

wavec=zeros(NumSnaps,1);
wavec(1)=demandscan.WaveLmin;
wavec(end)=demandscan.WaveLmax;

if (deltaL>6)
   for i=2:(NumSnaps-1)
      wavec(i)=round(centerscalc(i));
   end
end
if (deltaL<=6 && deltaL>=2)
    for i=2:(NumSnaps-1)
      wavec(i)=round(centerscalc(i)*10)*0.1; 
   end
end
if (deltaL<2)
    for i=2:(NumSnaps-1)
      wavec(i)=round(centerscalc(i)*100)*0.01; 
   end
end
% this must correspond with the code being excecuted for each geometry

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% LELAND2 MODIFICATION (Intermediate Brilliance 1/2 %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rewrite brilliance_ifit.m as a function with inputs defining a brilliance
% window to sample over. Note that the original brilliance_ifit.m is still used
% by guide_bot to do the optimization. This functional version is used for
% post optimization analysis.
l_IM = 'function IntermediateBrilliance_ref_ifit(Hdiv, Vdiv, Hsize, Vsize, brillname)\n';
IM_Str = [l_IM iFitStr '\n'];
clear l_IM
l_IM{1} = '';
l_IM{end+1} = '%% Change to Intermediate Brilliance Window';
l_IM{end+1} = ['p.sizeX = num2str(Hsize);'];
l_IM{end+1} = ['p.sizeY = num2str(Vsize);'];
l_IM{end+1} = ['p.divreq_x = num2str(Hdiv) ;'];
l_IM{end+1} = ['p.divreq_y = num2str(Vdiv) ;'];
for i=1:length(l_IM)
    IM_Str=[IM_Str l_IM{i} '\n'];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% END OF LELAND2 MODIFICATION (Intermediate Brilliance 1/2 %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

l{1}='';
l{end+1}='%%options_single_home.compile=0;';
l{end+1}='options_single_home.gravitation=1;';
l{end+1}=['options_single_home.ncount=' num2str(options_general.ncount_multiplier_analyze) '*5e7;'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% LELAND MODIFICATION %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Include machine list and fix mispelling (note that there are a few other
% mispelling of single as signle both here and in mcstas_bot but I think it
% is best that the developer sort out how they want to handle this.

% Include machine list
if isfield(options_general, 'machines')
    l{end+1} = ['options_single_home.machines =' '''' options_general.machines '''' ';'];
end

% Modify MPI line
% Original line (Note it was also misspelled):
% l{end+1}='options_signle_home.mpi=7;';
% New line:
l{end+1}=['options_single_home.mpi=' num2str(options_general.mpi) ';'];

% Fix misspelling
% Original line:
% l{end+1}='options_signle_home.dir=[bpath ''/'' filenamebrill];';
% New line:
l{end+1}='options_single_home.dir=[bpath ''/'' filenamebrill];';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% END OF LELAND MODIFICATION %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%l{end+1}='options_single_home.mpi=7;';
%l{end+1}='options_single_home.dir=[bpath ''/'' filenamebrill];';
l{end+1}='';
l{end+1}='';
l{end+1}='%% rest of options for cluster';
l{end+1}='if (select==2); load NUMCORES.DAT; warning off; else NUMCORES=0; end;';
l{end+1}='';
l{end+1}='options_single_cluster.mpi=NUMCORES;';
l{end+1}='options_single_cluster.Display=[];';
l{end+1}='';
l{end+1}='%% general opptions';
l{end+1}='options_single_cluster.mpi=NUMCORES;';
l{end+1}='options_single_cluster.Display=[];';
l{end+1}='options_single_cluster.dir=filenamebrill;';
l{end+1}='options_single_cluster.gravitation=1;';
l{end+1}='options_single_cluster.compile=0;';
l{end+1}=['options_single_cluster.ncount=' num2str(options_general.ncount_multiplier_analyze) '*5e8;'];
l{end+1}='';
l{end+1}='options_single={options_single_home options_single_cluster};';
l{end+1}='';
l{end+1}='';
l{end+1}='%% Run for fom wavelengths';
l{end+1}='monitor_fom_ref=mcstas(''brilliance.instr'',p,options_single{select});';
l{end+1}='';
strings='wavecenters=[ ';
for i=1:NumSnaps
    strings=[strings num2str(wavec(i)) ' '];
end
strings=[strings '];'];
l{end+1}=strings;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% LELAND MODIFICATION %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Need to vectorize snapwidth because of Monochromator case.
% See my other modifications further down for more info.

% Original Line:
% l{end+1}='snapwidth=0.01;';

% New Line:
l{end+1}=['snapwidth=0.01*ones(1,' num2str(NumSnaps) ');'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% END OF LELAND MODIFICATION %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

l{end+1}='';
l{end+1}='p.WaveMin=''0.1'';';
%if (demandscan.WaveLmax<3)
%l{end+1}=['p.WaveMax=''' num2str(demandscan.WaveLmax*3) ''';'];
%else
%l{end+1}=['p.WaveMax=''' num2str(demandscan.WaveLmax*2) ''';'];
%end
if max(ismember(fieldnames(options_general),'max_wavelength_investigated_multiplier')) == 0
    if (demands.WaveLmax<3)
    l{end+1}=['p.WaveMax=''' num2str(demandscan.WaveLmax*3) ''';'];
    else
    l{end+1}=['p.WaveMax=''' num2str(demandscan.WaveLmax*2) ''';'];
    end
else 
    l{end+1}=['p.WaveMax=''' num2str(demandscan.WaveLmax*options_general.max_wavelength_investigated_multiplier) ''';'];
end
l{end+1}='MaxWB=str2num(p.WaveMax);';
l{end+1}='';
l{end+1}='%% Run for all wavelengths';
l{end+1}=['options_single{select}.dir=[filenamebrill ''waveLarge''];'];
l{end+1}='monitor_ALLW_ref=mcstas(''brilliance.instr'',p,options_single{select});';
l{end+1}='';
l{end+1}='%% Run for individual wavelength snapshots';
l{end+1}='options_single_cluster.ncount=1e8;';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% LELAND MODIFICATION   %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% When no monochromator is present mcstas_bot sets the wavelength range of
% the wave1 through wave5 simulations to size snapwidth (which is extremely
% small). However, if you have a monochromator, then the wave1 through
% wave5 simulations need to be sampled over the wavelength range 
% corresponding to the acceptance of the monochromator.

% Determine if a Monochromator in reflection mode is present
for i = 1:numel(modulelist)
    if strcmp(modules{modulelist(i)},'M') % Check for monochromator
        monoopts = Parse_options(globalinfo.options{i});
        if strcmp(monoopts.BeamDir, 'reflect') % Check if monochromator is in reflection mode
            % Conditions met, need to calculate new snapwidths due to monochromator
            
            % Caluclation requires guide m-value
            if isfield(defaults, 'm')
                mval_temp = defaults.m;
            else
                disp('WARNING: Energy binning uses defaults.m value which is currently not defined, using m=2 for determining bin sizes')
                mval_temp = 2;
            end
            
            % Calculation requires monochromator mosaic
            if isfield(monoopts, 'mos')
                mos_temp = monoopts.mos;
            elseif isfield(monoopts, 'Hmos')
                mos_temp = monoopts.Hmos;
            elseif isfield(defaults, 'mos')
                mos_temp = defaults.mos;
            elseif isfield(defaults, 'Hmos')
                mos_temp = defaults.Hmos;
            else
                disp('WARNING: Monochromator mosaic is not defined, get ready for a crash')
            end
            
            % Calculation requires monochromator D-spacing
            if isfield(monoopts, 'dM')
                dM_temp = monoopts.dM;
            elseif isfield(defaults, 'dM')
                dM_temp = defaults.dM;
            else
                disp('WARNING: Monocromator D-spacing is not defined, get ready for a crash')
            end
                
            % Calculate snapwidth for each wavelength center (wavec)
            l{end+1} = '';
            l{end+1} = '%% Changing snapwidths to monochromator integration width';
            for i = 1:NumSnaps
                % Note that I hardwire in a 3*FWHM for the snapwidth.
                % This is different than the width used during
                % optimization, which is set using the Ebin monochromator option.
                snapwidth_calc = 3*2*dM_temp*(mos_temp/60*pi/180 + mval_temp*wavec(i)*0.0024)/2*cos(asin(wavec(i)/2/dM_temp));
                l{end+1} = ['snapwidth(' num2str(i) ') = ' num2str(snapwidth_calc) ';'];
            end
        end
    end
    clear monoopts dM_temp mval_temp mos_temp
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% END OF LELAND MODIFICATION   %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:NumSnaps
l{end+1}='';    
l{end+1}=['options_single{select}.dir=[filenamebrill ''wave' num2str(i) '''];'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% LELAND MODIFICATION %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Need to vectorize snapwidth because of Monochromator case.
% See my other modification above this for more info.

% Original 2 Lines:
% l{end+1}=['p.WaveMin=num2str(wavecenters(' num2str(i) ') - snapwidth*0.5);'];
% l{end+1}=['p.WaveMax=num2str(wavecenters(' num2str(i) ') + snapwidth*0.5);'];

% New Lines
l{end+1}=['p.WaveMin=num2str(wavecenters(' num2str(i) ') - snapwidth(' num2str(i) ')*0.5);'];
l{end+1}=['p.WaveMax=num2str(wavecenters(' num2str(i) ') + snapwidth(' num2str(i) ')*0.5);'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% END OF LELAND MODIFICATION %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
l{end+1}=['monitor_W_ref.wave' num2str(i) '=mcstas(''brilliance.instr'',p,options_single{select});'];
l{end+1}='';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% LELAND2 MODIFICATION (Intermediate Brilliance 2/2 %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:length(l)
    IM_Str=[IM_Str l{i} '\n'];
end

l_IM{1} = 'save([''../output/brill_ref/'' brillname] ,''monitor_W_ref'',''monitor_ALLW_ref'',''monitor_fom_ref'')';
IM_Str=[IM_Str l_IM{1} '\n'];
if isfield(options_general,'Intermediate_Brilliance')
    fid = fopen(['./' Project_name '/brilliance_refference/IntermediateBrilliance_ref_ifit.m'], 'w');
    fprintf(fid,IM_Str);
    fclose(fid);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% END OF LELAND2 MODIFICATION (Intermediate Brilliance 2/2 %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

l{end+1}='save([filenamebrill ''.mat'']);';
l{end+1}=['save(''../output/brill_ref/brilliance_ref' scanname '.mat'',''monitor_W_ref'',''monitor_ALLW_ref'',''monitor_fom_ref'')'];
% should only save the nessecary data.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% LELAND2 MODIFICATION (Project List 1)%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First check if this is part of a project list or a single run.
if isfield(options_general, 'projectlist')
    if strcmp(options_general.projectlist, 'single')
        check_single_run = 1;
    else
        check_single_run = 0;
    end
else
    check_single_run = 1;
end

l{end+1} = ['singlerun = ' num2str(check_single_run) ';'];
l{end+1} = 'if singlerun == 0';
l{end+1} = '    load([fileparts(bpath) ''/runs_left.mat'']);';
l{end+1} = '    if runs_left == 1';
l{end+1} = '        runs_left = runs_left - 1;';
l{end+1} = '        save([fileparts(bpath) ''/runs_left.mat''], ''runs_left'')';
l{end+1} = '        cd ..';
l{end+1} = '        cd ..';
l{end+1} = '        addpath(''guide_bot_source'')';
l{end+1} = '        RunProjectList()';
l{end+1} = '    else';
l{end+1} = '        runs_left = runs_left - 1;';
l{end+1} = '        save([fileparts(bpath) ''/runs_left.mat''], ''runs_left'')';
l{end+1} = '    end';
l{end+1} = 'end';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% END OF LELAND2 MODIFICATION (Project List 1)%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for i=1:length(l)
    iFitStr=[iFitStr l{i} '\n'];
end
clear l


fid = fopen(['./' Project_name '/brilliance_refference/brilliance_ifit' scanname '.m'], 'w');
write=iFitStr;
fprintf(fid,write);
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% LELAND2 MODIFICATION (clustering 1/2) %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Include NCNR version of batch script

if strcmp(options_general.cluster, 'NCNR')
        % This is how the NCNR does it on the rocks cluster
    gb_loc = fileparts(mfilename('fullpath'));
    gb_loc = gb_loc(1:end-length('guide_bot_source'));
    ifit_loc = which('mcstas.m');
    ifit_loc = ifit_loc(1:end-length('/Applications/McStas/mcstas.m'));
    
    % fprintf does not like pathsep \ since it is reserved as an escape
    % sequence. Thus, need to change to \\ pathsep.
    if ispc
        gb_loc_new = gb_loc(1);
        ifit_loc_new = ifit_loc(1);
        for i = 2:numel(gb_loc)
            if strcmp(gb_loc(i),'\')
                gb_loc_new = [gb_loc_new '\\'];
            else
                gb_loc_new = [gb_loc_new gb_loc(i)];
            end
        end
        for i = 2:numel(ifit_loc)
            if strcmp(ifit_loc(i),'\')
                ifit_loc_new = [ifit_loc_new '\\'];
            else
                ifit_loc_new = [ifit_loc_new ifit_loc(i)];
            end
        end
        gb_loc = gb_loc_new;
        ifit_loc = ifit_loc_new;
    end
        
    
    l{1}='#!/bin/bash';
    l{end+1}='#';
    l{end+1}='# Made by mcstas_bot';
    l{end+1}='';
    l{end+1}='# Gives a set of matlab commands to iFit. NUMCORES.DAT read from matlab script';
    l{end+1}=['nohup matlab -nosplash -nodisplay -nodesktop -r "addpath(genpath(''' ifit_loc ''')); cd ' gb_loc Project_name '/brilliance_refference/; tic; try; run ' gb_loc Project_name '/brilliance_refference/run_brilliance_ifit' scanname '.m; catch; end; toc; quit" > matlab_output.log &'];
else
    % This is Mads code for doing it elsewhere
    l{1}='#!/bin/bash';
    l{end+1}='#';
    l{end+1}='# Made by mcstas_bot';
    l{end+1}='';
    l{end+1}=['#SBATCH --job-name=brilliance'];
    l{end+1}=['#SBATCH --error=err_brilliance.txt'];
    l{end+1}=['#SBATCH --output=out_brilliance.txt'];
    l{end+1}='#SBATCH --nodes=1-1';
    if strcmp(options_general.queue,'verylong')
        l{end+1}='#SBATCH --exclude r3n9b7';
    elseif strcmp(options_general.queue,'long')
        l{end+1}='#SBATCH --exclude r1n26,r1n19';
    end
    l{end+1}=['#SBATCH --partition ' options_general.queue ];
    l{end+1}=['#SBATCH --time ' options_general.time ];
    l{end+1}='# the --exclusive is needed when running OpenMPI';
    l{end+1}='# it will allocate 1x12 core per node';
    l{end+1}='#SBATCH --exclusive';
    l{end+1}='';
    l{end+1}='NUMCORES=`echo "$SLURM_NNODES 12 * p "| dc`';
    l{end+1}='echo $NUMCORES > NUMCORES.DAT';
    l{end+1}='';
    for jj = 1:length(options_general.modules_node)
        l{end+1}=['module load ' options_general.modules_node{jj}];
    end
    l{end+1}='';
    l{end+1}='# Gives a set of matlab commands to iFit. NUMCORES.DAT read from matlab script';
    l{end+1}=['cat run_brilliance_ifit' scanname '.m | ' options_general.cluster_path 'run_ifit.sh ' options_general.cluster_path 'MATLAB_Compiler_Runtime/v716/'];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% END OF LELAND2 MODIFICATION (clustering 1/2) %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


BatchStr='';
for i=1:length(l)
  BatchStr=[BatchStr l{i} '\n'];
end
clear l

fid = fopen(['./' Project_name '/brilliance_refference/brilliance' scanname '.batch'], 'w');
write=BatchStr;
fprintf(fid,write);
fclose(fid);

if ~strcmp(requirements.source,'Virtual_in')
fid = fopen(['./' Project_name '/launch_all.sh' ],'a');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% MODIFIED BY LELAND2 (clustering 2/2) %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Added NCNR modification which does not use SLURM
switch options_general.cluster
    case 'NCNR'
        fprintf(fid,['cd brilliance_refference\nsh brilliance' scanname '.batch\ncd ..\n']);
    otherwise
        fprintf(fid,['cd brilliance_refference\nsbatch brilliance' scanname '.batch\ncd ..\n']);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% END OF MODIFIED BY LELAND2 (clustering 2/2) %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fclose(fid);
end

fid = fopen(['./' Project_name '/brilliance_refference/run_brilliance_ifit' scanname '.m'], 'w');
fprintf(fid,['run brilliance_ifit' scanname '.m \nexit\n']);
fclose(fid);

if scani==1 && scanj==1 % This will still give to many compiles of brill in case of 2d scan
fid = fopen(['./' Project_name '/compile_all.sh'],'a');
fprintf(fid,['cd brilliance_refference\nmcrun -c -n 0 --mpi brilliance.instr \ncd ..\n']);
fclose(fid);

fid = fopen(['./' Project_name '/compile_all_py.sh' ],'a');
fprintf(fid,['cd brilliance_refference\nmcrun-py -c -n 0 --mpi=1 brilliance.instr \ncd ..\n']);
fclose(fid);
end

end 

end % ending scan loop

end

