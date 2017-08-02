IM_l{1}=['Project_name=''' Project_name ''';'];
IM_l{end+1}=['instrument_name = ''' filename '_' fnamesIM{j} ''';'];
IM_l{end+1}=['filename=''' fnamesIM{j} '_' filename ''';'];
IM_l{end+1}=['inputstring=''' input ''';'];
if scan.mode==0
  IM_l{end+1}=['scanname='''';'];
else
  IM_l{end+1}=['scanname=''' scanname ''';'];
end
IM_l{end+1}='';
IM_l{end+1}=['scani=' num2str(scani) ';'];
IM_l{end+1}=['scanj=' num2str(scanj) ';'];
IM_l{end+1}='cpath=pwd;';
IM_l{end+1}='i=0;done=''1'';';
IM_l{end+1}='%% creates a new filename';
IM_l{end+1}='while(str2num(done)==1); i=1+i; [a,done]=unix([''test -d \'' cpath ''/'' filename scanname num2str(i) '' && echo "1" || echo "0"'']); end;';
IM_l{end+1}='filename=[filename scanname num2str(i)];';
IM_l{end+1}='';
IM_l{end+1}='if ispc==1';
IM_l{end+1}=' select=1;';
IM_l{end+1}='else';
IM_l{end+1}=' %% checks if the script is on the ESS cluster';
IM_l{end+1}=['[a,check]=unix(''test -f ' options_general.cluster_path 'clusterid && echo "2" || echo "1"'')'];
IM_l{end+1}=' select=str2num(check);';
IM_l{end+1}='end';
IM_l{end+1}='';
IM_l{end+1}='%% primary parameters';

%%% Params here %%%
IM_l{end+1}=['alldata = load(''' filename scanname '1_all.mat'');'];
IM_l{end+1} = 'monitor = alldata.monitor;';
IM_l{end+1} = 'p = monitor.Data.Parameters;';                         

% Change Brilliance Window
if isfield(options_general,'Brilliance_Window')
    % Determine which component is being set up
    icheck = 1;
    if strcmp(fnamesIM{j},'FullInstrument')
        IM_module = 'FullInstrument';
        IM_index = '0';
    else
        while isempty(str2num(fnamesIM{j}(icheck))); icheck = icheck+1; end
        IM_module = fnamesIM{j}(1:icheck-1);
        IM_index = fnamesIM{j}(icheck:end);
    end
   
    % Set BrillHsize
    if isfield(options_general.Brilliance_Window,'Hsize')
        IM_Hsize = options_general.Brilliance_Window.Hsize;        
        if ischar(IM_Hsize)
            switch IM_Hsize
                % Need to set Hsize to startx of component
                case 'component'
                    % If M_module then startx is a stored calculation
                    if strcmp(IM_module, 'M_module')
                        IM_l{end+1} = ['IM_fh = fopen(''' filename '1_geometry.dat'');'];
                        IM_l{end+1} = 'while ~strcmp(fgetl(IM_fh), ''M''); end';
                        IM_l{end+1} = 'line = fgetl(IM_fh);';
                        IM_l{end+1} = 'fclose(IM_fh);';
                        IM_l{end+1} = 'Mono_Dims = str2num(line(5:end)); %%[Full Mono Width, Projected Mono Width, Mono Height]';
                        IM_l{end+1}= 'BrillHsize = Mono_Dims(2);';
                    elseif strcmp(IM_module, 'FullInstrument')
                        disp('options.Intermediate_Brilliance = 0 with a custom window of dimensions ''component'' will be set to FOM')
                        IM_l{end+1}= 'BrillHsize = p.sizeX;';
                    else
                        IM_l{end+1}= ['BrillHsize = p.startx' IM_index ';']; 
                    end
                case 'FOM'
                    IM_l{end+1}= 'BrillHsize = p.sizeX;';
                otherwise
                    disp('options.Brilliance_Window.Hsize is set to invalid name, setting to FOM by default');
                    IM_l{end+1}= 'BrillHsize = p.sizeX;';
            end
        else
            IM_l{end+1}= ['BrillHsize = ' num2str(IM_Hsize) '/100;'];
        end
    else
        disp('options.Brilliance_Window.Hsize is not defined, setting to FOM by default');
        IM_l{end+1}= 'BrillHsize = p.sizeX;';
    end
    
    % Set BrillVsize
    if isfield(options_general.Brilliance_Window,'Vsize')
        IM_Vsize = options_general.Brilliance_Window.Vsize;        
        if ischar(IM_Vsize)
            switch IM_Vsize
                % Need to set Vsize to starty of component
                case 'component'
                    % If M_module then starty is a stored calculation
                    if strcmp(IM_module, 'M_module')
                        IM_l{end+1} = ['IM_fh = fopen(''' filename '1_geometry.dat'');'];
                        IM_l{end+1} = 'while ~strcmp(fgetl(IM_fh), ''M''); end';
                        IM_l{end+1} = 'line = fgetl(IM_fh);';
                        IM_l{end+1} = 'fclose(IM_fh);';
                        IM_l{end+1} = 'Mono_Dims = str2num(line(5:end)); %%[Full Mono Width, Projected Mono Width, Mono Height]';
                        IM_l{end+1}= 'BrillVsize = Mono_Dims(3);';
                    elseif strcmp(IM_module, 'FullInstrument')
                        disp('options.Intermediate_Brilliance = 0 with a custom window of dimensions ''component'' will be set to FOM')
                        IM_l{end+1}= 'BrillVsize = p.sizeY;';
                    else
                        IM_l{end+1}= ['BrillVsize = p.starty' IM_index ';']; 
                    end
                case 'FOM'
                % Need to set Vsize to FOM of component
                    IM_l{end+1}= 'BrillVsize = p.sizeY;';
                otherwise
                    disp('options.Brilliance_Window.Vsize is set to invalid name, setting to FOM by default');
                    IM_l{end+1}= 'BrillVsize = p.sizeY;';
            end
        else
            IM_l{end+1}= ['BrillVsize = ' num2str(IM_Vsize) '/100;'];
        end
    else
        disp('options.Brilliance_Window.Vsize is not defined, setting to FOM by default');
        IM_l{end+1}= 'BrillVsize = p.sizeY;';
    end
    
    % Set BrillHdiv
    if isfield(options_general.Brilliance_Window,'Hdiv')
        IM_Hdiv = options_general.Brilliance_Window.Hdiv;        
        if ischar(IM_Hdiv)
            switch IM_Hdiv
                % Not a valid option
                case 'component'
                    disp('options.Brilliance_Window.Hdiv cannot be set to component, resetting to FOM by default');
                    IM_l{end+1}= 'BrillHdiv = p.divreq_x;'; 
                case 'FOM'
                % Need to set Hdiv to FOM                    
                    IM_l{end+1}= 'BrillHdiv = p.divreq_x;';
                otherwise
                    disp('options.Brilliance_Window.Hdiv is set to invalid name, setting to FOM by default');
                    IM_l{end+1}= 'BrillHdiv = p.divreq_x;';
            end
        else
            IM_l{end+1}= ['BrillHdiv = ' num2str(IM_Hdiv) ';'];
        end
    else
        disp('options.Brilliance_Window.Hdiv is not defined, setting to FOM by default');
        IM_l{end+1}= 'BrillHdiv = p.divreq_x;';
    end
    
    % Set BrillVdiv
    if isfield(options_general.Brilliance_Window,'Vdiv')
        IM_Vdiv = options_general.Brilliance_Window.Vdiv;        
        if ischar(IM_Vdiv)
            switch IM_Vdiv
                % Not a valid option
                case 'component'
                    disp('options.Brilliance_Window.Vdiv cannot be set to component, resetting to FOM by default');
                    IM_l{end+1}= 'BrillVdiv = p.divreq_y;'; 
                case 'FOM'
                % Need to set Vdiv to FOM                    
                    IM_l{end+1}= 'BrillVdiv = p.divreq_y;';
                otherwise
                    disp('options.Brilliance_Window.Vdiv is set to invalid name, setting to FOM by default');
                    IM_l{end+1}= 'BrillVdiv = p.divreq_y;';
            end
        else
            IM_l{end+1}= ['BrillVdiv = ' num2str(IM_Vdiv) ';'];
        end
    else
        disp('options.Brilliance_Window.Vdiv is not defined, setting to FOM by default');
        IM_l{end+1}= 'BrillVdiv = p.divreq_y;';
    end
    
else
    disp('options.Brilliance_Window is not defined, setting all fields to FOM by default')
    IM_l{end+1}= 'BrillHsize = p.sizeX;';
    IM_l{end+1}= 'BrillVsize = p.sizeY;';
    IM_l{end+1}= 'BrillHdiv = p.divreq_x;';
    IM_l{end+1}= 'BrillVdiv = p.divreq_y;';
end

% Define Brillendx and Brellendy
if strcmp(IM_module, 'M_module')
    IM_l{end+1} = ['IM_fh = fopen(''' filename '1_geometry.dat'');'];
    IM_l{end+1} = 'while ~strcmp(fgetl(IM_fh), ''M''); end';
    IM_l{end+1} = 'line = fgetl(IM_fh);';
    IM_l{end+1} = 'fclose(IM_fh);';
    IM_l{end+1} = 'Mono_Dims = str2num(line(5:end)); %%[Full Mono Width, Projected Mono Width, Mono Height]';
    IM_l{end+1}= 'Brillendx = Mono_Dims(2);';
    IM_l{end+1}= 'Brillendy = Mono_Dims(3);';
elseif strcmp(IM_module, 'FullInstrument')
    % Full with dummy variables, these values will be replaced in the
    % Initialize section of the .inst file
    IM_l{end+1}= 'Brillendx = 1;';
    IM_l{end+1}= 'Brillendy = 1;';
else
    IM_l{end+1}= ['Brillendx = p.startx' IM_index ';'];
    IM_l{end+1}= ['Brillendy = p.starty' IM_index ';'];
end

% Define Brillsample_dist
if strcmp(IM_module, 'FullInstrument')
    IM_l{end+1}= 'Brillsample_dist = p.sample_dist;';
else
    IM_l{end+1}= 'Brillsample_dist = 0.0001;';
end

IM_l{end+1}='';
IM_l{end+1} = '%% Write corresponding ifit brilliance file for normalization';
IM_l{end+1}='fid = fopen([cpath ''/../brilliance_refference/run_brilliance_ifit'' filename ''.m''], ''a'');';
IM_l{end+1}='fprintf(fid,[''clear all;clc;close all;\\nIntermediateBrilliance_ref_ifit('' num2str(BrillHdiv) '', '' num2str(BrillVdiv) '', '' num2str(BrillHsize) '', '' num2str(BrillVsize) '', ''''brilliance_ref_'' filename ''.mat'''')\\nexit'']);';
IM_l{end+1}='fclose(fid);';

ifit_loc = which('mcstas.m');
ifit_loc = ifit_loc(1:end-length('/Applications/McStas/mcstas.m'));
if ispc % Fix escape sequence issue with Windows
    ifit_loc_new = ifit_loc(1);
    for i = 2:numel(ifit_loc)
        if strcmp(ifit_loc(i),'\')
            ifit_loc_new = [ifit_loc_new '\\'];
        else
            ifit_loc_new = [ifit_loc_new ifit_loc(i)];
        end
    end
    ifit_loc = ifit_loc_new;
end

IM_l{end+1} = '%% Make sure brilliance file has time to be fully written before executing';
IM_l{end+1} = 'pause(30)';    
IM_l{end+1} = '%% Execute brilliance file in background';    
IM_l{end+1} = ['system([''matlab -nosplash -nodisplay -nodesktop -r "addpath(genpath(''''' ifit_loc ''''')); cd '' cpath ''/../brilliance_refference; tic; try; run_brilliance_ifit'' filename ''; catch; end; toc; quit" > matlab_output.log &''])'];
IM_l{end+1} = '';
IM_l{end+1} = '%% Redefining Brilliance Window';
IM_l{end+1} = 'monitor_temp = monitor(1);';
IM_l{end+1} = 'Param_temp = monitor(1).Data.Parameters;';
IM_l{end+1} = 'Param_temp.BrillsizeX = BrillHsize;';
IM_l{end+1} = 'Param_temp.BrillsizeY = BrillVsize;';
IM_l{end+1} = 'Param_temp.Brilldivreq_x = BrillHdiv;';
IM_l{end+1} = 'Param_temp.Brilldivreq_y = BrillVdiv;';
IM_l{end+1} = 'Param_temp.Brillendx = Brillendx;';
IM_l{end+1} = 'Param_temp.Brillendy = Brillendy;';
IM_l{end+1} = 'Param_temp.Brillsample_dist = Brillsample_dist;';
IM_l{end+1} = 'set(monitor_temp,''Data.Parameters'',Param_temp);';
IM_l{end+1} = 'monitor = monitor_temp;';
IM_l{end+1} = '';
IM_l{end+1} = '%% Setting scale factors for plotting';
IM_l{end+1} = ['div1ds = ''' num2str(div1d) ''';'];
IM_l{end+1} = ['div2ds = ''' num2str(div2d) ''';'];
IM_l{end+1} = ['psd1ds = ''' num2str(psd1d) ''';'];
IM_l{end+1} = ['psd2ds = ''' num2str(psd2d) ''';'];
IM_l{end+1} = ['apsds = ''' num2str(apsd) ''';'];
IM_l{end+1} = ['adivs = ''' num2str(adiv) ''';'];
IM_l{end+1} = '';

IM_l{end+1}='';
IM_l{end+1}='%%options_home.compile=0;';
if options_general.gravity == 1
IM_l{end+1}='options_home.gravitation=1;';
else
IM_l{end+1}='options_home.gravitation=0;';    
end
IM_l{end+1}=['options_home.ncount=' num2str(options_general.ncount_multiplier_optimize) '*1e6;'];
IM_l{end+1}=['options_home.mpi=' num2str(options_general.mpi) ';'];
if ispc
    IM_l{end+1}='options_home.dir=[cpath ''\\'' filename];';
else
    IM_l{end+1}='options_home.dir=[cpath ''/'' filename];';
end
IM_l{end+1}='options_home.OutputFcn=''fminplot'';';
IM_l{end+1}='options_home.MaxFunEvals=750;';
IM_l{end+1}='options_home.monitors=''Div2d_sample_B'';';
IM_l{end+1}='options_home.type=''maximize'';';
IM_l{end+1}='options_home.optimizer=''fminpso'';';
IM_l{end+1}='options_home.TolFun =''0.08%%'';';
IM_l{end+1}='';
IM_l{end+1}='';
if options_general.gravity == 1
IM_l{end+1}='options_home.gravitation=1;';
else
IM_l{end+1}='options_home.gravitation=0;';
end
if isfield(options_general,'machines')
    IM_l{end+1} = ['options_home.machines=' '''' options_general.machines '''' ';'];
end
IM_l{end+1}='%% rest of options for cluster';
IM_l{end+1}='if (select==2); load NUMCORES.DAT; warning off; else NUMCORES=0; end;';
IM_l{end+1}='';
IM_l{end+1}='options_cluster.mpi=NUMCORES;';
IM_l{end+1}='options_cluster.Display=[];';
IM_l{end+1}='options_cluster.OutputFcn=''''';
IM_l{end+1}='';
IM_l{end+1}='%% optimization options';
IM_l{end+1}='options_cluster.type=''maximize'';';
IM_l{end+1}='options_cluster.optimizer=''fminpso'';';
IM_l{end+1}='options_cluster.TolFun =''0.03%%'';';
IM_l{end+1}='options_cluster.TolX =''0.07%%'';';
IM_l{end+1}='options_cluster.mode=''optimize'';';
IM_l{end+1}='options_cluster.monitors=''Div2d_sample_B'';';
IM_l{end+1}='options_cluster.MaxFunEvals=10000;';
IM_l{end+1}='options_cluster.MaxIter=10000;';
IM_l{end+1}='';
IM_l{end+1}='%% general opptions';
IM_l{end+1}='options_cluster.dir=filename;';
if options_general.gravity == 1
IM_l{end+1}='options_cluster.gravitation=1;';
else
IM_l{end+1}='options_cluster.gravitation=0;';
end
if isfield(options_general,'machines')
    IM_l{end+1} = ['options_cluster.machines=' '''' options_general.machines '''' ';'];
end
IM_l{end+1}='options_cluster.overwrite=1;';
IM_l{end+1}='options_home.overwrite=1;';
IM_l{end+1}='options_cluster.compile=0;';
IM_l{end+1}=['options_cluster.ncount=' num2str(options_general.ncount_multiplier_optimize) '*9e6;']; %9e6 standard
IM_l{end+1}='';
IM_l{end+1}='options={options_home options_cluster};';
IM_l{end+1}='';
IM_l{end+1}='options{select}';
IM_l{end+1}='';
IM_l{end+1}='%%------------------------ Analyzing the resulting optimum ------';
IM_l{end+1}='';
IM_l{end+1}='%%options_single_home.compile=0;';
if options_general.gravity == 1
IM_l{end+1}='options_single_home.gravitation=1;';
else
IM_l{end+1}='options_single_home.gravitation=0;';
end
if isfield(options_general,'machines')
    IM_l{end+1} = ['options_single_home.machines=' '''' options_general.machines '''' ';'];
end
IM_l{end+1}=['options_single_home.ncount=' num2str(options_general.ncount_multiplier_analyze) '*1e7;'];
IM_l{end+1}=['options_single_home.mpi=' num2str(options_general.mpi) ';'];
if ispc
    IM_l{end+1}='options_single_home.dir=[cpath ''\\'' filename ''waveALL''];';
else
    IM_l{end+1}='options_single_home.dir=[cpath ''/'' filename ''waveALL''];';
end
IM_l{end+1}='';
IM_l{end+1}='';
IM_l{end+1}='%% rest of options for cluster';
IM_l{end+1}='if (select==2); load NUMCORES.DAT; warning off; else NUMCORES=0; end;';
IM_l{end+1}='';
IM_l{end+1}='';
IM_l{end+1}='%% general opptions';
IM_l{end+1}='options_single_cluster.mpi=NUMCORES;';
IM_l{end+1}='options_single_cluster.Display=[];';
IM_l{end+1}='options_single_cluster.dir=[filename ''waveALL''];';
if options_general.gravity == 1
IM_l{end+1}='options_single_cluster.gravitation=1;';
else
IM_l{end+1}='options_single_cluster.gravitation=0;';
end
if isfield(options_general,'machines')
    IM_l{end+1} = ['options_single_cluster.machines=' '''' options_general.machines '''' ';'];
end
IM_l{end+1}='options_single_cluster.compile=0;';
IM_l{end+1}=['options_single_cluster.ncount=' num2str(options_general.ncount_multiplier_analyze) '*5e8;']; % 5e8 standard
IM_l{end+1}='';
IM_l{end+1}='options_single={options_single_home options_single_cluster};';
IM_l{end+1}='';
strings='wavecenters=[ ';
for i=1:NumSnaps
    strings=[strings num2str(wavec(i)) ' '];
end
strings=[strings '];'];
IM_l{end+1}=strings;
IM_l{end+1}=['snapwidth=0.01*ones(1,' num2str(NumSnaps) ');'];
if strcmp(options_general.optimizer_mode,'realistic_source')
    IM_l{end+1}='optimal_ess=monitor(1).Data.Parameters;';
    IM_l{end+1}='names_ess=fieldnames(optimal_ess);';
    if strcmp(text_mode,'''')
        IM_l{end+1}='for i=1:length(names_ess)';
        IM_l{end+1}='optimal_ess.(names_ess{i})=num2str(optimal_ess.(names_ess{i}));';
        IM_l{end+1}='end';
    end
    
    IM_l{end+1}='optimal=monitor(1).Data.Parameters;';
    for index = 1:length(McStasStr.input_ess)
        IM_l{end+1}=['optimal=rmfield(optimal,''' McStasStr.input_ess{index} ''');'];
    end
    IM_l{end+1}='names=fieldnames(optimal);';
    if strcmp(text_mode,'''')
        IM_l{end+1}='for i=1:length(names)';
        IM_l{end+1}='optimal.(names{i})=num2str(optimal.(names{i}));';
        IM_l{end+1}='end';
    end
elseif strcmp(options_general.optimizer_mode,'ideal_source')
    IM_l{end+1}='optimal=monitor(1).Data.Parameters;';
    IM_l{end+1}='names=fieldnames(optimal);';
    if strcmp(text_mode,'''')
        IM_l{end+1}='for i=1:length(names)';
        IM_l{end+1}='optimal.(names{i})=num2str(optimal.(names{i}));';
        IM_l{end+1}='end';
    end
    IM_l{end+1}='optimal_ess = optimal;'; % uses the fact that all the ess values to be used are the defualt values of the mcstas file
    IM_l{end+1}='names_ess=fieldnames(optimal_ess);';
elseif strcmp(options_general.optimizer_mode,'combined')
    IM_l{end+1}='optimal=monitor(1).Data.Parameters;';
    IM_l{end+1}='names=fieldnames(optimal);';
    if strcmp(text_mode,'''')
        IM_l{end+1}='for i=1:length(names)';
        IM_l{end+1}='optimal.(names{i})=num2str(optimal.(names{i}));';
        IM_l{end+1}='end';
    end
    IM_l{end+1}='optimal_ess=monitor(1).Data.Parameters;';
    IM_l{end+1}='names_ess=fieldnames(optimal_ess);';
    if strcmp(text_mode,'''')
        IM_l{end+1}='for i=1:length(names)';
        IM_l{end+1}='optimal_ess.(names_ess{i})=num2str(optimal_ess.(names_ess{i}));';
        IM_l{end+1}='end';
    end
end
IM_l{end+1}='';
IM_l{end+1}='%% Run for fom wavelengths';
IM_l{end+1}='optimal_visualizer = optimal;';
IM_l{end+1}='optimal_visualizer.file_name = [filename ''_geometry.dat'']';
IM_l{end+1}='monitor_fom=mcstas([instrument_name ''_analyze.instr''],optimal_visualizer,options_single{select});';
IM_l{end+1}='';
IM_l{end+1}='optimal.WaveMin=''0.1'';';
IM_l{end+1}='optimal_ess.WaveMin=''0.1'';';
% Determine if a Monochromator in reflection mode is present
for i = 1:numel(modulelist)
    if strcmp(modules{modulelist(i)},'M') % Check for monochromator
        monoopts = Parse_options(globalinfo.options{i});
        if strcmp(monoopts.BeamDir, 'reflect') || monoopts.Ebin <= 0 % Check if monochromator is in reflection mode or dynamic
            monoindex = num2str(1+ length(modulelist) - i);
            IM_l{end+1} = ['optimal.Ebin' monoindex ' = ''-3'';'];
            IM_l{end+1} = ['optimal_ess.Ebin' monoindex ' = ''-3'';'];
        end
    end
end
if max(ismember(fieldnames(options_general),'max_wavelength_investigated_multiplier')) == 0
    if (demands.WaveLmax<3)
    IM_l{end+1}=['optimal.WaveMax=''' num2str(demandscan.WaveLmax*3) ''';'];
    IM_l{end+1}=['optimal_ess.WaveMax=''' num2str(demandscan.WaveLmax*3) ''';'];
    else
    IM_l{end+1}=['optimal.WaveMax=''' num2str(demandscan.WaveLmax*2) ''';'];
    IM_l{end+1}=['optimal_ess.WaveMax=''' num2str(demandscan.WaveLmax*2) ''';'];
    end
else 
    IM_l{end+1}=['optimal.WaveMax=''' num2str(demandscan.WaveLmax*options_general.max_wavelength_investigated_multiplier) ''';'];
    IM_l{end+1}=['optimal_ess.WaveMax=''' num2str(demandscan.WaveLmax*options_general.max_wavelength_investigated_multiplier) ''';'];
end

IM_l{end+1}='MaxWB=str2num(optimal.WaveMax);';
IM_l{end+1}='';
IM_l{end+1}='%% Run for all wavelengths';
IM_l{end+1}=['options_single{select}.dir=[filename ''waveLarge''];'];
IM_l{end+1}='optimal_visualizer = optimal;';
IM_l{end+1}='optimal_visualizer.file_name = [filename ''_geometry.dat''];';
IM_l{end+1}='monitor_ALLW=mcstas([instrument_name ''_analyze.instr''],optimal_visualizer,options_single{select});';
IM_l{end+1}='';
IM_l{end+1}='%% Run for rough ESS source spectrum';
IM_l{end+1}=['options_single{select}.dir=[filename ''waveESS''];'];
%l{end+1}='optimal_ess = optimal;'; % these two are destinct in latest update
if strcmp(source_component,'ESS_pancake')
IM_l{end+1}='optimal_ess.scan1=num2str(scani);';
IM_l{end+1}='optimal_ess.scan2=num2str(scanj);';
end
IM_l{end+1}='monitor_ESSW=mcstas([instrument_name ''_analyze_ess.instr''],optimal_ess,options_single{select});';
IM_l{end+1}='';
IM_l{end+1}='%% Run for robustness check';
IM_l{end+1}='';
IM_l{end+1}='degraded=optimal;';
IM_l{end+1}='';
IM_l{end+1}='for i=1:length(names)';
IM_l{end+1}='tmpname=names{i}';
IM_l{end+1}=' if strcmp(tmpname(1:1),''m'') && length(tmpname)<4';
if strcmp(text_mode,'''')
    IM_l{end+1}='   if ischar(optimal.(names{i}))';
    IM_l{end+1}='     degraded.(names{i})=num2str(0.8*str2num(optimal.(names{i})))';
    IM_l{end+1}='   else';
    IM_l{end+1}='     degraded.(names{i})=num2str(0.8*optimal.(names{i}))';
    IM_l{end+1}='   end';
else
    IM_l{end+1}='   if ischar(optimal.(names{i}))';
    IM_l{end+1}='    degraded.(names{i})=0.8*str2num(optimal.(names{i}))';
    IM_l{end+1}='   else';
    IM_l{end+1}='    degraded.(names{i})=0.8*optimal.(names{i})';
    IM_l{end+1}='   end';
end
IM_l{end+1}=' end';
IM_l{end+1}=' if length(tmpname)>5';
IM_l{end+1}='  if strcmp(tmpname(1:5),''alpha'')';
if strcmp(text_mode,'''')
     IM_l{end+1}='   if ischar(optimal.(names{i}))';
     IM_l{end+1}='     degraded.(names{i})=num2str(1.4*str2num(optimal.(names{i})))';
     IM_l{end+1}='   else';
     IM_l{end+1}='     degraded.(names{i})=num2str(1.4*optimal.(names{i}))';
     IM_l{end+1}='   end';
else
     IM_l{end+1}='   if ischar(optimal.(names{i}))';
     IM_l{end+1}='     degraded.(names{i})=1.4*str2num(optimal.(names{i}))';
     IM_l{end+1}='   else';
     IM_l{end+1}='     degraded.(names{i})=1.4*optimal.(names{i})';
     IM_l{end+1}='   end';
end
IM_l{end+1}='  end';
IM_l{end+1}=' end';
IM_l{end+1}='end';
IM_l{end+1}='degraded_ess=optimal_ess;';
IM_l{end+1}='for i=1:length(names_ess)';
IM_l{end+1}='tmpname=names_ess{i}';
IM_l{end+1}=' if strcmp(tmpname(1:1),''m'') && length(tmpname)<4';
if strcmp(text_mode,'''')
    IM_l{end+1}='   if ischar(optimal_ess.(names_ess{i}))';
    IM_l{end+1}='     degraded_ess.(names_ess{i})=num2str(0.8*str2num(optimal_ess.(names_ess{i})))';
    IM_l{end+1}='   else';
    IM_l{end+1}='     degraded_ess.(names_ess{i})=num2str(0.8*optimal_ess.(names_ess{i}))';
    IM_l{end+1}='   end';
else
    IM_l{end+1}='   if ischar(optimal_ess.(names_ess{i}))';
    IM_l{end+1}='    degraded_ess.(names_ess{i})=0.8*str2num(optimal_ess.(names_ess{i}))';
    IM_l{end+1}='   else';
    IM_l{end+1}='    degraded_ess.(names_ess{i})=0.8*optimal_ess.(names_ess{i})';
    IM_l{end+1}='   end';
end
IM_l{end+1}=' end';
IM_l{end+1}=' if length(tmpname)>5';
IM_l{end+1}='  if strcmp(tmpname(1:5),''alpha'')';
if strcmp(text_mode,'''')
     IM_l{end+1}='   if ischar(optimal.(names_ess{i}))';
     IM_l{end+1}='     degraded_ess.(names_ess{i})=num2str(1.4*str2num(optimal_ess.(names_ess{i})))';
     IM_l{end+1}='   else';
     IM_l{end+1}='     degraded_ess.(names{i})=num2str(1.4*optimal_ess.(names_ess{i}))';
     IM_l{end+1}='   end';
else
     IM_l{end+1}='   if ischar(optimal_ess.(names_ess{i}))';
     IM_l{end+1}='     degraded_ess.(names_ess{i})=1.4*str2num(optimal_ess.(names_ess{i}))';
     IM_l{end+1}='   else';
     IM_l{end+1}='     degraded_ess.(names_ess{i})=1.4*optimal_ess.(names_ess{i})';
     IM_l{end+1}='   end';
end
IM_l{end+1}='  end';
IM_l{end+1}=' end';
IM_l{end+1}='end';
IM_l{end+1}='';
IM_l{end+1}='%% Run for all wavelengths';
IM_l{end+1}=['options_single{select}.dir=[filename ''Large_degraded''];'];
IM_l{end+1}='degraded_visualizer = degraded;';
IM_l{end+1}='degraded_visualizer.file_name = [filename ''_geometry.dat''];';
IM_l{end+1}='monitor_ALLW_degraded=mcstas([instrument_name ''_analyze.instr''],degraded_visualizer,options_single{select});';
IM_l{end+1}='';
IM_l{end+1}='%% Run for rough ESS source spectrum';
IM_l{end+1}=['options_single{select}.dir=[filename ''ESS_degraded''];'];
%l{end+1}='degraded_ess = degraded;';
if strcmp(source_component,'ESS_pancake')
IM_l{end+1}='degraded_ess.scan1=num2str(scani);';
IM_l{end+1}='degraded_ess.scan2=num2str(scanj);';
end
IM_l{end+1}='monitor_ESSW_degraded=mcstas([instrument_name ''_analyze_ess.instr''],degraded_ess,options_single{select});';
IM_l{end+1}='';
IM_l{end+1}='%% Run for individual wavelength snapshots';
IM_l{end+1}=['options_single_cluster.ncount=' num2str(options_general.ncount_multiplier_analyze) '*1e8;'];
IM_l{end+1}='options_single={options_single_home options_single_cluster};';

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
            if isfield(defaults, 'm') % MADS: Could have set m to something else in the module before mono
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
            IM_l{end+1} = '';
            IM_l{end+1} = '%% Changing snapwidths to monochromator integration width';
            for i = 1:NumSnaps
                % Note that I hardwire in a 3*FWHM for the snapwidth.
                % This is different than the width used during
                % optimization, which is set using the Ebin monochromator option.
                snapwidth_calc = 3*2*dM_temp*(mos_temp/60*pi/180 + mval_temp*wavec(i)*0.0024)/2*cos(asin(wavec(i)/2/dM_temp));
                IM_l{end+1} = ['snapwidth(' num2str(i) ') = ' num2str(snapwidth_calc) ';'];
            end
        end
    end
    clear dM_temp mval_temp mos_temp
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% END OF LELAND MODIFICATION   %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:NumSnaps
IM_l{end+1}='';    
IM_l{end+1}=['options_single{select}.dir=[filename ''wave' num2str(i) '''];'];
% guide_bot was updated to allow two different text_modes dependent on the
% iFit version, and this needs to be taken into account for the edits by
% Leland.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% LELAND MODIFICATION %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Need to vectorize snapwidth because of Monochromator case.
% See my other modification above this for more info.

% Original 2 Lines:
% l{end+1}=['optimal.WaveMin=num2str(wavecenters(' num2str(i) ') - snapwidth*0.5);'];
% l{end+1}=['optimal.WaveMax=num2str(wavecenters(' num2str(i) ') + snapwidth*0.5);'];

% Lelands suggested change:
%l{end+1}=['optimal.WaveMin=num2str(wavecenters(' num2str(i) ') - snapwidth(' num2str(i) ')*0.5);'];
%l{end+1}=['optimal.WaveMax=num2str(wavecenters(' num2str(i) ') + snapwidth(' num2str(i) ')*0.5);'];

% Implemented change
if strcmp(text_mode,'''')
    IM_l{end+1}=['optimal.WaveMin=num2str(wavecenters(' num2str(i) ') - snapwidth(' num2str(i) ')*0.5);'];
    IM_l{end+1}=['optimal.WaveMax=num2str(wavecenters(' num2str(i) ') + snapwidth(' num2str(i) ')*0.5);'];
else
    IM_l{end+1}=['optimal.WaveMin=wavecenters(' num2str(i) ') - snapwidth(' num2str(i) ')*0.5;'];
    IM_l{end+1}=['optimal.WaveMax=wavecenters(' num2str(i) ') + snapwidth(' num2str(i) ')*0.5;'];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% END OF LELAND MODIFICATION %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%if strcmp(text_mode,'''')
%l{end+1}=['optimal.WaveMin=num2str(wavecenters(' num2str(i) ') - snapwidth*0.5);'];
%l{end+1}=['optimal.WaveMax=num2str(wavecenters(' num2str(i) ') + snapwidth*0.5);'];
%else
%l{end+1}=['optimal.WaveMin=wavecenters(' num2str(i) ') - snapwidth*0.5;'];
%l{end+1}=['optimal.WaveMax=wavecenters(' num2str(i) ') + snapwidth*0.5;'];    
%end
IM_l{end+1}='optimal_visualizer = optimal;';
IM_l{end+1}='optimal_visualizer.file_name = [filename ''_geometry.dat''];';
IM_l{end+1}=['monitor_W.wave' num2str(i) '=mcstas([instrument_name ''_analyze.instr''],optimal_visualizer,options_single{select});'];
if detailed_absolute_run
IM_l{end+1}=['options_single{select}.dir=[filename ''wave_ess' num2str(i) '''];'];
if strcmp(text_mode,'''')
IM_l{end+1}=['optimal_ess.WaveMin=num2str(wavecenters(' num2str(i) ') - snapwidth(' num2str(i) ')*0.5);'];
IM_l{end+1}=['optimal_ess.WaveMax=num2str(wavecenters(' num2str(i) ') + snapwidth(' num2str(i) ')*0.5);'];
else
IM_l{end+1}=['optimal_ess.WaveMin=wavecenters(' num2str(i) ') - snapwidth(' num2str(i) ')*0.5;'];
IM_l{end+1}=['optimal_ess.WaveMax=wavecenters(' num2str(i) ') + snapwidth(' num2str(i) ')*0.5;'];    
end
if strcmp(source_component,'ESS_pancake')
IM_l{end+1}='optimal_ess.scan1=num2str(scani);';
IM_l{end+1}='optimal_ess.scan2=num2str(scanj);';
end
IM_l{end+1}=['monitor_W_ess.wave' num2str(i) '=mcstas([instrument_name ''_analyze_ess.instr''],optimal_ess,options_single{select});'];
end
IM_l{end+1}='';
end
IM_l{end+1}='';

IM_l{end+1}='save([filename ''_all.mat'']);';
IM_l{end+1}='save([cpath ''/../output/analysis/'' filename scanname ''_all.mat'']);';

IM_l{end+1}='';
IM_l{end+1}='fid = fopen([cpath ''/../output/analysis/'' filename ''_analyze_all_ifit.m''], ''a'');';
IM_l{end+1}='fprintf(fid,[''clear all;clc;close all;\\nIntermediateBrilliance_Analysis('' num2str(BrillHdiv) '', '' num2str(BrillVdiv) '', 100*'' num2str(BrillHsize) '', 100*'' num2str(BrillVsize) '', '' div1ds '', '' div2ds '', '' psd1ds '', '' psd2ds '', '' apsds '', '' adivs '', '''''' filename scanname ''_all.mat'''', ''''brilliance_ref_'' filename ''.mat'''', ''''default'''')\\n'']);';
IM_l{end+1}='fclose(fid);';


IM_iFitStr='';

for k=1:length(IM_l)
    IM_iFitStr=[IM_iFitStr IM_l{k} '\n'];
end
%clear mon_l

fid = fopen(['./' Project_name '/' filename '/' filename scanname '_' fnamesIM{j} '_ifit.m'], 'w');
fprintf(fid,IM_iFitStr);
fclose(fid);

clear IM_iFitStr IM_l



