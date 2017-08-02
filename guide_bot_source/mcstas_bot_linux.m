function mcstas_bot(input,demands,requirements,defaults,Project_name,options_general)
% guide bot release 1


% Variables to select which cluster to use:
% ESSS DMSC
% PSI MCC

% Relevant changes:
% module names
% queue names
% paths
%   - to clusterid file
%   - to ifit file
%   - in brilliance transfer
if ~isfield(options_general,'cluster')
    options_general.cluster='ESSS';
    disp('Creating batch files for ESSS DMSC cluster as default.')
end

if strcmp(options_general.cluster,'ESS')
    options_general.cluster='ESSS';
end

switch options_general.cluster
    case 'ESSS'
        options_general.cluster_path = '/users/mbertelsen/';
        if ~isfield(options_general,'queue')
           options_general.queue = 'long';
           disp('ESS DMSC long queue selected as default.')
        end
        switch options_general.queue
            case 'express'
                options_general.time = '1:59:00';
            case 'long'
                options_general.time = '23:59:00';
            case 'verylong'
                options_general.time = '2-23:59:00';
            otherwise
                disp(['ERROR, options.queue = ' options_general.queue ' is invalid when using ESSS DMSC cluster']);
        end
        % list of nessecary modules
        options_general.modules_node{1} = 'mcstas/2.0';
        options_general.modules_node{2} = 'intel';
        options_general.modules_node{3} = 'openmpi/intel.qlc';
        
        options_general.modules_compile{1} = 'mcstas/2.0';
        options_general.modules_compile{2} = 'openmpi/intel.qlc';
    case 'PSI'
        options_general.cluster_path = '/home/l_bertelsen/';
        if ~isfield(options_general,'queue')
            options_general.queue = 'medium';
            disp('PSI medium queue selected as default.')
        end
        switch options_general.queue
            case 'test'
                options_general.time = '1:59:00';
            case 'short'
                options_general.time = '0:59:00';
            case 'medium'
                options_general.time = '23:59:00';
            case 'long'
                options_general.time = '2-23:59:00';
            otherwise
                disp(['ERROR, options.queue = ' options_general.queue ' is invalid when using PSI mcc cluster']);
        end
        
        
        % list of nessecary modules
        % does not seem to work at the PSI cluster
        options_general.modules_node{1} = 'mcstas/mcstas_2.0';
        options_general.modules_node{2} = 'mpi/openmpi-1.4.3-intel-12.0';
        
        options_general.modules_compile{1} = 'mcstas/mcstas_2.0';
        options_general.modules_compile{2} = 'mpi/openmpi-1.4.3-intel-12.0';
        
    otherwise
        disp(['ERROR, options.cluster = ' options_general.cluster ' not valid. Use ESSS or PSI.'])
end

if strcmp(requirements.source,'cold') || strcmp(requirements.source,'thermal') || strcmp(requirements.source,'bispectral')
    source_component = 'ESS_moderator_long';
elseif strcmp(requirements.source,'cold_pancake') || strcmp(requirements.source,'thermal_pancake') || strcmp(requirements.source,'bispectral_pancake')
    source_component = 'ESS_pancake';
    if requirements.moderator_size_y > 0.10
        disp(['ERROR, moderator higher than 10cm, no data for this case!']);
    end
elseif strcmp(requirements.source,'3cm_pancake_cold') || strcmp(requirements.source,'3cm_pancake_thermal') || strcmp(requirements.source,'TDR_thermal_fixed') || strcmp(requirements.source,'TDR_cold_fixed') || strcmp(requirements.source,'OT_6cm_thermal') || strcmp(requirements.source,'OT_6cm_cold')  || strcmp(requirements.source,'OT_3cm_thermal') || strcmp(requirements.source,'OT_3cm_narrow_cold') || strcmp(requirements.source,'OT_3cm_wide_cold') || strcmp(requirements.source,'TDR_like_thermal') || strcmp(requirements.source,'TDR_like_cold') || strcmp(requirements.source,'Butterfly_heimdal_6cm_thermal')
    source_component = 'locked_source_gen';
    
    if strcmp(requirements.source,'3cm_pancake_cold')
        requirements.moderator_size_x = 0.12;
        requirements.moderator_size_y = 0.03;
        locked_source_gen_filename = 'FirstModerator_cold.dat'; 
    elseif strcmp(requirements.source,'3cm_pancake_thermal')
        requirements.moderator_size_x = 0.10;
        requirements.moderator_size_y = 0.03;
        locked_source_gen_filename = 'FirstModerator_thermal.dat'; 
    elseif strcmp(requirements.source,'TDR_thermal_fixed')
        requirements.moderator_size_x = 0.12;
        requirements.moderator_size_y = 0.12;
        locked_source_gen_filename = 'TDR_thermal.dat'; 
    elseif strcmp(requirements.source,'TDR_cold_fixed')
        requirements.moderator_size_x = 0.12;
        requirements.moderator_size_y = 0.12;
        locked_source_gen_filename = 'TDR_cold.dat'; 
    elseif strcmp(requirements.source,'OT_6cm_cold')
        requirements.moderator_size_x = 0.06;
        requirements.moderator_size_y = 0.06;
        locked_source_gen_filename = 'SecondModerator_6cm_cold.dat'; 
    elseif strcmp(requirements.source,'OT_6cm_thermal')
        requirements.moderator_size_x = 0.06;
        requirements.moderator_size_y = 0.06;
        locked_source_gen_filename = 'SecondModerator_6cm_thermal.dat'; 
    elseif strcmp(requirements.source,'Butterfly_heimdal_6cm_thermal')
        requirements.moderator_size_x = 0.126;
        requirements.moderator_size_y = 0.06;
        locked_source_gen_filename = 'SecondModerator_6cm_thermal.dat';         
    elseif strcmp(requirements.source,'OT_3cm_wide_cold')
        requirements.moderator_size_x = 0.06;
        requirements.moderator_size_y = 0.03;
        locked_source_gen_filename = 'SecondModerator_3cm_wide_cold.dat'; 
    elseif strcmp(requirements.source,'OT_3cm_narrow_cold')
        requirements.moderator_size_x = 0.03;
        requirements.moderator_size_y = 0.03;
        locked_source_gen_filename = 'SecondModerator_3cm_narrow_cold.dat'; 
    elseif strcmp(requirements.source,'OT_3cm_thermal')
        requirements.moderator_size_x = 0.06;
        requirements.moderator_size_y = 0.03;
        locked_source_gen_filename = 'SecondModerator_3cm_thermal.dat'; 
    elseif strcmp(requirements.source,'TDR_like_cold')
        requirements.moderator_size_x = 0.12;
        requirements.moderator_size_y = 0.06;
        locked_source_gen_filename = 'SecondModerator_TDR_like_cold.dat'; 
    elseif strcmp(requirements.source,'TDR_like_thermal')
        requirements.moderator_size_x = 0.10;
        requirements.moderator_size_y = 0.06;
        locked_source_gen_filename = 'SecondModerator_TDR_like_thermal.dat'; 
    else
        disp('ERROR, Mistake in source selecter')
    end
    
elseif strcmp(requirements.source,'Virtual_in')
    source_component = 'Virtual_in';
    options_general_names = fieldnames(options_general);
    if max(ismember(options_general_names,'virtual_in_repeat_analyze')) < 0.5 || max(ismember(options_general_names,'virtual_in_repeat_optimize')) < 0.5
       disp(['ERROR, when using Virtual_in, one needs to set option.virtual_in_repeat_analyze and option.virtual_in_repeat_optimize, now they are 10 and 1 as default.'])
       options_general.virtual_in_repeat_analyze = 10;
       options_general.virtual_in_repeat_optimize = 1;
    end
else
    disp(['ERROR, unknown requirements.source! ' requirements.source]) 
end

if ~isfield(options_general,'weighted_optimization')
    options_general.weighted_optimization=0;
else
    if options_general.weighted_optimization ~= 0
       options_general.weighted_optimization = 1; 
    end
end

if ~isfield(options_general,'mpi')
    options_general.mpi=2;
end

McStasStr.minimalist=options_general.minimalist;
% 1 Enable extraction calculations (minimalist)
% 0 Enable gap calculation which assumes same divergence req.
% -1 Only the end of the guide is calculated, every width/height free

if ~isfield(options_general,'minimalist_factor')
   requirements.minimalist_factor = 1;
else
   requirements.minimalist_factor = options_general.minimalist_factor; 
end

if ~isfield(options_general,'weighted_optimization')
    options_general.weighted_optimization=0;
else
    if options_general.weighted_optimization ~= 0
       options_general.weighted_optimization = 1; 
    end
end

NumSnaps=5;

% variables: Hdiv,Vdiv,Hsize,Vsize,WaveLmin,WaveLmax,Dist,Mod_sample
scan.possiblenames={'Hdiv' 'Vdiv' 'Hsize' 'Vsize' 'WaveLmin' 'WaveLmax' 'moderator_size_x' 'moderator_size_y' 'minimalist_factor'};
scan.possiblenames_mcstas={'divreq_x' 'divreq_y' 'sizeX' 'sizeY' 'WaveMin' 'WaveMax' 'mod_x' 'mod_y' 'minimalist_factor'};
scan.demands_to_requirements = 6;

% alternate way of finding the scaned variables:
scan.mode=0;

for ii = 1:scan.demands_to_requirements
    if length(demands.(scan.possiblenames{ii}))>1.5
        scan.mode=1;
        scan.demands(ii)=1;
        scan.tmplength(ii)=length(demands.(scan.possiblenames{ii}));
    else
        scan.demands(ii)=0;
    end
end
for ii = scan.demands_to_requirements+1:length(scan.possiblenames)
    if length(requirements.(scan.possiblenames{ii}))>1.5
        scan.mode=1;
        scan.demands(ii)=1;
        scan.tmplength(ii)=length(requirements.(scan.possiblenames{ii}));
    else
        scan.demands(ii)=0;
    end    
end

options_active = fieldnames(options_general);
if sum(ismember(options_active,'locked'))>0.5
   locked = options_general.locked; 
   scan.locked_mode = 1;
   disp(['scan.locked_mode = ' num2str(scan.locked_mode)])
else
   scan.locked_mode = 0;
   %disp(['scan.locked_mode = ' num2str(scan.locked_mode)])
end

if ~isfield(options_general,'gravity')
    options_general.gravity=1;
    disp('Simulating with gravity as default.')
else
    if ~(options_general.gravity == 0)
        options_general.gravity = 1;
    end 
end

if max(ismember(fieldnames(options_general),'max_wavelength_investigated_multiplier')) == 1
    disp(['Analyzing performance from 0 to ' num2str(demands.WaveLmax*options_general.max_wavelength_investigated_multiplier) ' AA'])
end

if ~isfield(options_general,'minimalist_direction')
        globalinfo.minimalist_direction_horizontal = 1;
        globalinfo.minimalist_direction_vertical = 1;
    if options_general.minimalist > -0.5
    disp('Using the minimalist principle in both directions, control with options.minimalist_direction = ''both'' ''horizontal'' or ''vertical''')
    end
else
    if strcmp(options_general.minimalist_direction,'horizontal')
        globalinfo.minimalist_direction_horizontal = 1;
        globalinfo.minimalist_direction_vertical = 0;
    elseif strcmp(options_general.minimalist_direction,'vertical')
        globalinfo.minimalist_direction_horizontal = 0;
        globalinfo.minimalist_direction_vertical = 1;
    elseif strcmp(options_general.minimalist_direction,'both')
        globalinfo.minimalist_direction_horizontal = 1;
        globalinfo.minimalist_direction_vertical = 1;
    else
        disp('ERROR: Unknown minimalist_direction value, can be vertical, horizontal or both. defaulting to both')
        globalinfo.minimalist_direction_horizontal = 1;
        globalinfo.minimalist_direction_vertical = 1;
    end 
end

if ~isfield(options_general,'phase_space_requirement')
   globalinfo.PS_req = 'homogen';
   if McStasStr.minimalist == 1
      disp('Requireing homogeneous phase-space from source as default, change with options.phase_space_requirement = ''total''') 
   end
else
   switch options_general.phase_space_requirement
       case 'total'
           globalinfo.PS_req = 'total';
       case 'homogen'
           globalinfo.PS_req = 'homogen';
       case 'homogeneous'
           globalinfo.PS_req = 'homogen';
       otherwise
           disp(['Unkonwn PS_req ' options_general.phase_space_requirement '. Use total or homogen'])
   end
end

if ~isfield(options_general,'focusing')
    globalinfo.focusing = 0;
    if McStasStr.minimalist == 1
        disp('Assuming a non focusing guide as default. Change with options.focusing = 1')
    end
else
    switch options_general.focusing
        case 'yes'
            globalinfo.focusing = 1;
        case 'no'
            globalinfo.focusing = 0;
        case 1
            globalinfo.focusing = 1;
        case 0
            globalinfo.focusing = 0;
        otherwise
            disp(['Unknown focusing option, ' options_general.focusing '. Use 1 or 0.'])
    end
end




scan.locked.list = zeros(length(scan.possiblenames),1);
if scan.locked_mode == 1
% if a field names Hdiv exists, lock with the string contained in that field
% example locked.Hdiv = 'Vdiv'; will lock these.
% if length(Vdiv) = 1, use Hdiv's values
% if length(Vdiv) = length(Hdiv) use both vectors
% if length(Vdiv) ~= length(Hdiv) || 1 Error!
locked_names = fieldnames(locked);
for i=1:length(locked_names)
   for j = 1:length(scan.possiblenames) 
       if strcmp(locked_names{i},scan.possiblenames{j})
            % This (j) is the parameter which should be locked down.
            % Now check what it should be locked to.
            locked_child = j;
            for t = 1:length(scan.possiblenames)
               if strcmp(locked.(locked_names{i}),scan.possiblenames{t}) 
                    locked_parent = t;
               end
            end
            % syntax in code: locked.Hdiv = 'Vdiv'
            % Hdiv is child Vdiv is parent.
            if scan.tmplength(locked_child) == 1 || scan.tmplength(locked_child) == scan.tmplength(locked_parent)
                scan.demands(locked_child)=0; % kill combination scanning for child
                scan.locked.names{i} = scan.possiblenames{locked_child};
                scan.locked.names_parent{i} = scan.possiblenames{locked_parent};
                scan.locked.namesmcstas{i} = scan.possiblenames_mcstas{locked_child};
                scan.locked.list(locked_child) = true;
                
                if scan.tmplength(locked_child) == 1
                    % In this case, just use the same values as the parent.
                    % This if is needed because scans can be demands or requirements
                    if locked_parent < scan.demands_to_requirements+0.5;
                        scan.locked.values{i} = demands.(scan.possiblenames{locked_parent}); 
                    else
                        scan.locked.values{i} = requirements.(scan.possiblenames{locked_parent}); 
                    end
                elseif scan.tmplength(locked_child) == scan.tmplength(locked_parent)
                    % In this case, use the vector given at the child
                    if locked_parent < scan.demands_to_requirements+0.5;
                        scan.locked.values{i} = demands.(scan.possiblenames{locked_child}); 
                    else
                        scan.locked.values{i} = requirements.(scan.possiblenames{locked_child}); 
                    end
                end
            else
               disp('ERROR, locked is used wrong. Please refer to manual for examples.') 
            end
       end
   end
end

if length(scan.locked.names) < 0.5
    disp('ERROR, Something went wrong with your options.locked assignment')
end
end


% Manual way below, automatic remade above.
% % Checking for scan mode (not automized as the order of demands is unknown)
% scan.mode=0;
% if length(demands.Hdiv)>1.5
%     scan.mode=1;
%     scan.demands(1)=1;
%     scan.tmplength(1)=length(demands.Hdiv);
% else
%     scan.demands(1)=0;
% end
% if length(demands.Vdiv)>1.5
%     scan.mode=1;
%     scan.demands(2)=1;
%     scan.tmplength(2)=length(demands.Vdiv);
% else
%     scan.demands(2)=0;
% end
% if length(demands.Hsize)>1.5
%     scan.mode=1;
%     scan.demands(3)=1;
%     scan.tmplength(3)=length(demands.Hsize);
% else
%     scan.demands(3)=0;
% end
% if length(demands.Vsize)>1.5
%     scan.mode=1;
%     scan.demands(4)=1;
%     scan.tmplength(4)=length(demands.Vsize);
% else
%     scan.demands(4)=0;
% end
% if length(demands.WaveLmin)>1.5
%     scan.mode=1;
%     scan.demands(5)=1;
%     scan.tmplength(5)=length(demands.WaveLmin);
% else
%     scan.demands(5)=0;
% end
% if length(demands.WaveLmax)>1.5
%     scan.mode=1;
%     scan.demands(6)=1;
%     scan.tmplength(6)=length(demands.WaveLmax);
% else
%     scan.demands(6)=0;
% end % AFTER THIS: REQUIREMENTS, below 6.5 demand, above requirement
% if length(requirements.moderator_size_x)>1.5
%     scan.mode=1;
%     scan.demands(7)=1;
%     scan.tmplength(7)=length(requirements.moderator_size_x);
% else
%     scan.demands(7)=0;
% end
% if length(requirements.moderator_size_y)>1.5
%     scan.mode=1;
%     scan.demands(8)=1;
%     scan.tmplength(8)=length(requirements.moderator_size_y);
% else
%     scan.demands(8)=0;
% end
% % Scanning over these last to might prove difficult as they are used
% % outside of the ifit files. Not yet.
% % if length(demands.Dist)>1.5
% %     scan.mode=1;
% %     scan.demands(7)=1;
% % else
% %     scan.demands(7)=0;
% % end
% % if length(demands.Mod_sample)>1.5
% %     scan.mode=1;
% %     scan.demands(8)=1;
% % else
% %     scan.demands(8)=0;
% % end

for i=find(scan.demands)
    if sum(ismember('names',fieldnames(scan)))>0.5 %If the names field is not made
       scan.names{end+1}=scan.possiblenames{i};
       scan.namesmcstas{end+1}=scan.possiblenames_mcstas{i};
    else
       scan.names{1}=scan.possiblenames{i};
       scan.namesmcstas{1}=scan.possiblenames_mcstas{i};
    end
end

if scan.mode==1
   scan.dimension=sum(scan.demands);
   scan.identifiers=find(scan.demands);
   scan.length=scan.tmplength(scan.identifiers);

   if scan.dimension==1
       scan.length(2)=1;
   end
   if scan.dimension>2.5
       disp('ERROR, cannot do more than a 2d scan!')
   end
% old error checking, should not be needed.   
%for i=find(scan.demands)
%   if scan.length ~= scan.tmplength(i)
%        disp('ERROR, scan dimension mismatch!')
%   end
%end
end

if scan.mode==0
    scan.length(1)=1;
    scan.length(2)=1;
    scan.names={''};
    scan.dimension=0;
end

% Handling the input string
sl=length(input);

% code for finding the different patterns embeded in the string:

% The modules which are checked for: 
modules={'S' 'G' 'K' 'E' 'P' 'C' 'Cg' 'Selene' 'Slit'};
% Add user modules to the above list.
% Planned additions: C B Es 

% List of modules which break line of sight.
losbreakers = {'K' 'C' 'Selene' 'Cg'};
% Planned additions: B Es

% List of modules which require ray tracing to calculate any los break
raytracereq = {'K' 'C' 'Cg'};
% Planned additions: B



% On to unraveling the input strings
% find spaces in input string (Bug if no spaces present)
j=0;
for i=1:sl
    if strcmp(input(i),' '); 
        j=j+1;
        spaces(j)=i;    
    end
end

if exist('spaces')<0.5
   part{1}=input;
else

for i=1:length(spaces)+1
       if i==1
            part{i}=input(1:(spaces(i)-1));
       elseif (i==length(spaces)+1)
            part{i}=input((spaces(i-1)+1):end);
       else
            part{i}=input((spaces(i-1)+1):(spaces(i)-1));
       end
end
end

% Check for easy options
if isfield(options_general,'absolute_intensity_run')
    if options_general.absolute_intensity_run == 1
        detailed_absolute_run = 1;
    else 
        detailed_absolute_run = 0;
    end
else
    detailed_absolute_run = 0;
end


% Initialize variables
last=length(part);
globalinfo.nummodules=length(part);
globalinfo.fixedlength(1:length(part))=0;
globalinfo.fixedlengthdata(1:length(part))=0;
globalinfo.maxlength(1:length(part))=0;
globalinfo.maxlengthdata(1:length(part))=0;
globalinfo.minlength(1:length(part))=0;
globalinfo.minlengthdata(1:length(part))=0;
globalinfo.fixedstart(1:length(part))=0;
globalinfo.fixedstartpoint(1:length(part))=0;
globalinfo.maxstart(1:length(part))=0;
globalinfo.maxstartpoint(1:length(part))=0;
globalinfo.minstart(1:length(part))=0;
globalinfo.minstartpoint(1:length(part))=0;
globalinfo.rotlogic(1:length(part))=0;
globalinfo.rotsign(1:length(part))=0;
globalinfo.kinkoverwrite(1:length(part))=0;
globalinfo.losbreaks=0;
globalinfo.raytracereq=0;
globalinfo.selene_optimize_first_slit=0;
globalinfo.selene_moderator_focus=0;
extraoptions{length(part)+1}='dummy';
McStasStr.initialize_seg = '';
start_first(1:length(part)) = -1;

parlist=-(-length(part):-1);
for par=1:length(part)
    
    % Variable cor_index (corresponding_index) is the index which will be
    % used in the main loop of guide bot to index the modules. It is 1 for
    % the last component and rises from the sample towards the moderator.
    cor_index = 1+length(part)-par;
    
%for par=parlist
    % find ( in part(i) string
    j=0;
    leftp=-1;rightp=-1;
    for i=1:length(part{par})
        if strcmp(part{par}(i),'('); 
            j=j+1;
            leftp=i;    
        end
    end

    % find ) in part(i) string
    j=0;
    for i=1:length(part{par})
        if strcmp(part{par}(i),')'); 
            j=j+1;
            rightp=i;    
        end
    end

    %if (length(leftp)>1); disp('ERROR, input have more than one ('); end;
    %if (length(rightp)>1); disp('ERROR, input have more than one )'); end;
    
    if (leftp==-1 || rightp==-1) %Case: no additional commands
        letter{par}=part{par};
        options{par}='';
    else % Case: Additional commands given
        letter{par}=part{par}(1:leftp-1);
        optionstr{par}=part{par}(leftp+1:rightp-1);
        
        commap=-1;j=0;
        for i=1:length(optionstr{par})
            if strcmp(optionstr{par}(i),','); 
                j=j+1;
                commap(j)=i;    
            end
        end
        
        if (commap ~= -1) % If there is commas, there is several commands
            for i=1:length(commap)+1
                   if i==1
                        options{par}{i}=optionstr{par}(1:(commap(i)-1));
                   elseif (i==length(commap)+1)
                        options{par}{i}=optionstr{par}((commap(i-1)+1):end);
                   else
                        options{par}{i}=optionstr{par}((commap(i-1)+1):(commap(i)-1));
                   end
            end
        else % No commas, means just one additional command was given
            options{par}{1}=optionstr{par};
        end
    end
    
    % Add the information to the modulelist variable
    for j=1:length(modules)
       if strcmp(letter{par},modules{j})
           modulelist(par)=j;
       end
    end
    
    
    % Counting the number of line of sight breakers
    for losb = 1:length(losbreakers)
       if strcmp(letter{par},losbreakers{losb})
           globalinfo.losbreaks = globalinfo.losbreaks + 1;
       end
    end
    
    % Counting the number of line of sight breakers which require ray trace
    for losb = 1:length(raytracereq)
       if strcmp(letter{par},losbreakers{losb})
           globalinfo.raytracereq = globalinfo.raytracereq + 1;
       end
    end
    
    % How to store information on the different los breaks?
    
    % Collecting the globalinfo needed for all modules
    % These values are needed for all, but also delivered for each in the
    % options cell which is not strictly nessecary.
    
    % The index numbers of the kinks
    if strcmp(letter{par},'K')
        globalinfo.kink_logic(cor_index)=true;
        if isfield(globalinfo,'kinks')
            globalinfo.kinks(end+1)=cor_index;
        else
            globalinfo.kinks=cor_index;
        end
        for i=1:length(options{par})
           tmp=options{par}(i);
           keyword='rotd=';
            klength=length(keyword);
            if length(tmp{1})>klength
               if strcmp(tmp{1}(1:klength),keyword)
                   if strcmp(tmp{1}(klength+1:end),'h')
                       globalinfo.rotlogic(cor_index)=1;
                   elseif strcmp(tmp{1}(klength+1:end),'v')
                       globalinfo.rotlogic(cor_index)=-1;
                   end
               end
            end
            tmp=options{par}(i);
            keyword='rots=';
            klength=length(keyword);
            if length(tmp{1})>klength
               if strcmp(tmp{1}(1:klength),keyword)
                   if str2num(tmp{1}(klength+1:end))>0
                       globalinfo.rotsign(cor_index)=1;
                   else
                       globalinfo.rotsign(cor_index)=-1;
                   end
               end
            end
            tmp=options{par}(i);
            keyword='rot=';
            klength=length(keyword);
            if length(tmp{1})>klength
               if strcmp(tmp{1}(1:klength),keyword)
                   globalinfo.kinkoverwrite(cor_index)=1;
                   % actually value of rot is read within K module
               end
            end
        end
    end
    
    % The index numbers of the curved sections
    if strcmp(letter{par},'C') || strcmp(letter{par},'Cg')
        globalinfo.curve_logic(cor_index)=true;
        if isfield(globalinfo,'curves')
            globalinfo.curves(end+1)=cor_index;
        else
            globalinfo.curves=cor_index;
        end
        for i=1:length(options{par})
           tmp=options{par}(i);
           keyword='rotd=';
            klength=length(keyword);
            if length(tmp{1})>klength
               if strcmp(tmp{1}(1:klength),keyword)
                   if strcmp(tmp{1}(klength+1:end),'h')
                       globalinfo.rotlogic(cor_index)=1;
                   elseif strcmp(tmp{1}(klength+1:end),'v')
                       globalinfo.rotlogic(cor_index)=-1;
                   end
               end
            end
            tmp=options{par}(i);
            keyword='rots=';
            klength=length(keyword);
            if length(tmp{1})>klength
               if strcmp(tmp{1}(1:klength),keyword)
                   if str2num(tmp{1}(klength+1:end))>0
                       globalinfo.rotsign(cor_index)=1;
                   else
                       globalinfo.rotsign(cor_index)=-1;
                   end
               end
            end
            tmp=options{par}(i);
            keyword='rot=';
            klength=length(keyword);
            if length(tmp{1})>klength
               if strcmp(tmp{1}(1:klength),keyword)
                   globalinfo.curveoverwrite(cor_index)=1;
                   % actually value of rot is read within C module
               end
            end
        end
        
        if length(options{par})>0.5
          for i=1:length(options{par})
            tmp=options{par}(i);
            keyword='minStartWidth=';
            klength=length(keyword);
            if length(tmp{1})>klength;
               if strcmp(tmp{1}(1:klength),keyword)
                   % send this option to the previous non C module module
                   if par < length(part)-0.5
                       if length(extraoptions{par+1})>0.5
                        extraoptions{par+1}{end+1}=['minStartWidth=' tmp{1}(klength+1:end)];
                       else
                        extraoptions{par+1}{1}=['minStartWidth=' tmp{1}(klength+1:end)];
                       end
                        options{par}{i}='changed=';
                   else
                      disp('ERROR, C should not be the last module')
                   end
               end
            end
            keyword='minStartHeight=';
            klength=length(keyword);
            if length(tmp{1})>klength
               if strcmp(tmp{1}(1:klength),keyword)
                   % send this option to the previous module
                   if par < length(part)-0.5
                       if length(extraoptions{par+1})>0.5
                        extraoptions{par+1}{end+1}=['minStartHeight=' tmp{1}(klength+1:end)];
                       else
                        extraoptions{par+1}{1}=['minStartHeight=' tmp{1}(klength+1:end)];   
                       end
                        options{par}{i}='changed=';
                   else
                      disp('ERROR, C should not be the last module')
                   end
               end
            end
            keyword='maxStartWidth=';
            klength=length(keyword);
            if length(tmp{1})>klength;
               if strcmp(tmp{1}(1:klength),keyword)
                   % send this option to the previous module
                   if par < length(part)-0.5
                       if length(extraoptions{par+1})>0.5
                        extraoptions{par+1}{end+1}=['maxStartWidth=' tmp{1}(klength+1:end)];
                       else
                        extraoptions{par+1}{1}=['maxStartWidth=' tmp{1}(klength+1:end)];
                       end
                        options{par}{i}='changed=';
                   else
                      disp('ERROR, C should not be the last module')
                   end
               end
            end
            keyword='maxStartHeight=';
            klength=length(keyword);
            if length(tmp{1})>klength
               if strcmp(tmp{1}(1:klength),keyword)
                   % send this option to the previous module
                   if par < length(part)-0.5
                       if length(extraoptions{par+1})>0.5
                        extraoptions{par+1}{end+1}=['maxStartHeight=' tmp{1}(klength+1:end)];
                       else
                        extraoptions{par+1}{1}=['maxStartHeight=' tmp{1}(klength+1:end)];   
                       end
                        options{par}{i}='changed=';
                   else
                      disp('ERROR, C should not be the last module')
                   end
               end
            end
            keyword='StartHeight=';
            klength=length(keyword);
            if length(tmp{1})>klength
               if strcmp(tmp{1}(1:klength),keyword)
                   % send this option to the previous module
                   if par < length(part)-0.5
                       if length(extraoptions{par+1})>0.5
                        extraoptions{par+1}{end+1}=['StartHeight=' tmp{1}(klength+1:end)];
                       else
                        extraoptions{par+1}{1}=['StartHeight=' tmp{1}(klength+1:end)];   
                       end
                        options{par}{i}='changed=';
                   else
                      disp('ERROR, C should not be the last module')
                   end
               end
            end
            keyword='StartWidth=';
            klength=length(keyword);
            if length(tmp{1})>klength
               if strcmp(tmp{1}(1:klength),keyword)
                   % send this option to the previous module
                   if par < length(part)-0.5
                       if length(extraoptions{par+1})>0.5
                        extraoptions{par+1}{end+1}=['StartWidth=' tmp{1}(klength+1:end)];
                       else
                        extraoptions{par+1}{1}=['StartWidth=' tmp{1}(klength+1:end)];   
                       end
                        options{par}{i}='changed=';
                   else
                      disp('ERROR, C should not be the last module')
                   end
               end
            end
          end
        end
        
        % In case of a series of C modules
        if length(extraoptions{par})>0.5
          for i=1:length(extraoptions{par})
            tmp=extraoptions{par}(i);
            keyword='minStartWidth=';
            klength=length(keyword);
            if length(tmp{1})>klength;
               if strcmp(tmp{1}(1:klength),keyword)
                   % send this option to the previous non C module module
                   if par < length(part)-0.5
                       if length(extraoptions{par+1})>0.5
                        extraoptions{par+1}{end+1}=['minStartWidth=' tmp{1}(klength+1:end)];
                       else
                        extraoptions{par+1}{1}=['minStartWidth=' tmp{1}(klength+1:end)];
                       end
                        extraoptions{par}{i}='changed=';
                   else
                      disp('ERROR, C should not be the last module')
                   end
               end
            end
            keyword='minStartHeight=';
            klength=length(keyword);
            if length(tmp{1})>klength
               if strcmp(tmp{1}(1:klength),keyword)
                   % send this option to the previous module
                   if par < length(part)-0.5
                       if length(extraoptions{par+1})>0.5
                        extraoptions{par+1}{end+1}=['minStartHeight=' tmp{1}(klength+1:end)];
                       else
                        extraoptions{par+1}{1}=['minStartHeight=' tmp{1}(klength+1:end)];   
                       end
                        extraoptions{par}{i}='changed=';
                   else
                      disp('ERROR, C should not be the last module')
                   end
               end
            end
            keyword='maxStartWidth=';
            klength=length(keyword);
            if length(tmp{1})>klength;
               if strcmp(tmp{1}(1:klength),keyword)
                   % send this option to the previous module
                   if par < length(part)-0.5
                       if length(extraoptions{par+1})>0.5
                        extraoptions{par+1}{end+1}=['maxStartWidth=' tmp{1}(klength+1:end)];
                       else
                        extraoptions{par+1}{1}=['maxStartWidth=' tmp{1}(klength+1:end)];
                       end
                        extraoptions{par}{i}='changed=';
                   else
                      disp('ERROR, C should not be the last module')
                   end
               end
            end
            keyword='maxStartHeight=';
            klength=length(keyword);
            if length(tmp{1})>klength
               if strcmp(tmp{1}(1:klength),keyword)
                   % send this option to the previous module
                   if par < length(part)-0.5
                       if length(extraoptions{par+1})>0.5
                        extraoptions{par+1}{end+1}=['maxStartHeight=' tmp{1}(klength+1:end)];
                       else
                        extraoptions{par+1}{1}=['maxStartHeight=' tmp{1}(klength+1:end)];   
                       end
                        extraoptions{par}{i}='changed=';
                   else
                      disp('ERROR, C should not be the last module')
                   end
               end
            end
            keyword='StartHeight=';
            klength=length(keyword);
            if length(tmp{1})>klength
               if strcmp(tmp{1}(1:klength),keyword)
                   % send this option to the previous module
                   if par < length(part)-0.5
                       if length(extraoptions{par+1})>0.5
                        extraoptions{par+1}{end+1}=['StartHeight=' tmp{1}(klength+1:end)];
                       else
                        extraoptions{par+1}{1}=['StartHeight=' tmp{1}(klength+1:end)];   
                       end
                        extraoptions{par}{i}='changed=';
                   else
                      disp('ERROR, C should not be the last module')
                   end
               end
            end
            keyword='StartWidth=';
            klength=length(keyword);
            if length(tmp{1})>klength
               if strcmp(tmp{1}(1:klength),keyword)
                   % send this option to the previous module
                   if par < length(part)-0.5
                       if length(extraoptions{par+1})>0.5
                        extraoptions{par+1}{end+1}=['StartWidth=' tmp{1}(klength+1:end)];
                       else
                        extraoptions{par+1}{1}=['StartWidth=' tmp{1}(klength+1:end)];   
                       end
                        extraoptions{par}{i}='changed=';
                   else
                      disp('ERROR, C should not be the last module')
                   end
               end
            end
          end
        end
    end    
    
    if strcmp(letter{par},'Selene')
        globalinfo.Selene_logic(cor_index)=true;
        if isfield(globalinfo,'Selenes')
            globalinfo.Selenes(end+1)=cor_index;
            else21
            globalinfo.Selenes=cor_index;
        end
    end
    
    if strcmp(letter{par},'Selene') && par == 1
        for i=1:length(options{par})
               tmp=options{par}(i);
               keyword='optimize_first_slit';
               klength=length(keyword);
               if length(tmp{1})>=klength
                    if strcmp(tmp{1}(1:klength),keyword)
                        globalinfo.selene_optimize_first_slit = 1;
                        if strcmp(tmp{1}(end),'v')
                           globalinfo.selene_optimize_first_slit_v=1; 
                        else
                           globalinfo.selene_optimize_first_slit_v=0;  
                        end
                        if strcmp(tmp{1}(end),'h')
                           globalinfo.selene_optimize_first_slit_h=1;
                        else
                           globalinfo.selene_optimize_first_slit_h=0;
                        end
                        options{par}{i}='changed=';
                    end
               end
               tmp=options{par}(i);
               keyword='moderator_focus';
               klength=length(keyword);
               if length(tmp{1})>=klength
                    if strcmp(tmp{1}(1:klength),keyword)
                        globalinfo.selene_moderator_focus = 1;
                        if strcmp(tmp{1}(end),'v')
                           globalinfo.selene_moderator_focus_v=1; 
                        else
                           globalinfo.selene_moderator_focus_v=0; 
                        end 
                        if strcmp(tmp{1}(end),'h')
                           globalinfo.selene_moderator_focus_h=1; 
                        else
                           globalinfo.selene_moderator_focus_h=0;  
                        end
                        options{par}{i}='changed=';
                    end
               end 
        end
            
    end

    if strcmp(letter{par},'Slit') && McStasStr.minimalist > -0.5
         disp('Warning, the Slit module is not designed to work unless minimalist is set to -1.')
         disp('It does give meaningfull results, but the slit might as well not be there.')
    end
    
    if strcmp(letter{par},'Slit')
    % Force a very small fixedlength to tell guide_bot what length this
    % module takes.
    globalinfo.fixedlength(cor_index)=1;
    globalinfo.fixedlengthdata(cor_index)=0.00002;
    % In principle the above value should be 2*1E-6, but I am uncertain of
    % the matlab acuraccy when writing files. It is a insigicicant distance
    % uncertaincy in any case.
      if length(options{par})>0.5
        for i=1:length(options{par})
            tmp=options{par}(i);
            keyword='minStartWidth=';
            klength=length(keyword);
            if length(tmp{1})>klength;
               if strcmp(tmp{1}(1:klength),keyword)
                   % send this option to the previous module
                   if par < length(part)-0.5
                       if length(extraoptions{par+1})>0.5
                        extraoptions{par+1}{end+1}=['minStartWidth=' tmp{1}(klength+1:end)];
                       else
                        extraoptions{par+1}{1}=['minStartWidth=' tmp{1}(klength+1:end)];
                       end
                        options{par}{i}='changed=';
                   else
                      disp('ERROR, minStartWidth specified for last slit')
                   end
               end
            end
            keyword='minStartHeight=';
            klength=length(keyword);
            if length(tmp{1})>klength
               if strcmp(tmp{1}(1:klength),keyword)
                   % send this option to the previous module
                   if par < length(part)-0.5
                       if length(extraoptions{par+1})>0.5
                        extraoptions{par+1}{end+1}=['minStartHeight=' tmp{1}(klength+1:end)];
                       else
                        extraoptions{par+1}{1}=['minStartHeight=' tmp{1}(klength+1:end)];   
                       end
                        options{par}{i}='changed=';
                   else
                      disp('ERROR, minStartHeight specified for last slit')
                   end
               end
            end
            keyword='maxStartWidth=';
            klength=length(keyword);
            if length(tmp{1})>klength;
               if strcmp(tmp{1}(1:klength),keyword)
                   % send this option to the previous module
                   if par < length(part)-0.5
                       if length(extraoptions{par+1})>0.5
                        extraoptions{par+1}{end+1}=['maxStartWidth=' tmp{1}(klength+1:end)];
                       else
                        extraoptions{par+1}{1}=['maxStartWidth=' tmp{1}(klength+1:end)];
                       end
                        options{par}{i}='changed=';
                   else
                      disp('ERROR, maxStartWidth specified for last slit')
                   end
               end
            end
            keyword='maxStartHeight=';
            klength=length(keyword);
            if length(tmp{1})>klength
               if strcmp(tmp{1}(1:klength),keyword)
                   % send this option to the previous module
                   if par < length(part)-0.5
                       if length(extraoptions{par+1})>0.5
                        extraoptions{par+1}{end+1}=['maxStartHeight=' tmp{1}(klength+1:end)];
                       else
                        extraoptions{par+1}{1}=['maxStartHeight=' tmp{1}(klength+1:end)];   
                       end
                        options{par}{i}='changed=';
                   else
                      disp('ERROR, maxStartHeight specified for last slit')
                   end
               end
            end
            keyword='StartHeight=';
            klength=length(keyword);
            if length(tmp{1})>klength
               if strcmp(tmp{1}(1:klength),keyword)
                   % send this option to the previous module
                   if par < length(part)-0.5
                       if length(extraoptions{par+1})>0.5
                        extraoptions{par+1}{end+1}=['StartHeight=' tmp{1}(klength+1:end)];
                       else
                        extraoptions{par+1}{1}=['StartHeight=' tmp{1}(klength+1:end)];   
                       end
                        options{par}{i}='changed=';
                   else
                      disp('ERROR, StartHeight specified for last slit')
                   end
               end
            end
            keyword='StartWidth=';
            klength=length(keyword);
            if length(tmp{1})>klength
               if strcmp(tmp{1}(1:klength),keyword)
                   % send this option to the previous module
                   if par < length(part)-0.5
                       if length(extraoptions{par+1})>0.5
                        extraoptions{par+1}{end+1}=['StartWidth=' tmp{1}(klength+1:end)];
                       else
                        extraoptions{par+1}{1}=['StartWidth=' tmp{1}(klength+1:end)];   
                       end
                        options{par}{i}='changed=';
                   else
                      disp('ERROR, StartWidth specified for last slit')
                   end
               end
            end
        end
      end
    end
    
    % The index numbers for the fixedstart modules and their startpoint
    if length(options{par})>0.5
        for i=1:length(options{par})
            tmp=options{par}(i);
            keyword='start=';
            klength=length(keyword);
            if length(tmp{1})>klength
               if strcmp(tmp{1}(1:klength),keyword)
                   startp=str2num(tmp{1}(klength+1:end));  
                   if (startp < max(globalinfo.fixedstartpoint))
                            disp('ERROR, a later component starts before an earlier!')
                   end
                      globalinfo.fixedstart(cor_index)=1;
                      globalinfo.fixedstartpoint(cor_index)=startp;  
               end
            end
            keyword='maxstart=';
            klength=length(keyword);
            if length(tmp{1})>klength
               if strcmp(tmp{1}(1:klength),keyword)
                   startp=str2num(tmp{1}(klength+1:end));  
                   %if (startp < max(globalinfo.fixedstartpoint))
                   %         disp('ERROR, a latter component starts before an earlier!')
                   %end
                      globalinfo.maxstart(cor_index)=1;
                      globalinfo.maxstartpoint(cor_index)=startp;
               end
            end
            keyword='minstart=';
            klength=length(keyword);
            if length(tmp{1})>klength
               if strcmp(tmp{1}(1:klength),keyword)
                   startp=str2num(tmp{1}(klength+1:end));  
                      globalinfo.minstart(cor_index)=1;
                      globalinfo.minstartpoint(cor_index)=startp;
               end
            end            
            keyword='length=';
            klength=length(keyword);
            if length(tmp{1})>klength
               if strcmp(tmp{1}(1:klength),keyword) 
                   lengthdata=str2num(tmp{1}(klength+1:end));
                             globalinfo.fixedlength(cor_index)=1;
                             globalinfo.fixedlengthdata(cor_index)=lengthdata;
               end
            end
            keyword='maxlength=';
            klength=length(keyword);
            if length(tmp{1})>klength
               if strcmp(tmp{1}(1:klength),keyword) 
                   lengthdata=str2num(tmp{1}(klength+1:end));
                             globalinfo.maxlength(cor_index)=1;
                             globalinfo.maxlengthdata(cor_index)=lengthdata;
               end
            end
            keyword='minlength=';
            klength=length(keyword);
            if length(tmp{1})>klength
               if strcmp(tmp{1}(1:klength),keyword) 
                   lengthdata=str2num(tmp{1}(klength+1:end));
                   if globalinfo.minlength(cor_index) == 1
                       if globalinfo.minlengthdata(cor_index) < lengthdata
                           globalinfo.minlengthdata(cor_index)=lengthdata;
                       end
                   else
                       globalinfo.minlengthdata(cor_index)=lengthdata;
                   end
                   globalinfo.minlength(cor_index)=1;
               end
            end
            keyword='los_divide=';
            klength=length(keyword);
            if length(tmp{1})>klength
               if strcmp(tmp{1}(1:klength),keyword)
                   if strcmp(tmp{1}(end),'m')
                    data=str2num(tmp{1}(klength+1:end-1));
                    if data == 0
                    globalinfo.los_divide_mode(cor_index)=1;    
                    else
                    globalinfo.los_divide_mode(cor_index)=2;
                    globalinfo.minlength(cor_index)=1;
                    end
                   elseif strcmp(tmp{1}(end),'s')
                    if length(tmp{1}) <= 3; disp('ERROR,los_divide=Xm_s not used correctly!'); end;
                    data=str2num(tmp{1}(klength+1:end-3));
                    if data == 0
                    globalinfo.los_divide_mode(cor_index)=1;     
                    else
                    globalinfo.los_divide_mode(cor_index)=2;
                    globalinfo.minlength(cor_index)=1;
                    end
                   elseif strcmp(tmp{1}(end),'e') % los_divide=5m_e
                    if length(tmp{1}) <= 3; disp('ERROR,los_divide=Xm_e not used correctly!'); end;
                    data=str2num(tmp{1}(klength+1:end-3));
                    if data == 0
                    globalinfo.los_divide_mode(cor_index)=1;
                    data = 1;
                    else
                    globalinfo.los_divide_mode(cor_index)=3;
                    globalinfo.minlength(cor_index)=1;
                    end
                   else
                    data=str2num(tmp{1}(klength+1:end));
                    globalinfo.los_divide_mode(cor_index)=1;   
                   end
                   if globalinfo.los_divide_mode(cor_index) >= 2
                    if globalinfo.minlength(cor_index) == 1
                        if globalinfo.minlengthdata(cor_index) < data
                            globalinfo.minlengthdata(cor_index)=data;
                        end
                    else
                        globalinfo.minlength(cor_index)=1;
                        globalinfo.minlengthdata(cor_index)=data;
                    end  
                   end
                   globalinfo.los_divide_logic(cor_index)=true;
                   globalinfo.los_divide_data(cor_index)=data;
                   globalinfo.los_divide_option_num(cor_index) = i;
                   % could add lines to divert this to changed
               end
            end
            keyword='los_start=';
            klength=length(keyword);
            if length(tmp{1})>klength
               if strcmp(tmp{1}(1:klength),keyword)
                   if strcmp(tmp{1}(end),'m')
                    data=str2num(tmp{1}(klength+1:end-1));
                    if data == 0
                    globalinfo.los_start_mode(cor_index)=1;
                    else
                    globalinfo.los_start_mode(cor_index)=2;
                    globalinfo.minlength(cor_index)=1;
                    end
                   elseif strcmp(tmp{1}(end),'s')
                    if length(tmp{1}) <= 3; disp('ERROR,los_start=Xm_s not used correctly!'); end;
                    data=str2num(tmp{1}(klength+1:end-3));
                    if data == 0
                    globalinfo.los_start_mode(cor_index)=1;
                    else
                    globalinfo.los_start_mode(cor_index)=2;
                    globalinfo.minlength(cor_index)=1;
                    end
                   elseif strcmp(tmp{1}(end),'e')
                    if length(tmp{1}) <= 3; disp('ERROR,los_start=Xm_e not used correctly!'); end;
                    data=str2num(tmp{1}(klength+1:end-3));
                    if data == 0
                    globalinfo.los_start_mode(cor_index)=1;
                    data = 1;
                    else
                    globalinfo.los_start_mode(cor_index)=3;
                    globalinfo.minlength(cor_index)=1;
                    end
                   else
                    data=str2num(tmp{1}(klength+1:end));
                    globalinfo.los_start_logic(cor_index)=true;
                    globalinfo.los_start_data(cor_index)=data;
                    globalinfo.los_start_mode(cor_index)=1;   
                   end
                   % could add lines to divert this to changed
                   globalinfo.los_start_logic(cor_index)=true;
                   globalinfo.los_start_data(cor_index)=data;
                   globalinfo.los_start_option_num(cor_index) = i;
                   
                   if globalinfo.los_start_mode(cor_index) >= 2
                    if globalinfo.minlength(cor_index) == 1
                        if globalinfo.minlengthdata(cor_index) < data
                            globalinfo.minlengthdata(cor_index)=data;
                        end
                    else
                        globalinfo.minlength(cor_index)=1;
                        globalinfo.minlengthdata(cor_index)=data;
                    end  
                   end
                   if start_first(cor_index) == -1
                        start_first(cor_index) = 1;
                   end
               end
            end
            keyword='los_end=';
            klength=length(keyword);
            if length(tmp{1})>klength
               if strcmp(tmp{1}(1:klength),keyword)
                   if strcmp(tmp{1}(end),'m')
                    data=str2num(tmp{1}(klength+1:end-1));
                    if data == 0
                    globalinfo.los_end_mode(cor_index)=1;
                    else
                    globalinfo.los_end_mode(cor_index)=2;
                    globalinfo.minlength(cor_index)=1;
                    end
                   elseif strcmp(tmp{1}(end),'s')
                    if length(tmp{1}) <= 3; disp('ERROR,los_end=Xm_s not used correctly!'); end;
                    data=str2num(tmp{1}(klength+1:end-3));
                    if data == 0
                    globalinfo.los_end_mode(cor_index)=1;
                    else
                    globalinfo.los_end_mode(cor_index)=2;
                    globalinfo.minlength(cor_index)=1;
                    end
                   elseif strcmp(tmp{1}(end),'e')
                    if length(tmp{1}) <= 3; disp('ERROR,los_end=Xm_e not used correctly!'); end;
                    data=str2num(tmp{1}(klength+1:end-3));
                    if data == 0;
                    globalinfo.los_end_mode(cor_index)=1;
                    data = 1;
                    else
                    globalinfo.los_end_mode(cor_index)=3;
                    globalinfo.minlength(cor_index)=1;
                    end
                   else
                    data=str2num(tmp{1}(klength+1:end));
                    globalinfo.los_end_mode(cor_index)=1;   
                   end
                   globalinfo.los_end_logic(cor_index)=true;
                   globalinfo.los_end_data(cor_index)=data;
                   globalinfo.los_end_option_num(cor_index) = i;
                   
                   if globalinfo.los_end_mode(cor_index) >= 2
                    if globalinfo.minlength(cor_index) == 1
                        if globalinfo.minlengthdata(cor_index) < data
                            globalinfo.minlengthdata(cor_index)=data;
                        end
                    else
                        globalinfo.minlength(cor_index)=1;
                        globalinfo.minlengthdata(cor_index)=data;
                    end
                   end
                   % could add lines to divert this to changed
                   if start_first(cor_index) == -1
                       start_first(cor_index) = 0;
                   end
               end
            end
        end
    end
    % In principle all these can be rewritte as one
    % Check for min/max end because these need to go to another module
    if length(options{par})>0.5
        for i=1:length(options{par})
            tmp=options{par}(i);
            keyword='minEndWidth=';
            klength=length(keyword);
            if length(tmp{1})>klength;
               if strcmp(tmp{1}(1:klength),keyword)
                   % send this option to the previous module
                   if par < length(part)-0.5
                       if length(extraoptions{par+1})>0.5
                        extraoptions{par+1}{end+1}=['minStartWidth=' tmp{1}(klength+1:end)];
                       else
                        extraoptions{par+1}{1}=['minStartWidth=' tmp{1}(klength+1:end)];
                       end
                        options{par}{i}='changed=';
                   else
                      disp('BETA FEATURE, minEndWidth specified for last guide')
                   end
               end
            end
            keyword='minEndHeight=';
            klength=length(keyword);
            if length(tmp{1})>klength
               if strcmp(tmp{1}(1:klength),keyword)
                   % send this option to the previous module
                   if par < length(part)-0.5
                       if length(extraoptions{par+1})>0.5
                        extraoptions{par+1}{end+1}=['minStartHeight=' tmp{1}(klength+1:end)];
                       else
                        extraoptions{par+1}{1}=['minStartHeight=' tmp{1}(klength+1:end)];   
                       end
                        options{par}{i}='changed=';
                   else
                      disp('BETA FEATURE, minEndHeight specified for last guide')
                   end
               end
            end
            keyword='maxEndWidth=';
            klength=length(keyword);
            if length(tmp{1})>klength;
               if strcmp(tmp{1}(1:klength),keyword)
                   % send this option to the previous module
                   if par < length(part)-0.5
                       if length(extraoptions{par+1})>0.5
                        extraoptions{par+1}{end+1}=['maxStartWidth=' tmp{1}(klength+1:end)];
                       else
                        extraoptions{par+1}{1}=['maxStartWidth=' tmp{1}(klength+1:end)];
                       end
                        options{par}{i}='changed=';
                   else
                      disp('BETA FEATURE, maxEndWidth specified for last guide')
                   end
               end
            end
            keyword='maxEndHeight=';
            klength=length(keyword);
            if length(tmp{1})>klength
               if strcmp(tmp{1}(1:klength),keyword)
                   % send this option to the previous module
                   if par < length(part)-0.5
                       if length(extraoptions{par+1})>0.5
                        extraoptions{par+1}{end+1}=['maxStartHeight=' tmp{1}(klength+1:end)];
                       else
                        extraoptions{par+1}{1}=['maxStartHeight=' tmp{1}(klength+1:end)];   
                       end
                        options{par}{i}='changed=';
                   else
                      disp('BETA FEATURE, maxEndHeight specified for last guide')
                   end
               end
            end
            keyword='EndHeight=';
            klength=length(keyword);
            if length(tmp{1})>klength
               if strcmp(tmp{1}(1:klength),keyword)
                   % send this option to the previous module
                   if par < length(part)-0.5
                       if length(extraoptions{par+1})>0.5
                        extraoptions{par+1}{end+1}=['StartHeight=' tmp{1}(klength+1:end)];
                       else
                        extraoptions{par+1}{1}=['StartHeight=' tmp{1}(klength+1:end)];   
                       end
                        options{par}{i}='changed=';
                   else
                      disp('BETA FEATURE, EndHeight specified for last guide')
                   end
               end
            end
            keyword='EndWidth=';
            klength=length(keyword);
            if length(tmp{1})>klength
               if strcmp(tmp{1}(1:klength),keyword)
                   % send this option to the previous module
                   if par < length(part)-0.5
                       if length(extraoptions{par+1})>0.5
                        extraoptions{par+1}{end+1}=['StartWidth=' tmp{1}(klength+1:end)];
                       else
                        extraoptions{par+1}{1}=['StartWidth=' tmp{1}(klength+1:end)];   
                       end
                        options{par}{i}='changed=';
                   else
                      disp('BETA FEATURE, EndWidth specified for last guide')
                   end
               end
            end
            keyword='max_smallaxis_x=';
            klength=length(keyword);
            if length(tmp{1})>klength
               if strcmp(tmp{1}(1:klength),keyword)
                   % send this option to the previous module
                   if par < length(part)-0.5
                       if str2num(tmp{1}(klength+1:end)) < 0.20 % standard width limit
                           if length(extraoptions{par+1})>0.5
                            extraoptions{par+1}{end+1}=['maxStartWidth=' num2str(0.95*str2num(tmp{1}(klength+1:end)))];
                           else
                            extraoptions{par+1}{1}=['maxStartWidth=' num2str(0.95*str2num(tmp{1}(klength+1:end)))];   
                           end
                       end
                   %else
                      %disp('ERROR, maxEndWidth specified for last guide')
                      % Will be handled by the E_module itself.
                   end
               end
            end
            keyword='max_smallaxis_y=';
            klength=length(keyword);
            if length(tmp{1})>klength
               if strcmp(tmp{1}(1:klength),keyword)
                   % send this option to the previous module
                   if par < length(part)-0.5
                       if str2num(tmp{1}(klength+1:end)) < 0.20 % standard width limit
                           if length(extraoptions{par+1})>0.5
                            extraoptions{par+1}{end+1}=['maxStartHeight=' num2str(0.95*str2num(tmp{1}(klength+1:end)))];
                           else
                            extraoptions{par+1}{1}=['maxStartHeight=' num2str(0.95*str2num(tmp{1}(klength+1:end)))];   
                           end
                       end
                   %else
                      %disp('ERROR, maxEndWidth specified for last guide')
                      % Will be handled by the E_module itself.
                   end
               end
            end
            keyword='max_smallaxis=';
            klength=length(keyword);
            if length(tmp{1})>klength
               if strcmp(tmp{1}(1:klength),keyword)
                   % send this option to the previous module
                   if par < length(part)-0.5
                       if str2num(tmp{1}(klength+1:end)) < 0.20 % standard width limit
                           if length(extraoptions{par+1})>0.5
                            extraoptions{par+1}{end+1}=['maxStartHeight=' num2str(0.95*str2num(tmp{1}(klength+1:end)))];
                            extraoptions{par+1}{end+1}=['maxStartWidth=' num2str(0.95*str2num(tmp{1}(klength+1:end)))];
                           else
                            extraoptions{par+1}{1}=['maxStartHeight=' num2str(0.95*str2num(tmp{1}(klength+1:end)))];   
                            extraoptions{par+1}{end+1}=['maxStartWidth=' num2str(0.95*str2num(tmp{1}(klength+1:end)))];   
                           end
                       end
                   %else
                      %disp('ERROR, maxEndWidth specified for last guide')
                      % Will be handled by the E_module itself.
                   end
               end
            end
        end
    end
end

for i=1:length(part)
   if ~isempty(extraoptions{i})
       for j=1:length(extraoptions{i})
           options{i}{end+1}=extraoptions{i}{j};
       end
   end    
end

ERRORlogic = modulelist == 0;
if sum(ERRORlogic)>0; disp('ERROR, unknown module used!'); end;

% DEFAULTS

% checking kink options and applying defaults if missing
if isfield(globalinfo,'kinks')
  for i=1:length(globalinfo.kinks)
    if globalinfo.rotlogic(globalinfo.kinks(i))==0
        globalinfo.rotlogic(globalinfo.kinks(i))=1;
    end
    if globalinfo.rotsign(globalinfo.kinks(i))==0
        globalinfo.rotsign(globalinfo.kinks(i))=1;
    end
  end
end
%if ~sum(ismember(fieldnames(globalinfo),'kinkoverwrite')); globalinfo.kinkoverwrite=0; end;

% checking curve options and applying defaults if missing
if isfield(globalinfo,'curves')
  for i=1:length(globalinfo.curves)
    if globalinfo.rotlogic(globalinfo.curves(i))==0
        globalinfo.rotlogic(globalinfo.curves(i))=1;
    end
    if globalinfo.rotsign(globalinfo.curves(i))==0
        globalinfo.rotsign(globalinfo.curves(i))=1;
    end
  end
end
%if ~sum(ismember(fieldnames(globalinfo),'kinkoverwrite')); globalinfo.kinkoverwrite=0; end;
% Each part is now split into letter and options (options in paranthesis)

% Pad the logic strings for K,C and Selene to avoid exceeding matrix dimensions
max_cor_index = length(part);
if isfield(globalinfo,'kink_logic')
   if length(globalinfo.kink_logic) < max_cor_index
       globalinfo.kink_logic(max_cor_index)=false;
   end
else
    globalinfo.kink_logic(max_cor_index)=false;
end
if isfield(globalinfo,'curve_logic')
   if length(globalinfo.curve_logic) < max_cor_index
       globalinfo.curve_logic(max_cor_index)=false;
   end
else
    globalinfo.curve_logic(max_cor_index)=false;
end
if isfield(globalinfo,'Selene_logic')
   if length(globalinfo.Selene_logic) < max_cor_index
       globalinfo.Selene_logic(max_cor_index)=false;
   end
else
    globalinfo.Selene_logic(max_cor_index)=false;
end
if isfield(globalinfo,'kinkoverwrite')
   if length(globalinfo.kinkoverwrite) < max_cor_index
       globalinfo.kinkoverwrite(max_cor_index)=false;
   end
else
    globalinfo.kinkoverwrite(max_cor_index)=false;
end

if isfield(globalinfo,'curveoverwrite')
   if length(globalinfo.curveoverwrite) < max_cor_index
       globalinfo.curveoverwrite(max_cor_index)=false;
   end
else
    globalinfo.curveoverwrite(max_cor_index)=false;
end

if ~isfield(McStasStr,'minimalist')
   McStasStr.minimalist=1;
   disp('options.minimalist not set, defaulting to 1')
end

if ~isfield(requirements,'source')
    requirements.source='cold';
    disp('requirements.source not set, defaulting to cold source')
end

globalinfo.modules = modules;
for ii = 1:length(modulelist)
globalinfo.modulelist(ii) = modulelist(length(modulelist)-ii+1);
end

% Add minlength and maxlength defaults for modules which have such.
module_defaults.E.minlength=2;
module_defaults.P.minlength=1;
module_defaults.K.maxlength=2.5;

module_defaults_names = fieldnames(module_defaults);

for ii = 1:length(globalinfo.modulelist)
   for jj = 1:length(module_defaults_names) 
        if strcmp(globalinfo.modules{globalinfo.modulelist(ii)},module_defaults_names{jj})
            default_properties = fieldnames(module_defaults.(module_defaults_names{jj}));
            for tt = 1:length(default_properties)
                switch default_properties{tt}
                    case 'minlength'
                        if globalinfo.minlength(ii) == 0
                           globalinfo.minlength(ii) = 1;
                           globalinfo.minlengthdata(ii) = module_defaults.(module_defaults_names{jj}).minlength;
                        end
                    case 'maxlength'
                        if globalinfo.maxlength(ii) == 0
                           globalinfo.maxlength(ii) = 1;
                           globalinfo.maxlengthdata(ii) = module_defaults.(module_defaults_names{jj}).maxlength;
                        end
                    otherwise
                        disp('ERROR, default code malfunctioning!')
                end
            end
        end
   end
end

% Making a list of raytracers in general.
% They all start a while loop if they are the first to be called, and add a
% optimized ratio between the two angles if they are after.
globalinfo.raytracers = false(1,max_cor_index);
if isfield(globalinfo,'kink_logic')
    globalinfo.raytracers = max([globalinfo.kink_logic ; globalinfo.raytracers]);
end
if isfield(globalinfo,'kinkoverwrite')
overwriten_logic = globalinfo.raytracers + globalinfo.kinkoverwrite == 2;
globalinfo.raytracers(overwriten_logic)=0;
end
if isfield(globalinfo,'curve_logic')
    globalinfo.raytracers = max([globalinfo.curve_logic ; globalinfo.raytracers]);
end
if isfield(globalinfo,'curveoverwrite')
overwriten_logic = globalinfo.raytracers + globalinfo.curveoverwrite == 2;
globalinfo.raytracers(overwriten_logic)=0;
end

filename='';
for i=1:length(modulelist)
    filename=[filename modules{modulelist(i)}];
end

disp(filename)

origfilename=filename;
unsure=1;j=0;
while (unsure==1)
    j=j+1;
    if exist(['./' Project_name '/' filename ])>0.5
       if exist(['./' Project_name '/' filename '/input.txt'])>0.5
           fid=fopen(['./' Project_name '/' filename '/input.txt']);
           %test=textscan(fid,'%s');
           test=fgetl(fid);
           fclose(fid)
           if strcmp(input,test)
               unsure = 0; % Just overwrite, same input!
               disp('overwriting old instrument folder without deleting')
           end
       else
           unsure=0;
       end
    else
        unsure=0;
    end
    
    if unsure==1; filename=[origfilename '_alt' num2str(j)]; end;
end


[status,message,messageid]=mkdir([Project_name '/' filename]);
[status,message,messageid]=copyfile([options_general.guide_bot_path '/components/*comp'],['./' Project_name '/' filename ]);
[status,message,messageid]=copyfile([options_general.guide_bot_path '/components/*c'],['./' Project_name '/' filename ]);
[status,message,messageid]=copyfile([options_general.guide_bot_path '/components/*h'],['./' Project_name '/' filename ]);
[status,message,messageid]=copyfile([options_general.guide_bot_path '/components/*txt'],['./' Project_name '/' filename ]);
[status,message,messageid]=copyfile([options_general.guide_bot_path '/components/*dat'],['./' Project_name '/' filename ]);

% Check if the first module restrict the entrance to exit.
fixedend_modules={'C' 'Cg'};


globalinfo.first_fixedend = 0;
for check = 1:length(fixedend_modules)
    if strcmp(modules{modulelist(1)},fixedend_modules{check})
        globalinfo.first_fixedend = 1;
        disp(['Module with end dimensions fixed by start dimensions used first (' fixedend_modules{check} ').'])
    end
end

% Time to create ghost_globalinfo.

% Pseudo code

% get list of los_divide los_start and los_end
% If there is no balance between los_start and los_end, add los_start=0 at last or a los_end=1 at 1, in some cases both.
% examples:
% S K S(los_divide=0.5) K S S : add both los_start and los_end
% S K S(los_start=0.5) K S(los_end=0.5) S: add los_start
% S K S(los_start=0.25,los_end=0.75) K S S: add los_start and los_end
% if the first los_ option from start is divide or end: add los_start
% if the first los_ option from end is divide or start: add los_end

% santitize:
% One module can not have more than one los_divide
% If there is no los_divide, there can be a los_start, los_end or both.
% Only S,P,E can have these options. Selene should be added at some point.
% There must be at least one los breaker between in every los interval.

% Starting from the sample, count the amout of ghost_indexes needed for
% each module.
% no los markings:    0
% los_divide:         1
% los_start:          1
% los_end:            1
% both end and start: 2 (should be in the order end then start).

% globalinfo.ghost_indexes should be written

% ghost_globalinfo should be a subset of the real globalinfo, but with the
% added indexes to make it the correct length.
% Additionally it should have the los_start and los_end list where there
% can be maximum one los_ option per index. The los_divide option disappears

% ghost_globalinfo.los_start_logic: logic vector (is there a los_start here)
% ghost_globalinfo.los_start_data: value vector (is there a los_start here)
% ghost_globalinfo.los_start_mode: value vector unit [1 = percentage, 2 = m]
% ghost_globalinfo.los_end_logic: logic vector (is there a los_start here)
% ghost_globalinfo.los_end_data: value vector (is there a los_start here)
% ghost_globalinfo.los_end_mode: value vector unit [1 = percentage, 2 = m]
% ghost_globalinfo.original_index: vector with the original index for the current module

% remember to pad the logics before this loop.
if isfield(globalinfo,'los_divide_logic')
   if length(globalinfo.los_divide_logic) < globalinfo.nummodules
       globalinfo.los_divide_logic(max_cor_index)=false;
   end
else
    globalinfo.los_divide_logic(max_cor_index)=false;
end
if isfield(globalinfo,'los_start_logic')
   if length(globalinfo.los_start_logic) < globalinfo.nummodules
       globalinfo.los_start_logic(max_cor_index)=false;
   end
else
    globalinfo.los_start_logic(max_cor_index)=false;
end
if isfield(globalinfo,'los_end_logic')
   if length(globalinfo.los_end_logic) < globalinfo.nummodules
       globalinfo.los_end_logic(max_cor_index)=false;
   end
else
    globalinfo.los_end_logic(max_cor_index)=false;
end


% This code ensure the correct min length and max length of modules using
% los_options.
for ii = 1:max_cor_index 
   rel_mode = 0; abs_mode = 0;
   rel_data = 0; abs_data = 0;
   abs_mode_start = 0; abs_mode_end = 0;
   abs_data_start = 0; abs_data_end = 0;
   rel_pos = 0; abs_start_pos = 0; abs_end_pos = 0;
      if globalinfo.los_start_logic(ii)
         if globalinfo.los_start_mode(ii) == 1
             rel_data(end+1) = globalinfo.los_start_data(ii);
             rel_pos(end+1) = globalinfo.los_start_option_num(ii);
         elseif globalinfo.los_start_mode(ii) == 2
             abs_data_start(end+1) = globalinfo.los_start_data(ii);
             abs_start_pos(end+1) = globalinfo.los_start_option_num(ii);
         elseif globalinfo.los_start_mode(ii) == 3
             abs_data_end(end+1) = globalinfo.los_start_data(ii);
             abs_end_pos(end+1) = globalinfo.los_start_option_num(ii);
         end
      end
      if globalinfo.los_end_logic(ii)
         if globalinfo.los_end_mode(ii) == 1 
             rel_data(end+1) = globalinfo.los_end_data(ii);
             rel_pos(end+1) = globalinfo.los_end_option_num(ii);
         elseif globalinfo.los_end_mode(ii) == 2
             abs_data_start(end+1) = globalinfo.los_end_data(ii);
             abs_start_pos(end+1) = globalinfo.los_end_option_num(ii);
         elseif globalinfo.los_end_mode(ii) == 3
             abs_data_end(end+1) = globalinfo.los_end_data(ii);
             abs_end_pos(end+1) = globalinfo.los_end_option_num(ii);
         end
      end
      if globalinfo.los_divide_logic(ii)
         if globalinfo.los_divide_mode(ii) == 1
             rel_data(end+1) = globalinfo.los_divide_data(ii);
             rel_pos(end+1) = globalinfo.los_divide_option_num(ii);
         elseif globalinfo.los_divide_mode(ii) == 2
             abs_data_start(end+1) = globalinfo.los_divide_data(ii);
             abs_start_pos(end+1) = globalinfo.los_divide_option_num(ii);
         elseif globalinfo.los_divide_mode(ii) == 3
             abs_data_end(end+1) = globalinfo.los_divide_data(ii);
             abs_end_pos(end+1) = globalinfo.los_divide_option_num(ii);
         end
      end
      
      if length(rel_data) > 1
         rel_mode = 1;
         rel_data = rel_data(2:end);
         rel_pos = rel_pos(2:end);
      else
         rel_mode = 0;
      end
      
      if length(abs_data_start) > 1
          abs_mode = 1; abs_mode_start = 1;
          abs_data_start = abs_data_start(2:end);
          abs_start_pos = abs_start_pos(2:end);
      end
      
      if length(abs_data_end) > 1
          abs_mode = 1; abs_mode_end = 1;
          abs_data_end = abs_data_end(2:end);
          abs_end_pos = abs_end_pos(2:end);
      end
          
      if rel_mode == 1 && abs_mode == 1
          temp_minlength = 0;temp_maxlength=0; % first index is bogus
          for j = 1:length(rel_data)
               if abs_mode_start
                for t = 1:length(abs_data_start)
                    if rel_pos(j) > abs_start_pos(t)
                        temp_minlength(end+1) = abs_data_start(t)/rel_data(j);
                    else
                        temp_maxlength(end+1) = abs_data_start(t)/rel_data(j);
                    end
                end
               end
               if abs_mode_end
                for t = 1:length(abs_data_end)
                    if rel_pos(j) < abs_end_pos(t)
                        temp_minlength(end+1) = abs_data_end(t)/(1-rel_data(j));
                    else
                        temp_maxlength(end+1) = abs_data_end(t)/(1-rel_data(j));
                    end
                end
               end
          end
          
          
          
          if max(temp_minlength) > globalinfo.minlengthdata(ii)
              globalinfo.minlength(ii)=1;
              globalinfo.minlengthdata(ii) = max(temp_minlength);
              %disp(['changed globalinfo.minlengthdata(' num2str(ii) ') to ' num2str(max(temp_minlength)) ' (rel).'])
          end
          
          if length(temp_maxlength)>1
           temp_maxlength = temp_maxlength(2:end); % remove bogus
           % doesnt matter for minlength, as it is only used in max.
           if globalinfo.maxlengthdata(ii) == 0;
              globalinfo.maxlength(ii)=1;
              globalinfo.maxlengthdata(ii) = min(temp_maxlength);
              %disp(['changed globalinfo.maxlengthdata(' num2str(ii) ') to ' num2str(min(temp_maxlength)) ' (rel).'])
           else    
            if min(temp_maxlength) < globalinfo.maxlengthdata(ii)
              globalinfo.maxlength(ii)=1;
              globalinfo.maxlengthdata(ii) = min(temp_maxlength);
              %disp(['changed globalinfo.maxlengthdata(' num2str(ii) ') to ' num2str(min(temp_maxlength)) ' (rel).'])
            end
           end
          end
      end
      
      if abs_mode_start && abs_mode_end
         % In this case there is more than one abs_data in this module.
         % if a m_s is before a m_e, there needs to be the sum of these lengths in minlength! 
         temp_minlength = 0;
         temp_maxlength = 0;
         for jj = 1:length(abs_data_start)
             for tt = 1:length(abs_data_end)
                 if abs_start_pos(jj) < abs_end_pos(tt)
                     temp_minlength = max([temp_minlength abs_data_start(jj) + abs_data_end(tt)]);
                 else
                 % Unlikely situation, but should be accounted for
                     if temp_maxlength == 0
                         temp_maxlength = abs_data_start(jj) + abs_data_end(tt);
                     else
                         temp_maxlength = min([temp_maxlength abs_data_start(jj) + abs_data_end(tt)]);
                     end
                 end
             end
         end
         
         if temp_minlength > globalinfo.minlengthdata(ii)
              globalinfo.minlengthdata(ii) = temp_minlength;
              %disp(['changed globalinfo.minlengthdata(' num2str(ii) ') to ' num2str(temp_minlength) ' (abs).'])
          end
          
          if temp_maxlength ~= 0
           if globalinfo.maxlengthdata(ii) == 0;
              globalinfo.maxlength(ii)=1;
              globalinfo.maxlengthdata(ii) = temp_maxlength;
              %disp(['changed globalinfo.maxlengthdata(' num2str(ii) ') to ' num2str(temp_maxlength) ' (abs).'])
           else    
            if temp_maxlength < globalinfo.maxlengthdata(ii)
              globalinfo.maxlength(ii)=1;
              globalinfo.maxlengthdata(ii) = temp_maxlength;
              %disp(['changed globalinfo.maxlengthdata(' num2str(ii) ') to ' num2str(temp_maxlength) ' (abs).'])
            end
           end
          end    
         
      end
      
      for j = 1:length(rel_data)
         for t = 1:length(rel_data)
             if t~=j
                 if rel_pos(t) > rel_pos(j)
                    if rel_data(t) < rel_data(j)
                        disp('ERROR, los option with relative position in wrong order!')
                    end
                 elseif rel_pos(t) < rel_pos(j)
                    if rel_data(t) > rel_data(j)
                        disp('ERROR, los option with relative position in wrong order!')
                    end
                 end
             end
         end         
      end
end

ghost_index = 0;
for index=1:globalinfo.nummodules
   los_start_found = 0;
   if globalinfo.los_divide_logic(index)
       ghost_index = ghost_index + 1;
       ghost_globalinfo.real_index(ghost_index) = index;
       ghost_globalinfo.los_start_logic(ghost_index)=true;
       ghost_globalinfo.los_start_data(ghost_index)=globalinfo.los_divide_data(index);
       ghost_globalinfo.los_start_mode(ghost_index)=globalinfo.los_divide_mode(index);
       
       ghost_index = ghost_index + 1;
       % Add here as well.
       ghost_globalinfo.real_index(ghost_index) = index;
       ghost_globalinfo.los_end_logic(ghost_index)=true;
       ghost_globalinfo.los_end_data(ghost_index)=globalinfo.los_divide_data(index);
       ghost_globalinfo.los_end_mode(ghost_index)=globalinfo.los_divide_mode(index);
   else
       % check if start or end should be first here.
       
       % This can be wrong, and cause start start end end instead of start end start end
       if start_first(index) == 0
           % with end first, a single module will have end then start when seen from the source.
           if globalinfo.los_start_logic(index)
               ghost_index = ghost_index + 1;
               % Add here as well.
               ghost_globalinfo.real_index(ghost_index) = index;
               ghost_globalinfo.los_start_logic(ghost_index)=true;
               ghost_globalinfo.los_start_data(ghost_index)=globalinfo.los_start_data(index);
               ghost_globalinfo.los_start_mode(ghost_index)=globalinfo.los_start_mode(index);
           end        
           if globalinfo.los_end_logic(index)
               ghost_index = ghost_index + 1;
               % Add here as well.
               ghost_globalinfo.real_index(ghost_index) = index;
               ghost_globalinfo.los_end_logic(ghost_index)=true;
               ghost_globalinfo.los_end_data(ghost_index)=globalinfo.los_end_data(index);
               ghost_globalinfo.los_end_mode(ghost_index)=globalinfo.los_end_mode(index);
           end
       else
           % with end first, a single module will have start then end when seen from the source.
           if globalinfo.los_end_logic(index)
               ghost_index = ghost_index + 1;
               % Add here as well.
               ghost_globalinfo.real_index(ghost_index) = index;
               ghost_globalinfo.los_end_logic(ghost_index)=true;
               ghost_globalinfo.los_end_data(ghost_index)=globalinfo.los_end_data(index);
               ghost_globalinfo.los_end_mode(ghost_index)=globalinfo.los_end_mode(index);
           end
           if globalinfo.los_start_logic(index)
               ghost_index = ghost_index + 1;
               % Add here as well.
               ghost_globalinfo.real_index(ghost_index) = index;
               ghost_globalinfo.los_start_logic(ghost_index)=true;
               ghost_globalinfo.los_start_data(ghost_index)=globalinfo.los_start_data(index);
               ghost_globalinfo.los_start_mode(ghost_index)=globalinfo.los_start_mode(index);
           end        
       end
            
   end
   ghost_index = ghost_index + 1;
   % Add code here to transfer from globalinfo to ghost_globalinfo
   % ghost_globalinfo.rotlogic(ghost_index) = globalinfo.rotlogic(index);
   ghost_globalinfo.real_index(ghost_index) = index;
end
ghost_globalinfo.length=length(ghost_globalinfo.real_index);



if sum(globalinfo.raytracers) == 0 % wrong.
    % DISABLE LOS SYSTEM
    ghost_globalinfo.los_logic(1:ghost_globalinfo.length) = 0;
else
    % If the los system should not be disabled, then all these variables
    % are meaningful and allocated.
    if ~isfield(ghost_globalinfo,'los_start_logic')
        ghost_globalinfo.los_start_logic(ghost_globalinfo.length) = 1;
        ghost_globalinfo.los_start_data(ghost_globalinfo.length) = 0;
        ghost_globalinfo.los_start_mode(ghost_globalinfo.length) = 1;
    end
    if ~isfield(ghost_globalinfo,'los_end_logic')
        ghost_globalinfo.los_end_logic(1) = 1;
        ghost_globalinfo.los_end_data(1) = 1;
        ghost_globalinfo.los_end_mode(1) = 1;
    end
        
    if find(ghost_globalinfo.los_start_logic,1) < find(ghost_globalinfo.los_end_logic,1)
       % In case of unbalanced 
        ghost_globalinfo.los_end_logic(1) = 1;
        ghost_globalinfo.los_end_data(1) = 1;
        ghost_globalinfo.los_end_mode(1) = 1;
    end
    ghost_globalinfo.los_end_indexs = find(ghost_globalinfo.los_end_logic);
    ghost_globalinfo.los_start_indexs = find(ghost_globalinfo.los_start_logic);
    
    if ghost_globalinfo.los_end_indexs(end) > ghost_globalinfo.los_start_indexs(end)
        % In this case, there is a missing los_start, add it at the start of the guide
        % It is used that there can not be a los_end or los_divide in the first module
        ghost_globalinfo.los_start_logic(ghost_globalinfo.length) = 1;
        ghost_globalinfo.los_start_data(ghost_globalinfo.length) = 0;
        ghost_globalinfo.los_start_mode(ghost_globalinfo.length) = 1;
        % Needs to be performed again to add the start which was just added
        ghost_globalinfo.los_start_indexs = find(ghost_globalinfo.los_start_logic);
    end
%else
%        % los_end at sample
%        ghost_globalinfo.los_end_logic(1) = 1;
%        ghost_globalinfo.los_end_data(1) = 1;
%        ghost_globalinfo.los_end_mode(1) = 1;
%        ghost_globalinfo.los_end_indexs = find(ghost_globalinfo.los_end_logic);
        
        % los_start at guide start
%        ghost_globalinfo.los_start_logic(ghost_globalinfo.length) = 1;
%        ghost_globalinfo.los_start_data(ghost_globalinfo.length) = 0;
%        ghost_globalinfo.los_start_mode(ghost_globalinfo.length) = 1;
%        ghost_globalinfo.los_start_indexs = find(ghost_globalinfo.los_start_logic);
    

% Pad ghost_globalinfo logicals
if isfield(ghost_globalinfo,'los_start_logic')
   if length(ghost_globalinfo.los_start_logic) < ghost_globalinfo.length
       ghost_globalinfo.los_start_logic(ghost_globalinfo.length)=false;
   end
else
    globalinfo.los_start_logic(ghost_globalinfo.length)=false;
end
if isfield(ghost_globalinfo,'los_end_logic')
   if length(ghost_globalinfo.los_end_logic) < ghost_globalinfo.length
       ghost_globalinfo.los_end_logic(ghost_globalinfo.length)=false;
   end
else
    globalinfo.los_end_logic(ghost_globalinfo.length)=false;
end

% if mode 2 is used, remember to add a minlength option to that module

% Code for determining the sections of the guides where line of sight
% should be broken. Basicly all section from los_start_logic to
% los_end_logic.

% Notice the first l_s is not l_s = 0.
% If the guide is S(l_s) K S(l_e,l_s) K S(l_e) S
% Index:          6      5 4          3 2      1
% Ghost_index     10 9   8 7 6   5    4 3 2    1
% los1:              9   8 7 6   
% los2:                          5    4 3 2

for los_index = 1:length(ghost_globalinfo.los_start_indexs)
    los_index_control = length(ghost_globalinfo.los_start_indexs) - los_index+1;
    start_value =  ghost_globalinfo.los_start_indexs(los_index);
    end_value = ghost_globalinfo.los_end_indexs(los_index);
    
    if ghost_globalinfo.los_start_mode(start_value) == 1
        if ghost_globalinfo.los_start_data(start_value) == 1
            % Do not take the first ghost element into consideration
            start_value = start_value - 1;
            %disp(' --------------------------------- * ----- Changed active_los start')
        end
    end
    if ghost_globalinfo.los_end_mode(end_value) == 1
        if ghost_globalinfo.los_end_data(end_value) == 0
            % No not take the last ghost element into consideration
            end_value = end_value + 1;
            %disp(' --------------------------------- * ----- Changed active_los end')
        end
    end
    % Old version before fix
    %ghost_globalinfo.active_los{los_index_control} = ghost_globalinfo.los_start_indexs(los_index):-1:ghost_globalinfo.los_end_indexs(los_index);
    ghost_globalinfo.active_los{los_index_control} = start_value:-1:end_value;
end

% Add sanatizing here.
% Check if there is a los breaker in each active_los section.

% Easy way to avoid running the raytracer if there is no need for it and
% avoid accesing undeclared variables.
% Remember length(active_los) is now meaningless, which is okay.
% ghost_globalinfo.active_los{end+1} = 0;

% Tunrs out that it is not needed to have seperate los_start_data and
% los_end_data together with los_end_mode and los_start_mode.

% They are collected to the unified los_data and los_mode here:

for index = 1:ghost_globalinfo.length;
    if ghost_globalinfo.los_start_logic(index)
        ghost_globalinfo.los_logic(index) = ghost_globalinfo.los_start_logic(index);
        ghost_globalinfo.los_data(index) = ghost_globalinfo.los_start_data(index);
        ghost_globalinfo.los_mode(index) = ghost_globalinfo.los_start_mode(index);
    elseif ghost_globalinfo.los_end_logic(index)
        ghost_globalinfo.los_logic(index) = ghost_globalinfo.los_end_logic(index);
        ghost_globalinfo.los_data(index) = ghost_globalinfo.los_end_data(index);
        ghost_globalinfo.los_mode(index) = ghost_globalinfo.los_end_mode(index);
    else
        ghost_globalinfo.los_logic(index) = 0;
        ghost_globalinfo.los_data(index) = 0;
        ghost_globalinfo.los_mode(index) = 0;
    end
end

end % ends los enabled section


if isfield(ghost_globalinfo,'active_los')
  for ii = 1:length(ghost_globalinfo.active_los)
   if sum(globalinfo.Selene_logic(ghost_globalinfo.real_index(ghost_globalinfo.active_los{ii}))) > 0
       disp('ERROR, Selene can no be part of a los section!')
   end
  end
end


% initialize the McStasStr structure which contains the info writen to
% the instr file.
McStasStr.trace='';
McStasStr.trace_seg='';
McStasStr.declare{1}='u';
%McStasStr.declare{2}='mod_x';
%McStasStr.declare{3}='mod_y';
McStasStr.declare{2}='Lambda0';
McStasStr.declare{3}='dLambda';
McStasStr.declare{4}='var_divreq_x';
McStasStr.declare{5}='var_divreq_y';
McStasStr.declare{6}='x_div';
McStasStr.declare{7}='y_div';
if isfield(options_general,'fom_weight')
    McStasStr.declare{8}='this_lambda';
end 
McStasStr.declareint{1}='flag';
McStasStr.declareint{2}='i';
McStasStr.declareint{3}='part';
enddeclare=7; % For brilliance refference (get declare 1:7)
McStasStr.initialize='';
% This way of adding input to McStasStr is not optimal.
% Should be a loop over the input instead.
% Or should be done inside the scan loop to enable use of demandscan.
McStasStr.input{1}='sizeX';
if scan.demands(3)==0 && scan.locked.list(3) == 0
McStasStr.inputvalue(1)=demands.Hsize/100;
else
McStasStr.inputvalue(1)=demands.Hsize(1)/100;
end
McStasStr.input{2}='sizeY';
if scan.demands(4)==0 && scan.locked.list(4) == 0
McStasStr.inputvalue(2)=demands.Vsize/100;
else
McStasStr.inputvalue(2)=demands.Vsize(1)/100;
end
McStasStr.input{3}='divreq_x';
if scan.demands(1)==0 && scan.locked.list(1) == 0
McStasStr.inputvalue(3)=demands.Hdiv;
else
McStasStr.inputvalue(3)=demands.Hdiv(1);
end
McStasStr.input{4}='divreq_y';
if scan.demands(2)==0 && scan.locked.list(2) == 0
McStasStr.inputvalue(4)=demands.Vdiv;
else
McStasStr.inputvalue(4)=demands.Vdiv(2);
end
% The Selene guides includes the sample space, set to sample_dist to zero.
% This is done inside the Selene module, by using that McStasStr.input{5}
% is the string which controls sample_dist. Hence this number is stored in
% globalinfo.
globalinfo.sample_dist_index = 5;
McStasStr.input{globalinfo.sample_dist_index}='sample_dist';
McStasStr.inputvalue(globalinfo.sample_dist_index)=demands.Dist;
McStasStr.input{6}='WaveMin';
if scan.demands(5)==0 && scan.locked.list(5) == 0
McStasStr.inputvalue(6)=demands.WaveLmin;
else
McStasStr.inputvalue(6)=demands.WaveLmin(1);
end
McStasStr.input{7}='WaveMax';
if scan.demands(6)==0 && scan.locked.list(6) == 0
McStasStr.inputvalue(7)=demands.WaveLmax;
else
McStasStr.inputvalue(7)=demands.WaveLmax(1);
end
McStasStr.input{8}='mod_x';
if scan.demands(7)==0 && scan.locked.list(7) == 0
McStasStr.inputvalue(8)=requirements.moderator_size_x;
else
McStasStr.inputvalue(8)=requirements.moderator_size_x(1);
end
McStasStr.input{9}='mod_y';
if scan.demands(8)==0 && scan.locked.list(8) == 0
McStasStr.inputvalue(9)=requirements.moderator_size_y;
else
McStasStr.inputvalue(9)=requirements.moderator_size_y(1);
end
McStasStr.input{10}='minimalist_factor';
if scan.demands(9)==0 && scan.locked.list(9) == 0
McStasStr.inputvalue(10)=requirements.minimalist_factor;
else
McStasStr.inputvalue(10)=requirements.minimalist_factor(1);
end
McStasStr.input{11}='Mod_sample';
McStasStr.inputvalue(11)=demands.Mod_sample;
McStasStr.input{12}='closest_element';
McStasStr.inputvalue(12)=requirements.closest_element;


endinput=9; % For brilliance refference (get input 1:7)


if strcmp(source_component,'ESS_moderator_long')
% Declarations and input for ess instrument files
McStasStr.declare_ess{1}='ESS_cyl_radius'; % Added recently
McStasStr.declare_ess{end+1}='ESS_angle';
McStasStr.declare_ess{end+1}='target_point[2][6]';
McStasStr.declare_ess{end+1}='guide_point[2]';
McStasStr.declare_ess{end+1}='mc_cold_frac';
McStasStr.declareint_ess{1}='mod_target_index';

McStasStr.input_ess{1} = 'beamport_angle_zero_center';

if isfield(options_general,'beamport')
    McStasStr.inputvalue_ess(1)=options_general.beamport; 
else
    McStasStr.inputvalue_ess(1)=0;
end
McStasStr.input_ess{2} = 'manual_guide_angle';
McStasStr.inputvalue_ess(2)=0;
McStasStr.input_ess{3} = 'thermal_wing_length';
McStasStr.inputvalue_ess(3)=0.12;
McStasStr.input_ess{4} = 'thermal_wing_angle'; % not active yet
McStasStr.inputvalue_ess(4)=30;
McStasStr.input_ess{5} = 'mod_target_index_input';
McStasStr.optimize_ess(5) = 0; % Nothing will be left without inputvalue

if strcmp(requirements.source,'cold')
McStasStr.inputvalue_ess(5)=1;
elseif strcmp(requirements.source,'thermal')
McStasStr.inputvalue_ess(5)=2;
elseif strcmp(requirements.source,'bispectral')
McStasStr.inputvalue_ess(5)=3;
else
McStasStr.inputvalue_ess(5)=3;
disp('Source needs to be cold, thermal or bispectral! Using defaults.')
end
elseif strcmp(requirements.source,'ESS_pancake')
    % Nothing to declare in this case.
    % char command_string is declared manually later
elseif strcmp(requirements.source,'virtual_in')
    % Nothing to declare in this case.
end

% input for fixedstart 'num' og fixed 'length'

if length(McStasStr.input)>length(McStasStr.inputvalue)
    disp('ERROR, not alll permanant inputs have standard values')
end

% The next logical section determines how guide_start should be declared
% and optimized including the ranges possible
linp=length(McStasStr.input);
% TO DO: Calculate the max guide_start based on phase-space requirements
% The above could take the place of requirement.closest_element

if globalinfo.fixedstart(end)==1
    McStasStr.input{linp+1}='guide_start';
    McStasStr.inputvalue(linp+1)=globalinfo.fixedstartpoint(end);
elseif (strcmp(modules{modulelist(1)},'Selene') && globalinfo.selene_optimize_first_slit == 0 && globalinfo.selene_moderator_focus == 0)
    % If total standard selene (no options)
    McStasStr.input{linp+1}='guide_start';
    McStasStr.optimize(linp+1)=0;
    McStasStr.inputvalue(linp+1)=requirements.closest_element;
elseif (strcmp(modules{modulelist(1)},'Selene') && globalinfo.selene_moderator_focus == 1)
    % If anything does moderator focus
    McStasStr.input{linp+1}='guide_start';
    McStasStr.optimize(linp+1)=0;
    McStasStr.inputvalue(linp+1)=0.001;
    % Selene with moderator_focus = 0 and optimize_first_slit = 1 will go
    % to normal guide_start calculations.
else
    
    if sum(globalinfo.fixedstart(1:end))==0 % There is no fixed start
        % Allocate guide_start as optimized input as normal, but
        % check if the cummulative fixed length interferes with the
        % range of guide_start
        
        McStasStr.input{linp+1}='guide_start';
        McStasStr.inputvalue(linp+1)=requirements.closest_element;
        McStasStr.optimize(linp+1)=1;
        McStasStr.optimvals.min(linp+1)=requirements.closest_element;
        if demands.Mod_sample-demands.Dist-sum(globalinfo.fixedlengthdata(1:last))-sum(1-globalinfo.fixedlength(1:last))*0.1<requirements.latest_start
            % modify range
            McStasStr.optimvals.max(linp+1)=demands.Mod_sample-demands.Dist-sum(globalinfo.fixedlengthdata(1:last))-sum(1-globalinfo.fixedlength(1:last))*0.1;
            % above should include minlength!
        else
            McStasStr.optimvals.max(linp+1)=requirements.latest_start;
        end
        McStasStr.optimvals.guess(linp+1)=mean([McStasStr.optimvals.max(linp+1) McStasStr.optimvals.min(linp+1)]);
    else % There is a fixed start module
        fppos=find(globalinfo.fixedstart(1:end));
        lastfp=fppos(end);
        
        %check if all modules between last and lastfp have fixed length
        if sum(globalinfo.fixedlength(lastfp:last-1))>last-lastfp-0.5
            if globalinfo.fixedlength(end)==1
                % In this case guide_start is actually fixed
                % guide_start=lastfixed-sum length (from lastfp to last)
                McStasStr.input{linp+1}='guide_start';
                McStasStr.inputvalue(linp+1)=globalinfo.fixedstartpoint(lastfp)-sum(globalinfo.fixedlengthdata(lastfp:last));
            else
                tmp_check=globalinfo.fixedstartpoint(lastfp)-sum(globalinfo.fixedlengthdata(lastfp+1:last))-sum(globalinfo.minlengthdata(lastfp+1:last)); % lastfp + 1 because the length is in the other direction
                % Allocate as optimized input and check if there is another
                % max
                McStasStr.input{linp+1}='guide_start';
                McStasStr.inputvalue(linp+1)=requirements.closest_element;
                McStasStr.optimize(linp+1)=1;
                McStasStr.optimvals.min(linp+1)=requirements.closest_element;
                if tmp_check<requirements.latest_start
                    % In this case guide_start should not be above tmp_check
                    % remember warning here
                    McStasStr.optimvals.max(linp+1)=tmp_check;
                else
                    McStasStr.optimvals.max(linp+1)=requirements.latest_start;

                end
                McStasStr.optimvals.guess(linp+1)=mean([McStasStr.optimvals.max(linp+1) McStasStr.optimvals.min(linp+1)]);
            end
        else
            % Even though there is modules with free length left, it might
            % not be enough to allow the normal range of guide_start!

            % Need to add the possibility for minlength here.
            
            McStasStr.input{linp+1}='guide_start';
            McStasStr.inputvalue(linp+1)=requirements.closest_element;
            McStasStr.optimize(linp+1)=1;
            McStasStr.optimvals.min(linp+1)=requirements.closest_element;
            %if globalinfo.fixedstartpoint(lastfp)-sum(globalinfo.fixedlengthdata(lastfp:last))-sum(1-globalinfo.fixedlength(lastfp+1:last))*0.1<requirements.latest_start
            if globalinfo.fixedstartpoint(lastfp)-sum(globalinfo.fixedlengthdata(lastfp:last))-sum(globalinfo.minlengthdata(lastfp:last))-sum(1-globalinfo.fixedlength(lastfp+1:last))*0.1<requirements.latest_start
                % modify range
                %McStasStr.optimvals.max(linp+1)=globalinfo.fixedstartpoint(lastfp)-sum(globalinfo.fixedlengthdata(lastfp:last))-sum(1-globalinfo.fixedlength(lastfp+1:last))*0.1;
                McStasStr.optimvals.max(linp+1)=globalinfo.fixedstartpoint(lastfp)-sum(globalinfo.fixedlengthdata(lastfp:last))-sum(globalinfo.minlengthdata(lastfp:last))-sum(1-globalinfo.fixedlength(lastfp+1:last))*0.1;
            else
                McStasStr.optimvals.max(linp+1)=requirements.latest_start;
            end
            McStasStr.optimvals.guess(linp+1)=mean([McStasStr.optimvals.max(linp+1) McStasStr.optimvals.min(linp+1)]);
        end
    end
end
% The above calculations of the guide start range assumes no problems with
% min/max startpoints, which is not exactly the case. Should be expanded to
% take care of these situations.

if sum(ismember(McStasStr.input,'guide_start'))>0.5 && sum(ismember(fieldnames(McStasStr),'optimize')) > 0.5
   guide_start_index=find(ismember(McStasStr.input,'guide_start'),1);
   if length(McStasStr.optimize) >= guide_start_index
    if McStasStr.optimize(guide_start_index) == 1
     if globalinfo.maxstart(last)==1
       if globalinfo.maxstartpoint(last)<McStasStr.optimvals.max(guide_start_index)
            McStasStr.optimvals.max(guide_start_index)=globalinfo.maxstartpoint(last);
       end
     end
     if globalinfo.minstart(last)==1
       if globalinfo.minstartpoint(last)>McStasStr.optimvals.min(guide_start_index)
            McStasStr.optimvals.min(guide_start_index)=globalinfo.mintartpoint(last);
       end
     end
    end
   end
end



% ADD: if fixedstart is chosen for the module closest to the moderator,
% change that fixed start to guide_start

for i=find(globalinfo.fixedstart)
McStasStr.input{end+1}=['fixedStart' num2str(i)];
McStasStr.inputvalue(end+1)=globalinfo.fixedstartpoint(i);
end
for i=find(globalinfo.maxstart)
McStasStr.input{end+1}=['maxStart' num2str(i)];
McStasStr.inputvalue(end+1)=globalinfo.maxstartpoint(i);
end
for i=find(globalinfo.minstart)
McStasStr.input{end+1}=['minStart' num2str(i)];
McStasStr.inputvalue(end+1)=globalinfo.minstartpoint(i);
end


%length can be declare in other places, maybe not here
for i=find(globalinfo.fixedlength)
McStasStr.input{end+1}=['length' num2str(i)];
McStasStr.inputvalue(end+1)=globalinfo.fixedlengthdata(i);
end
for i=find(globalinfo.maxlength)
McStasStr.input{end+1}=['maxlength' num2str(i)];
McStasStr.inputvalue(end+1)=globalinfo.maxlengthdata(i);
end
for i=find(globalinfo.minlength)
McStasStr.input{end+1}=['minlength' num2str(i)];
McStasStr.inputvalue(end+1)=globalinfo.minlengthdata(i);
end

% INSERT FUNCTION FOR GUESSTIMATING
% This function should take the globalinfo string and if needed the input
% string in order to calculate how a decent guide could look and use this
% as guesses. The most essential part of this is the length of each module,
% as guide kink guide, should have a guess og two guides of equal length
% and a small kink length.


% SMALL CODE FOR VISUALIZING ACTIVE SECTION
if isfield(ghost_globalinfo,'active_los')
disp('Showing line of sight ========')
string = '';
for ii = 1:length(globalinfo.modulelist)
   ii_used = length(globalinfo.modulelist) - ii + 1;
   string = [string ' ' globalinfo.modules{globalinfo.modulelist(ii_used)}];
end  
disp(string)
for ii = 1:length(ghost_globalinfo.active_los)
    los_range = ghost_globalinfo.real_index(ghost_globalinfo.active_los{ii});
    globalinfo.nummodules;
    white_spaces = globalinfo.nummodules - los_range(1);
    string = '';
    for jj = 1:white_spaces; string = [string '  ']; end
    for jj = 1:length(los_range)
        if jj == 1
         string = [string ' -'];    
        else
        if los_range(jj) ~= los_range(jj-1)
         string = [string '--'];
        end
        end
    end
    disp(string)
end
disp('==============================')
end

% Can be changed to 1 to go back to old raytracer.
globalinfo.old_point_calc = 0;

% This could give problems in the end.
% It real_index(index+1) is used to identify if the last element was a line
% of sight breaker which needs special attention in point_calc.
% It shoould only be used to check that is not the case. Alternertive wanted.
%ghost_globalinfo.real_index(end+1)=last+1;

%guesstimate = guesstimator(

% init pcalcstring
McStasStr.pcalcstring='';

% ------------------------------------------------------------------------
% main loop writing the optical components to the strings
% Initialize the needed 
loopend=length(modulelist)-1;
last=length(modulelist);
ghost_last = ghost_globalinfo.length;
ghost_globalinfo.los_section = 0;
point_calc_runs = 0;
if isfield(ghost_globalinfo,'active_los')
   point_calc_run_limit = length(ghost_globalinfo.active_los);
else
   point_calc_run_limit = 0;
end
for loop=0:loopend
    %index=loop+1;
    index=length(modulelist)-loop;
    
    if globalinfo.old_point_calc == 0
    if point_calc_runs < point_calc_run_limit
        if ghost_globalinfo.real_index(ghost_globalinfo.active_los{ghost_globalinfo.los_section+1}(1)) == index
          point_calc_runs = point_calc_runs + 1;
          ghost_globalinfo.los_section = ghost_globalinfo.los_section + 1;
          McStasStr.pcalcstring = '';
          % These two lines are handy for debugging. Might enable a verbose
          %disp(['running p_calc_loop' num2str(point_calc_runs) ' with los_section ' num2str(ghost_globalinfo.los_section)])
          %disp(num2str(ghost_globalinfo.active_los{ghost_globalinfo.los_section}))
          for ghost_index = ghost_globalinfo.active_los{ghost_globalinfo.los_section}
              [McStasStr] = point_calc_writer(McStasStr,ghost_index,ghost_last,globalinfo,ghost_globalinfo,defaults,1);
               % Put this directly into the initialize string
          end
          McStasStr.initialize = [McStasStr.pcalcstring '\n\n' McStasStr.initialize] ;
        end
    end
    end
    
    % call a specific module
    % it is named something like E_module.m
    run=[modules{modulelist(loop+1)} '_module']; 
    %disp([run options{loop+1}])
    disp(run)
    disp(options{loop+1})
    [McStasStr] = eval([run '(McStasStr,index,last,requirements,demands,options{loop+1},globalinfo,ghost_globalinfo,defaults)']);
    if globalinfo.old_point_calc %1==1 %globalinfo.kinkoverwrite==0 && sum(ismember(fieldnames(globalinfo),'kinks'))>0.5
       [McStasStr] = old_point_calc_writer_old(McStasStr,index,last,options{loop+1},globalinfo,defaults);
    end
end
% ------------------------------------------------------------------------


% The above line of sight raytracers used parts of the overall. 
% For debugging and plotting purposes.

McStasStr.pcalcstring='';
for ghost_index = ghost_globalinfo.length:-1:1
   [McStasStr] = point_calc_writer(McStasStr,ghost_index,ghost_last,globalinfo,ghost_globalinfo,defaults,0); % the last zero disables the raytracer parts
end
McStasStr.initialize = [McStasStr.initialize '\n\n' McStasStr.pcalcstring];


% Build up McStasStr.visualizer_str.
McStasStr.visualizer_str='';
for ghost_index = ghost_globalinfo.length:-1:1
    [McStasStr] = visualizer_output(McStasStr,ghost_index,ghost_last,globalinfo,ghost_globalinfo,defaults); 
end



% Modify the initialize section
l{1}='dLambda = 0.5*(WaveMax - WaveMin);';
l{end+1}='Lambda0 = dLambda+WaveMin;';
l{end+1}='';
l{end+1}='var_divreq_x = divreq_x;';
l{end+1}='var_divreq_y = divreq_y;';
l{end+1}='';
l{end+1}='u=1e-4;';

if globalinfo.fixedstart(end)==1;
l{end+1}='';
l{end+1}='// the following line is needed in order to allow the control of the guide start from the guide_start parameter, instead of the forced start used in the input';
l{end+1}=['fixedStart' num2str(globalinfo.nummodules) ' = guide_start;'];    
end

initstring='';
for i=1:length(l)
    initstring=[initstring l{i} '\n'];
end
clear l;
saved_inject_point = McStasStr.initialize;
McStasStr.initialize=[initstring '\n\n' McStasStr.initialize];
% needs temp name in order to not overwrite data needed to initialize_ess_seg.
%McStasStr.initialize_seg_final=['\nINITIALIZE\n%%{\n' McStasStr.initialize '\n' McStasStr.initialize_seg '\n' McStasStr.pcalcstring '\n%%}\n'];
McStasStr.initialize_seg_final=['\nINITIALIZE\n%%{\n' McStasStr.initialize '\n' McStasStr.initialize_seg];
%McStasStr.initialize=['\nINITIALIZE\n%%{\n' McStasStr.initialize '\n' McStasStr.pcalcstring '\n%%}\n'];
McStasStr.initialize=['\nINITIALIZE\n%%{\n' McStasStr.initialize];

% inistring is used in the following    

if strcmp(source_component,'ESS_moderator_long')
% implementation using mod_target_index and modern ESS source
l{1}='ESS_cyl_radius = mod_x;';
l{end+1} = 'mod_target_index = round(mod_target_index_input);';
l{end+1} = '// calculations nessecary for aiming the guide at the moderator';
l{end+1} = 'if (mod_target_index == 0) {';
l{end+1} = '   ESS_angle = manual_guide_angle;'; % default = 0.
l{end+1} = '}';
l{end+1} = 'else {';
l{end+1} = '  // calculate target points / guide point';
l{end+1} = '  // target for cold moderator';
l{end+1} = '  target_point[0][1] = ESS_cyl_radius;';
l{end+1} = '  target_point[1][1] = 0;';
l{end+1} = '  // targets for thermal moderator';
l{end+1} = '  target_point[0][2] = (ESS_cyl_radius+0.5*thermal_wing_length)*cos(DEG2RAD*thermal_wing_angle);';
l{end+1} = '  target_point[1][2] = (ESS_cyl_radius+0.5*thermal_wing_length)*sin(DEG2RAD*thermal_wing_angle);';
l{end+1} = '  target_point[0][3] = (ESS_cyl_radius+0.5*thermal_wing_length)*cos(DEG2RAD*thermal_wing_angle);';
l{end+1} = '  target_point[1][3] = -(ESS_cyl_radius+0.5*thermal_wing_length)*sin(DEG2RAD*thermal_wing_angle);';
l{end+1} = '  // targets for the point between cold and thermal moderator';
l{end+1} = '  target_point[0][4] = ESS_cyl_radius*cos(DEG2RAD*thermal_wing_angle);';
l{end+1} = '  target_point[1][4] = ESS_cyl_radius*sin(DEG2RAD*thermal_wing_angle);';
l{end+1} = '  target_point[0][5] = ESS_cyl_radius*cos(DEG2RAD*thermal_wing_angle);';
l{end+1} = '  target_point[1][5] = -ESS_cyl_radius*sin(DEG2RAD*thermal_wing_angle);';
l{end+1} = '  // guide_point';
l{end+1} = '  guide_point[0] = (ESS_cyl_radius + guide_start)*cos(DEG2RAD*beamport_angle_zero_center);';
l{end+1} = '  guide_point[1] = (ESS_cyl_radius + guide_start)*sin(DEG2RAD*beamport_angle_zero_center);';
l{end+1} = '';
l{end+1} = '  ESS_angle = - RAD2DEG*(atan(guide_point[1]/guide_point[0]) - atan((guide_point[1]-target_point[1][mod_target_index])/(guide_point[0]-target_point[0][mod_target_index])) );';
l{end+1} = '}';
l{end+1} = 'if (mod_target_index == 0)';
if strcmp(requirements.source,'cold')
l{end+1} = '   mc_cold_frac = 1;'; % normally 0.7 because of bug
elseif strcmp(requirements.source,'thermal')
l{end+1} = '   mc_cold_frac = 0;'; % normally 0.3 because of bug
elseif strcmp(requirements.source,'bispectral')
l{end+1} = '   mc_cold_frac = 0.5;'; % Factor 2 wrong
else
l{end+1} = '   mc_cold_frac = 0.5;';
disp('Source needs to be cold, thermal or bispectral! Using defaults.')
end
l{end+1} = 'else if (mod_target_index == 1)';
l{end+1} = '   mc_cold_frac = 1;'; % needs to be 0.75
l{end+1} = 'else if (mod_target_index == 2 || mod_target_index == 3)';
l{end+1} = '   mc_cold_frac = 0;'; % needs to be 0.35
l{end+1} = 'else if (mod_target_index == 4 || mod_target_index == 5)';
l{end+1} = '   mc_cold_frac = 0.5;'; % factor 2 error
l{end+1} = 'else printf("\\n - Invalid mod_target_index! use int between 0 and 5.\\n");';
elseif strcmp(source_component,'ESS_pancake') 
l{1}='';
l{end+1} = '#if defined (USE_MPI)';
l{end+1} = 'MPI_MASTER(';
l{end+1} = '#endif';
l{end+1} = 'if (mod_y>0.1) printf("G_ERROR moderator height outside of boundary.\\n");';
if strcmp(requirements.source,'cold_pancake')
l{end+1} = 'sprintf(command_string,"./write_cold %%f",mod_y*100);';
elseif strcmp(requirements.source,'thermal_pancake')
l{end+1} = 'sprintf(command_string,"./write_thermal %%f",mod_y*100);';
elseif strcmp(requirements.source,'bispectral_pancake')
l{end+1} = 'sprintf(command_string,"./write_bispectral %%f",mod_y*100);';    
end
if max(ismember(fieldnames(options_general),'kaspar_hack'))
l{end+1} = '// system(command_string);';
l{end+1} = '// sprintf(command_string,"mv pancake_spectrum.dat pancake_spectrum_%%i_%%i.dat",scan1,scan2);';
l{end+1} = '// system(command_string);';       
else
l{end+1} = 'system(command_string);';
l{end+1} = 'sprintf(command_string,"mv pancake_spectrum.dat pancake_spectrum_%%i_%%i.dat",scan1,scan2);';
l{end+1} = 'system(command_string);';
end

%if max(ismember(fieldnames(options_general),'kaspar_hack')) Error in copy
%paste. Checked with micro version 2.

l{end+1} = '#if defined (USE_MPI) || defined(USE_THREADS)';
l{end+1} = ');';
l{end+1} = '#endif';
l{end+1} = 'sleep(1);';
if max(ismember(fieldnames(options_general),'kaspar_hack'))
l{end+1} = 'sprintf(name_string,"custom_spectrum.dat");';    
else
l{end+1} = 'sprintf(name_string,"pancake_spectrum_%%i_%%i.dat",scan1,scan2);';
end
l{end+1} = '';

elseif strcmp(source_component,'Virtual_in')
l{1}='';

elseif strcmp(source_component,'locked_source_gen')
l{1}='';
l{end+1} = ['sprintf(name_string,"' locked_source_gen_filename '");'];
% add necessary initialize code for Virtual_in source here

end

initstring_ess=[initstring '\n'];
for i=1:length(l)
    initstring_ess=[initstring_ess l{i} '\n'];
end
clear l;
McStasStr.initialize_ess=[initstring_ess '\n\n' saved_inject_point];
%McStasStr.initialize_ess_seg = ['\nINITIALIZE\n%%{\n' McStasStr.initialize_ess '\n' McStasStr.initialize_seg '\n' McStasStr.pcalcstring '\n%%}\n'];
McStasStr.initialize_ess_seg = ['\nINITIALIZE\n%%{\n' McStasStr.initialize_ess '\n' McStasStr.initialize_seg];
%McStasStr.initialize_ess=['\nINITIALIZE\n%%{\n' McStasStr.initialize_ess '\n' McStasStr.pcalcstring '\n%%}\n'];
McStasStr.initialize_ess=['\nINITIALIZE\n%%{\n' McStasStr.initialize_ess];

% Overwrites the data needed in the above lines with a better name.
McStasStr.initialize_seg = McStasStr.initialize_seg_final;



% Modify the trace section to make it a complete instrument
        % add a source before the guide
        % add the detectors after the guide
if strcmp(source_component,'Virtual_in')
l{1}='COMPONENT Origin = Progress_bar()';
l{end+1}=' AT (0,0,0) ABSOLUTE';
l{end+1}='';
l{end+1}='COMPONENT source = Virtual_input(';
l{end+1}='   filename = "guide_bot_virtual_in.txt", type = "text", verbose = 1,';
l{end+1}=['   repeat_count = ' num2str(options_general.virtual_in_repeat_optimize) ', smooth = 1)'];
l{end+1}='   AT (0, 0, 0) RELATIVE Origin';
if isfield(options_general,'fom_weight')
l_opt{1}='  EXTEND';
l_opt{end+1}='    %%{';
l_opt{end+1}='    this_lambda = 2*PI/(sqrt(vx*vx+vy*vy+vz*vz)*V2K);';
l_opt{end+1}=['    if (this_lambda>' num2str(options_general.fom_weight) ') {'];
% can't use WaveMax in extend for some reason, one option would be to solve
% the problem in McStas, but since guide_bot already know the WaveMax it is
% just substituted directly.
%l{end+1}=['        p*=this_lambda/(' num2str(options_general.fom_weight) '-WaveMax)-WaveMax/(' num2str(options_general.fom_weight) '-WaveMax);'];
l_opt{end+1}=['        p*=this_lambda/(' num2str(options_general.fom_weight) '-' num2str(demands.WaveLmax) ')-' num2str(demands.WaveLmax) '/(' num2str(options_general.fom_weight) '-' num2str(demands.WaveLmax) ');'];
l_opt{end+1}='    }';
l_opt{end+1}='    %%}';
l_opt{end+1}='';
else
l_opt{1}='';    
end

l_opt{end+1}='COMPONENT StartOfGuide = Arm()';
l_opt{end+1}='AT (0,0,guide_start) RELATIVE source';

l_ana{1}='';
l_ana{end+1}='COMPONENT StartOfGuide = Arm()';
l_ana{end+1}='AT (0,0,guide_start) RELATIVE source';

else
l{1}='COMPONENT Origin = Progress_bar()';
l{end+1}=' AT (0,0,0) ABSOLUTE';
l{end+1}='';
l{end+1}='COMPONENT source = Source_simple(';
if ~strcmp(modules{modulelist(1)},'Selene') || globalinfo.selene_moderator_focus == 1
    % The first module is not a selene guide, standard procedure.
    l{end+1}='    yheight = mod_y, xwidth = mod_x,';
else
    % The first module is a Selene guide, special attention to focus and
    % moderator size is needed.
    % The calculations and declarations are in the selene module.
    % turned out only size was nessecary
    l{end+1}='    yheight = selene_mod_y, xwidth = selene_mod_x,';
end
if ~(strcmp(modules{modulelist(1)},'Selene') && globalinfo.selene_moderator_focus == 1) % Bug 22 jan 2014
    l{end+1}=['    dist = guide_start, focus_xw = startx' num2str(last) ', focus_yh = starty' num2str(last) ','];
else 
    % This might need modification in case of moderator focus in one
    % direction and not in the other.
    l{end+1}=['    dist = selene_distance1, focus_xw = selene_entry_x' num2str(last) '*3/4, focus_yh = selene_entry_y' num2str(last) '*3/4,'];
end
l{end+1}='    lambda0 = Lambda0, dlambda = dLambda, flux = 100, gauss = 0)';
l{end+1}='AT (0, 0, 0) RELATIVE Origin';
l{end+1}='ROTATED (0,0,0) RELATIVE Origin';
if isfield(options_general,'fom_weight')
l_opt{1}='  EXTEND';
l_opt{end+1}='    %%{';
l_opt{end+1}='    this_lambda = 2*PI/(sqrt(vx*vx+vy*vy+vz*vz)*V2K);';
l_opt{end+1}=['    if (this_lambda>' num2str(options_general.fom_weight) ') {'];
% can't use WaveMax in extend for some reason, one option would be to solve
% the problem in McStas, but since guide_bot already know the WaveMax it is
% just substituted directly.
%l{end+1}=['        p*=this_lambda/(' num2str(options_general.fom_weight) '-WaveMax)-WaveMax/(' num2str(options_general.fom_weight) '-WaveMax);'];
l_opt{end+1}=['        p*=this_lambda/(' num2str(options_general.fom_weight) '-' num2str(demands.WaveLmax) ')-' num2str(demands.WaveLmax) '/(' num2str(options_general.fom_weight) '-' num2str(demands.WaveLmax) ');'];
l_opt{end+1}='    }';
l_opt{end+1}='    %%}';
l_opt{end+1}='';
else
l_opt{1}='';   
end

l_opt{end+1}='COMPONENT StartOfGuide = Arm()';
l_opt{end+1}='AT (0,0,guide_start) RELATIVE source';

l_ana{1}='';
l_ana{end+1}='COMPONENT StartOfGuide = Arm()';
l_ana{end+1}='AT (0,0,guide_start) RELATIVE source';
end


clear s;

% Running the official ESS moderator
if strcmp(source_component,'ESS_moderator_long')
s{1}='COMPONENT Origin = Progress_bar()';
s{end+1}=' AT (0,0,0) ABSOLUTE';
s{end+1}='';
s{end+1}='COMPONENT ESS_source = ESS_moderator_long(';
s{end+1}='          yheight=mod_y,cyl_radius=ESS_cyl_radius,width_t=thermal_wing_length,beamport_angle=beamport_angle_zero_center+30,';
if ~(strcmp(modules{modulelist(1)},'Selene') && globalinfo.selene_moderator_focus == 1)
    s{end+1}=['          dist=guide_start, focus_xw=1.3*startx' num2str(last) ', focus_yh=1.05*starty' num2str(last) ','];
else
    s{end+1}=['          dist=selene_distance' num2str(last) ', focus_xw=1.3*selene_entry_x' num2str(last) '*3/4, focus_yh=1.05*selene_entry_y' num2str(last) '*3/4,'];
end
s{end+1}='          Lmin=WaveMin, Lmax=WaveMax, cold_frac=mc_cold_frac,';
s{end+1}='          T=50, tau=287e-6, tau1=0, tau2=20e-6, d=2.857e-3, n=20, width_c=0, nu=14,';
s{end+1}='          n2=5, chi2=0.9, I0=8.21e11, I2=3.29e11,branch1=1, branch2=0.5, branch_tail=0.14350,';
s{end+1}='          T_t=325, tau_t=80e-6, tau1_t=400e-6, tau2_t=12e-6, chi2_t=2.5, I0_t=13.5e11,';
s{end+1}='          I2_t=27.6e10, branch1_t=0.5,branch2_t=0.5)';
s{end+1}='  AT (0, 0, 0) RELATIVE Origin';
s{end+1}='';
s{end+1}='  COMPONENT StartOfGuide = Arm()';
s{end+1}='AT (0,0,guide_start+ESS_cyl_radius) RELATIVE Origin';
s{end+1}='ROTATED (0,ESS_angle,0) RELATIVE Origin';
elseif strcmp(source_component,'ESS_pancake') || strcmp(source_component,'locked_source_gen')
s{1}='COMPONENT Origin = Progress_bar()';
s{end+1}=' AT (0,0,0) ABSOLUTE';
s{end+1}='';
s{end+1}='COMPONENT source = Source_gen(';
if ~strcmp(modules{modulelist(1)},'Selene') || globalinfo.selene_moderator_focus == 1
    % The first module is not a selene guide, standard procedure.
    s{end+1}='    yheight = mod_y, xwidth = mod_x,';
else
    % The first module is a Selene guide, special attention to focus and
    % moderator size is needed.
    % The calculations and declarations are in the selene module.
    % turned out only size was nessecary
    s{end+1}='    yheight = selene_mod_y, xwidth = selene_mod_x,';
end
if ~(strcmp(modules{modulelist(1)},'Selene') && globalinfo.selene_moderator_focus == 1) % Recent change
    s{end+1}=['    dist = guide_start, focus_xw = startx' num2str(last) ', focus_yh = starty' num2str(last) ','];
else 
    % This might need modification in case of moderator focus in one
    % direction and not in the other.
    s{end+1}=['    dist = selene_distance1, focus_xw = selene_entry_x' num2str(last) '*3/4, focus_yh = selene_entry_y' num2str(last) '*3/4,'];
end
s{end+1}='    flux_file=name_string,flux_file_perAA=1,';
s{end+1}='    lambda0 = Lambda0, dlambda = dLambda)';
s{end+1}='AT (0, 0, 0) RELATIVE Origin';
s{end+1}='ROTATED (0,0,0) RELATIVE Origin';
s{end+1}='';
s{end+1}='COMPONENT StartOfGuide = Arm()';
s{end+1}='AT (0,0,guide_start) RELATIVE source';
elseif strcmp(source_component,'Virtual_in')
s{1}='COMPONENT Origin = Progress_bar()';
s{end+1}=' AT (0,0,0) ABSOLUTE';
s{end+1}='';    
s{end+1} = 'COMPONENT source = Virtual_input(';
s{end+1} = '   filename = "guide_bot_virtual_in.txt", type = "text", verbose = 1,';
s{end+1}=['   repeat_count = ' num2str(options_general.virtual_in_repeat_analyze) ', smooth = 1)'];
s{end+1} = '   AT (0, 0, 0) RELATIVE Origin';
s{end+1}='';
s{end+1}='COMPONENT StartOfGuide = Arm()';
s{end+1}='AT (0,0,guide_start) RELATIVE source';
end
    
source_for_bt_refference=l; % The source used in refference 

tracestring='';
for i=1:length(s)
    tracestring=[tracestring s{i} '\n'];
end
McStasStr.trace_ess_source=[tracestring '\n\n' McStasStr.trace '\n\n'];

tracestring='';
for i=1:length(s)
    tracestring=[tracestring s{i} '\n'];
end
McStasStr.trace_ess_seg=[tracestring '\n\n' McStasStr.trace_seg '\n\n'];

%tracestring='';
%for i=1:length(l)
%    tracestring=[tracestring l{i} '\n'];
%end
%McStasStr.trace=[tracestring '\n\n' McStasStr.trace '\n\n'];
%McStasStr.trace_seg=[tracestring '\n\n' McStasStr.trace_seg '\n\n']; 

tracestring='';
for i=1:length(l)
    tracestring=[tracestring l{i} '\n'];
end
for i=1:length(l_opt)
    tracestring=[tracestring l_opt{i} '\n'];
end
McStasStr.trace_opt=[tracestring '\n\n' McStasStr.trace '\n\n'];
McStasStr.trace_seg_opt=[tracestring '\n\n' McStasStr.trace_seg '\n\n'];

tracestring='';
for i=1:length(l)
    tracestring=[tracestring l{i} '\n'];
end
for i=1:length(l_ana)
    tracestring=[tracestring l_ana{i} '\n'];
end
McStasStr.trace_ana=[tracestring '\n\n' McStasStr.trace '\n\n'];
McStasStr.trace_seg_ana=[tracestring '\n\n' McStasStr.trace_seg '\n\n'];

clear l;





% These could be moved to a seperate text file and read
l{1}='';
if options_general.weighted_optimization == 1
    if ~isfield(demands,'weighted_Hdiv')
       disp('ERROR, need to specify demand.weighted_Hdiv for weighted optimization, or set options.weighted_optimization = 0');
    end 
    if ~isfield(demands,'weighted_Vdiv')
       disp('ERROR, need to specify demand.weighted_Vdiv for weighted optimization, or set options.weighted_optimization = 0');
    end 
    if ~isfield(demands,'weighted_Hsize')
       disp('ERROR, need to specify demand.weighted_Hsize for weighted optimization, or set options.weighted_optimization = 0');
    end
    if ~isfield(demands,'weighted_Vsize')
       disp('ERROR, need to specify demand.weighted_Vsize for weighted optimization, or set options.weighted_optimization = 0');
    end
    declare_index = length(McStasStr.declare);
    McStasStr.declare{declare_index+1} = 'optimize_weight_factor';
    McStasStr.initialize = [McStasStr.initialize '\n optimize_weight_factor = (sizeX*sizeY*divreq_x*divreq_y - ' num2str(demands.weighted_Hdiv) '*' num2str(demands.weighted_Vdiv) '*' num2str(demands.weighted_Hsize*0.01) '*' num2str(demands.weighted_Vsize*0.01) ')/(' num2str(demands.weighted_Hdiv) '*' num2str(demands.weighted_Vdiv) '*' num2str(demands.weighted_Hsize*0.01) '*' num2str(demands.weighted_Vsize*0.01) ');'];
    
    l{end+1}='COMPONENT Weighted_optimization = Arm()';
    l{end+1}='  AT (0, 0, sample_dist) RELATIVE PREVIOUS';    
    l{end+1}='  EXTEND';
    l{end+1}='    %%{';
    l{end+1}='    x_div = RAD2DEG*atan2(vx,vz);';
    l{end+1}='    y_div = RAD2DEG*atan2(vy,vz);';
    l{end+1}=['    if (x_div > -' num2str(demands.weighted_Hdiv) ' && x_div < ' num2str(demands.weighted_Hdiv) ' && y_div > -' num2str(demands.weighted_Vdiv) ' && y_div < ' num2str(demands.weighted_Vdiv) ') {']; 
    l{end+1}=['      if (x > -' num2str(demands.weighted_Hsize*0.01) ' && x < ' num2str(demands.weighted_Hsize*0.01) ' && y > -' num2str(demands.weighted_Vsize*0.01) ' && y < ' num2str(demands.weighted_Vsize*0.01) ') ']; 
    l{end+1}='         p *= optimize_weight_factor;';
    l{end+1}='    }';
    l{end+1}='    %%}';
    l{end+1}='';
    l{end+1}='COMPONENT Div2d_sample_B = Divergence_monitor(';
    l{end+1}='    nh = 20, nv = 20, filename = "Div2d_sample_B", xwidth = sizeX,';
    l{end+1}='    yheight = sizeY, maxdiv_h = divreq_x, maxdiv_v = divreq_y)';
    l{end+1}='  AT (0, 0, u) RELATIVE PREVIOUS';
else
l{end+1}='COMPONENT Div2d_sample_B = Divergence_monitor(';
l{end+1}='    nh = 20, nv = 20, filename = "Div2d_sample_B", xwidth = sizeX,';
l{end+1}='    yheight = sizeY, maxdiv_h = divreq_x, maxdiv_v = divreq_y)';
l{end+1}='  AT (0, 0,sample_dist) RELATIVE PREVIOUS';    
end
l{end+1}='';
l{end+1}='FINALLY';
l{end+1}='%%{';
%l{end+1}='#if defined (USE_MPI)';
%l{end+1}='MPI_MASTER(';
%l{end+1}='#endif';
%l{end+1}='';
%l{end+1}='printf("The var_divreq_x variable is=%%lf \\n",var_divreq_x);';
%l{end+1}='printf("The var_divreq_y variable is=%%lf \\n",var_divreq_y);';
%l{end+1}='printf("The guide_start variable is=%%lf \\n",guide_start);';
%l{end+1}=['printf("The startx' num2str(last) ' variable is=%%lf \\n",startx' num2str(last) ');'];
%l{end+1}=['printf("The endx' num2str(last) ' variable is=%%lf \\n",endx' num2str(last) ');'];
%l{end+1}=['printf("The starty' num2str(last) ' variable is=%%lf \\n",starty' num2str(last) ');'];
%l{end+1}=['printf("The endy' num2str(last) ' variable is=%%lf \\n",endy' num2str(last) ');'];
%l{end+1}='printf("tan(var_divreq_x*DEG2RAD) = %%lf \\n",tan(var_divreq_x*DEG2RAD));';
%l{end+1}='printf("tan(var_divreq_y*DEG2RAD) = %%lf \\n",tan(var_divreq_y*DEG2RAD));';
% l{end+1}=['for (i=0;i<' num2str(last+1) ';i++) {'];
% l{end+1}=['printf("startxpoint[%%d][1][1]=%%lf \\n",i,startxpoint[i][1][1]);'];
% l{end+1}=['printf("startxpoint[%%d][2][1]=%%lf \\n",i,startxpoint[i][2][1]);'];
% l{end+1}=['printf("startxpoint[%%d][1][2]=%%lf \\n",i,startxpoint[i][1][2]);'];
% l{end+1}=['printf("startxpoint[%%d][2][2]=%%lf \\n",i,startxpoint[i][2][2]);'];
% l{end+1}='}';
%l{end+1}='#if defined (USE_MPI) || defined(USE_THREADS)';
%l{end+1}='   );';
%l{end+1}='#endif';
l{end+1}='%%}';

% The following commands controll how much larger than the sample the
% different monitors are.
% If these are increased, it does have consequences for the brilliance scan
% that are not automatically takken into account. At present the moderator
% is at least 4 times larger than the sample, and have a divergence 6 times
% larger than needed.

if isfield(options_general,'psd2d_multiplier')
   psd2d = options_general.psd2d_multiplier;
   disp(['Showing performance on 2d psd monitors to ' num2str(psd2d) ' times the figure of merit'])
else
   psd2d=3;
end

if isfield(options_general,'psd1d_multiplier')
   psd1d = options_general.psd1d_multiplier;
   disp(['Showing performance on 1d psd monitors to ' num2str(psd1d) ' times the figure of merit'])
else
   psd1d=4;
end

if isfield(options_general,'div2d_multiplier')
   div2d = options_general.div2d_multiplier;
   disp(['Showing performance on 2d divergence monitors to ' num2str(div2d) ' times the figure of merit'])
else
   div2d=3;
end

if isfield(options_general,'div1d_multiplier')
   div1d = options_general.div1d_multiplier;
   disp(['Showing performance on 1d divergence monitors to ' num2str(div1d) ' times the figure of merit'])
else
   div1d=3;
end

if isfield(options_general,'apsd_multiplier')
   apsd = options_general.apsd_multiplier;
   disp(['Showing performance on the spatial part of acceptance monitors to ' num2str(apsd) ' times the figure of merit'])
else
   apsd=2;
end

if isfield(options_general,'adiv_multiplier')
   adiv = options_general.adiv_multiplier;
   disp(['Showing performance on the divergence part of acceptance monitors to ' num2str(adiv) ' times the figure of merit'])
else
   adiv=2;
end

if isfield(options_general,'ncount_multiplier')
   options_general.ncount_multiplier_optimize = options_general.ncount_multiplier;
   options_general.ncount_multiplier_analyze = options_general.ncount_multiplier;
   disp(['multiplying normal ncount with a factor of' num2str(options_general.ncount_multiplier) 'for both optimization and analysis'])
else
   if isfield(options_general,'ncount_multiplier_analyze')
    %options_general.ncount_multiplier_analyze = options_general.ncount_multiplier_analyze;
    disp(['multiplying normal ncount with a factor of ' num2str(options_general.ncount_multiplier_analyze) ' for both analysis'])
   else
       options_general.ncount_multiplier_analyze=1;
   end
   if isfield(options_general,'ncount_multiplier_optimize')
    %options_general.ncount_multiplier_optimize = options_general.ncount_multiplier_optimize;
    disp(['multiplying normal ncount with a factor of ' num2str(options_general.ncount_multiplier_optimize) ' for optimization'])
   else 
       options_general.ncount_multiplier_optimize=1;
   end
end


p2d=num2str(psd2d);
p1d=num2str(psd1d);
d2d=num2str(div2d);
d1d=num2str(div1d);
ap=num2str(apsd);
ad=num2str(adiv);

clear t;
t{1}=    'COMPONENT Lmon_guide_end = L_monitor(';
t{end+1}='    nL = 300, filename = "Lmon_guide_end", xwidth = endx1*1.1, restore_neutron=1,';
t{end+1}='    yheight = endy1*1.1, Lmin = WaveMin, Lmax = WaveMax)';
t{end+1}='  AT (0, 0,u) RELATIVE PREVIOUS';
t{end+1}='';
t{end+1}='COMPONENT Div2d_sample_B = Divergence_monitor(';
t{end+1}='    nh = 200, nv = 200, filename = "Div2d_sample_B", xwidth = sizeX,';
t{end+1}='    yheight = sizeY, maxdiv_h = divreq_x, maxdiv_v = divreq_y,restore_neutron=1)';
t{end+1}='  AT (0, 0,sample_dist) RELATIVE PREVIOUS';    
t{end+1}='  EXTEND';
t{end+1}='    %%{';
t{end+1}='    x_div = RAD2DEG*atan2(vx,vz);';
t{end+1}='    y_div = RAD2DEG*atan2(vy,vz);';
t{end+1}='    if (SCATTERED) flag=1; else flag=0;';
t{end+1}='    %%}';
t{end+1}='';
t{end+1}='COMPONENT Div2d_sample = Divergence_monitor(';
t{end+1}='    nh = 200, nv = 200, filename = "Div2d_sample", xwidth = sizeX, restore_neutron=1,';
t{end+1}=['    yheight = sizeY, maxdiv_h = divreq_x*' d2d ', maxdiv_v = divreq_y*' d2d ')'];
t{end+1}='  AT (0, 0,u) RELATIVE PREVIOUS';    
t{end+1}='';
t{end+1}='COMPONENT PSD_sample = PSD_monitor(';
t{end+1}='    nx = 200, ny = 200, filename = "PSD_sample", restore_neutron=1,';
t{end+1}=['    xwidth = sizeX*' p2d ', yheight = sizeY*' p2d ')'];
t{end+1}='  AT (0, 0,u) RELATIVE PREVIOUS';
t{end+1}='';
t{end+1}='COMPONENT HPSD_sample = PSDlin_monitor(';
t{end+1}=['    nx = 100, filename="HPSD_sample",xwidth = sizeX*' p1d ', yheight = sizeY, restore_neutron=1)'];
t{end+1}='  AT (0, 0,u) RELATIVE PREVIOUS';
t{end+1}='';
t{end+1}='COMPONENT VPSD_sample = PSDlin_y_monitor(';
t{end+1}=['    ny = 100, filename="VPSD_sample",xwidth = sizeX, yheight = sizeY*' p1d ', restore_neutron=1)'];
t{end+1}='  AT (0, 0,u) RELATIVE PREVIOUS';
t{end+1}='';
t{end+1}='COMPONENT Hdiv_sample = Hdiv_monitor(';
t{end+1}='    nh = 200, filename = "Hdiv_sample", xwidth = sizeX,';
t{end+1}=['    yheight = sizeY, h_maxdiv = divreq_x*' d1d ', restore_neutron=1)'];
t{end+1}='  AT (0, 0,u) RELATIVE PREVIOUS';
t{end+1}='';
t{end+1}='COMPONENT Vdiv_sample = Vdiv_monitor(';
t{end+1}='    nv = 200, filename = "Vdiv_sample", xwidth = sizeX, restore_neutron=1,';
t{end+1}=['    yheight = sizeY, v_maxdiv = divreq_y*' d1d ')'];
t{end+1}='  AT (0, 0,u) RELATIVE PREVIOUS';
t{end+1}='';
t{end+1}='COMPONENT acceptance_x_divx = DivPos_monitor(';
t{end+1}='   nh = 200, ndiv = 200, filename = "acceptance_x_divx",';
t{end+1}=['   xwidth = sizeX*' ap ', yheight = sizeY, maxdiv_h = divreq_x*' ad ', restore_neutron=1)'];
t{end+1}=' AT (0, 0, u) RELATIVE PREVIOUS';
t{end+1}='';
t{end+1}='COMPONENT acceptance_y_divy = DivPos_y_monitor(';
t{end+1}='   npos = 200, ndiv = 200, filename = "acceptance_y_divy",';
t{end+1}=['   xwidth = sizeX, yheight = sizeY*' ap ', maxdiv = divreq_y*' ad ', restore_neutron=1)'];
t{end+1}=' AT (0, 0, u) RELATIVE PREVIOUS';
t{end+1}='';
% From this point, everything is within maxdiv and can be expressed in BT
t{end+1}='COMPONENT Lmon_sample_B = L_monitor(';
t{end+1}='    nL = 300, filename = "Lmon_sample_B", xwidth = sizeX, restore_neutron=1,';
t{end+1}='    yheight = sizeY, Lmin = WaveMin, Lmax = WaveMax) WHEN flag';
t{end+1}='  AT (0, 0,u) RELATIVE PREVIOUS';
t{end+1}='';
t{end+1}='COMPONENT Div2d_sample_maxdiv = Divergence_monitor(';
t{end+1}='    nh = 200, nv = 200, filename = "Div2d_sample_maxdiv", xwidth = sizeX, restore_neutron=1,';
t{end+1}=['    yheight = sizeY, maxdiv_h = divreq_x*' d2d ', maxdiv_v = divreq_y*' d2d ')'];
t{end+1}='  AT (0, 0,u) RELATIVE PREVIOUS';    
t{end+1}='';
t{end+1}='COMPONENT PSD_sample_maxdiv = PSD_monitor(';
t{end+1}='    nx = 200, ny = 200, filename = "PSD_sample_maxdiv", restore_neutron=1,';
t{end+1}=['    xwidth = sizeX*' p2d ', yheight = sizeY*' p2d ') WHEN (x_div<divreq_x && x_div>-divreq_x && y_div<divreq_y && y_div>-divreq_y)'];
t{end+1}='  AT (0, 0,u) RELATIVE PREVIOUS';
t{end+1}='';
t{end+1}='COMPONENT HPSD_sample_maxdiv = PSDlin_monitor(';
t{end+1}=['    nx = 100, filename="HPSD_sample_maxdiv",xwidth = sizeX*' p1d ', yheight = sizeY, restore_neutron=1)'];
t{end+1}='  WHEN (x_div<divreq_x && x_div>-divreq_x && y_div<divreq_y && y_div>-divreq_y)';
t{end+1}='  AT (0, 0,u) RELATIVE PREVIOUS';
t{end+1}='';
t{end+1}='COMPONENT VPSD_sample_maxdiv = PSDlin_y_monitor(';
t{end+1}=['    ny = 100, filename="VPSD_sample_maxdiv",xwidth = sizeX, yheight = sizeY*' p1d ', restore_neutron=1)'];
t{end+1}='  WHEN (x_div<divreq_x && x_div>-divreq_x && y_div<divreq_y && y_div>-divreq_y)';
t{end+1}='  AT (0, 0,u) RELATIVE PREVIOUS';
t{end+1}='';
t{end+1}='COMPONENT Hdiv_sample_maxdiv = Hdiv_monitor(';
t{end+1}='    nh = 200, filename = "Hdiv_sample_maxdiv", xwidth = sizeX, restore_neutron=1,';
t{end+1}=['    yheight = sizeY, h_maxdiv = divreq_x*' d1d ') WHEN (y_div<divreq_y && y_div>-divreq_y)'];
t{end+1}='  AT (0, 0,u) RELATIVE PREVIOUS';
t{end+1}='';
t{end+1}='COMPONENT Vdiv_sample_maxdiv = Vdiv_monitor(';
t{end+1}='    nv = 200, filename = "Vdiv_sample_maxdiv", xwidth = sizeX, restore_neutron=1,';
t{end+1}=['    yheight = sizeY, v_maxdiv = divreq_y*' d1d ') WHEN (x_div<divreq_x && x_div>-divreq_x)'];
t{end+1}='  AT (0, 0,u) RELATIVE PREVIOUS';
t{end+1}='';
t{end+1}='COMPONENT acceptance_x_divx_maxdiv = DivPos_monitor(';
t{end+1}='   nh = 200, ndiv = 200, filename = "acceptance_x_divx_maxdiv", restore_neutron=1,';
t{end+1}=['   xwidth = sizeX*' ap ', yheight = sizeY, maxdiv_h = divreq_x*' ad ') WHEN (y_div<divreq_y && y_div>-divreq_y)'];
t{end+1}=' AT (0, 0, u) RELATIVE PREVIOUS';
t{end+1}='';
t{end+1}='COMPONENT acceptance_y_divy_maxdiv = DivPos_y_monitor(';
t{end+1}='   npos = 200, ndiv = 200, filename = "acceptance_y_divy_maxdiv", restore_neutron=1,';
t{end+1}=['   xwidth = sizeX, yheight = sizeY*' ap ', maxdiv = divreq_y*' ad ') WHEN (x_div<divreq_x && x_div>-divreq_x)'];
t{end+1}=' AT (0, 0, u) RELATIVE PREVIOUS';
t{end+1}='';
t{end+1}='COMPONENT Lmon_sample = L_monitor(';
t{end+1}='    nL = 300, filename = "Lmon_sample", xwidth = sizeX, restore_neutron=1,';
t{end+1}='    yheight = sizeY, Lmin = WaveMin, Lmax = WaveMax)';
t{end+1}='  AT (0, 0,u) RELATIVE PREVIOUS';

if isfield(options_general,'monitor_file')
    fid = fopen(options_general.monitor_file);
    tline = fgetl(fid);
    t{end+1}='';
    while ischar(tline)
        t{end+1} = tline;
        tline = fgetl(fid);
    end
    t{end+1}='';
    fclose(fid);
end

% This function is called now as it does not need the latter part of t
brilliance_ref(McStasStr,source_for_bt_refference,t,input,filename,Project_name,requirements,demands,endinput,enddeclare,NumSnaps,scan,options_general)


t{end+1}='FINALLY';
t{end+1}='%%{';
%t{end+1}='';
%t{end+1}='#if defined (USE_MPI)';
%t{end+1}='MPI_MASTER(';
%t{end+1}='#endif';
% t{end+1}='';
% t{end+1}='printf("The var_divreq_x variable is=%%lf \\n",var_divreq_x);';
% t{end+1}='printf("The var_divreq_y variable is=%%lf \\n",var_divreq_y);';
% t{end+1}='printf("The guide_start variable is=%%lf \\n",guide_start);';
% t{end+1}=['printf("The startx' num2str(last) ' variable is=%%lf \\n",startx' num2str(last) ');'];
% t{end+1}=['printf("The endx' num2str(last) ' variable is=%%lf \\n",endx' num2str(last) ');'];
% t{end+1}=['printf("The starty' num2str(last) ' variable is=%%lf \\n",starty' num2str(last) ');'];
% t{end+1}=['printf("The endy' num2str(last) ' variable is=%%lf \\n\\n",endy' num2str(last) ');'];
% t{end+1}=['for (i=1;i<' num2str(last+1) ';i++) {'];
% t{end+1}=['printf("startxpoint[%%d][1][1]=%%lf \\n",i,startxpoint[i][1][1]);'];
% t{end+1}=['printf("startxpoint[%%d][2][1]=%%lf \\n",i,startxpoint[i][2][1]);'];
% t{end+1}=['printf("startxpoint[%%d][1][2]=%%lf \\n",i,startxpoint[i][1][2]);'];
% t{end+1}=['printf("startxpoint[%%d][2][2]=%%lf \\n",i,startxpoint[i][2][2]);'];
% t{end+1}='}';
%t{end+1}='';
%t{end+1}='#if defined (USE_MPI) || defined(USE_THREADS)';
%t{end+1}='   );';
%t{end+1}='#endif';
t{end+1}='%%}';

% Tilføj PSD monitors til overstående og måske endda forskellige versioner
% af monitore, f.eks med og uden maxdiv begrænsning.

tracestring_optimize='';
tracestring_analyze='';
for i=1:length(l)
    tracestring_optimize=[tracestring_optimize l{i} '\n'];
end
for i=1:length(t)
    tracestring_analyze=[tracestring_analyze t{i} '\n'];
end
McStasStr.trace_optimize='';McStasStr.trace_analyze='';
McStasStr.trace_optimize=['\n%%}\nTRACE\n' McStasStr.trace_opt tracestring_optimize '\n\n'];
McStasStr.trace_analyze=['\n%%}\nTRACE\n' McStasStr.trace_ana tracestring_analyze '\n\n'];
McStasStr.trace_optimize_seg=['\n%%}\nTRACE\n' McStasStr.trace_seg_opt tracestring_optimize '\n\n'];
McStasStr.trace_analyze_seg=['\n%%}\nTRACE\n' McStasStr.trace_seg_ana tracestring_analyze '\n\n'];
McStasStr.trace_analyze_ess=['\n%%}\nTRACE\n' McStasStr.trace_ess_source tracestring_analyze '\n\n'];
McStasStr.trace_analyze_ess_seg=['\n%%}\nTRACE\n' McStasStr.trace_ess_seg tracestring_analyze '\n\n'];
clear l;clear t;


% start declare section
McStasStr.declarestring='\n\nDECLARE\n%%{\n';
for i=1:length(McStasStr.declare) 
    McStasStr.declarestring = [McStasStr.declarestring 'double ' McStasStr.declare{i} ';\n' ];
end
for i=1:length(McStasStr.declareint) 
    McStasStr.declarestring = [McStasStr.declarestring 'int ' McStasStr.declareint{i} ';\n' ];
end
McStasStr.declarestring = [McStasStr.declarestring 'FILE *fp;\n' ]; % declare pointer for writing files

if strcmp(source_component,'ESS_moderator_long')
% Just for the added ESS declarations
McStasStr.declarestring_ess='';
for i=1:length(McStasStr.declare_ess) 
    McStasStr.declarestring_ess = [McStasStr.declarestring_ess 'double ' McStasStr.declare_ess{i} ';\n' ];
end
for i=1:length(McStasStr.declareint_ess) 
    McStasStr.declarestring_ess = [McStasStr.declarestring_ess 'int ' McStasStr.declareint_ess{i} ';\n' ];
end
elseif strcmp(source_component,'ESS_pancake') || strcmp(source_component,'locked_source_gen')
% This is a temporary way to do this.
McStasStr.declarestring_ess=['char command_string[100];\nchar name_string[40];\n'];
elseif strcmp(source_component,'Virtual_in')
McStasStr.declarestring_ess=['double virtual_dummy;\n'];    
end

% Just for the added seg declarations
McStasStr.declarestring_seg='';
if isfield(McStasStr,'declare_seg')
 for i=1:length(McStasStr.declare_seg) 
    McStasStr.declarestring_seg = [McStasStr.declarestring_seg 'double ' McStasStr.declare_seg{i} ';\n' ];
 end
end
if isfield(McStasStr,'declareint_seg')
 for i=1:length(McStasStr.declareint_seg) 
    McStasStr.declarestring_seg = [McStasStr.declarestring_seg 'int ' McStasStr.declareint_seg{i} ';\n' ];
 end
end
   
    
    

% Combine, thins needs to be done in the right order as they overwrite
% their own variables names (in order to keep the amount of variables down)
McStasStr.declarestring_ess_seg = [McStasStr.declarestring McStasStr.declarestring_ess McStasStr.declarestring_seg '%%}\n'];
McStasStr.declarestring_ess = [McStasStr.declarestring McStasStr.declarestring_ess '%%}\n'];
McStasStr.declarestring_seg = [McStasStr.declarestring McStasStr.declarestring_seg '%%}\n'];
McStasStr.declarestring = [McStasStr.declarestring '%%}\n'];
% end declare section

% the things done to avoid accesing inputvalue(i) where i is to large are
% stupid. Instead, check if inputvalue(imax) exists, if not, create it with
% a zero to indicate that it does not have an input value.
last_element = length(McStasStr.input);
if length(McStasStr.optimize) < last_element
    McStasStr.optimize(last_element)=0;
end

% input string handling for raw input.
McStasStr.inputstring=['DEFINE INSTRUMENT ' filename '(\n'];
for i=1:length(McStasStr.input)%-1 
    if i == length(McStasStr.input); delimiter = ''; else; delimiter = ',\n'; end;
    %delimiter = ',';
    if (i<=numel(McStasStr.inputvalue))
        if McStasStr.optimize(i) ~= 1%(McStasStr.inputvalue(i)~=0) %This change needs validation
          McStasStr.inputstring = [McStasStr.inputstring McStasStr.input{i} '=' num2str(McStasStr.inputvalue(i)) delimiter ];
        else
          McStasStr.inputstring = [McStasStr.inputstring McStasStr.input{i} delimiter ];
        end
    else
    McStasStr.inputstring = [McStasStr.inputstring McStasStr.input{i}  delimiter ];
    end
end

% seg input string addition
McStasStr.inputstring_seg=',\n';
if sum(ismember(fieldnames(McStasStr),'inputvalue_seg')) < 0.5
    McStasStr.inputvalue_seg = 0;
end

if isfield(McStasStr,'input_seg')
for i=1:length(McStasStr.input_seg)
    if i == length(McStasStr.input_seg); delimiter = ''; else; delimiter = ',\n'; end;
    if (i<=numel(McStasStr.inputvalue_seg))
        if McStasStr.optimize_seg(i) ~= 1%(McStasStr.inputvalue(i)~=0) %This change needs validation
          McStasStr.inputstring_seg = [McStasStr.inputstring_seg McStasStr.input_seg{i} '=' num2str(McStasStr.inputvalue_seg(i)) delimiter ];
        else
          McStasStr.inputstring_seg = [McStasStr.inputstring_seg McStasStr.input_seg{i} delimiter ];
        end
    else
    McStasStr.inputstring_seg = [McStasStr.inputstring_seg McStasStr.input_seg{i}  delimiter ];
    end
end
else
    McStasStr.inputstring_seg = '';
end


if strcmp(source_component,'ESS_moderator_long')
% ESS input string addition
McStasStr.inputstring_ess=',\n';
if sum(ismember(fieldnames(McStasStr),'inputvalue_ess')) < 0.5
    % not really nessecary as it will never happen, but better safe.
    McStasStr.inputvalue_ess = 0;
end
    
for i=1:length(McStasStr.input_ess)%-1 
    if i == length(McStasStr.input_ess); delimiter = ''; else; delimiter = ',\n'; end;
    if (i<=numel(McStasStr.inputvalue_ess))
        if McStasStr.optimize_ess(i) ~= 1%(McStasStr.inputvalue(i)~=0) %This change needs validation
          McStasStr.inputstring_ess = [McStasStr.inputstring_ess McStasStr.input_ess{i} '=' num2str(McStasStr.inputvalue_ess(i)) delimiter ];
        else
          McStasStr.inputstring_ess = [McStasStr.inputstring_ess McStasStr.input_ess{i} delimiter ];
        end
    else
    McStasStr.inputstring_ess = [McStasStr.inputstring_ess McStasStr.input_ess{i} delimiter ];
    end
end

elseif strcmp(source_component,'ESS_pancake')
    McStasStr.inputstring_ess = ',\nint scan1=1,\nint scan2=1';
elseif strcmp(source_component,'locked_source_gen')    
    McStasStr.inputstring_ess = ',\ndummy_param1=0';
elseif strcmp(source_component,'Virtual_in')    
    McStasStr.inputstring_ess = ',\nVirtual_dummy2=0';
end



% Write the mcstas file:
        % write the instrument name using the input string (remove spaces)
        % write the input section using the input list (standard values?)
        % write the declare section using the declare list
        % write the trace section using the trace string
        % write the finally string, which should write critical information
        %       about the parameters calculated to std output.
fid = fopen(['./' Project_name '/' filename '/' filename '_optimize.instr'], 'w');
write=[McStasStr.inputstring ')' McStasStr.declarestring McStasStr.initialize McStasStr.trace_optimize '\nEND'];
fprintf(fid,write);
fclose(fid);
        
fid = fopen(['./' Project_name '/' filename '/' filename '_analyze.instr'], 'w');
write=[McStasStr.inputstring ',\nstring file_name="test.dat"\n)' McStasStr.declarestring McStasStr.initialize McStasStr.visualizer_str McStasStr.trace_analyze '\nEND'];
fprintf(fid,write);
fclose(fid);

fid = fopen(['./' Project_name '/' filename '/' filename '_analyze_ess.instr'], 'w');
write=[McStasStr.inputstring McStasStr.inputstring_ess ')' McStasStr.declarestring_ess McStasStr.initialize_ess McStasStr.trace_analyze_ess '\nEND'];
fprintf(fid,write);
fclose(fid);

fid = fopen(['./' Project_name '/' filename '/' filename '_optimize_coating.instr'], 'w');
write=[McStasStr.inputstring McStasStr.inputstring_seg ')' McStasStr.declarestring_seg McStasStr.initialize_seg McStasStr.trace_optimize_seg '\nEND'];
fprintf(fid,write);
fclose(fid);

fid = fopen(['./' Project_name '/' filename '/' filename '_analyze_coating.instr'], 'w');
write=[McStasStr.inputstring McStasStr.inputstring_seg ')' McStasStr.declarestring_seg McStasStr.initialize_seg McStasStr.trace_analyze_seg '\nEND'];
fprintf(fid,write);
fclose(fid);

fid = fopen(['./' Project_name '/' filename '/' filename '_analyze_ess_coating.instr'], 'w');
write=[McStasStr.inputstring McStasStr.inputstring_ess McStasStr.inputstring_seg ')' McStasStr.declarestring_ess_seg McStasStr.initialize_ess_seg McStasStr.trace_analyze_ess_seg '\nEND'];
fprintf(fid,write);
fclose(fid);


% McStas instrument files writen!

% Now I have to write the iFit file!

% Hardcode for 0,1,2 dimension runs :(
% times to run ifit generator:
% scan.dimension = 0 => 1
% scan.dimension = 1 => 1
% scan.dimension = 2 => 2
% scan.dimension = 3 => 3 ...

%if scan.dimension == 0
%    genruns = 1;
%else
%   genruns = scan.dimension;
%end

if exist('McStasStr.optimvals.max(linp+1)') > 0
protected_guide_start_max = McStasStr.optimvals.max(linp+1);
protected_guide_start_min = McStasStr.optimvals.min(linp+1);
protected_guide_start_guess = McStasStr.optimvals.guess(linp+1);
end


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
    
    % what if a 2d scan is one demand and one requirement? <- Handled.
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
                if strcmp(scan.locked.names_parent{num_locked},scan.names{1}) % NEED PARENT HERE
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
    
    
    
       
    % If the minimalist principle is used, the range og guide_start could
    % be affected. This is calculated and used to avoid wasting the
    % optimizers time and possibly confusing ifit.
    if McStasStr.optimize(linp+1)==1 && McStasStr.minimalist==1 
        if strcmp(globalinfo.PS_req,'homogen')
        % Maximum reasonable guide_start value calculation:
        Deg2rad = 3.14159/180;
        if exist('McStasStr.optimvals.max(linp+1)') > 0
        McStasStr.optimvals.max(linp+1) = protected_guide_start_max;
        McStasStr.optimvals.min(linp+1) = protected_guide_start_min;
        McStasStr.optimvals.guess(linp+1) = protected_guide_start_guess;
        end
        
        
        if globalinfo.focusing == 1
            PS_requirement_x = sqrt(reqscan.minimalist_factor)*2*demandscan.Hsize*demandscan.Hdiv*Deg2rad/100;
            PS_requirement_y = sqrt(reqscan.minimalist_factor)*2*demandscan.Vsize*demandscan.Vdiv*Deg2rad/100;
        else
            PS_requirement_x = sqrt(reqscan.minimalist_factor)*2*(demandscan.Hsize/100+2*demandscan.Dist*tand(demandscan.Hdiv))*demandscan.Hdiv*Deg2rad;
            PS_requirement_y = sqrt(reqscan.minimalist_factor)*2*(demandscan.Vsize/100+2*demandscan.Dist*tand(demandscan.Vdiv))*demandscan.Vdiv*Deg2rad;
        end
        
        %L_max_x=reqscan.moderator_size_x^2/(8*tand(demandscan.Hdiv)*(demandscan.Hsize/100+2*demandscan.Dist*tand(demandscan.Hdiv)));
        %L_max_y=reqscan.moderator_size_y^2/(8*tand(demandscan.Vdiv)*(demandscan.Vsize/100+2*demandscan.Dist*tand(demandscan.Vdiv)));

        L_max_x=reqscan.moderator_size_x^2/(4*PS_requirement_x);
        L_max_y=reqscan.moderator_size_y^2/(4*PS_requirement_y);
        
        % any value larger than the minimum of the two will give some
        % complex lengths.
        L_max = min([L_max_x L_max_y]);
              
        %requirements.guide_min_opening_x=requirements.moderator_size_x/6;
        %requirements.guide_min_opening_y=requirements.moderator_size_y/6;
        %reqscan.guide_min_opening_x=reqscan.moderator_size_x/6;
        %reqscan.guide_min_opening_y=reqscan.moderator_size_y/6;
        
        requirements.guide_min_opening_x=0.015;
        requirements.guide_min_opening_y=0.015;
        reqscan.guide_min_opening_x=0.015;
        reqscan.guide_min_opening_y=0.015;        
        
        % I have no idear why this needs to be here, can be removed.
        %if reqscan.guide_min_opening_x>0.059
        %    reqscan.guide_min_opening_x=0.059;
        %end
        %if reqscan.guide_min_opening_y>0.059
        %    reqscan.guide_min_opening_y=0.059;
        %end
        
        % simple factor error for a long time
        %L_min_x=(reqscan.moderator_size_x^2-(reqscan.moderator_size_x/2-reqscan.guide_min_opening_x)^2)/(4*PS_requirement_x)
        %L_min_y=(reqscan.moderator_size_y^2-(reqscan.moderator_size_y/2-reqscan.guide_min_opening_y)^2)/(4*PS_requirement_y)
        
        % correct expressions
        %L_min_x=(reqscan.moderator_size_x^2/4-(reqscan.moderator_size_x/2-reqscan.guide_min_opening_x)^2)/(PS_requirement_x)
        %L_min_y=(reqscan.moderator_size_y^2/4-(reqscan.moderator_size_y/2-reqscan.guide_min_opening_y)^2)/(PS_requirement_y)
        
        % correct and simple expressions
        L_min_x=reqscan.guide_min_opening_x*(reqscan.moderator_size_x-reqscan.guide_min_opening_x)/PS_requirement_x;
        L_min_y=reqscan.guide_min_opening_y*(reqscan.moderator_size_y-reqscan.guide_min_opening_y)/PS_requirement_y;
        
        
        L_min = max([L_min_y L_min_x]);

        % TEMPORARY CODE, THIS CAN BE IMPROVED VASTLY 
        % ( But should not be nessecary, just addedd for safety )
        if L_min > L_max
            L_min = L_max - 0.01;
        end

        if McStasStr.optimvals.max(linp+1)>L_max
            if L_max>reqscan.closest_element
                McStasStr.optimvals.max(linp+1)=L_max;
            else
                McStasStr.optimvals.max(linp+1)=reqscan.closest_element+0.01;
                disp('Should use options.phase_space_requirement=''total'', as there is not enough homogeneous!')
            end
            
        end
        % The last part is needed in some cases which is disturbing, should
        % really never happen that L_min > requirements.latest_start
        if L_min > reqscan.closest_element && L_min < reqscan.latest_start
            McStasStr.optimvals.min(linp+1)=L_min;
        end
        McStasStr.optimvals.guess(linp+1)=mean([McStasStr.optimvals.max(linp+1) McStasStr.optimvals.min(linp+1)]);
        
        elseif strcmp(globalinfo.PS_req,'total')
        % only a new L_min is relevant, as there is no max length where the
        % equation can not be solved.
        Deg2rad = 3.14159/180;
        if exist('McStasStr.optimvals.max(linp+1)') > 0
        McStasStr.optimvals.max(linp+1) = protected_guide_start_max;
        McStasStr.optimvals.min(linp+1) = protected_guide_start_min;
        McStasStr.optimvals.guess(linp+1) = protected_guide_start_guess;
        end
        
        requirements.guide_min_opening_x=0.015;
        requirements.guide_min_opening_y=0.015;
        reqscan.guide_min_opening_x=0.015;
        reqscan.guide_min_opening_y=0.015;
        
        
        if globalinfo.focusing == 1
            PS_requirement_x = sqrt(reqscan.minimalist_factor)*2*demandscan.Hsize*demandscan.Hdiv*Deg2rad/100;
            PS_requirement_y = sqrt(reqscan.minimalist_factor)*2*demandscan.Vsize*demandscan.Vdiv*Deg2rad/100;
        else
            PS_requirement_x = sqrt(reqscan.minimalist_factor)*2*(demandscan.Hsize/100+2*demandscan.Dist*tand(demandscan.Hdiv))*demandscan.Hdiv*Deg2rad;
            PS_requirement_y = sqrt(reqscan.minimalist_factor)*2*(demandscan.Vsize/100+2*demandscan.Dist*tand(demandscan.Vdiv))*demandscan.Vdiv*Deg2rad;
        end
        
        % W = PS*guide_start / mod = > guide_start = W*mod/PS
        
        
        %disp('min lengths for')
        %disp(['demandscan.Hsize = ' num2str(demandscan.Hsize) ', demandscan.Hdiv = ' num2str(demandscan.Hdiv) ])
        L_min_x = reqscan.guide_min_opening_x*reqscan.moderator_size_x/PS_requirement_x;
        L_min_y = reqscan.guide_min_opening_y*reqscan.moderator_size_y/PS_requirement_y;
       
        L_min = max([L_min_y L_min_x]);
        
        if L_min > reqscan.closest_element 
            if L_min < reqscan.latest_start
            McStasStr.optimvals.min(linp+1)=L_min;
            else
               disp('Consider increasing latest start, as the guide opening will be very small in one dimension') 
               disp(['demandscan.Hsize = ' num2str(demandscan.Hsize) ', demandscan.Hdiv = ' num2str(demandscan.Hdiv) ])
               disp(['demandscan.Vsize = ' num2str(demandscan.Vsize) ', demandscan.Vdiv = ' num2str(demandscan.Vdiv) ])
               % it would be better to fix guide_start to the maximum
               % value, but difficult under the scan environment
               if McStasStr.optimvals.max(linp+1)-0.01 > McStasStr.optimvals.min(linp+1)
                McStasStr.optimvals.min(linp+1)=McStasStr.optimvals.max(linp+1)-0.01;
                disp('Only allowing long distances for guide_start')
               else
                disp('ERROR, the narrow band for guide_start means the simulation will probably fail unless the first module is S.')
               end
            end
        end
        
        McStasStr.optimvals.guess(linp+1)=mean([McStasStr.optimvals.max(linp+1) McStasStr.optimvals.min(linp+1)]);
        end
    end

clear l;
l{1}=['Project_name=''' Project_name ''';'];
l{end+1}=['instrument_name=''' filename ''';'];
l{end+1}=['filename=''' filename ''';'];
l{end+1}=['inputstring=''' input ''';'];
if scan.mode==0
  l{end+1}=['scanname='''';'];
else
  l{end+1}=['scanname=''' scanname ''';'];
end
l{end+1}='';
l{end+1}=['scani=' num2str(scani) ';'];
l{end+1}=['scanj=' num2str(scanj) ';'];
l{end+1}='cpath=pwd;';
l{end+1}='i=0;done=''1''';
l{end+1}='%% creates a new filename';
l{end+1}='while(str2num(done)==1); i=1+i; [a,done]=unix([''test -d \'' cpath ''/'' filename scanname num2str(i) '' && echo "1" || echo "0"'']); end;';
l{end+1}='filename=[filename scanname num2str(i)];';
l{end+1}='';
l{end+1}='if ispc==1';
l{end+1}=' select=1;';
l{end+1}='else';
l{end+1}=' %% checks if the script is on the ESS cluster';
l{end+1}=[ '[a,check]=unix(''test -f ' options_general.cluster_path 'clusterid && echo "2" || echo "1"'')'];
l{end+1}=' select=str2num(check);';
l{end+1}='end';
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
for i=1:length(McStasStr.input)
    if (optimize(i)==1)
        stringmin=num2str(McStasStr.optimvals.min(i));
        stringmax=num2str(McStasStr.optimvals.max(i));
        stringguess=num2str(McStasStr.optimvals.guess(i));
        iFitStr=[iFitStr 'p.' McStasStr.input{i} '= [' stringmin ' ' stringguess ' ' stringmax ']; \n']; %? Do i need a value for everything?
    else
        if scan.mode==0
        iFitStr=[iFitStr 'p.' McStasStr.input{i} '=''' num2str(inputvalue(i)) ''';\n']; %? Do i need a value for everything?
        else
            
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
                elseif (strcmp(variable_name,'moderator_size_x') || strcmp(variable_name,'moderator_size_y') || strcmp(variable_name,'minimalist_factor'))
                  iFitStr=[iFitStr 'p.' McStasStr.input{i} '=''' num2str(reqscan.(variable_name)) ''';\n'];  
                else
                  iFitStr=[iFitStr 'p.' McStasStr.input{i} '=''' num2str(demandscan.(variable_name)) ''';\n'];
                end
            else
               % write normally
               iFitStr=[iFitStr 'p.' McStasStr.input{i} '=''' num2str(inputvalue(i)) ''';\n'];                
            end
            
            % This code is a result of the suboptimal way the scanned
            % variables are loaded into McStasStr.input.
            % Legacy code, the above is a rewrite but not a fix.
%             if strcmp(scan.namesmcstas{1},McStasStr.input{i})
%                % Add to the ifit string in the three cases
%                 identifier = 1;
%                 if (strcmp(scan.names{identifier},'Hsize') || strcmp(scan.names{identifier},'Vsize'))
%                   iFitStr=[iFitStr 'p.' McStasStr.input{i} '=''' num2str(demands.(scan.names{identifier})(scani)/100) ''';\n'];    
%                 elseif (strcmp(scan.names{identifier},'moderator_size_x') || strcmp(scan.names{identifier},'moderator_size_y'))
%                   iFitStr=[iFitStr 'p.' McStasStr.input{i} '=''' num2str(requirements.(scan.names{identifier})(scani)) ''';\n'];  
%                 else
%                   iFitStr=[iFitStr 'p.' McStasStr.input{i} '=''' num2str(demands.(scan.names{identifier})(scani)) ''';\n'];
%                 end
%             elseif scan.dimension == 2 && strcmp(scan.namesmcstas{2},McStasStr.input{i})
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
%             else % Need to add another endif to check if the name is on the lock list. If it is on the lock list, add it using scani/scanj.
%                % write normally
%                iFitStr=[iFitStr 'p.' McStasStr.input{i} '=''' num2str(inputvalue(i)) ''';\n'];
%             end
            
            %if match==0
            %    iFitStr=[iFitStr 'p.' McStasStr.input{i} '=''' num2str(inputvalue(i)) ''';\n'];
            %else
            %    if (strcmp(scan.names{stest},'Hsize') || strcmp(scan.names{stest},'Vsize'))
            %      iFitStr=[iFitStr 'p.' McStasStr.input{i} '=''' num2str(demands.(scan.names{stest})(scani)/100) ''';\n'];    
            %    elseif (strcmp(scan.names{stest},'moderator_size_x') || strcmp(scan.names{stest},'moderator_size_y'))
            %      iFitStr=[iFitStr 'p.' McStasStr.input{i} '=''' num2str(requirements.(scan.names{stest})(scani)) ''';\n'];  
            %    else
            %      iFitStr=[iFitStr 'p.' McStasStr.input{i} '=''' num2str(demands.(scan.names{stest})(scani)) ''';\n'];
            %    end
            %end
        end
    end
end


List=0:(NumSnaps-1);
wavec=zeros(NumSnaps,1);
deltaL=demandscan.WaveLmax-demandscan.WaveLmin;
centerscalc=demandscan.WaveLmin+deltaL/(NumSnaps-1)*List;
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

% check if it is nessecary to run an optimizer
if sum(optimize) == 0
    need_optimizer = 0;
else
    need_optimizer = 1;
end


l{1}='';
l{end+1}='%%options_home.compile=0;';
if options_general.gravity == 1
l{end+1}='options_home.gravitation=1;';
else
l{end+1}='options_home.gravitation=0;';    
end
l{end+1}=['options_home.ncount=' num2str(options_general.ncount_multiplier_optimize) '*9e6;'];
l{end+1}=['options_home.mpi=' num2str(options_general.mpi) ';'];
l{end+1}='options_home.dir=[cpath ''/'' filename];';
l{end+1}='options_home.OutputFcn=''fminplot'';';
l{end+1}='options_home.MaxFunEvals=2500;';
l{end+1}='options_home.MaxIter=2500;';
l{end+1}='options_home.monitors=''Div2d_sample_B'';';
l{end+1}='options_home.type=''maximize'';';
l{end+1}='options_home.optimizer=''fminpso'';';
l{end+1}='options_home.TolFun =''0.03%%'';';
l{end+1}='options_home.TolX =''0.07%%'';';
if options_general.gravity == 1
l{end+1}='options_home.gravitation=1;';
else
l{end+1}='options_home.gravitation=0;';
end
%%%%LELAND MODIFICATION TO INCLUDE MACHINE LIST%%%%%
if isfield(options_general,'machines')
    l{end+1} = ['options_home.machines=' '''' options_general.machines '''' ';'];
end
%%%%%END OF LELAND MADIFICATION%%%%%
l{end+1}='';
l{end+1}='';
l{end+1}='%% rest of options for cluster';
l{end+1}='if (select==2); load NUMCORES.DAT; warning off; else NUMCORES=0; end;';
l{end+1}='';
l{end+1}='options_cluster.mpi=NUMCORES;';
l{end+1}='options_cluster.Display=[];';
l{end+1}='';
l{end+1}='%% optimization options';
l{end+1}='options_cluster.type=''maximize'';';
l{end+1}='options_cluster.optimizer=''fminpso'';';
l{end+1}='options_cluster.TolFun =''0.03%%'';';
l{end+1}='options_cluster.TolX =''0.07%%'';';
l{end+1}='options_cluster.mode=''optimize'';';
l{end+1}='options_cluster.monitors=''Div2d_sample_B'';';
l{end+1}='options_cluster.MaxFunEvals=10000;';
l{end+1}='options_cluster.MaxIter=10000;';
l{end+1}='';
l{end+1}='%% general opptions';
l{end+1}='options_cluster.dir=filename;';
if options_general.gravity == 1
l{end+1}='options_cluster.gravitation=1;';
else
l{end+1}='options_cluster.gravitation=0;';
end
%%%%LELAND MODIFICATION TO INCLUDE MACHINE LIST%%%%%
if isfield(options_general,'machines')
    l{end+1} = ['options_cluster.machines=' '''' options_general.machines '''' ';'];
end
%%%%%END OF LELAND MADIFICATION%%%%%
l{end+1}='options_cluster.overwrite=1;';
l{end+1}='options_cluster.compile=0;';
l{end+1}=['options_cluster.ncount=' num2str(options_general.ncount_multiplier_optimize) '*9e6;']; %9e6 standard
l{end+1}='';
l{end+1}='options={options_home options_cluster};';
l{end+1}='';
l{end+1}='options{select}';
l{end+1}='';
if need_optimizer
l{end+1}='[pars,monitor,m,o]=mcstas([instrument_name ''_optimize.instr''],p,options{select});';
l{end+1}=['optimized = 1;'];
else
    l{end+1}=['optimized = 0;'];
end 
l{end+1}='';
%l{end+1}='save([filename ''optim_res.mat'']);'; Debugging
% Insert code here for doing a coating optimization on the resulting
% optimal geometrical parameters. Use _optimize_seg.instr
l{end+1}='';
l{end+1}='%%------------------------ Analyzing the resulting optimum ------';
l{end+1}='';
l{end+1}='%%options_single_home.compile=0;';
if options_general.gravity == 1
l{end+1}='options_single_home.gravitation=1;';
else
l{end+1}='options_single_home.gravitation=0;';
end
%%%%LELAND MODIFICATION TO INCLUDE MACHINE LIST%%%%%
if isfield(options_general,'machines')
    l{end+1} = ['options_single_home.machines=' '''' options_general.machines '''' ';'];
end
%%%%%END OF LELAND MADIFICATION%%%%%
l{end+1}=['options_single_home.ncount=' num2str(options_general.ncount_multiplier_analyze) '*5e8;'];
l{end+1}=['options_single_home.mpi=' num2str(options_general.mpi) ';'];
l{end+1}='options_single_home.dir=[cpath ''/'' filename ''waveALL''];';
l{end+1}='';
l{end+1}='';
l{end+1}='%% rest of options for cluster';
l{end+1}='if (select==2); load NUMCORES.DAT; warning off; else NUMCORES=0; end;';
l{end+1}='';
l{end+1}='%% general opptions';
l{end+1}='options_single_cluster.mpi=NUMCORES;';
l{end+1}='options_single_cluster.Display=[];';
l{end+1}='options_single_cluster.dir=[filename ''waveALL''];';
if options_general.gravity == 1
l{end+1}='options_single_cluster.gravitation=1;';
else
l{end+1}='options_single_cluster.gravitation=0;';
end
%%%%LELAND MODIFICATION TO INCLUDE MACHINE LIST%%%%%
if isfield(options_general,'machines')
    l{end+1} = ['options_single_cluster.machines=' '''' options_general.machines '''' ';'];
end
%%%%%END OF LELAND MADIFICATION%%%%%
l{end+1}='options_single_cluster.compile=0;';
l{end+1}=['options_single_cluster.ncount=' num2str(options_general.ncount_multiplier_analyze) '*5e8;']; % 5e8 standard
l{end+1}='';
l{end+1}='options_single={options_single_home options_single_cluster};';
l{end+1}='';
strings='wavecenters=[ ';
for i=1:NumSnaps
    strings=[strings num2str(wavec(i)) ' '];
end
strings=[strings '];'];
l{end+1}=strings;
l{end+1}='snapwidth=0.01;';
l{end+1}='';
if need_optimizer
l{end+1}='optimal=monitor.Data.Parameters;';
l{end+1}='names=fieldnames(optimal);';
l{end+1}='for i=1:length(names)';
l{end+1}='optimal.(names{i})=num2str(optimal.(names{i}));';
l{end+1}='end';
else
l{end+1}='optimal=p;';
l{end+1}='names=fieldnames(optimal);';
end
l{end+1}='';
l{end+1}='%% Run for fom wavelengths';
l{end+1}='optimal_visualizer = optimal;';
l{end+1}='optimal_visualizer.file_name = [filename ''_geometry.dat'']';
l{end+1}='monitor_fom=mcstas([instrument_name ''_analyze.instr''],optimal_visualizer,options_single{select});';
l{end+1}='';
l{end+1}='optimal.WaveMin=''0.1'';';

if max(ismember(fieldnames(options_general),'max_wavelength_investigated_multiplier')) == 0
    if (demands.WaveLmax<3)
    l{end+1}=['optimal.WaveMax=''' num2str(demandscan.WaveLmax*3) ''';'];
    else
    l{end+1}=['optimal.WaveMax=''' num2str(demandscan.WaveLmax*2) ''';'];
    end
else 
    l{end+1}=['optimal.WaveMax=''' num2str(demandscan.WaveLmax*options_general.max_wavelength_investigated_multiplier) ''';'];
end

l{end+1}='MaxWB=str2num(optimal.WaveMax);';
l{end+1}='';
l{end+1}='%% Run for all wavelengths';
l{end+1}=['options_single{select}.dir=[filename ''waveLarge''];'];
l{end+1}='optimal_visualizer = optimal;';
l{end+1}='optimal_visualizer.file_name = [filename ''_geometry.dat''];';
l{end+1}='monitor_ALLW=mcstas([instrument_name ''_analyze.instr''],optimal_visualizer,options_single{select});';
l{end+1}='';
l{end+1}='%% Run for rough ESS source spectrum';
l{end+1}=['options_single{select}.dir=[filename ''waveESS''];'];
l{end+1}='optimal_ess = optimal;';
if strcmp(source_component,'ESS_pancake')
l{end+1}='optimal_ess.scan1=num2str(scani);';
l{end+1}='optimal_ess.scan2=num2str(scanj);';
end
l{end+1}='monitor_ESSW=mcstas([instrument_name ''_analyze_ess.instr''],optimal_ess,options_single{select});';
l{end+1}='';
l{end+1}='%% Run for robustness check';
l{end+1}='';
l{end+1}='degraded=optimal;';
l{end+1}='';
l{end+1}='for i=1:length(names)';
l{end+1}='tmpname=names{i}';
l{end+1}=' if strcmp(tmpname(1:1),''m'') && length(tmpname)<4';
l{end+1}='  degraded.(names{i})=num2str(0.8*str2num(optimal.(names{i})))';
l{end+1}=' end';
l{end+1}=' if length(tmpname)>5';
l{end+1}='  if strcmp(tmpname(1:5),''alpha'')';
l{end+1}='   degraded.(names{i})=num2str(1.4*str2num(optimal.(names{i})))';
l{end+1}='  end';
l{end+1}=' end';
l{end+1}='end';
l{end+1}='';
l{end+1}='%% Run for all wavelengths';
l{end+1}=['options_single{select}.dir=[filename ''Large_degraded''];'];
l{end+1}='degraded_visualizer = degraded;';
l{end+1}='degraded_visualizer.file_name = [filename ''_geometry.dat''];';
l{end+1}='monitor_ALLW_degraded=mcstas([instrument_name ''_analyze.instr''],degraded_visualizer,options_single{select});';
l{end+1}='';
l{end+1}='%% Run for rough ESS source spectrum';
l{end+1}=['options_single{select}.dir=[filename ''ESS_degraded''];'];
l{end+1}='degraded_ess = degraded;';
if strcmp(source_component,'ESS_pancake')
l{end+1}='degraded_ess.scan1=num2str(scani);';
l{end+1}='degraded_ess.scan2=num2str(scanj);';
end
l{end+1}='monitor_ESSW_degraded=mcstas([instrument_name ''_analyze_ess.instr''],degraded_ess,options_single{select});';
l{end+1}='';
l{end+1}='%% Run for individual wavelength snapshots';
l{end+1}=['options_signle_cluster.ncount=' num2str(options_general.ncount_multiplier_analyze) '*1e8;'];
l{end+1}='options_single={options_single_home options_single_cluster};';
for i=1:NumSnaps
l{end+1}='';    
l{end+1}=['options_single{select}.dir=[filename ''wave' num2str(i) '''];'];
l{end+1}=['optimal.WaveMin=num2str(wavecenters(' num2str(i) ') - snapwidth*0.5);'];
l{end+1}=['optimal.WaveMax=num2str(wavecenters(' num2str(i) ') + snapwidth*0.5);'];
l{end+1}='optimal_visualizer = optimal;';
l{end+1}='optimal_visualizer.file_name = [filename ''_geometry.dat''];';
l{end+1}=['monitor_W.wave' num2str(i) '=mcstas([instrument_name ''_analyze.instr''],optimal_visualizer,options_single{select});'];
if detailed_absolute_run
l{end+1}=['options_single{select}.dir=[filename ''wave_ess' num2str(i) '''];'];
l{end+1}='optimal_ess = optimal;';
if strcmp(source_component,'ESS_pancake')
l{end+1}='optimal_ess.scan1=num2str(scani);';
l{end+1}='optimal_ess.scan2=num2str(scanj);';
end
l{end+1}=['monitor_W_ess.wave' num2str(i) '=mcstas([instrument_name ''_analyze_ess.instr''],optimal_ess,options_single{select});'];
end
l{end+1}='';
end
l{end+1}='';
l{end+1}='check1 = size(monitor_ALLW);';
l{end+1}='check2 = size(monitor_ALLW_degraded);';
l{end+1}='check3 = size(monitor_ESSW);';
l{end+1}='check4 = size(monitor_ESSW_degraded);';
l{end+1}='';
l{end+1}='check_logic5 = true;';
l{end+1}='for i=1:length(fieldnames(monitor_W))';
l{end+1}='  field_name_vector = fieldnames(monitor_W);';
l{end+1}='  check5 = size(monitor_W.(field_name_vector{i}));';
l{end+1}='  if check5(1) < 3';
l{end+1}='    check_logic5 = false;';
l{end+1}='  end';
l{end+1}='end';
l{end+1}='';
l{end+1}='input_length = length(fieldnames(p));';
if need_optimizer
l{end+1}='output_length = length(fieldnames(monitor.Data.Parameters));';
end
l{end+1}='';
l{end+1}='check_logic1 = check1(1) > 3;';
l{end+1}='check_logic2 = check2(1) > 3;';
l{end+1}='check_logic3 = check3(1) > 3;';
l{end+1}='check_logic4 = check4(1) > 3;';
if need_optimizer
l{end+1}='input_output_logic = input_length == output_length;';
end
l{end+1}='';
if need_optimizer
l{end+1}='if check_logic1 + check_logic2 + check_logic3 + check_logic4 + check_logic5 + input_output_logic == 6';
else
l{end+1}='if check_logic1 + check_logic2 + check_logic3 + check_logic4 + check_logic5 == 5';
end
l{end+1}='%% Everything went well! Save the data and report the succes to the logbook.';
l{end+1}='';
l{end+1}='save([filename ''_all.mat'']);';
l{end+1}='save([cpath ''/../output/analysis/'' instrument_name scanname ''_all.mat'']);';
l{end+1}='copyfile([filename ''_geometry.dat''],[cpath ''/../output/analysis/'' filename ''_geometry.dat''])';
l{end+1}='';
l{end+1}='fid = fopen([cpath ''/../output/analysis/analyze_all_ifit.m''], ''a'');';
l{end+1}='fprintf(fid,[''clear all;clc;close all;\\n'' instrument_name scanname ''_ifit_analyse\\n'']);';
l{end+1}='fclose(fid);';
l{end+1}='';
l{end+1}='string='''';';
l{end+1}='for i=1:length(names)';
l{end+1}='string=[string  names{i} ''='' optimal.(names{i}) ''\\n''];';
l{end+1}='end';
l{end+1}='';
l{end+1}='fid = fopen([cpath ''/'' filename ''.par''], ''w'');';
l{end+1}='fprintf(fid,string);';
l{end+1}='fclose(fid);';
l{end+1}='';
l{end+1}='flux=monitor_fom(11).Data.values(1);';
l{end+1}='';
l{end+1}='fid = fopen([cpath ''/../master_record-done'' scanname ''.txt''],''a'');';
l{end+1}='fprintf(fid,[num2str(flux) '' = '' filename '' - '' inputstring ''\\n''])'; % Add a value describing the flux and BT
l{end+1}='fclose(fid);';
l{end+1}='';
l{end+1}='else';
l{end+1}='%% Well, something went wrong, not all needed data is present. Add to suggested reruns!';
l{end+1}='';
l{end+1}='i=0;done=''1'';';
l{end+1}='while(str2num(done)==1);';
l{end+1}=' i=1+i;';
l{end+1}=' if i==1';
l{end+1}='  fail_name=''suggested_reruns-fail.sh'';';
l{end+1}=' else';
l{end+1}='  fail_name=[''suggested_reruns-fail'' num2str(i) ''.sh''];';
l{end+1}=' end;';
l{end+1}=' [a,done]=unix([''test -d \'' cpath ''/../'' fail_name '' && echo "1" || echo "0"'']);';
l{end+1}='end;';
l{end+1}='';
l{end+1}='fid = fopen([cpath ''/../'' fail_name],''a'');';
l{end+1}='fprintf(fid,[''cd ./'' instrument_name ''\\n sbatch'' instrument_name scanname ''.batch \\n cd .. \\n'']);';
l{end+1}='fclose(fid);';
l{end+1}='';
l{end+1}='end';
% Want to make basic plots directly by the output script.
% Was not possible duo to cluster/MATLAB constraints

analysis_start=length(l);
l{end+1}='';
l{end+1}='MaxIndex=length(wavecenters)+1;';
l{end+1}='Foverall=figure';
l{end+1}='set(Foverall, ''Position'', [0 0 700 1000])';
l{end+1}='subplot(4,2,1)';
l{end+1}='if select==1; set(gcf, ''Renderer'', ''painters''); end;';
l{end+1}='Fposdiv=figure';
l{end+1}='set(Fposdiv, ''Position'', [0 0 700 1000])';
l{end+1}='subplot(MaxIndex,2,1)';
l{end+1}='if select==1; set(gcf, ''Renderer'', ''painters''); end;';
l{end+1}='Facceptance=figure';
l{end+1}='set(Facceptance, ''Position'', [0 0 700 1000])';
l{end+1}='subplot(MaxIndex,2,1)';
l{end+1}='if select==1; set(gcf, ''Renderer'', ''painters''); end;';
l{end+1}='background=figure';
l{end+1}='set(background, ''Position'', [0 0 700 1000])';
l{end+1}='subplot(4,2,1)';
l{end+1}='if select==1; set(gcf, ''Renderer'', ''painters''); end;';
l{end+1}='';
l{end+1}='colors={''r'' ''g'' ''b'' ''k'' ''m''};';
l{end+1}='linecolor=''k'';';
l{end+1}='linethick=2;';
l{end+1}='';
% This code is unecessary as the brilliance files are made so that the
% moderator size is always sufficient.
%l{end+1}='brill_logic1 = str2num(p.mod_x) >= str2num(p.sizeX);';
%l{end+1}='brill_logic2 = str2num(p.mod_y) >= str2num(p.sizeY);';
l{end+1}='reflogic=exist([''../brill_ref/brilliance_ref'' scanname ''.mat'' ])>0.5;';
%l{end+1}='';
%l{end+1}='reflogic = reffile == 1 & brill_logic1 + brill_logic2 == 2;';
l{end+1}='';
l{end+1}='if reflogic';
l{end+1}='load([''../brill_ref/brilliance_ref'' scanname ''.mat'' ]);';
l{end+1}='end';
l{end+1}='';
l{end+1}='for i=1:MaxIndex';
l{end+1}='   if i==MaxIndex';
%l{end+1}='    monitors=monitor_ALLW;';
%l{end+1}='    if reflogic; monitors_ref=monitor_ALLW_ref; end;';
l{end+1}='    monitors=monitor_fom;';
l{end+1}='    if reflogic; monitors_ref=monitor_fom_ref; end;';
l{end+1}='   else';
l{end+1}='    Wnames=fieldnames(monitor_W);';
l{end+1}='    monitors=monitor_W.(Wnames{i});';
l{end+1}='    if reflogic; monitors_ref=monitor_W_ref.(Wnames{i}); end;';
l{end+1}='   end';
l{end+1}='   ';
l{end+1}='   if reflogic';
l{end+1}='   guide_end_lambda = monitors(8);';
l{end+1}='   DIV2D=monitors(1)/monitors_ref(1).Data.Mean;';
l{end+1}='   PSD2D=monitors(11)/monitors_ref(11).Data.Mean;';
l{end+1}='   HPSD=monitors(4)/monitors_ref(4).Data.Mean;';
l{end+1}='   VPSD=monitors(13)/monitors_ref(13).Data.Mean;';
l{end+1}='   HDIV=monitors(6)/monitors_ref(6).Data.Mean;';
l{end+1}='   VDIV=monitors(15)/monitors_ref(15).Data.Mean;';
l{end+1}='   ACCPX=monitors(17)/monitors_ref(17).Data.Mean;';
l{end+1}='   ACCPY=monitors(19)/monitors_ref(19).Data.Mean;';
%l{end+1}='   LAMBDA_B_RAW=monitors(10).Data.Mean;';
l{end+1}='   LAMBDA_B_RAW=monitors(10).Data.values(1);';
l{end+1}='   LAMBDA_B_RAW_monitor=monitors(10);';
l{end+1}='   LAMBDA_B=monitors(10)/monitors_ref(10).Data.Mean;';
l{end+1}='   DIV2D_B=monitors(3)/monitors_ref(3).Data.Mean;';
l{end+1}='   PSD2D_B=monitors(12)/monitors_ref(12).Data.Mean;';
% Allways use mean, the brilliance measurement is always made on a
% moderator of sufficient size for even phase-space illumination.
%if reqscan.moderator_size_x*100 > 4*demandscan.Hsize; %cm and m
l{end+1}='   HPSD_B=monitors(5)/monitors_ref(5).Data.Mean;';
%else
%l{end+1}='   HPSD_B=monitors(13)/monitors_ref(13).Data.Max;';
%end
%if reqscan.moderator_size_y*100 > 4*demandscan.Vsize; %cm and m
l{end+1}='   VPSD_B=monitors(14)/monitors_ref(14).Data.Mean;';
%else
%l{end+1}='   VPSD_B=monitors(14)/monitors_ref(14).Data.Max;';
%end
l{end+1}='   HDIV_B=monitors(7)/monitors_ref(7).Data.Mean;';
l{end+1}='   VDIV_B=monitors(16)/monitors_ref(16).Data.Mean;';
l{end+1}='   HACCP_B=monitors(18)/monitors_ref(18).Data.Mean;';
l{end+1}='   VACCP_B=monitors(20)/monitors_ref(20).Data.Mean;';
l{end+1}='   LAMBDA=monitors(9)/monitors_ref(9).Data.Mean;';
l{end+1}='   LAMBDA_RAW_monitor=monitors(9);';
l{end+1}='   AROUND_SAMPLE_BT=monitors(2).Data.values(1)/monitors_ref(2).Data.values(1);';
l{end+1}='   else';
l{end+1}='   guide_end_lambda = monitors(8);';
l{end+1}='   DIV2D=monitors(2);';
l{end+1}='   PSD2D=monitors(11);';
l{end+1}='   HPSD=monitors(4);';
l{end+1}='   VPSD=monitors(13);';
l{end+1}='   HDIV=monitors(6);';
l{end+1}='   VDIV=monitors(15);';
l{end+1}='   HACCP=monitors(17);';
l{end+1}='   VACCP=monitors(19);';
l{end+1}='   LAMBDA_B_RAW=monitors(10).Data.values(1);';
l{end+1}='   LAMBDA_B_RAW_monitor=monitors(10);';
l{end+1}='   LAMBDA_B=monitors(10);';
l{end+1}='   DIV2D_B=monitors(3);';
l{end+1}='   PSD2D_B=monitors(12);';
l{end+1}='   HPSD_B=monitors(5);';
l{end+1}='   VPSD_B=monitors(14);';
l{end+1}='   HDIV_B=monitors(7);';
l{end+1}='   VDIV_B=monitors(16);';
l{end+1}='   HACCP_B=monitors(18);';
l{end+1}='   VACCP_B=monitors(20);';
l{end+1}='   LAMBDA=monitors(9);';
l{end+1}='   LAMBDA_RAW_monitor=monitors(9);';
l{end+1}='   end';
l{end+1}='   ';
l{end+1}='   if i==MaxIndex';
l{end+1}='      figure(Foverall);';
l{end+1}='      subplot(4,2,7:8)';
l{end+1}='      axis off;';
l{end+1}='      text(0.5,1,Project_name,''interpreter'',''none'')';
l{end+1}='      text(0,0.83,[instrument_name '' - '' filename '' - ''  inputstring],''interpreter'',''none'')';
l{end+1}='      text(0,0.66,[''Sample size: Horizontal='' p.sizeX ''m, Vertical='' p.sizeY ''m, Divergence requirement: Horizontal='' p.divreq_x ''deg, Vertical= '' p.divreq_y ''deg, Sample distance ='' p.sample_dist ''m''],''interpreter'',''none'');';
l{end+1}='      text(0,0.5,[''Intensity on sample of 100 emitted = '' num2str(LAMBDA.Data.values(1)) '' (no divergence limit) '' num2str(LAMBDA_B_RAW) '' (width divergence limits)'' ],''interpreter'',''none'');';
l{end+1}='      text(0,0.33,[''Intensity at guide end of 100 emitted = '' num2str(guide_end_lambda.Data.values(1)) '' (no divergence limit)'' ],''interpreter'',''none'');';
l{end+1}='   if reflogic';
l{end+1}='      text(0,0.16,[''BT on sample = '' num2str(LAMBDA_B.Data.values(1)) '' (width divergence limits)'' ],''interpreter'',''none'');';
l{end+1}='      text(0,0,[''BT near sample = '' num2str(AROUND_SAMPLE_BT) '' (width divergence limits)'' ],''interpreter'',''none'');';
l{end+1}='   end';
l{end+1}='   ';
l{end+1}='      figure(Foverall)';
l{end+1}='      subplot(4,2,1:2) ';
l{end+1}='      ';
l{end+1}='      %%plot(LAMBDA_B,''k'')';
l{end+1}='    if reflogic';
l{end+1}='      axis([0 MaxWB 0 1])';
l{end+1}='    else';
l{end+1}='      axis([0 MaxWB 0 LAMBDA_B.Data.Max*1.1])';
l{end+1}='    end';
l{end+1}='      %%set(gca,''XTick'',[0 0.5 1.0 1.5 2.0 2.5 3.0])';
l{end+1}='      title([''Lambda dependence, + is wavelengths for 1d graphs. I='' num2str(LAMBDA_B_RAW)])';
l{end+1}='      ylabel(''Brilliance transfer'')';
l{end+1}='    if reflogic';
l{end+1}='      markerheight=0.1;';
l{end+1}='    else';
l{end+1}='      markerheight=LAMBDA_B.Data.Max*0.2;';
l{end+1}='    end';
l{end+1}='        for j=1:length(wavecenters)';
l{end+1}='                hold on';
l{end+1}='                plot(wavecenters(j),markerheight,[''+'' colors{j}],''MarkerSize'',12)';
l{end+1}='                hold off';
l{end+1}='        end';
l{end+1}='   else';
l{end+1}='   ';
l{end+1}='    figure(Foverall)';
l{end+1}='    subplot(4,2,3)';
l{end+1}='    hold on';
l{end+1}='    plot(HDIV_B,colors{i})';
l{end+1}='    box';
l{end+1}='    if reflogic';
l{end+1}=['    axis([' num2str(-div1d*demandscan.Hdiv) ' ' num2str(div1d*demandscan.Hdiv) ' 0 1])'];
l{end+1}='    else';
l{end+1}=['    if i==MaxIndex-1; axis([' num2str(-3*demandscan.Hdiv) ' ' num2str(3*demandscan.Hdiv) ' 0 HDIV_B.Data.Max*1.1]); end;'];
l{end+1}='    end';
l{end+1}='    %%set(gca,''XTick'',[-1 -0.5 -0.3 0 0.3 0.5 1])';
l{end+1}='    title(''Horizontal div'')';
l{end+1}='    ylabel(''Brilliance transfer'')';
l{end+1}='';
l{end+1}='    subplot(4,2,4)';
l{end+1}='    hold on';
l{end+1}='    plot(VDIV_B,colors{i})';
l{end+1}='    box';
l{end+1}='    if reflogic';
l{end+1}=['    axis([' num2str(-div1d*demandscan.Vdiv) ' ' num2str(div1d*demandscan.Vdiv) ' 0 1])'];
l{end+1}='    else';
l{end+1}=['    if i==MaxIndex-1; axis([' num2str(-3*demandscan.Vdiv) ' ' num2str(3*demandscan.Vdiv) ' 0 VDIV_B.Data.Max*1.1]); end;'];
l{end+1}='    end';
l{end+1}='    %%set(gca,''XTick'',[-1 -0.5 -0.3 0 0.3 0.5 1])';
l{end+1}='    title(''Vertical div'')';
l{end+1}='    ylabel(''Brilliance transfer'')';
l{end+1}='';
l{end+1}='    subplot(4,2,5)';
l{end+1}='    hold on';
l{end+1}='    plot(HPSD_B,colors{i})';
l{end+1}='    box';
l{end+1}='    %%axis([-0.75 0.75 0 1])';
l{end+1}='    if reflogic';
l{end+1}=['    axis([' num2str(-0.5*psd1d*demandscan.Hsize/100) ' ' num2str(0.5*psd1d*demandscan.Hsize/100) ' 0 1])'];
l{end+1}='    else';
l{end+1}=['    if i==MaxIndex-1; axis([' num2str(-3*demandscan.Hsize/100) ' ' num2str(3*demandscan.Hsize/100) ' 0 HPSD_B.Data.Max*1.1]); end;'];
l{end+1}='    end';
l{end+1}='    %%set(gca,''XTick'',[-0.75 -0.5 -0.25 0 0.25 0.5 0.75])';
l{end+1}='    title(''Horizontal psd'')';
l{end+1}='    xlabel(''Horizontal position'')';
l{end+1}='    ylabel(''Brilliance transfer'')';
l{end+1}='';
l{end+1}='    subplot(4,2,6)';
l{end+1}='    hold on';
l{end+1}='    plot(VPSD_B,colors{i})';
l{end+1}='    box';
l{end+1}='    %%axis([-1 1 0 1])';
l{end+1}='    if reflogic';
l{end+1}=['    axis([' num2str(-0.5*psd1d*demandscan.Vsize/100) ' ' num2str(0.5*psd1d*demandscan.Vsize/100) ' 0 1])'];
l{end+1}='    else';
l{end+1}=['    if i==MaxIndex-1; axis([' num2str(-3*demandscan.Vsize/100) ' ' num2str(3*demandscan.Vsize/100) ' 0 VPSD_B.Data.Max*1.1]); end;'];
l{end+1}='    end';
l{end+1}='    %%set(gca,''XTick'',[-1 -0.5 -0.25 0 0.25 0.5 1])';
l{end+1}='    title(''Vertical psd'') ';
l{end+1}='    xlabel(''Vertical position'')';
l{end+1}='    ylabel(''Brilliance transfer'') ';
l{end+1}='   ';
l{end+1}='   end';
l{end+1}='   ';
l{end+1}='   ';
l{end+1}='    figure(Fposdiv)';
l{end+1}='    subplot(MaxIndex,2,(i-1)*2+1)';
l{end+1}='    plot(DIV2D_B)';
l{end+1}='    box';
l{end+1}='    view([-0.5 90]);';
l{end+1}='    %%axis([-1 1 -1 1])';
l{end+1}=['   axis([' num2str(-div2d*demandscan.Hdiv) ' ' num2str(div2d*demandscan.Hdiv) ' ' num2str(-div2d*demandscan.Vdiv) ' ' num2str(div2d*demandscan.Vdiv) '])'];
l{end+1}='    maxi=max(DIV2D_B)+100;';
l{end+1}='    hold on';
l{end+1}=['    x=' num2str(demandscan.Hdiv) ';'];
l{end+1}=['    y=' num2str(demandscan.Vdiv) ';'];
l{end+1}='    plot3([-x x], [-y -y], [maxi maxi],''color'',linecolor,''LineWidth'',linethick)';
l{end+1}='    plot3([-x x], [y y], [maxi maxi],''color'',linecolor,''LineWidth'',linethick)';
l{end+1}='    plot3([-x -x], [-y y], [maxi maxi],''color'',linecolor,''LineWidth'',linethick)';
l{end+1}='    plot3([x x], [-y y], [maxi maxi],''color'',linecolor,''LineWidth'',linethick)';
l{end+1}='    hold off';
l{end+1}='    %%set(gca,''XTick'',[-1 -0.5 -0.3 0 0.3 0.5 1])';
l{end+1}='    %%set(gca,''YTick'',[-1 -0.5 -0.3 0 0.3 0.5 1])';
l{end+1}='    xlabel(''Horizontal divergence [deg]'');';
l{end+1}='    ylabel(''Vertical divergence [deg]'');';
l{end+1}='    if i==MaxIndex';
l{end+1}='        title(''2d div (wavelength interval)'')';
l{end+1}='    else';
l{end+1}='        title([''2d div Lambda='' num2str(wavecenters(i)) ''A''])';
l{end+1}='    end';
l{end+1}='';
l{end+1}='';
l{end+1}='    subplot(MaxIndex,2,(i-1)*2+2)';
l{end+1}='    plot(PSD2D_B)';
l{end+1}='    box';
l{end+1}='    view([-0.5 90]);';
l{end+1}='    %%axis([-1 1 -1 1])';
l{end+1}=['   axis([' num2str(-0.5*psd2d*demandscan.Hsize) ' ' num2str(0.5*psd2d*demandscan.Hsize) ' ' num2str(-0.5*psd2d*demandscan.Vsize) ' ' num2str(0.5*psd2d*demandscan.Vsize) '])'];
l{end+1}='    maxi=max(PSD2D_B)+100;';
l{end+1}='    hold on';
l{end+1}=['    x=' num2str(demandscan.Hsize*0.5) ';'];
l{end+1}=['    y=' num2str(demandscan.Vsize*0.5) ';'];
l{end+1}='    plot3([-x x], [-y -y], [maxi maxi],''color'',linecolor,''LineWidth'',linethick)';
l{end+1}='    plot3([-x x], [y y], [maxi maxi],''color'',linecolor,''LineWidth'',linethick)';
l{end+1}='    plot3([-x -x], [-y y], [maxi maxi],''color'',linecolor,''LineWidth'',linethick)';
l{end+1}='    plot3([x x], [-y y], [maxi maxi],''color'',linecolor,''LineWidth'',linethick)';
l{end+1}='    hold off';
l{end+1}='    %%set(gca,''XTick'',[-1 -0.5 -0.3 0 0.3 0.5 1])';
l{end+1}='    %%set(gca,''YTick'',[-1 -0.5 -0.3 0 0.3 0.5 1])';
l{end+1}='    xlabel(''Horizontal position [cm]'');';
l{end+1}='    ylabel(''Vertical position [cm]'');';
l{end+1}='    if i==MaxIndex';
l{end+1}='        title(''2d psd (wavelength interval)'')';
l{end+1}='    else';
l{end+1}='        title([''2d psd Lambda='' num2str(wavecenters(i)) ''A''])';
l{end+1}='    end';
l{end+1}='   ';
l{end+1}='    ';
l{end+1}='    figure(Facceptance)';
l{end+1}='    subplot(MaxIndex,2,(i-1)*2+1)';
l{end+1}='    box';
l{end+1}='    plot(HACCP_B)';
l{end+1}='    view([-0.5 90]);';
l{end+1}='    %%axis([-1 1 -1 1])';
l{end+1}=['   axis([' num2str(-0.5*apsd*demandscan.Hsize/100) ' ' num2str(0.5*apsd*demandscan.Hsize/100) ' ' num2str(-adiv*demandscan.Hdiv) ' ' num2str(adiv*demandscan.Hdiv) '])'];
l{end+1}='    maxi=max(HACCP_B)+100;';
l{end+1}='    hold on';
l{end+1}=['    x=' num2str(demandscan.Hsize*0.01*0.5) ';'];
l{end+1}=['    y=' num2str(demandscan.Hdiv) ';'];
l{end+1}='    plot3([-x x], [-y -y], [maxi maxi],''color'',linecolor,''LineWidth'',linethick)';
l{end+1}='    plot3([-x x], [y y], [maxi maxi],''color'',linecolor,''LineWidth'',linethick)';
l{end+1}='    plot3([-x -x], [-y y], [maxi maxi],''color'',linecolor,''LineWidth'',linethick)';
l{end+1}='    plot3([x x], [-y y], [maxi maxi],''color'',linecolor,''LineWidth'',linethick)';
l{end+1}='    hold off';
l{end+1}='    %%set(gca,''XTick'',[-1 -0.5 -0.3 0 0.3 0.5 1])';
l{end+1}='    %%set(gca,''YTick'',[-1 -0.5 -0.3 0 0.3 0.5 1])';
l{end+1}='    xlabel(''Horizontal position [m]'')';
l{end+1}='    ylabel(''Horizontal divergence [deg]'')';
l{end+1}='    if i==MaxIndex';
l{end+1}='    title(''Horizontal acceptance diagram (wavelength interval)'')';
l{end+1}='    else';
l{end+1}='    title([''Horizontal accep. dia. Wavelength='' num2str(wavecenters(i)) ''A'']);';
l{end+1}='    end    ';
l{end+1}='';
l{end+1}='    subplot(MaxIndex,2,(i-1)*2+2)';
l{end+1}='    plot(VACCP_B)';
l{end+1}='    box';
l{end+1}='    view([-0.5 90]);';
l{end+1}='    %%axis([-1 1 -1 1])';
l{end+1}=['   axis([' num2str(-0.5*apsd*demandscan.Vsize/100) ' ' num2str(0.5*apsd*demandscan.Vsize/100) ' ' num2str(-adiv*demandscan.Vdiv) ' ' num2str(adiv*demandscan.Vdiv) '])'];
l{end+1}='    maxi=max(VACCP_B)+100;';
l{end+1}='    hold on';
l{end+1}=['    x=' num2str(demandscan.Vsize*0.01*0.5) ';'];
l{end+1}=['    y=' num2str(demandscan.Vdiv) ';'];
l{end+1}='    plot3([-x x], [-y -y], [maxi maxi],''color'',linecolor,''LineWidth'',linethick)';
l{end+1}='    plot3([-x x], [y y], [maxi maxi],''color'',linecolor,''LineWidth'',linethick)';
l{end+1}='    plot3([-x -x], [-y y], [maxi maxi],''color'',linecolor,''LineWidth'',linethick)';
l{end+1}='    plot3([x x], [-y y], [maxi maxi],''color'',linecolor,''LineWidth'',linethick)';
l{end+1}='    hold off';
l{end+1}='    %%set(gca,''XTick'',[-1 -0.5 -0.3 0 0.3 0.5 1])';
l{end+1}='    %%set(gca,''YTick'',[-1 -0.5 -0.3 0 0.3 0.5 1])';
l{end+1}='    xlabel(''Vertical position [m]'')';
l{end+1}='    ylabel(''Vertical divergence [deg]'')';
l{end+1}='    if i==MaxIndex';
l{end+1}='    title(''Vertical acceptance diagram (wavelength interval)'')';
l{end+1}='    else';
l{end+1}='    title([''Vertical accep. dia. Wavelength='' num2str(wavecenters(i)) ''A'']);';
l{end+1}='    end';
l{end+1}='    ';
l{end+1}='end';
l{end+1}='';
l{end+1}='if reflogic';
l{end+1}='flux=monitor_fom(10).Data.values(1)/monitor_fom_ref(10).Data.values(1);';
l{end+1}='fluxtext='' BT '';';
l{end+1}='else';
l{end+1}='flux=monitor_fom(10).Data.values(1);';
l{end+1}='fluxtext='' Absolute '';';
l{end+1}='end';
l{end+1}='';
l{end+1}='fid = fopen([cpath ''/master_record-analyzed'' scanname ''.txt''],''a'');';
l{end+1}='fprintf(fid,[num2str(flux) fluxtext ''= '' filename '' - '' inputstring ''\\n''])'; % Add a value describing the flux and BT
l{end+1}='fclose(fid);';
l{end+1}='';
l{end+1}='set(Foverall, ''paperpositionmode'', ''auto'');';
l{end+1}='set(Fposdiv, ''paperpositionmode'', ''auto'');';
l{end+1}='set(Facceptance, ''paperpositionmode'', ''auto'');';
l{end+1}='';
l{end+1}='if reflogic';
l{end+1}='hold on';
l{end+1}='figure(Foverall)';
l{end+1}='subplot(4,2,1:2)';
l{end+1}='hold on';
l{end+1}='plot(monitor_ALLW_degraded(10)/monitor_ALLW_ref(10),''r'')';
l{end+1}='plot(monitor_ALLW(10)/monitor_ALLW_ref(10),''k'')';
l{end+1}='ylabel(''Brilliance transfer'')';
l{end+1}='title([''Lambda dependence, + is wavelengths for 1d graphs. I='' num2str(LAMBDA_B_RAW)])';
l{end+1}='hold off';
l{end+1}='box';
l{end+1}='else';
l{end+1}='figure(Foverall)';
l{end+1}='subplot(4,2,1:2)';
l{end+1}='hold on';
l{end+1}='plot(monitor_ALLW_degraded(10),''r'')';
l{end+1}='plot(monitor_ALLW(10),''k'')';
l{end+1}='ylim([0 monitor_ALLW(10).Data.Max*1.1+1e-10])';
l{end+1}='hold off';
l{end+1}='end';
l{end+1}='';
l{end+1}='figure(background)';
%l{end+1}='      subplot(4,2,5:6)';
%l{end+1}='      axis off;';
%l{end+1}='      text(0.5,1,Project_name,''interpreter'',''none'')';
%l{end+1}='      text(0,0.83,[instrument_name '' - '' filename '' - ''  inputstring],''interpreter'',''none'')';
%l{end+1}='      text(0,0.66,[''Sample size: Horizontal='' p.sizeX ''m, Vertical='' p.sizeY ''m, Divergence requirement: Horizontal='' p.divreq_x ''deg, Vertical= '' p.divreq_y ''deg, Sample distance ='' p.sample_dist ''m''],''interpreter'',''none'');';
%l{end+1}='      text(0,0.5,[''Intensity on sample of 100 emitted = '' num2str(LAMBDA.Data.values(1)) '' (no divergence limit) '' num2str(LAMBDA_B_RAW) '' (width divergence limits)'' ],''interpreter'',''none'');';
%l{end+1}='      text(0,0.33,[''Intensity near sample of 100 emitted = '' num2str(PSD2D.Data.values(1)) '' (no divergence limit)'' ],''interpreter'',''none'');';
l{end+1}='subplot(4,2,1:4)';
l{end+1}='plot(monitor_ALLW(9)/monitor_ALLW(8),''k'');';
l{end+1}='hold on';
l{end+1}='plot(monitor_ALLW(10)/monitor_ALLW(8),''b'');';
l{end+1}='hold off';
l{end+1}='box on';
l{end+1}='title(''Signal compared to total guide output'') ';
l{end+1}='xlabel(''Wavelength [AA]'')';
l{end+1}='ylabel(''Signal fraction [Unitless]'') ';
l{end+1}='';
l{end+1}='print(Foverall,''-dpng'',''-r300'',[cpath ''/'' filename ''_overall_pure.png''])';
l{end+1}='print(Fposdiv,''-dpng'',''-r300'',[cpath ''/'' filename ''_posdiv_pure.png''])';
l{end+1}='print(Facceptance,''-dpng'',''-r300'',[cpath ''/'' filename ''_acceptance_pure.png''])';
l{end+1}='';
%l{end+1}='print(Foverall,''-depsc'',[cpath ''/'' filename ''_overall.eps''])';
%l{end+1}='print(Fposdiv,''-depsc'',[cpath ''/'' filename ''_posdiv.eps''])';
%l{end+1}='print(Facceptance,''-depsc'',[cpath ''/'' filename ''_acceptance.eps''])';

% Experimental section which makes the same three graphs for the selected
% source
l{end+1}='';
l{end+1}='if exist(''OnlyPureFigs'') < 0.5';
l{end+1}='close(Foverall);close(Fposdiv);close(Facceptance);';
l{end+1}='MaxIndex=length(wavecenters)+1;';
l{end+1}='Foverall_ess=figure';
l{end+1}='set(Foverall_ess, ''Position'', [0 0 700 1000])';
l{end+1}='subplot(4,2,1)';
l{end+1}='if select==1; set(gcf, ''Renderer'', ''painters''); end;';
l{end+1}='Fposdiv_ess=figure';
l{end+1}='set(Fposdiv_ess, ''Position'', [0 0 700 1000])';
l{end+1}='subplot(MaxIndex,2,1)';
l{end+1}='if select==1; set(gcf, ''Renderer'', ''painters''); end;';
l{end+1}='Facceptance_ess=figure';
l{end+1}='set(Facceptance_ess, ''Position'', [0 0 700 1000])';
l{end+1}='subplot(MaxIndex,2,1)';
l{end+1}='if select==1; set(gcf, ''Renderer'', ''painters''); end;';
l{end+1}='';
l{end+1}='colors={''r'' ''g'' ''b'' ''k'' ''m''};';
l{end+1}='linecolor=''k'';';
l{end+1}='linethick=2;';
l{end+1}='';
%l{end+1}='brill_logic1 = str2num(p.mod_x) >= str2num(p.sizeX);';
%l{end+1}='brill_logic2 = str2num(p.mod_y) >= str2num(p.sizeY);';
%l{end+1}='reffile=exist([''../brill_ref/brilliance_ref'' scanname ''.mat'' ])>0.5;';
%l{end+1}='';
%l{end+1}='reflogic = reffile == 1 & brill_logic1 + brill_logic2 == 2;';
%l{end+1}='';
%l{end+1}='if reflogic';
%l{end+1}='load([''../brill_ref/brilliance_ref'' scanname ''.mat'' ]);';
%l{end+1}='end';
l{end+1}='Wnames=fieldnames(monitor_W_ess);';
l{end+1}='for i=1:length(wavecenters)';
l{end+1}='  current_monitor=monitor_W_ess.(Wnames{i});';
l{end+1}='  for j=2:length(current_monitor)';
l{end+1}='     max_val(i,j) = current_monitor(j).Data.Max;';  
l{end+1}='  end';
l{end+1}='end';
l{end+1}='';
l{end+1}='for j=2:length(current_monitor)';
l{end+1}='     overall_max_val(j) = max(max_val(:,j));';  
l{end+1}='end';
l{end+1}='';
l{end+1}='   HPSD_num = 4;';
l{end+1}='   VPSD_num = 13;';
l{end+1}='   HDIV_num = 6;';
l{end+1}='   VDIV_num = 15;';
l{end+1}='';
l{end+1}='   max_HPSD = overall_max_val(HPSD_num);';
l{end+1}='   max_VPSD = overall_max_val(VPSD_num);';
l{end+1}='   max_HDIV = overall_max_val(HDIV_num);';
l{end+1}='   max_VDIV = overall_max_val(VDIV_num);';


% didn't work
%l{end+1}='for j=2:length(current_monitor)';
%l{end+1}=' for i=1:length(wavecenters)';
%l{end+1}='     monitor_W_ess.(Wnames{i}).Data.guide_bot_overall_maxval = overall_max_val(j); ';  
%l{end+1}=' end';
%l{end+1}='end';

if detailed_absolute_run
l{end+1}='for i=1:MaxIndex';
l{end+1}='   if i==MaxIndex';
%l{end+1}='    monitors=monitor_ALLW;';
%l{end+1}='    if reflogic; monitors_ref=monitor_ALLW_ref; end;';
l{end+1}='    monitors=monitor_ESSW;';
%l{end+1}='    if reflogic; monitors_ref=monitor_fom_ref; end;';
l{end+1}='   else';
l{end+1}='    monitors=monitor_W_ess.(Wnames{i});';
%l{end+1}='    if reflogic; monitors_ref=monitor_W_ref.(Wnames{i}); end;';
l{end+1}='   end';
l{end+1}='   ';
% l{end+1}='   if reflogic';
% l{end+1}='   DIV2D=monitors(2)/monitors_ref(2).Data.Mean;';
% l{end+1}='   PSD2D=monitors(3)/monitors_ref(3).Data.Mean;';
% l{end+1}='   HPSD=monitors(4)/monitors_ref(4).Data.Mean;';
% l{end+1}='   VPSD=monitors(5)/monitors_ref(5).Data.Mean;';
% l{end+1}='   HDIV=monitors(6)/monitors_ref(6).Data.Mean;';
% l{end+1}='   VDIV=monitors(7)/monitors_ref(7).Data.Mean;';
% l{end+1}='   HACCP=monitors(8)/monitors_ref(8).Data.Mean;';
% l{end+1}='   VACCP=monitors(9)/monitors_ref(9).Data.Mean;';
% l{end+1}='   LAMBDA_B_RAW=monitors(10).Data.Mean;';
% l{end+1}='   LAMBDA_B=monitors(10)/monitors_ref(10).Data.Mean;';
% l{end+1}='   DIV2D_B=monitors(11)/monitors_ref(11).Data.Mean;';
% l{end+1}='   PSD2D_B=monitors(12)/monitors_ref(12).Data.Mean;';
% if reqscan.moderator_size_x*100 > 4*demandscan.Hsize; %cm and m
% l{end+1}='   HPSD_B=monitors(13)/monitors_ref(13).Data.Mean;';
% else
% l{end+1}='   HPSD_B=monitors(13)/monitors_ref(13).Data.max;';
% end
% if reqscan.moderator_size_y*100 > 4*demandscan.Vsize; %cm and m
% l{end+1}='   VPSD_B=monitors(14)/monitors_ref(14).Data.Mean;';
% else
% l{end+1}='   VPSD_B=monitors(14)/monitors_ref(14).Data.max;';
% end
% l{end+1}='   HDIV_B=monitors(15)/monitors_ref(15).Data.Mean;';
% l{end+1}='   VDIV_B=monitors(16)/monitors_ref(16).Data.Mean;';
% l{end+1}='   HACCP_B=monitors(17)/monitors_ref(17).Data.Mean;';
% l{end+1}='   VACCP_B=monitors(18)/monitors_ref(18).Data.Mean;';
% l{end+1}='   LAMBDA=monitors(19)/monitors_ref(19).Data.Mean;';
% l{end+1}='   AROUND_SAMPLE_BT=monitors(12).Data.values(1)/monitors_ref(12).Data.values(1);';
% l{end+1}='   else';
l{end+1}='   guide_end_lambda = monitors(8);';
l{end+1}='   DIV2D=monitors(1);';
l{end+1}='   PSD2D=monitors(11);';
l{end+1}='   HPSD=monitors(4);';
l{end+1}='   VPSD=monitors(13);';
l{end+1}='   HDIV=monitors(6);';
l{end+1}='   VDIV=monitors(15);';
l{end+1}='   HACCP=monitors(17);';
l{end+1}='   VACCP=monitors(19);';
l{end+1}='   LAMBDA_B_RAW=monitors(10).Data.values(1);';
l{end+1}='   LAMBDA_B_RAW_monitor=monitors(10);';
l{end+1}='   LAMBDA_B=monitors(10);';
l{end+1}='   DIV2D_B=monitors(3);';
l{end+1}='   PSD2D_B=monitors(12);';
l{end+1}='   HPSD_B=monitors(5);';
l{end+1}='   VPSD_B=monitors(14);';
l{end+1}='   HDIV_B=monitors(7);';
l{end+1}='   VDIV_B=monitors(16);';
l{end+1}='   HACCP_B=monitors(18);';
l{end+1}='   VACCP_B=monitors(20);';
l{end+1}='   LAMBDA=monitors(9);';
l{end+1}='   LAMBDA_RAW_monitor=monitors(9);';
%l{end+1}='   end';
l{end+1}='   ';
l{end+1}='   if i==MaxIndex';
l{end+1}='      figure(Foverall_ess);';
l{end+1}='      subplot(4,2,7:8)';
l{end+1}='      axis off;';
l{end+1}='      text(0.5,1,Project_name,''interpreter'',''none'')';
l{end+1}='      text(0,0.83,[instrument_name '' - '' filename '' - ''  inputstring],''interpreter'',''none'')';
l{end+1}='      text(0,0.66,[''Sample size: Horizontal='' p.sizeX ''m, Vertical='' p.sizeY ''m, Divergence requirement: Horizontal='' p.divreq_x ''deg, Vertical= '' p.divreq_y ''deg, Sample distance ='' p.sample_dist ''m''],''interpreter'',''none'');';
l{end+1}='      text(0,0.5,[''Intensity on sample of 100 emitted = '' num2str(LAMBDA.Data.values(1)) '' (no divergence limit) '' num2str(LAMBDA_B_RAW) '' (width divergence limits)'' ],''interpreter'',''none'');';
l{end+1}='      text(0,0.33,[''Intensity near sample of 100 emitted = '' num2str(PSD2D.Data.values(1)) '' (no divergence limit)'' ],''interpreter'',''none'');';
% l{end+1}='   if reflogic';
% l{end+1}='      text(0,0.16,[''BT on sample = '' num2str(LAMBDA_B.Data.values(1)) '' (width divergence limits)'' ]);';
% l{end+1}='      text(0,0,[''BT near sample = '' num2str(AROUND_SAMPLE_BT) '' (width divergence limits)'' ]);';
% l{end+1}='   end';
l{end+1}='   ';
l{end+1}='      figure(Foverall_ess)';
l{end+1}='      subplot(4,2,1:2) ';
%l{end+1}='      ';
%l{end+1}='      %%plot(LAMBDA_B,''k'')';
%l{end+1}='    if reflogic';
%l{end+1}='      axis([0 MaxWB 0 1])';
%l{end+1}='    else';
%l{end+1}='      axis([0 MaxWB 0 LAMBDA_B.Data.Max*1.1])';
%l{end+1}='    end';
l{end+1}='     figure(Foverall_ess)';
l{end+1}='     subplot(4,2,1:2)';
l{end+1}='     hold on';
%monitor_ESSW(10).Data.Parameters.sizeX
l{end+1}='     sample_area_cm2 = 10000*LAMBDA.Data.Parameters.sizeX*LAMBDA.Data.Parameters.sizeY;';
l{end+1}='     norm_factor = LAMBDA.Data.array_1d(1)/(LAMBDA.Data.xlimits(2)-LAMBDA.Data.xlimits(1))/sample_area_cm2;';
l{end+1}='';
l{end+1}='     plot(LAMBDA*norm_factor,''k'')';
l{end+1}='     ylim([0 LAMBDA.Data.Max*norm_factor*1.1+1e-10]);';
l{end+1}='     hold off';
% Should do something cool here to ensure ticks at wavelength snapshot positions
l{end+1}='      %%set(gca,''XTick'',[0 0.5 1.0 1.5 2.0 2.5 3.0])';
l{end+1}='      title([''Lambda dependence, + is wavelengths for 1d graphs. I='' num2str(LAMBDA_B_RAW)])';
l{end+1}='      ylabel(''Flux n/s/cm^2/AA'')';
l{end+1}='      box';
%l{end+1}='    if reflogic';
%l{end+1}='      markerheight=0.1;';
%l{end+1}='    else';
l{end+1}='      markerheight=LAMBDA.Data.Max*norm_factor*0.1;';
%l{end+1}='    end';
l{end+1}='        for j=1:length(wavecenters)';
l{end+1}='                hold on';
l{end+1}='                plot(wavecenters(j),markerheight,[''+'' colors{j}],''MarkerSize'',12)';
l{end+1}='                hold off';
l{end+1}='        end';
l{end+1}='   else';
l{end+1}='   ';
l{end+1}='    figure(Foverall_ess)';
l{end+1}='    subplot(4,2,3)';
l{end+1}='    hold on';
l{end+1}='    plot(HDIV,colors{i})';
l{end+1}='    box';
%l{end+1}='    if reflogic';
%l{end+1}=['    axis([' num2str(-div1d*demandscan.Hdiv) ' ' num2str(div1d*demandscan.Hdiv) ' 0 1])'];
%l{end+1}='    else';
l{end+1}=['    if i==MaxIndex-1; axis([' num2str(-3*demandscan.Hdiv) ' ' num2str(3*demandscan.Hdiv) ' 0 max_HDIV*1.1+1e-10]); end;'];
%l{end+1}='    end';
l{end+1}='    %%set(gca,''XTick'',[-1 -0.5 -0.3 0 0.3 0.5 1])';
l{end+1}='    title(''Horizontal div'')';
l{end+1}='    ylabel(''Arb intensity unit'')';
l{end+1}='';
l{end+1}='    subplot(4,2,4)';
l{end+1}='    hold on';
l{end+1}='    plot(VDIV,colors{i})';
l{end+1}='    box';
% l{end+1}='    if reflogic';
% l{end+1}=['    axis([' num2str(-div1d*demandscan.Vdiv) ' ' num2str(div1d*demandscan.Vdiv) ' 0 1])'];
% l{end+1}='    else';
l{end+1}=['    if i==MaxIndex-1; axis([' num2str(-3*demandscan.Vdiv) ' ' num2str(3*demandscan.Vdiv) ' 0 max_VDIV*1.1+1e-10]); end;'];
% l{end+1}='    end';
l{end+1}='    %%set(gca,''XTick'',[-1 -0.5 -0.3 0 0.3 0.5 1])';
l{end+1}='    title(''Vertical div'')';
l{end+1}='    ylabel(''Arb intensity unit'')';
l{end+1}='';
l{end+1}='    subplot(4,2,5)';
l{end+1}='    hold on';
l{end+1}='    monitor_area_cm2 = 10000*HPSD.Data.Parameters.sizeY*(HPSD.Data.xlimits(2)-HPSD.Data.xlimits(1));';
l{end+1}='    number_of_bins = HPSD.Data.array_1d;';
l{end+1}='    wavelength_band = HPSD.Data.Parameters.WaveMax-HPSD.Data.Parameters.WaveMin;';
l{end+1}='    norm_factor = number_of_bins/monitor_area_cm2/wavelength_band;';
l{end+1}='    plot(HPSD*norm_factor,colors{i})';
l{end+1}='    box';
l{end+1}='    %%axis([-0.75 0.75 0 1])';
% l{end+1}='    if reflogic';
% l{end+1}=['    axis([' num2str(-0.5*psd1d*demandscan.Hsize/100) ' ' num2str(0.5*psd1d*demandscan.Hsize/100) ' 0 1])'];
% l{end+1}='    else';
l{end+1}=['    if i==MaxIndex-1; axis([' num2str(-3*demandscan.Hsize/100) ' ' num2str(3*demandscan.Hsize/100) ' 0 max_HPSD*norm_factor*1.1+1e-10]); end;'];
%l{end+1}='    end';
l{end+1}='    %%set(gca,''XTick'',[-0.75 -0.5 -0.25 0 0.25 0.5 0.75])';
l{end+1}='    title(''Horizontal psd'')';
l{end+1}='    xlabel(''Horizontal position'')';
l{end+1}='    ylabel(''Flux n/s/cm^2/AA'')';
l{end+1}='';
l{end+1}='    subplot(4,2,6)';
l{end+1}='    hold on';
l{end+1}='    monitor_area_cm2 = 10000*VPSD.Data.Parameters.sizeX*(VPSD.Data.xlimits(2)-VPSD.Data.xlimits(1));';
l{end+1}='    number_of_bins = VPSD.Data.array_1d;';
l{end+1}='    wavelength_band = VPSD.Data.Parameters.WaveMax-VPSD.Data.Parameters.WaveMin;';
l{end+1}='    norm_factor = number_of_bins/monitor_area_cm2/wavelength_band;';
l{end+1}='    plot(VPSD*norm_factor,colors{i})';
l{end+1}='    box';
l{end+1}='    %%axis([-1 1 0 1])';
% l{end+1}='    if reflogic';
% l{end+1}=['    axis([' num2str(-0.5*psd1d*demandscan.Vsize/100) ' ' num2str(0.5*psd1d*demandscan.Vsize/100) ' 0 1])'];
% l{end+1}='    else';
l{end+1}=['    if i==MaxIndex-1; axis([' num2str(-3*demandscan.Vsize/100) ' ' num2str(3*demandscan.Vsize/100) ' 0 max_VPSD*norm_factor*1.1+1e-10]); end;'];
%l{end+1}='    end';
l{end+1}='    %%set(gca,''XTick'',[-1 -0.5 -0.25 0 0.25 0.5 1])';
l{end+1}='    title(''Vertical psd'') ';
l{end+1}='    xlabel(''Vertical position'')';
l{end+1}='    ylabel(''Flux n/s/cm^2/AA'') ';
l{end+1}='   ';
l{end+1}='   end';
l{end+1}='   ';
l{end+1}='   ';
l{end+1}='    figure(Fposdiv_ess)';
l{end+1}='    subplot(MaxIndex,2,(i-1)*2+1)';
l{end+1}='    plot(DIV2D)';
l{end+1}='    box';
l{end+1}='    view([-0.5 90]);';
l{end+1}='    %%axis([-1 1 -1 1])';
l{end+1}=['   axis([' num2str(-div2d*demandscan.Hdiv) ' ' num2str(div2d*demandscan.Hdiv) ' ' num2str(-div2d*demandscan.Vdiv) ' ' num2str(div2d*demandscan.Vdiv) '])'];
l{end+1}='    maxi=max(DIV2D)+100;';
l{end+1}='    hold on';
l{end+1}=['    x=' num2str(demandscan.Hdiv) ';'];
l{end+1}=['    y=' num2str(demandscan.Vdiv) ';'];
l{end+1}='    plot3([-x x], [-y -y], [maxi maxi],''color'',linecolor,''LineWidth'',linethick)';
l{end+1}='    plot3([-x x], [y y], [maxi maxi],''color'',linecolor,''LineWidth'',linethick)';
l{end+1}='    plot3([-x -x], [-y y], [maxi maxi],''color'',linecolor,''LineWidth'',linethick)';
l{end+1}='    plot3([x x], [-y y], [maxi maxi],''color'',linecolor,''LineWidth'',linethick)';
l{end+1}='    hold off';
l{end+1}='    %%set(gca,''XTick'',[-1 -0.5 -0.3 0 0.3 0.5 1])';
l{end+1}='    %%set(gca,''YTick'',[-1 -0.5 -0.3 0 0.3 0.5 1])';
l{end+1}='    xlabel(''Horizontal divergence [deg]'');';
l{end+1}='    ylabel(''Vertical divergence [deg]'');';
l{end+1}='    if i==MaxIndex';
l{end+1}='        title(''2d div (wavelength interval)'')';
l{end+1}='    else';
l{end+1}='        title([''2d div Lambda='' num2str(wavecenters(i)) ''A''])';
l{end+1}='    end';
l{end+1}='';
l{end+1}='';
l{end+1}='    subplot(MaxIndex,2,(i-1)*2+2)';
l{end+1}='    plot(PSD2D)';
l{end+1}='    box';
l{end+1}='    view([-0.5 90]);';
l{end+1}='    %%axis([-1 1 -1 1])';
l{end+1}=['   axis([' num2str(-0.5*psd2d*demandscan.Hsize) ' ' num2str(0.5*psd2d*demandscan.Hsize) ' ' num2str(-0.5*psd2d*demandscan.Vsize) ' ' num2str(0.5*psd2d*demandscan.Vsize) '])'];
l{end+1}='    maxi=max(PSD2D)+100;';
l{end+1}='    hold on';
l{end+1}=['    x=' num2str(demandscan.Hsize*0.5) ';'];
l{end+1}=['    y=' num2str(demandscan.Vsize*0.5) ';'];
l{end+1}='    plot3([-x x], [-y -y], [maxi maxi],''color'',linecolor,''LineWidth'',linethick)';
l{end+1}='    plot3([-x x], [y y], [maxi maxi],''color'',linecolor,''LineWidth'',linethick)';
l{end+1}='    plot3([-x -x], [-y y], [maxi maxi],''color'',linecolor,''LineWidth'',linethick)';
l{end+1}='    plot3([x x], [-y y], [maxi maxi],''color'',linecolor,''LineWidth'',linethick)';
l{end+1}='    hold off';
l{end+1}='    %%set(gca,''XTick'',[-1 -0.5 -0.3 0 0.3 0.5 1])';
l{end+1}='    %%set(gca,''YTick'',[-1 -0.5 -0.3 0 0.3 0.5 1])';
l{end+1}='    xlabel(''Horizontal position [cm]'');';
l{end+1}='    ylabel(''Vertical position [cm]'');';
l{end+1}='    if i==MaxIndex';
l{end+1}='        title(''2d psd (wavelength interval)'')';
l{end+1}='    else';
l{end+1}='        title([''2d psd Lambda='' num2str(wavecenters(i)) ''A''])';
l{end+1}='    end';
l{end+1}='   ';
l{end+1}='    ';
l{end+1}='    figure(Facceptance_ess)';
l{end+1}='    subplot(MaxIndex,2,(i-1)*2+1)';
l{end+1}='    plot(HACCP)';
l{end+1}='    box';
l{end+1}='    view([-0.5 90]);';
l{end+1}='    %%axis([-1 1 -1 1])';
l{end+1}=['   axis([' num2str(-0.5*apsd*demandscan.Hsize/100) ' ' num2str(0.5*apsd*demandscan.Hsize/100) ' ' num2str(-adiv*demandscan.Hdiv) ' ' num2str(adiv*demandscan.Hdiv) '])'];
l{end+1}='    maxi=max(HACCP)+100;';
l{end+1}='    hold on';
l{end+1}=['    x=' num2str(demandscan.Hsize*0.01*0.5) ';'];
l{end+1}=['    y=' num2str(demandscan.Hdiv) ';'];
l{end+1}='    plot3([-x x], [-y -y], [maxi maxi],''color'',linecolor,''LineWidth'',linethick)';
l{end+1}='    plot3([-x x], [y y], [maxi maxi],''color'',linecolor,''LineWidth'',linethick)';
l{end+1}='    plot3([-x -x], [-y y], [maxi maxi],''color'',linecolor,''LineWidth'',linethick)';
l{end+1}='    plot3([x x], [-y y], [maxi maxi],''color'',linecolor,''LineWidth'',linethick)';
l{end+1}='    hold off';
l{end+1}='    %%set(gca,''XTick'',[-1 -0.5 -0.3 0 0.3 0.5 1])';
l{end+1}='    %%set(gca,''YTick'',[-1 -0.5 -0.3 0 0.3 0.5 1])';
l{end+1}='    xlabel(''Horizontal position [m]'')';
l{end+1}='    ylabel(''Horizontal divergence [deg]'')';
l{end+1}='    if i==MaxIndex';
l{end+1}='    title(''Horizontal acceptance diagram (wavelength interval)'')';
l{end+1}='    else';
l{end+1}='    title([''Horizontal accep. dia. Wavelength='' num2str(wavecenters(i)) ''A'']);';
l{end+1}='    end    ';
l{end+1}='';
l{end+1}='    subplot(MaxIndex,2,(i-1)*2+2)';
l{end+1}='    plot(VACCP)';
l{end+1}='    box';
l{end+1}='    view([-0.5 90]);';
l{end+1}='    %%axis([-1 1 -1 1])';
l{end+1}=['   axis([' num2str(-0.5*apsd*demandscan.Vsize/100) ' ' num2str(0.5*apsd*demandscan.Vsize/100) ' ' num2str(-adiv*demandscan.Vdiv) ' ' num2str(adiv*demandscan.Vdiv) '])'];
l{end+1}='    maxi=max(VACCP)+100;';
l{end+1}='    hold on';
l{end+1}=['    x=' num2str(demandscan.Vsize*0.01*0.5) ';'];
l{end+1}=['    y=' num2str(demandscan.Vdiv) ';'];
l{end+1}='    plot3([-x x], [-y -y], [maxi maxi],''color'',linecolor,''LineWidth'',linethick)';
l{end+1}='    plot3([-x x], [y y], [maxi maxi],''color'',linecolor,''LineWidth'',linethick)';
l{end+1}='    plot3([-x -x], [-y y], [maxi maxi],''color'',linecolor,''LineWidth'',linethick)';
l{end+1}='    plot3([x x], [-y y], [maxi maxi],''color'',linecolor,''LineWidth'',linethick)';
l{end+1}='    hold off';
l{end+1}='    %%set(gca,''XTick'',[-1 -0.5 -0.3 0 0.3 0.5 1])';
l{end+1}='    %%set(gca,''YTick'',[-1 -0.5 -0.3 0 0.3 0.5 1])';
l{end+1}='    xlabel(''Vertical position [m]'')';
l{end+1}='    ylabel(''Vertical divergence [deg]'')';
l{end+1}='    if i==MaxIndex';
l{end+1}='    title(''Vertical acceptance diagram (wavelength interval)'')';
l{end+1}='    else';
l{end+1}='    title([''Vertical accep. dia. Wavelength='' num2str(wavecenters(i)) ''A'']);';
l{end+1}='    end';
l{end+1}='    ';
l{end+1}='    ';
l{end+1}='   if i==MaxIndex';
l{end+1}='      figure(background)';
%l{end+1}='      subplot(4,2,5:6)';
%l{end+1}='      axis off;';
%l{end+1}='      text(0.5,1,Project_name,''interpreter'',''none'')';
%l{end+1}='      text(0,0.83,[instrument_name '' - '' filename '' - ''  inputstring],''interpreter'',''none'')';
%l{end+1}='      text(0,0.66,[''Sample size: Horizontal='' p.sizeX ''m, Vertical='' p.sizeY ''m, Divergence requirement: Horizontal='' p.divreq_x ''deg, Vertical= '' p.divreq_y ''deg, Sample distance ='' p.sample_dist ''m''],''interpreter'',''none'');';
%l{end+1}='      text(0,0.5,[''Intensity on sample of 100 emitted = '' num2str(LAMBDA.Data.values(1)) '' (no divergence limit) '' num2str(LAMBDA_B_RAW) '' (width divergence limits)'' ],''interpreter'',''none'');';
%l{end+1}='      text(0,0.33,[''Intensity near sample of 100 emitted = '' num2str(PSD2D.Data.values(1)) '' (no divergence limit)'' ],''interpreter'',''none'');';
l{end+1}='      subplot(4,2,5:8)';
l{end+1}='      plot(LAMBDA_RAW_monitor/guide_end_lambda,''k'');';
l{end+1}='      hold on';
l{end+1}='      plot(LAMBDA_B_RAW_monitor/guide_end_lambda,''b'');';
l{end+1}='      hold off';
l{end+1}='      box on';
l{end+1}='      title(''Signal compared to total guide output'') ';
l{end+1}='      xlabel(''Wavelength [AA]'')';
l{end+1}='      ylabel(''Signal fraction [Unitless]'') ';
l{end+1}='    end';
l{end+1}='    ';
l{end+1}='end';
l{end+1}='';
%l{end+1}='if reflogic';
%l{end+1}='flux=monitor_fom(10).Data.values(1)/monitor_fom_ref(10).Data.values(1);';
%l{end+1}='fluxtext='' BT '';';
%l{end+1}='else';
l{end+1}='flux=monitor_fom(11).Data.values(1);';
l{end+1}='fluxtext='' Absolute '';';
%l{end+1}='end';
l{end+1}='';
%l{end+1}='fid = fopen([cpath ''/master_record-analyzed'' scanname ''.txt''],''a'');';
%l{end+1}='fprintf(fid,[num2str(flux) fluxtext ''= '' filename '' - '' inputstring ''\\n''])'; % Add a value describing the flux and BT
%l{end+1}='fclose(fid);';
l{end+1}='';
l{end+1}='set(Foverall_ess, ''paperpositionmode'', ''auto'');';
l{end+1}='set(Fposdiv_ess, ''paperpositionmode'', ''auto'');';
l{end+1}='set(Facceptance_ess, ''paperpositionmode'', ''auto'');';
l{end+1}='set(background, ''paperpositionmode'', ''auto'');';
l{end+1}='';
%l{end+1}='if reflogic';
%l{end+1}='hold on';
%l{end+1}='figure(Foverall)';
%l{end+1}='subplot(4,2,1:2)';
%l{end+1}='hold on';
%l{end+1}='plot(monitor_ALLW_degraded(10)/monitor_ALLW_ref(10),''r'')';
%l{end+1}='plot(monitor_ALLW(10)/monitor_ALLW_ref(10),''k'')';
%l{end+1}='title([''Lambda dependence, + is wavelengths for 1d graphs. I='' num2str(LAMBDA_B_RAW)])';
%l{end+1}='hold off';
%l{end+1}='else';

%l{end+1}='end';
l{end+1}='';
l{end+1}='print(Foverall_ess,''-dpng'',''-r300'',[cpath ''/'' filename ''_overall_ess.png''])';
l{end+1}='print(Fposdiv_ess,''-dpng'',''-r300'',[cpath ''/'' filename ''_posdiv_ess.png''])';
l{end+1}='print(Facceptance_ess,''-dpng'',''-r300'',[cpath ''/'' filename ''_acceptance_ess.png''])';
l{end+1}='print(background,''-dpng'',''-r300'',[cpath ''/'' filename ''_background.png''])';
l{end+1}='end';
l{end+1}='';
l{end+1}='if exist(''visualizer.m'')';
l{end+1}='  close all;';
l{end+1}='  visualizer(filename);';
l{end+1}='end';
end



analysis=l;
% One could make some plots of the resulting monitors.
% plot1: Large subplot: Brilliance as lambda and 1d snapshots for PSD/DIV
% plot2: 2D Div for each snapshot and 2D Psd for each snapshot and overall
% plot4: 2D Acceptance for each snapshot and overall


% This makes it so that the run script will not plot the data, even on the
% home computer. This is done to avoid the cluster running code it is not
% designed to perform which could trick seg faults by accident.
%for i=1:length(l)
for i=1:analysis_start
    iFitStr=[iFitStr l{i} '\n'];
end
clear l;


fid = fopen(['./' Project_name '/' filename '/' filename scanname '_ifit.m'], 'w');
write=iFitStr;
fprintf(fid,write);
fclose(fid);


analysis_string=['load(''' filename scanname '_all.mat'')\ncpath=pwd;\n'];
for i=analysis_start+1:length(analysis)
    analysis_string=[analysis_string analysis{i} '\n'];
end
clear analysis;

[status,message,messageid]=mkdir(['./' Project_name '/output']);
[status,message,messageid]=mkdir(['./' Project_name '/output/brill_ref']);
[status,message,messageid]=mkdir(['./' Project_name '/output/analysis']);

fid = fopen(['./' Project_name '/output/analysis/' filename scanname '_ifit_analyse.m'], 'w');
write=iFitStr;
fprintf(fid,analysis_string);
fclose(fid);

% move the visualizer script to the output/analysis folder if it is not
% there allready.
not_found = 1;
namelist = dir(['./' Project_name '/output/analysis']);
for ii = 1:length(namelist)
    if strcmp(namelist(ii).name,'visualizer.m')
        not_found = 0;
    end
end
if not_found
[status,message,messageid]=copyfile([options_general.guide_bot_path '/guide_bot_source/visualizer.m'],['./' Project_name '/output/analysis/.']);
end

% Write batch script for the cluster

if strcmp(options_general.cluster,'PSI')
l{1}='#!/bin/bash -l';
else
l{1}='#!/bin/bash';
end
l{end+1}='#';
l{end+1}='# Made by mcstas_bot';
l{end+1}='';
l{end+1}=['#SBATCH --job-name=' filename scanname];
l{end+1}=['#SBATCH --error=err_' filename scanname '.txt'];
l{end+1}=['#SBATCH --output=out_' filename scanname '.txt'];
l{end+1}='#SBATCH --nodes=1-1';
if strcmp(options_general.queue,'verylong')
l{end+1}='#SBATCH --exclude r3n9b3';
elseif strcmp(options_general.queue,'long') && strcmp(options_general.cluster,'ESSS')
l{end+1}='#SBATCH --exclude r1n16,r2n18,r1n4,r1n6,r1n15,r1n38,r2n20';    
end
l{end+1}=['#SBATCH --partition ' options_general.queue];    
l{end+1}=['#SBATCH --time ' options_general.time];
l{end+1}='# the --exclusive is needed when running OpenMPI';
l{end+1}='# it will allocate 1x12 core per node';
l{end+1}='#SBATCH --exclusive';
l{end+1}='';
if strcmp(options_general.queue,'verylong')
l{end+1}='NUMCORES=`echo "$SLURM_NNODES 16 * p "| dc`';    
else
l{end+1}='NUMCORES=`echo "$SLURM_NNODES 12 * p "| dc`';
end
l{end+1}='echo $NUMCORES > NUMCORES.DAT';
l{end+1}='';
for jj = 1:length(options_general.modules_node)
l{end+1}=['module load ' options_general.modules_node{jj}]; 
end
l{end+1}='';
l{end+1}='# Gives a set of matlab commands to iFit. NUMCORES.DAT read from matlab script';
l{end+1}=['cat run_' filename scanname '_ifit.m | ' options_general.cluster_path 'run_ifit.sh ' options_general.cluster_path 'MATLAB_Compiler_Runtime/v716/'];

BatchStr='';
for i=1:length(l)
  BatchStr=[BatchStr l{i} '\n'];
end
clear l

fid = fopen(['./' Project_name '/' filename '/' filename scanname '.batch'], 'w');
write=BatchStr;
fprintf(fid,write);
fclose(fid);

fid = fopen(['./' Project_name '/' filename '/run_' filename scanname '_ifit.m'], 'w');
fprintf(fid,['run ' filename scanname '_ifit.m \nexit\n']);
fclose(fid);

%demandscan
%reqscan

% Add note of succesfull run to the project file
fid = fopen(['./' Project_name '/master_record-writen' scanname '.txt'],'a');
fprintf(fid,[filename ' ' input '\n']);
if scan.mode==1
   if scan.dimension == 1
      if scan.identifiers(1) < scan.demands_to_requirements+0.5
            string = [scan.names{1} '=' num2str(demandscan.(scan.names{1})) '\n'];
      else
            string = [scan.names{1} '=' num2str(reqscan.(scan.names{1})) '\n'];
      end
   else
      if scan.identifiers(1) < scan.demands_to_requirements+0.5
            string = [scan.names{1} '=' num2str(demandscan.(scan.names{1})) '\n'];
      else
            string = [scan.names{1} '=' num2str(reqscan.(scan.names{1})) '\n'];
      end
      if scan.identifiers(2) < scan.demands_to_requirements+0.5
            string = [string scan.names{2} '=' num2str(demandscan.(scan.names{2})) '\n'];
      else
            string = [string scan.names{2} '=' num2str(reqscan.(scan.names{2})) '\n'];
      end
   end
   fprintf(fid,string);
end
fclose(fid);

% For some versions of guide elliptical, some nested functions flags are
% needed. They need to be forced in this version.


if isfield(options_general,'enable_nested_compile')
    if options_general.enable_nested_compile == 1
        enable_nested_compile = 1;
    else
        enable_nested_compile = 0;
    end
else
    enable_nested_compile = 0;
end

% Lunch script for the cluster

fid = fopen(['./' Project_name '/launch_all.sh' ],'a');
fprintf(fid,['cd ' filename '\nsbatch ./' filename scanname '.batch\ncd .. \n']);
fclose(fid);

end %ending for loop over scani for the scan durration
end % end the scanj loop
% The following scripts only have to be done once, not for every scan point

% Quick sh script for compiling on the cluster

l{1}='#!/bin/bash';
l{end+1}='# Made by mcstas_bot';
l{end+1}='# Will compile the different .instr files ';
l{end+1}='';
if enable_nested_compile
l{end+1}='export MCSTAS_CFLAGS="-g -lm -O2 -fnested-functions"';
end
if options_general.gravity == 1
l{end+1}=['mcrun-2.0-py -n 0 -g -c --mpi=1 ' filename '_optimize.instr'];
l{end+1}=['mcrun-2.0-py -n 0 -g -c --mpi=1 ' filename '_optimize_coating.instr'];
l{end+1}=['mcrun-2.0-py -n 0 -g -c --mpi=1 ' filename '_analyze.instr'];
l{end+1}=['mcrun-2.0-py -n 0 -g -c --mpi=1 ' filename '_analyze_coating.instr'];
l{end+1}=['mcrun-2.0-py -n 0 -g -c --mpi=1 ' filename '_analyze_ess.instr'];
l{end+1}=['mcrun-2.0-py -n 0 -g -c --mpi=1 ' filename '_analyze_ess_coating.instr'];
else
l{end+1}=['mcrun-2.0-py -n 0 -c --mpi=1 ' filename '_optimize.instr'];
l{end+1}=['mcrun-2.0-py -n 0 -c --mpi=1 ' filename '_optimize_coating.instr'];
l{end+1}=['mcrun-2.0-py -n 0 -c --mpi=1 ' filename '_analyze.instr'];
l{end+1}=['mcrun-2.0-py -n 0 -c --mpi=1 ' filename '_analyze_coating.instr'];
l{end+1}=['mcrun-2.0-py -n 0 -c --mpi=1 ' filename '_analyze_ess.instr'];
l{end+1}=['mcrun-2.0-py -n 0 -c --mpi=1 ' filename '_analyze_ess_coating.instr'];    
end
if strcmp(source_component,'ESS_pancake')
l{end+1}=['gcc -lm write_cold.c -o write_cold'];
l{end+1}=['gcc -lm write_thermal.c -o write_thermal'];
l{end+1}=['gcc -lm write_bispectral.c -o write_bispectral'];
end

CompileStr='';
for i=1:length(l)
  CompileStr=[CompileStr l{i} '\n'];
end
clear l

fid = fopen(['./' Project_name '/' filename '/compile_' filename '_py.sh'], 'w');
write=CompileStr;
fprintf(fid,write);
fclose(fid);
% make the file executable
if ~ispc
unix(['chmod 744 ./' Project_name '/' filename '/compile_' filename '_py.sh']);
else
fid = fopen(['./' Project_name '/windows_support.sh'],'a');
fprintf(fid,['chmod 744 ./' filename '/compile_' filename '_py.sh\n']);
fclose(fid);    
end


l{1}='#!/bin/bash';
l{end+1}='# Made by mcstas_bot';
l{end+1}='# Will compile the different .instr files ';
l{end+1}='';
if enable_nested_compile
l{end+1}='export MCSTAS_CFLAGS="-g -lm -O2 -fnested-functions"';
end
if options_general.gravity == 1
l{end+1}=['mcrun-2.0 -n 0 -g -c --mpi ./' filename '_optimize.instr'];
l{end+1}=['mcrun-2.0 -n 0 -g -c --mpi ./' filename '_optimize_coating.instr'];
l{end+1}=['mcrun-2.0 -n 0 -g -c --mpi ./' filename '_analyze.instr'];
l{end+1}=['mcrun-2.0 -n 0 -g -c --mpi ./' filename '_analyze_coating.instr'];
l{end+1}=['mcrun-2.0 -n 0 -g -c --mpi ./' filename '_analyze_ess.instr'];
l{end+1}=['mcrun-2.0 -n 0 -g -c --mpi ./' filename '_analyze_ess_coating.instr'];
else
l{end+1}=['mcrun-2.0 -n 0 -c --mpi ./' filename '_optimize.instr'];
l{end+1}=['mcrun-2.0 -n 0 -c --mpi ./' filename '_optimize_coating.instr'];
l{end+1}=['mcrun-2.0 -n 0 -c --mpi ./' filename '_analyze.instr'];
l{end+1}=['mcrun-2.0 -n 0 -c --mpi ./' filename '_analyze_coating.instr'];
l{end+1}=['mcrun-2.0 -n 0 -c --mpi ./' filename '_analyze_ess.instr'];
l{end+1}=['mcrun-2.0 -n 0 -c --mpi ./' filename '_analyze_ess_coating.instr'];    
end

if strcmp(source_component,'ESS_pancake')
l{end+1}=['gcc -lm write_cold.c -o write_cold'];
l{end+1}=['gcc -lm write_thermal.c -o write_thermal'];
l{end+1}=['gcc -lm write_bispectral.c -o write_bispectral'];
end


CompileStr='';
for i=1:length(l)
  CompileStr=[CompileStr l{i} '\n'];
end
clear l

fid = fopen(['./' Project_name '/' filename '/compile_' filename '.sh'], 'w');
write=CompileStr;
fprintf(fid,write);
fclose(fid);
% make the file executable
if ~ispc
unix(['chmod 744 ./' Project_name '/' filename '/compile_' filename '.sh']);
else
fid = fopen(['./' Project_name '/windows_support.sh'],'a');
fprintf(fid,['chmod 744 ./' filename '/compile_' filename '.sh\n']);
fclose(fid);    
end

% Add a file in the directory telling what the input was

fid = fopen(['./' Project_name '/' filename '/input.txt'],'w');
fprintf(fid,[input '\n']);
fclose(fid);

% Add it to the compile script (needs to be initialized before this

fid = fopen(['./' Project_name '/compile_all.sh'],'a');
fprintf(fid,['cd ' filename '\n./compile_' filename '.sh\ncd ..\n']);
fclose(fid);

fid = fopen(['./' Project_name '/compile_all_py.sh' ],'a');
fprintf(fid,['cd ' filename '\n./compile_' filename '_py.sh\ncd ..\n']);
fclose(fid);

disp(' ')

end 
% The output also needs to consists of a list of which of the input
% parameters should be optimized to keep the minimalist concept, and
% ideally also limits for which the optimization would be sensible. 


% How will this fit in overall?
    % 1: loop over different input strings, either list or generated
    % 2: make McStas files for each (optimization and single analyse run)
    % 3: make sh script for compiling all these McStas simulations
    % 4: make ifit optimization files for each (cluster or single pc)
    % 5: make batch script for each (optimization & single) (cluster)
    % 6: make sh script for launching all these scripts (cluster)
    % 7: make sh script for halting all launched batchs (error!)
    % 8: make file describing what results are expected and where.
    
% Analyse ifit script
    % 1: read file from last step above
    % 2: check which simulations that are done
    % 3: read the ones that are done
    % 4: produce figures describing the performance for each solution
    % 5: produce mcdisplay that draws the guide if possible or latex?
    % 6: produce tex file with simple reports that include the figures
    % 7: write a list of the solutions in order of best performance


% Note:
    % It is actually possible to produce graphs on the cluster
    % It is not possible to compile latex on the cluster so far
    % Adding a collection of graphs in a large subplot should be done


