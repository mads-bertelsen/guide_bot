% This function performs two tasks.
% 1)  It allows a user to define the arguement of a guide_bot component as
%     a structure. The function parses the stucture into the long hand
%     specification that was originally implemented. For instance, you can
%     specify:

%     StraightStruct.minlength = 3;
%     StraightStruct.maxlength = 10;
%     StraightStruct.start = 20.5;

%     input{1} = S E S(StraightStruct)

%     This function will then allow guide_bot to read it nativly as:

%     input{1} = S E S(minlength=3,maxlength=10,start=20.5)

%     This functionality was primarily implemented for a monochromator since it
%     has so many values that must be defined in advance.


% 2)  It allows a developer to add non-user defined arguements to a
%     guide_bot component


% INPUTS:
%    part: A cell generated in mcstas_bot containing the seperate guide_bot component constructors that make up input.
%          In the above example:
%          part{1} = 'S', part{2} = 'E', and part{3} = 'S(StraightStruct)'
%
%    demands: The user defined demands structure.
%
%
% OUTPUTS:
%   parsed_part: The parsed version of part.
%                 In the above example:
%                 parsed_part{1} = 'S', parsed_part{2} = 'E', parsed_part{3} = 'S(minlength=3,maxlength=10,start=20.5)'


% Leland Harriger, NCNR September 2016

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% BEGINING OF MAIN FUNCTION %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function parsed_part = parse_gbstruct(part,demands)

% This first loop handles parsing into long hand any user defined structures.
% Loop through each element in part{i}
for par = 1:length(part)
    parsed_part{par} = Struct2Long(part{par});
end

% This second loop modifies, removes, or adds arguements not specified by users.
for par = 1:length(parsed_part)
    CompLetter = getfield(Long2Struct(parsed_part{par}),'CompLetter');
    switch(CompLetter)
        % For a Monochromator you must add 'length', and you may need to add 'start'.
        case 'M'
            index = length(parsed_part) - par + 1;
            parsed_part{par} = M_Add(parsed_part{par}, demands, index);
        otherwise
            % Do Nothing
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% END OF MAIN FUNCTION %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% BEGINING OF M_Add() SUPPORT FUNCTION %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Adds non-user defined arguements to the M Constructor.
function expanded_part = M_Add(Constructor, demands, index)

    % Convert long hand constructor to a structure
    StructArg = Long2Struct(Constructor);
    
    % First define length.
    monlength = StructArg.BladeN*StructArg.BladeW + (StructArg.BladeN - 1)*StructArg.BladeGap;
    if (isfield(StructArg, 'length') && StructArg.length ~= monlength)
        disp(['Monochromator length incorrectly defined, rewriting as length=' num2str(monlength)])
    end
    StructArg.length = monlength;
            
    % Now, if needed, define start.
    % start is required if M is the last component (index = 1). This is because the
    % length of M is always fixed. Thus, there is no way to 'stretch' it to
    % a distance other than demands.Dist + length1 from the sample.
    if index == 1 % add start
        monstart = demands.Mod_sample - demands.Dist - monlength;
        if (isfield(StructArg, 'start') && monstart ~= StructArg.start)
            disp(['Redefining the monochromator start position to be consistent with demands.Dist (start=' num2str(monstart) ').'])
        end
    StructArg.start = monstart;   
    end               
    
    
    % Convert back into long hand
    expanded_part = Struct2Long(StructArg);           
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% END OF M_Add() SUPPORT FUNCTION %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% BEGINING OF STRUCT2LONG() SUPPORT FUNCTION %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Converts a constructors arguement from a structure to the guide_bot long hand.
function LongArg = Struct2Long(Constructor)
    if isstruct(Constructor) % The Constructor is fully specified as a single structure. (Code defined)
        if length(fieldnames(Constructor)) > 2
            CompLetter = Constructor.CompLetter;
            Constructor = rmfield(Constructor, 'CompLetter');
            ArgFields = fieldnames(Constructor);
            parsed = [CompLetter '('];
            for i=1:numel(ArgFields)
                parsed = [parsed [ArgFields{i} '=' num2str(Constructor.(ArgFields{i})) ',']];
            end
            % removed extra comma, and add closing parenthesis
            LongArg = [parsed(1:end-1), ')'];
            
        else
            LongArg = Constructor.CompLetter;
        end        
    else % The Constructor is defined as a component letter with the structure as an arguement. (User defined)
        % find ( in part{par} string
        leftp=-1;rightp=-1;
        for i=1:length(Constructor)
            if strcmp(Constructor(i),'(');
                leftp=i;
            end
        end
        
        % find ) in part{par} string
        for i=1:length(Constructor)
            if strcmp(Constructor(i),')');
                rightp=i;
            end
        end
        
        % Determine the arguement of part{par}.
        if (leftp==-1 || rightp==-1)    %Case: no additional commands
            CompArg = '';
        else                            % Case: Additional commands given
            CompLetter = Constructor(1:leftp-1);
            CompArg = Constructor(leftp+1:rightp-1);
        end
        
        
        % Determine if the component arguement is a structure.
        struct_check = 0;
        try
            struct_check = isstruct((evalin('base',CompArg)));
        end
        
        
        % If it is a structure then parse it.
        if struct_check == 1
            % Get the structure from the base workspace
            CompArg = evalin('base',CompArg);
            ArgFields = fieldnames(CompArg);
            
            parsed = [CompLetter '('];
            for i=1:numel(ArgFields)
                parsed = [parsed [ArgFields{i} '=' num2str(CompArg.(ArgFields{i})) ',']];
            end
            % removed extra comma, and add closing parenthesis
            LongArg = [parsed(1:end-1), ')'];
        else
            LongArg = Constructor;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% END OF STRUCT2LONG() SUPPORT FUNCTION %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% BEGINING OF LONG2STRUCT() SUPPORT FUNCTION %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function StructArg = Long2Struct(Constructor)
    
    % find ( in part{par} string
    leftp=-1;rightp=-1;
    for i=1:length(Constructor)
        if strcmp(Constructor(i),'('); 
            leftp=i;    
        end
    end

    % find ) in part{par} string
    for i=1:length(Constructor)
        if strcmp(Constructor(i),')'); 
            rightp=i;    
        end
    end
    
    % Determine the arguement of part{par}.
    if (leftp==-1 || rightp==-1)    %Case: no additional commands
        CompLetter = Constructor;
        CompArg = '';
    else                            % Case: Additional commands given
        CompLetter = Constructor(1:leftp-1);
        CompArg = Constructor(leftp+1:rightp-1);
    end
    
    % Store letter and full arguement in output structure (StructArg)
    StructArg.CompLetter = CompLetter;
    
    % Split the arguement into individual commands
    commap=-1;j=0;
        for i=1:length(CompArg)
            if strcmp(CompArg(i),','); 
                j=j+1;
                commap(j)=i;    
            end
        end
                
        if (commap ~= -1) % If there are commas, there are several commands
            for i=1:length(commap)+1
                    if i==1
                        command = CompArg(1:(commap(i)-1));
                    elseif (i==length(commap)+1)
                        command = CompArg((commap(i-1)+1):end);
                    else
                        command = CompArg((commap(i-1)+1):(commap(i)-1));
                    end
                    try % Assume command is numeric
                        eval(['StructArg.' command ';']);
                    catch % If not numeric then store command as a string;
                        for n = 1:length(command)
                            if strcmp(command(n),'=');
                                equals=n;
                            end
                        end
                        eval(['StructArg.' command(1:equals-1) '=''' command(equals+1:end) ''';']);                                
                    end
                            
                
            end
        else % No commas, means just one additional command was given
            command = CompArg;
            try % Assume command is numeric
                eval(['StructArg.' command ';']);
            catch % If not numeric then store command as a string;
                for n = 1:length(command)
                    if strcmp(command(n),'=');
                        equals=n;
                    end
                end
                try
                    eval(['StructArg.' command(1:equals-1) '=''' command(equals+1:end) ''';']);
                end
            end
        end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% END OF LONG2STRUCT() SUPPORT FUNCTION %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%









    
        
        
        
        
        
        
        
        
        
        
        
    