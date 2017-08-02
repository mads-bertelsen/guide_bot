function McStasStr=guide_writer_length(McStasStr,index,last,demands,requirements,globalinfo,lengthfracpar)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

num=num2str(index);
numM1=num2str(index-1);
l{1}='';
strings = '';
if globalinfo.fixedlength(index)==1 % The module have a fixed length
    % Do nothing, the value length num is already defined and have
    % the right input as a mcstas input. (Done i main script as it is global)
else % The module does not have a fixed length
    if (globalinfo.fixedstart(index)==1) % The module have a fixed start point
        ldec=length(McStasStr.declare);
        McStasStr.declare{ldec+1}=['length' num];
        strings = ['length' num ' = endPoint' num ' - fixedStart' num ';'];
        McStasStr.initialize=[strings '\n\n' McStasStr.initialize];
    else % The module does not have a fixed start point
        if sum(globalinfo.fixedstart(index+1:end))==0 % There is no fixed start left
            nextfp=last;
            nextfpstr='guide_start';
            % val for index==1
            if sum(ismember(McStasStr.input,'guide_start'))>0.5
                guide_start_index=find(ismember(McStasStr.input,'guide_start'),1);
                if McStasStr.optimize(guide_start_index)==1
                    val=McStasStr.optimvals.max(guide_start_index);
                else
                    val=McStasStr.inputvalue(guide_start_index);
                end
            end
        else % There is another fixed start module left
            % fixed start declared as mcstas input in main script
            nextfp=find(globalinfo.fixedstart(index:last),1)+index-1;
            nextfpstr=['fixedStart' num2str(nextfp)];
            val=globalinfo.fixedstartpoint(nextfp);
        end
        
        if sum(globalinfo.minstart(index+1:end))==0
            % There is no maxstart left
            nextminfp=last; %Will this hurt anything?
            % or signal
            nextminfpexist=0;
            nextminfpstr='ERROR'; % should never be needed
        else
            nextminfp=find(globalinfo.minstart(index+1:last),1)+index;
            nextminfpexist=1;
            nextminfpstr=['minStart' num2str(nextminfp)];
        end
        
        
        [nextminorfp,fpormin]=min([nextfp nextminfp]);
        if nextminorfp == last
            nextminorfpstr=['guide_start'];
            valminorfp=val;
        elseif fpormin==1
            valminorfp=val;
            nextminorfpstr=['fixedStart' num2str(nextfp)];
        else
            valminorfp=globalinfo.minstartpoint(nextminorfp);
            nextminorfpstr=['minStart' num2str(nextminfp)];
        end
        
        if sum(globalinfo.maxstart(index+1:end))==0
            % There is no maxstart left
            nextmaxfp=last; % will this screw anything up?
            % or signal
            %nextmaxfpexist=0;
        else
            nextmaxfp=find(globalinfo.maxstart(index+1:last),1)+index;
            %nextmaxfpexist=1;
            nextmaxfpstr=['maxStart' num2str(nextmaxfp)];
        end
        
        
        if (index==1)
            % Only here if the module closest to the sample:
            % Does not have a fixed lenght
            % Does not have a fixed starting point
            % In these cases it is calculated as the rest
            if index==last % Module closest to the moderator and sample!
                % TODO
                % min/max/fixed start of the last module should be handled
                % through the guide_start parameter in main
                % min/max/fixed length of the last module likewise
                ldec=length(McStasStr.declare);
                McStasStr.declare{ldec+1}=['length' num];
                strings=['length' num ' = Mod_sample - sample_dist - guide_start;'];
                McStasStr.initialize=[strings '\n\n' McStasStr.initialize];
            else % First module in a series of more than one
                % if globalinfo.maxstart = 1, check min distance needed
                
                
                minimumneeded=0;
                for i = find(globalinfo.maxstart)
                    if i>index && sum(max(globalinfo.fixedlength(index+1:i),globalinfo.maxlength(index+1:i)))>i-index-0.5
                        minimumneeded = max([minimumneeded demands.Mod_sample-demands.Dist-globalinfo.maxstartpoint(i)-sum(globalinfo.fixedlengthdata(index+1:i)+globalinfo.maxlengthdata(index+1:i))]);
                    end
                end
                
                % Experimental code!
                for i = find(globalinfo.fixedstart)
                    if i>index && sum(max(globalinfo.fixedlength(index+1:i),globalinfo.maxlength(index+1:i)))>i-index-0.5
                        minimumneeded = max([minimumneeded demands.Mod_sample-demands.Dist-globalinfo.fixedstartpoint(i)-sum(globalinfo.fixedlengthdata(index+1:i)+globalinfo.maxlengthdata(index+1:i))]);
                    end
                end
                
                
                
                %                if nextmaxfp < nextminorfp && sum(globalinfo.fixedlengthdata(index+1:nextmaxfp)+globalinfo.maxlengthdata(index+1:nextmaxfp))<demands.Mod_sample-demands.Dist-globalinfo.maxstartpoint(nextmaxfp)% check if a min length1 is needed
                %                            minimumneeded=demands.Mod_sample-demands.Dist-globalinfo.maxstartpoint(nextmaxfp)-sum(globalinfo.fixedlengthdata+globalinfo.maxlengthdata);
                %                            if globalinfo.maxstart(index)==1 % this should be able to happen even if the nextmaxfp is not smaller than the nextminorfp?
                %                                minimumneeded=max([minimumneeded demands.Mod_sample-demands.Dist-globalinfo.maxstartpoint(index)]);
                %                            end
                %                else
                %                            % no maxfixp to think about
                %                            minimumneeded=0;
                %                end
                %if sum(max(globalinfo.fixedlength(index:nextminorfp),globalinfo.minlength(index:nextminorfp)))>nextminorfp-0.5
                %BUG? fixed back and forth a couple of times.. bugs exists
                % with each it seems?
                if sum(max(globalinfo.fixedlength(index+1:nextminorfp),globalinfo.minlength(index+1:nextminorfp)))>nextminorfp-1.5
                    % Every module between now and nextminorfp have a fixed
                    % or min length.
                    
                    % BUG?
                    if sum(globalinfo.fixedlength(index+1:nextfp))>nextfp-1.5
                    %if sum(globalinfo.fixedlength(index+1:nextfp))>nextminorfp-1.5
                    %if sum(globalinfo.fixedlength(index+1:nextminorfp))>nextminorfp-1.5    
                        % ONLY fixed lengths between now and nextminorfp
                        if nextmaxfp < nextminorfp
                            % need to take care of a maxfp
                            % not able to accomodate that!
                            disp('ERROR, the chosen startpoints/lengths are not possible')
                        else
                            % no maxfp to worry about
                            % do as always: calculate length1
                            %save(['debug' num])
                            ldec=length(McStasStr.declare);
                            McStasStr.declare{ldec+1}=['length' num];
                            % BUG?
                            strings = ['length' num ' = endPoint' num ' - ' nextfpstr ];
                            %strings = ['length' num ' = endPoint' num ' - ' nextminorfpstr ];
                            % BUG?
                            for i=find(globalinfo.fixedlength(index+1:nextfp))+index
                            %for i=find(globalinfo.fixedlength(index+1:nextminorfp))+index
                                strings = [ strings ' - length' num2str(i)];
                            end
                            strings = [ strings ';'];
                            McStasStr.initialize=[strings '\n\n' McStasStr.initialize];
                        end
                    else
                        % unbroken mix between fixed and min length between
                        % 1 and nextminorfp
                        % also need to check for maxlengths
                        
                        % calculate length1 range from the above info
                        linp=length(McStasStr.input);
                        McStasStr.input{linp+1}=['length' num];
                        McStasStr.optimize(linp+1)=1;
                        
                        %save(['debug' num '.mat'])
                        
                        %McStasStr.optimvals.max(linp+1)=demands.Mod_sample-demands.Dist-requirements.closest_element-last*0.1;
                        %McStasStr.optimvals.max(linp+1)=demands.Mod_sample-demands.Dist-val- ...
                        %    sum(globalinfo.fixedlength(index:nextfp).*globalinfo.fixedlengthdata(index:nextfp))-0.1*(nextfp-index);
                        McStasStr.optimvals.max(linp+1)=demands.Mod_sample-demands.Dist-valminorfp- ...
                            sum(globalinfo.fixedlengthdata(index+1:nextminorfp)+globalinfo.minlengthdata(index+1:nextminorfp))-0.1*(nextminorfp-index-1);
                        %sum(globalinfo.fixedlengthdata(index:nextminorfp)+globalinfo.minlengthdata(index:nextminorfp))-0.1*(nextminorfp-index-1);
                        % Was the above a bug?
                        if globalinfo.maxlength(index)
                            McStasStr.optimvals.max(linp+1)=min([McStasStr.optimvals.max(linp+1) globalinfo.maxlengthdata(index)]);
                        end
                        McStasStr.optimvals.min(linp+1)=max([globalinfo.minlengthdata(index) minimumneeded 0.1]);
                        %McStasStr.optimvals.guess(linp+1)=(McStasStr.optimvals.max(linp+1)-0.1)/(1+nextminorfp-index)+0.1;
                        McStasStr.optimvals.guess(linp+1)=mean([McStasStr.optimvals.min(linp+1) McStasStr.optimvals.max(linp+1)]);
                    end
                    
                    %elseif sum(globalinfo.fixedlength(index+1:nextfp))>nextfp-1.5
                else
                    % Here if there is a length between n and nextminorfp
                    % which is not restricted by a min or set value
                    %disp('debug message of the day')
                    
                    linp=length(McStasStr.input);
                    McStasStr.input{linp+1}=['length' num];
                    McStasStr.optimize(linp+1)=1;
                    
                    %save(['debug' num '.mat'])
                    
                    %McStasStr.optimvals.max(linp+1)=demands.Mod_sample-demands.Dist-requirements.closest_element-last*0.1;
                    %McStasStr.optimvals.max(linp+1)=demands.Mod_sample-demands.Dist-val- ...
                    %    sum(globalinfo.fixedlength(index:nextfp).*globalinfo.fixedlengthdata(index:nextfp))-0.1*(nextfp-index);
                    McStasStr.optimvals.max(linp+1)=demands.Mod_sample-demands.Dist-valminorfp- ...
                        sum(globalinfo.fixedlengthdata(index+1:nextminorfp)+globalinfo.minlengthdata(index+1:nextminorfp))-0.1*(nextminorfp-index-1);
                    %sum(globalinfo.fixedlengthdata(index:nextminorfp)+globalinfo.minlengthdata(index:nextminorfp))-0.1*(nextminorfp-index-1);
                    % Was the above a bug?
                    if globalinfo.maxlength(index)
                        McStasStr.optimvals.max(linp+1)=min([McStasStr.optimvals.max(linp+1) globalinfo.maxlengthdata(index)]);
                    end
                    McStasStr.optimvals.min(linp+1)=max([globalinfo.minlengthdata(index) minimumneeded 0.1]);
                    McStasStr.optimvals.guess(linp+1)=mean([McStasStr.optimvals.min(linp+1) McStasStr.optimvals.max(linp+1)]);
                    %McStasStr.optimvals.guess(linp+1)=(McStasStr.optimvals.max(linp+1)-0.1)/(1+nextminorfp-index)+0.1;
                    
                end
                %                if test==2
                %                % old code
                %                    % Above statement checks if there is a free length
                %                    % between the current module and the next fixed start
                %                    % If there is not, it can be calculated as follows:
                %                    ldec=length(McStasStr.declare);
                %                    McStasStr.declare{ldec+1}=['length' num];
                %                    strings = ['length' num ' = endPoint' num ' - ' nextfpstr ];
                %                    for i=find(globalinfo.fixedlength(index+1:nextfp))+index
                %                        strings = [ strings ' - length' num2str(i)];
                %                    end
                %                    strings = [ strings ';'];
                %                    McStasStr.initialize=[strings '\n\n' McStasStr.initialize];
                %                elseif sum(globalinfo.maxlength(index+1:nextfp))
                %                else % The length of the first module should be optimized:
                %                linp=length(McStasStr.input);
                %                McStasStr.input{linp+1}=['length' num];
                %                McStasStr.optimize(linp+1)=1;
                %                McStasStr.optimvals.min(linp+1)=0.1;
                %                %McStasStr.optimvals.max(linp+1)=demands.Mod_sample-demands.Dist-requirements.closest_element-last*0.1;
                %                McStasStr.optimvals.max(linp+1)=demands.Mod_sample-demands.Dist-val- ...
                %                    sum(globalinfo.fixedlength(index:nextfp).*globalinfo.fixedlengthdata(index:nextfp))-0.1*(nextfp-index);
                %                McStasStr.optimvals.guess(linp+1)=(McStasStr.optimvals.max(linp+1)-0.1)/(1+nextfp-index)+0.1;
                %                end
            end
        else
            % check if a minimum distance is needed, and adjust
            
            
            %save('debug')
            if sum(max(globalinfo.fixedlength(index+1:nextminorfp),globalinfo.minlength(index+1:nextminorfp)))>nextminorfp-index-0.5
                % There is an unbroken line between index and nextminorfp
                % consisting of minlength or fixedlength modules
                
                if sum(globalinfo.fixedlength(index+1:nextminorfp))>nextminorfp-index-0.5
                    % only fixed length between now and nextminorfp
                    % need to calculate length because there is no freedom
                    
                    %no free length left after this one
                    
                    % the variable containing the fixed point have been declared
                    % by the main script because it is global
                    McStasStr.declare{end+1}=['length' num];
                    strings = ['length' num ' = endPoint' num ' - ' nextminorfpstr ];
                    for i=find(globalinfo.fixedlength(index+1:nextfp))+index
                        strings = [ strings ' - length' num2str(i)];
                    end
                    strings = [ strings ';'];
                    McStasStr.initialize=[strings '\n\n' McStasStr.initialize];
                    
                    
                else
                    % mix of fixed and min length or only min, which reduces
                    % the amount of free length left
                    
                    % Right now i do not see how this is different from
                    % else, copy pasting code from there.
                    
                    % need to check for minstart1?
                    
                    
                    
                    % relativly normal
                    linp=length(McStasStr.input);
                    McStasStr.input{linp+1}=['lengthfrac' num];
                    McStasStr.optimize(linp+1)=1;
                    McStasStr.optimvals.min(linp+1)=lengthfracpar(1);
                    McStasStr.optimvals.max(linp+1)=lengthfracpar(2);
                    McStasStr.declare{end+1}=['length' num];
                    if length(lengthfracpar)>2.5
                        McStasStr.optimvals.guess(linp+1)=lengthfracpar(3);
                    else
                        McStasStr.optimvals.guess(linp+1)=1./(1+last-index); % needs to be updated
                    end
                    
                    if globalinfo.minstart(index)==1
                        nextminorfpstr=['minStart' num];
                    end
                    
                    if globalinfo.maxlength(index)<0.5
                        strings = ['length' num ' = minlengthneeded + (endPoint' num ' - ' nextminorfpstr ' - minlengthneeded'];
                        for i=find(globalinfo.fixedlength(index+1:nextfp))+index
                            strings = [ strings ' - length' num2str(i)];
                        end
                        for i=find(globalinfo.minlength(index+1:nextfp))+index
                            strings = [ strings ' - minlength' num2str(i)]; % Needs to be declared and set in main
                        end
                        strings = [ strings ')*lengthfrac' num ';'];
                        McStasStr.initialize=[strings '\n\n' McStasStr.initialize];
                    else
                        strings = ['length' num ' = minlengthneeded + (maxlength' num ' - minlengthneeded)*lengthfrac' num ';'];
                        if length(lengthfracpar)>2.5
                            McStasStr.optimvals.guess(linp+1)=lengthfracpar(3);
                        else
                            McStasStr.optimvals.guess(linp+1)=0.5;
                        end
                        McStasStr.initialize=[strings '\n\n' McStasStr.initialize];
                        % Newly added code to combat a specific bug
                        % might be symptom of a overall problem with the algorithm
                        strings = ['if (maxlength' num ' > endPoint' num ' - ' nextminorfpstr ') maxlength' num ' = endPoint' num ' - ' nextminorfpstr ';'];
                        McStasStr.initialize=[strings '\n' McStasStr.initialize];
                    end
                    
                    
                end
                
                
            else
                
                % relativly normal
                linp=length(McStasStr.input);
                McStasStr.input{linp+1}=['lengthfrac' num];
                McStasStr.optimize(linp+1)=1;
                if globalinfo.minlength(index)==1 % Testing
                McStasStr.optimvals.min(linp+1)=0;
                else
                McStasStr.optimvals.min(linp+1)=lengthfracpar(1);
                end

                McStasStr.optimvals.max(linp+1)=lengthfracpar(2);
                McStasStr.declare{end+1}=['length' num];
                if length(lengthfracpar)>2.5
                    McStasStr.optimvals.guess(linp+1)=lengthfracpar(3);
                else
                    McStasStr.optimvals.guess(linp+1)=1./(1+last-index); % needs to be updated
                end
                
                if globalinfo.minstart(index)==1
                     nextminorfpstr=['minStart' num];
                end
                
                if globalinfo.maxlength(index)<0.5
                    strings = ['length' num ' = minlengthneeded + (endPoint' num ' - ' nextminorfpstr ' - minlengthneeded'];
                    for i=find(globalinfo.fixedlength(index+1:nextfp))+index
                        strings = [ strings ' - length' num2str(i)];
                    end
                    for i=find(globalinfo.minlength(index+1:nextfp))+index
                        strings = [ strings ' - minlength' num2str(i)]; % Needs to be declared and set in main
                    end
                    strings = [ strings ')*lengthfrac' num ';'];
                    McStasStr.initialize=[strings '\n\n' McStasStr.initialize];
                else
                    
                    strings = ['length' num ' = minlengthneeded + (maxlength' num ' - minlengthneeded)*lengthfrac' num ';'];
                    McStasStr.optimvals.max(linp+1)=0.999;
                    if length(lengthfracpar)>2.5
                        if lengthfracpar(3)<0.05
                            McStasStr.optimvals.guess(linp+1)=0.5;
                        else
                            McStasStr.optimvals.guess(linp+1)=lengthfracpar(3);
                        end
                    else
                        McStasStr.optimvals.guess(linp+1)=0.5;
                    end
                    McStasStr.initialize=[strings '\n\n' McStasStr.initialize];
                    
                    % Newly added code to combat a specific bug
                    % might be symptom of a overall problem with the algorithm
                    strings = ['if (maxlength' num ' > endPoint' num ' - ' nextminorfpstr ') maxlength' num ' = endPoint' num ' - ' nextminorfpstr ';'];
                    McStasStr.initialize=[strings '\n' McStasStr.initialize];
                end
                
            end
            
            if sum(ismember(McStasStr.declare,'minlengthneeded'))<0.5 % not always needed
                McStasStr.declare{end+1}='minlengthneeded';
            end
            
            if sum(ismember(McStasStr.declare,'tmp_length'))<0.5 % not always needed
                McStasStr.declare{end+1}='tmp_length';
            end
            l{1} = 'minlengthneeded = 0;';
            l{end+1} = '';
            for i = find(globalinfo.maxstart)
                if i>index && sum(max(globalinfo.fixedlength(index+1:i),globalinfo.maxlength(index+1:i)))>i-index-0.5
                    % minneeded = max af ( minneeded og End(index)-nextmaxstart- sum)
                    l{end} = '';
                    l{end+1} = ['tmp_length = endPoint' num ' - maxStart' num2str(i)];
                        for j=find(globalinfo.fixedlength)
                            if j<=i && j>index
                                l{end}=[l{end} ' - length' num2str(j) ];
                            end
                        end
                        for j=find(globalinfo.maxlength)
                            if j<=i && j>index
                                l{end}=[l{end} ' - maxlength' num2str(j) ];
                            end
                        end
                        l{end} = [l{end} ';'];
                        l{end+1} = 'if (minlengthneeded < tmp_length) minlengthneeded = tmp_length;';
                        l{end+1} = '';
                end
            end
            
            if globalinfo.maxstart(index)==1
                l{end+1} = ['if (minlengthneeded < endPoint' num ' - maxStart' num ') minlengthneeded = endPoint' num ' - maxStart' num ';'];
                l{end+1} = '';
            end
            
            if globalinfo.minlength(index)==1
                l{end+1} = ['if (minlengthneeded < ' num2str(globalinfo.minlengthdata(index)) ') minlengthneeded = ' num2str(globalinfo.minlengthdata(index)) ';'];
                l{end+1} = '';
            end
            
            %            % the following logic just seems wrong. New and improved above:
            %            if nextmaxfp < nextminorfp
            %                 % this next logic decision should be taken inside mcstas
            %                 % actually, the above one too, sigh
            %                 minorfixedlength=sum(globalinfo.fixedlengthdata(index+1:nextmaxfp)+globalinfo.maxlengthdata(index+1:nextmaxfp));
            %                 % above line should be done in McStas, not matlab,
            %                 % otherwise inputs are circumvented OK for test though
            %
            %                 l{end+1} = '';
            %                 l{end+1} = '// The if statement below is nessecary because of complicated requirements in the input string';
            %                 l{end+1} = ['if ' num2str(minorfixedlength) '< endPoint' num '- maxStart' num2str(nextmaxfp) '{' ];
            %                 l{end+1} = ['  minlengthneeded = endPoint' num ' - maxStart' num2str(nextmaxfp) ' - ' num2str(minorfixedlength) ';'];  % maxStart num should be declared in main
            %                 l{end+1} = ['}'];
            %                 l{end+1} = ['else minlengthneeded = 0;'];
            %                 l{end+1} = '';
            %            else
            %                 l{end+1} = ['minlengthneeded = 0;'];
            %            end
            
            
            
            initstring='';
            for i=1:length(l)
                initstring=[initstring l{i} '\n'];
            end
            clear l;
            McStasStr.initialize=[initstring '\n\n' McStasStr.initialize];
            % old code
            %            if sum(globalinfo.fixedlength(index+1:nextfp))<nextfp-index-0.5 % Shouldn't it be more free length before the fixpoint!?
            %                % The above means checks if there is more modules left with
            %                % a free length
            %
            %                 linp=length(McStasStr.input);
            %                 McStasStr.input{linp+1}=['lengthfrac' num];
            %                 McStasStr.optimize(linp+1)=1;
            %                 McStasStr.optimvals.min(linp+1)=lengthfracpar(1);
            %                 McStasStr.optimvals.max(linp+1)=lengthfracpar(2);
            %                 McStasStr.declare{end+1}=['length' num];
            %                 if length(lengthfracpar)>2.5
            %                 McStasStr.optimvals.guess(linp+1)=lengthfracpar(3);
            %                 else
            %                 McStasStr.optimvals.guess(linp+1)=1./(1+last-index);
            %                 end
            %                strings = ['length' num ' = (endPoint' num ' - ' nextfpstr ];
            %                for i=find(globalinfo.fixedlength(index+1:nextfp))+index
            %                    strings = [ strings ' - length' num2str(i)];
            %                end
            %                strings = [ strings ')*lengthfrac' num ';'];
            %                McStasStr.initialize=[strings '\n\n' McStasStr.initialize];
            %            else
            %                %no free length left after this one
            %
            %                % the variable containing the fixed point have been declared
            %                % by the main script because it is global
            %                McStasStr.declare{end+1}=['length' num];
            %                strings = ['length' num ' = endPoint' num ' - ' nextfpstr ];
            %                for i=find(globalinfo.fixedlength(index+1:nextfp))+index
            %                    strings = [ strings ' - length' num2str(i)];
            %                end
            %                strings = [ strings ';'];
            %                McStasStr.initialize=[strings '\n\n' McStasStr.initialize];
            %            end
        end
    end
end

% EndPoint
ldec=length(McStasStr.declare);
McStasStr.declare{ldec+1}=['endPoint' num];
% Calculates the endPoint of the current module
if (index==1)
    strings=['endPoint' num ' = Mod_sample - sample_dist;'];
else
    strings=['endPoint' num ' = endPoint' numM1 ' - length' numM1 ' - u;'];
end
McStasStr.initialize=[strings '\n\n' McStasStr.initialize];

end

