function McStasStr=guide_writer_start(McStasStr,index,last,startxpar,startypar,locked,globalinfo)
%guide_writer_start works with the other guide_writers to make the basis of
%a guide module for the mcstas_bot program

num = num2str(index);

% Could maybe be moved out of the index == last if statement

if (index==last) && McStasStr.minimalist==1 %start can be calculated
    
    
    if globalinfo.minimalist_direction_horizontal == 1
        % declare start  
        ldec=length(McStasStr.declare);
        McStasStr.declare{ldec+1}=['startx' num];


        % var_divreq known
        % calculate dimensions of the guide opening so that it is illuminated
        % in phase-space

        if strcmp(globalinfo.PS_req,'homogen')
            if globalinfo.focusing == 0
                l{1}=['if (0<mod_x*mod_x - sqrt(minimalist_factor)*endx' num '*tan(var_divreq_x*DEG2RAD)*8*guide_start) {'];
                l{1+end}=['startx' num ' = mod_x/2 - sqrt(mod_x*mod_x - sqrt(minimalist_factor)*endx' num '*tan(var_divreq_x*DEG2RAD)*8*guide_start)/2;'];
                l{1+end}=['}'];
                l{1+end}=['else {'];
                l{1+end}=['startx' num ' = mod_x/2;'];
                l{1+end}=['}'];
            elseif globalinfo.focusing == 1
                l{1}=['if (0<mod_x*mod_x - sqrt(minimalist_factor)*(2*endx' num '*var_divreq_x*DEG2RAD-2*sample_dist*tan(divreq_x*DEG2RAD)*divreq_x*DEG2RAD)*4*guide_start) {'];
                l{1+end}=['startx' num ' = mod_x/2 - sqrt(mod_x*mod_x - sqrt(minimalist_factor)*(2*endx' num '*var_divreq_x*DEG2RAD-2*sample_dist*tan(divreq_x*DEG2RAD)*divreq_x*DEG2RAD)*4*guide_start)/2;'];
                l{1+end}=['}'];
                l{1+end}=['else {'];
                l{1+end}=['startx' num ' = mod_x/2;'];
                l{1+end}=['}'];
            end

        elseif strcmp(globalinfo.PS_req,'total')
            if globalinfo.focusing == 0
                l{1}=['startx' num ' = sqrt(minimalist_factor)*2*endx' num '*var_divreq_x*DEG2RAD*guide_start/mod_x;'];
            elseif globalinfo.focusing == 1
                l{1}=['startx' num ' = sqrt(minimalist_factor)*(2*endx' num '*var_divreq_x*DEG2RAD-2*sample_dist*tan(divreq_x*DEG2RAD)*divreq_x*DEG2RAD)*guide_start/mod_x;'];
            end
        else
            disp('ERROR, globalinfo.PS_req empty!')
        end

        % Add this calculation to the mcstas initialize section
        initstring='';
        for i=1:length(l)
            initstring=[initstring l{i} '\n'];
        end
        clear l;
        McStasStr.initialize=[initstring '\n\n' McStasStr.initialize]; 
    elseif globalinfo.minimalist_direction_horizontal == 0
        
        linp=length(McStasStr.input);
        McStasStr.input{linp+1}=['startx' num];
        if locked(1)<0
        McStasStr.optimize(linp+1)=1;    
        McStasStr.optimvals.min(linp+1)=startxpar(1);
        McStasStr.optimvals.max(linp+1)=startxpar(2);
        McStasStr.optimvals.guess(linp+1)=startxpar(3);
        else
        McStasStr.inputvalue(linp+1)=locked(1);
        end
        
    else
       disp('ERROR, (bug) globalinfo.minimalist_direction_horizontal not set correctly. Contact developer') 
    end
        
    
    if globalinfo.minimalist_direction_vertical == 1
        ldec=length(McStasStr.declare);
        McStasStr.declare{ldec+1}=['starty' num];

        if strcmp(globalinfo.PS_req,'homogen')
            if globalinfo.focusing == 0
                l{1}=['if (0<mod_y*mod_y - sqrt(minimalist_factor)*endy' num '*tan(var_divreq_y*DEG2RAD)*8*guide_start) {'];
                l{1+end}=['starty' num ' = mod_y/2 - sqrt(mod_y*mod_y - sqrt(minimalist_factor)*endy' num '*tan(var_divreq_y*DEG2RAD)*8*guide_start)/2;']; 
                l{1+end}=['}'];
                l{1+end}=['else {'];
                l{1+end}=['starty' num ' = mod_y/2;'];
                l{1+end}=['}'];
            elseif globalinfo.focusing == 1
                l{1}=['if (0<mod_y*mod_y - sqrt(minimalist_factor)*(2*endy' num '*var_divreq_y*DEG2RAD-2*sample_dist*tan(divreq_y*DEG2RAD)*divreq_y*DEG2RAD)*4*guide_start) {'];
                l{1+end}=['starty' num ' = mod_y/2 - sqrt(mod_y*mod_y - sqrt(minimalist_factor)*(2*endy' num '*var_divreq_y*DEG2RAD-2*sample_dist*tan(divreq_y*DEG2RAD)*divreq_y*DEG2RAD)*4*guide_start)/2;']; 
                l{1+end}=['}'];
                l{1+end}=['else {'];
                l{1+end}=['starty' num ' = mod_y/2;'];
                l{1+end}=['}'];
            end
        elseif strcmp(globalinfo.PS_req,'total')
            if globalinfo.focusing == 0
                l{1}=['starty' num ' = sqrt(minimalist_factor)*2*endy' num '*var_divreq_y*DEG2RAD*guide_start/mod_y;'];
            elseif globalinfo.focusing == 1
                l{1}=['starty' num ' = sqrt(minimalist_factor)*(2*endy' num '*var_divreq_y*DEG2RAD-2*sample_dist*tan(divreq_y*DEG2RAD)*divreq_y*DEG2RAD)*guide_start/mod_y;'];
            end
        else
            disp('ERROR, globalinfo.PS_req empty!')
        end

        % Add this calculation to the mcstas initialize section
        initstring='';
        for i=1:length(l)
            initstring=[initstring l{i} '\n'];
        end
        clear l;
        McStasStr.initialize=[initstring '\n\n' McStasStr.initialize]; 
    elseif globalinfo.minimalist_direction_vertical == 0
       
        linp=length(McStasStr.input);
        McStasStr.input{linp+1}=['starty' num];
        if locked(2)<0
        McStasStr.optimize(linp+1)=1;
        McStasStr.optimvals.min(linp+1)=startypar(1);
        McStasStr.optimvals.max(linp+1)=startypar(2);
        McStasStr.optimvals.guess(linp+1)=startypar(3);
        else
        McStasStr.inputvalue(linp+2)=locked(2);    
        end
    else
        disp('ERROR, (bug) globalinfo.minimalist_direction_vertical not set correctly. Contact developer')
    end
       
    
elseif (index == last - 1) && (globalinfo.first_fixedend == 1) && 2 == 1 %not in use! delete this code
    % This happens when the module closest to the moderator does not have
    % the possibility of having unequal start end end dimensions. Since
    % that module is first, the start width of this module is fixed.
    
    
    % declare startx num and starty num
    % assign these to endx num + 1 and endy num +1
    
    ldec=length(McStasStr.declare);
    McStasStr.declare{ldec+1}=['startx' num];
    McStasStr.declare{ldec+2}=['starty' num];
    
    % Maybe the modules which are the problem should add this decleration?
    l{1}=['startx' num ' = startx' num2str(index+1) ';'];
    l{2}=['starty' num ' = starty' num2str(index+1) ';'];
    
    initstring='';
    for i=1:length(l)
        initstring=[initstring l{i} '\n'];
    end
    clear l;
    McStasStr.initialize=[initstring '\n\n' McStasStr.initialize];
    
    
    % calculate required phasespace IF known (no gaps)
    % phase_x = endx num * tan(var_div_req_x)
    % phase_y = endy num * tan(var_div_req_y)
    
    % if gap, the length is known, and thus the phase space needed is:
    % phase_x = endx num * tan(var_div_req_x) + 2*gap_l*tan(var_div_req_x)
    % phase_y = endy num * tan(var_div_req_y) + 2*gap_l*tan(var_div_req_y)
    
    
    
    
else
    % can not be calculated and should be optimized
    % input start 
    % optimize start
    
    linp=length(McStasStr.input);
    McStasStr.input{linp+1}=['startx' num];
    if locked(1)<0
    McStasStr.optimize(linp+1)=1;    
    McStasStr.optimvals.min(linp+1)=startxpar(1);
    McStasStr.optimvals.max(linp+1)=startxpar(2);
    McStasStr.optimvals.guess(linp+1)=startxpar(3);
    else
    McStasStr.inputvalue(linp+1)=locked(1);
    end
    
    
    McStasStr.input{linp+2}=['starty' num];
    if locked(2)<0
    McStasStr.optimize(linp+2)=1;
    McStasStr.optimvals.min(linp+2)=startypar(1);
    McStasStr.optimvals.max(linp+2)=startypar(2);
    McStasStr.optimvals.guess(linp+2)=startypar(3);
    else
    McStasStr.inputvalue(linp+2)=locked(2);    
    end
    
    % think about what range start can be optimized within

    % update var_divreq_x/y to current value
    % var_divreq_x/y resembles the angle-space needed at this guide size
    l{1} = ['var_divreq_x = atan(endx' num '*tan(var_divreq_x*DEG2RAD)/startx' num ')*RAD2DEG;'];
    l{2} = ['var_divreq_y = atan(endy' num '*tan(var_divreq_y*DEG2RAD)/starty' num ')*RAD2DEG;'];
    
    
    initstring='';
    for i=1:length(l)
        initstring=[initstring l{i} '\n'];
    end
    clear l;
    McStasStr.initialize=[initstring '\n\n' McStasStr.initialize];    
end

end

