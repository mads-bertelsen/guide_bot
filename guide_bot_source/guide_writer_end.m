function McStasStr=guide_writer_end(McStasStr,index,last,startxpar,startypar,locked,globalinfo,optimize_end_logic)
%guide_writer_start works with the other guide_writers to make the basis of
%a guide module for the mcstas_bot program

num=num2str(index);

if (index==1) % first module after the sample
    
    
    if optimize_end_logic(1) == 0
        ldec=length(McStasStr.declare);
        McStasStr.declare{ldec+1}=['endx' num];
        l{1}=['endx' num ' = sizeX + 2*sample_dist*tan(divreq_x*DEG2RAD);'];

        initstring='';
        for i=1:length(l)
            initstring=[initstring l{i} '\n'];
        end
        clear l;
        McStasStr.initialize=[initstring '\n\n' McStasStr.initialize];
    else
        linp=length(McStasStr.input);
        McStasStr.input{linp+1}=['endx' num];
        if locked(1)<0
        McStasStr.optimize(linp+1)=1;    
        McStasStr.optimvals.min(linp+1)=startxpar(1);
        McStasStr.optimvals.max(linp+1)=startxpar(2);
        McStasStr.optimvals.guess(linp+1)=startxpar(3);
        else
        McStasStr.inputvalue(linp+1)=locked(1);
        % could introduce another option in locked to select an
        % parameterisation that optimizes endx based on the calculated
        % optimum.
        end
    end
    
    if optimize_end_logic(2) == 0
        ldec=length(McStasStr.declare);
        McStasStr.declare{ldec+1}=['endy' num];
        l{1}=['endy' num ' = sizeY + 2*sample_dist*tan(divreq_y*DEG2RAD);'];

        initstring='';
        for i=1:length(l)
            initstring=[initstring l{i} '\n'];
        end
        clear l;
        McStasStr.initialize=[initstring '\n\n' McStasStr.initialize];
    else
        linp=length(McStasStr.input);
        McStasStr.input{linp+1}=['endy' num];
        if locked(2)<0
        McStasStr.optimize(linp+1)=1;    
        McStasStr.optimvals.min(linp+1)=startypar(1);
        McStasStr.optimvals.max(linp+1)=startypar(2);
        McStasStr.optimvals.guess(linp+1)=startypar(3);
        else
        McStasStr.inputvalue(linp+1)=locked(2);
        % could introduce another option in locked to select an
        % parameterisation that optimizes endx based on the calculated
        % optimum.
        end
    end
    
else
    ldec=length(McStasStr.declare);
    McStasStr.declare{ldec+1}=['endx' num];
    McStasStr.declare{ldec+2}=['endy' num];
    
    l{1}=['endx' num ' = startx' num2str(index-1) ';'];
    l{2}=['endy' num ' = starty' num2str(index-1) ';'];
    
    initstring='';
    for i=1:length(l)
        initstring=[initstring l{i} '\n'];
    end
    clear l;
    McStasStr.initialize=[initstring '\n\n' McStasStr.initialize];    
end


end

