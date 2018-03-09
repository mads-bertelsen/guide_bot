function [McStasStr] = visualizer_output(McStasStr,index,last,globalinfo,ghost_globalinfo,defaults)

% printf("kink rot = %%lf.\n\n",rot2);

num = num2str(index);
numM1 = num2str(index-1);

real_index = ghost_globalinfo.real_index(index);
real_index_str = num2str(real_index);

real_num = num2str(ghost_globalinfo.real_index(index));

if index == last;
    l{1} = ['fp = fopen(file_name,"w");'];

    vis_string = '';
    for i=1:length(l)
            vis_string=[vis_string l{i} '\n'];
    end
    McStasStr.visualizer_str=[McStasStr.visualizer_str '\n' vis_string];
end

if index == last
   % print moderator info 
    
end

los_line = 0;
if ghost_globalinfo.los_logic(index)
   if ghost_globalinfo.los_start_logic(index)
       % start of this module is the start of a los boundary.
       los = 'fprintf(fp,"Los s\\n");';
       los_line = 1;
   else
       % start of this module is the end of a los boundary.
       if ghost_globalinfo.los_end_data(index) == 1 && ghost_globalinfo.los_end_mode(index) == 1 && index ==1
         los = 'fprintf(fp,"Los t\\n");';
         los_line = 1;
       else
         los = 'fprintf(fp,"Los e\\n");';
         los_line = 1;
       end
   end
end

if strcmp(globalinfo.modules(globalinfo.modulelist(ghost_globalinfo.real_index(index))),'S')%% corresponds to straight guide
    
    l{1} = 'fprintf(fp,"S\\n");';    
    l{end+1} = ['fprintf(fp,"x1: %%lf \\t %%lf \\t %%lf \\t %%lf\\n",startxpoint[' num '][1][1],startxpoint[' num '][2][1],startxpoint[' num '][1][2],startxpoint[' num '][2][2]);';];
    l{end+1} = ['fprintf(fp,"x2: %%lf \\t %%lf \\t %%lf \\t %%lf\\n",startxpoint[' numM1 '][1][1],startxpoint[' numM1 '][2][1],startxpoint[' numM1 '][1][2],startxpoint[' numM1 '][2][2]);';];
    l{end+1} = ['fprintf(fp,"y1: %%lf \\t %%lf \\t %%lf \\t %%lf\\n",startypoint[' num '][1][1],startypoint[' num '][2][1],startypoint[' num '][1][2],startypoint[' num '][2][2]);';];
    l{end+1} = ['fprintf(fp,"y2: %%lf \\t %%lf \\t %%lf \\t %%lf\\n",startypoint[' numM1 '][1][1],startypoint[' numM1 '][2][1],startypoint[' numM1 '][1][2],startypoint[' numM1 '][2][2]);';];
    if los_line; l{end+1} = los; end;
      
elseif strcmp(globalinfo.modules(globalinfo.modulelist(ghost_globalinfo.real_index(index))),'P')%% corresponds to parabolic guide

    l{1} = 'fprintf(fp,"P\\n");';
    l{end+1} = ['fprintf(fp,"x1: %%lf \\t %%lf \\t %%lf \\t %%lf\\n",startxpoint[' num '][1][1],startxpoint[' num '][2][1],startxpoint[' num '][1][2],startxpoint[' num '][2][2]);';];
    l{end+1} = ['fprintf(fp,"x2: %%lf \\t %%lf \\t %%lf \\t %%lf\\n",startxpoint[' numM1 '][1][1],startxpoint[' numM1 '][2][1],startxpoint[' numM1 '][1][2],startxpoint[' numM1 '][2][2]);';];
    l{end+1} = ['fprintf(fp,"y1: %%lf \\t %%lf \\t %%lf \\t %%lf\\n",startypoint[' num '][1][1],startypoint[' num '][2][1],startypoint[' num '][1][2],startypoint[' num '][2][2]);';];
    l{end+1} = ['fprintf(fp,"y2: %%lf \\t %%lf \\t %%lf \\t %%lf\\n",startypoint[' numM1 '][1][1],startypoint[' numM1 '][2][1],startypoint[' numM1 '][1][2],startypoint[' numM1 '][2][2]);';];
    l{end+1} = ['fprintf(fp,"xpos: %%lf \\t %%lf\\n",startXposition[' num '][1],startXposition[' num '][2]);'];
    l{end+1} = ['fprintf(fp,"ypos: %%lf \\t %%lf\\n",startYposition[' num '][1],startYposition[' num '][2]);'];
    l{end+1} = ['fprintf(fp,"xdir: %%lf \\t %%lf\\n",startXdirec[' num '][1],startXdirec[' num '][2]);'];
    l{end+1} = ['fprintf(fp,"ydir: %%lf \\t %%lf\\n",startYdirec[' num '][1],startYdirec[' num '][2]);'];
    l{end+1} = ['fprintf(fp,"xpars: %%lf \\t %%lf \\t %%lf \\t %%lf\\n",length' real_index_str ',smallaxis_parabolic_x' real_index_str ',Linx' real_index_str,',Loutx' real_index_str ');'];
    l{end+1} = ['fprintf(fp,"ypars: %%lf \\t %%lf \\t %%lf \\t %%lf\\n",length' real_index_str ',smallaxis_parabolic_y' real_index_str ',Liny' real_index_str,',Louty' real_index_str ');'];
    if los_line; l{end+1} = los; end;
    
elseif strcmp(globalinfo.modules(globalinfo.modulelist(ghost_globalinfo.real_index(index))),'E')%% corresponds to elliptic guide
    
    l{1} = 'fprintf(fp,"E\\n");';
    l{end+1} = ['fprintf(fp,"x1: %%lf \\t %%lf \\t %%lf \\t %%lf\\n",startxpoint[' num '][1][1],startxpoint[' num '][2][1],startxpoint[' num '][1][2],startxpoint[' num '][2][2]);';];
    l{end+1} = ['fprintf(fp,"x2: %%lf \\t %%lf \\t %%lf \\t %%lf\\n",startxpoint[' numM1 '][1][1],startxpoint[' numM1 '][2][1],startxpoint[' numM1 '][1][2],startxpoint[' numM1 '][2][2]);';];
    l{end+1} = ['fprintf(fp,"y1: %%lf \\t %%lf \\t %%lf \\t %%lf\\n",startypoint[' num '][1][1],startypoint[' num '][2][1],startypoint[' num '][1][2],startypoint[' num '][2][2]);';];
    l{end+1} = ['fprintf(fp,"y2: %%lf \\t %%lf \\t %%lf \\t %%lf\\n",startypoint[' numM1 '][1][1],startypoint[' numM1 '][2][1],startypoint[' numM1 '][1][2],startypoint[' numM1 '][2][2]);';];
    l{end+1} = ['fprintf(fp,"xpos: %%lf \\t %%lf\\n",startXposition[' num '][1],startXposition[' num '][2]);'];
    l{end+1} = ['fprintf(fp,"ypos: %%lf \\t %%lf\\n",startYposition[' num '][1],startYposition[' num '][2]);'];
    l{end+1} = ['fprintf(fp,"xdir: %%lf \\t %%lf\\n",startXdirec[' num '][1],startXdirec[' num '][2]);'];
    l{end+1} = ['fprintf(fp,"ydir: %%lf \\t %%lf\\n",startYdirec[' num '][1],startYdirec[' num '][2]);'];
    l{end+1} = ['fprintf(fp,"xpars: %%lf \\t %%lf \\t %%lf \\t %%lf\\n",length' real_index_str ',smallaxis_x' real_index_str ',Linx' real_index_str,',Loutx' real_index_str ');'];
    l{end+1} = ['fprintf(fp,"ypars: %%lf \\t %%lf \\t %%lf \\t %%lf\\n",length' real_index_str ',smallaxis_y' real_index_str ',Liny' real_index_str,',Louty' real_index_str ');'];
    if los_line; l{end+1} = los; end;
    
elseif strcmp(globalinfo.modules(globalinfo.modulelist(ghost_globalinfo.real_index(index))),'C') || strcmp(globalinfo.modules(globalinfo.modulelist(ghost_globalinfo.real_index(index))),'Cg')%% corresponds to curved guide%% corresponds to curved guide 
    
    l{1} = 'fprintf(fp,"C\\n");';
    l{end+1} = ['fprintf(fp,"x1: %%lf \\t %%lf \\t %%lf \\t %%lf\\n",startxpoint[' num '][1][1],startxpoint[' num '][2][1],startxpoint[' num '][1][2],startxpoint[' num '][2][2]);';];
    l{end+1} = ['fprintf(fp,"x2: %%lf \\t %%lf \\t %%lf \\t %%lf\\n",startxpoint[' numM1 '][1][1],startxpoint[' numM1 '][2][1],startxpoint[' numM1 '][1][2],startxpoint[' numM1 '][2][2]);';];
    l{end+1} = ['fprintf(fp,"y1: %%lf \\t %%lf \\t %%lf \\t %%lf\\n",startypoint[' num '][1][1],startypoint[' num '][2][1],startypoint[' num '][1][2],startypoint[' num '][2][2]);';];
    l{end+1} = ['fprintf(fp,"y2: %%lf \\t %%lf \\t %%lf \\t %%lf\\n",startypoint[' numM1 '][1][1],startypoint[' numM1 '][2][1],startypoint[' numM1 '][1][2],startypoint[' numM1 '][2][2]);';];
    if globalinfo.rotlogic(ghost_globalinfo.real_index(index)) > 0
    l{end+1} = ['fprintf(fp,"x: %%lf \\t %%lf \\t %%lf \\t %%lf \\t %%lf\\n",curve_radius' real_num ',curve_small_radius' real_num ',curveXcenter' real_num '[1],curveXcenter' real_num '[2],channels' real_num ');'];
    else
    l{end+1} = ['fprintf(fp,"y: %%lf \\t %%lf \\t %%lf \\t %%lf \\t %%lf\\n",curve_radius' real_num ',curve_small_radius' real_num ',curveYcenter' real_num '[1],curveYcenter' real_num '[2],channels' real_num ');'];
    end
    if los_line; l{end+1} = los; end;
    
elseif strcmp(globalinfo.modules(globalinfo.modulelist(ghost_globalinfo.real_index(index))),'K')%% corresponds to kink
    
    l{1} = 'fprintf(fp,"K\\n");';
    l{end+1} = ['fprintf(fp,"x1: %%lf \\t %%lf \\t %%lf \\t %%lf\\n",startxpoint[' num '][1][1],startxpoint[' num '][2][1],startxpoint[' num '][1][2],startxpoint[' num '][2][2]);';];
    l{end+1} = ['fprintf(fp,"x2: %%lf \\t %%lf \\t %%lf \\t %%lf\\n",startxpoint[' numM1 '][1][1],startxpoint[' numM1 '][2][1],startxpoint[' numM1 '][1][2],startxpoint[' numM1 '][2][2]);';];
    l{end+1} = ['fprintf(fp,"y1: %%lf \\t %%lf \\t %%lf \\t %%lf\\n",startypoint[' num '][1][1],startypoint[' num '][2][1],startypoint[' num '][1][2],startypoint[' num '][2][2]);';];
    l{end+1} = ['fprintf(fp,"y2: %%lf \\t %%lf \\t %%lf \\t %%lf\\n",startypoint[' numM1 '][1][1],startypoint[' numM1 '][2][1],startypoint[' numM1 '][1][2],startypoint[' numM1 '][2][2]);';];    

elseif strcmp(globalinfo.modules(globalinfo.modulelist(ghost_globalinfo.real_index(index))),'G')%% corresponds to Selene
    
    l{1} = 'fprintf(fp,"G\\n");';
    l{end+1} = ['fprintf(fp,"x1: %%lf \\t %%lf \\t %%lf \\t %%lf\\n",startxpoint[' num '][1][1],startxpoint[' num '][2][1],startxpoint[' num '][1][2],startxpoint[' num '][2][2]);';];
    l{end+1} = ['fprintf(fp,"x2: %%lf \\t %%lf \\t %%lf \\t %%lf\\n",startxpoint[' numM1 '][1][1],startxpoint[' numM1 '][2][1],startxpoint[' numM1 '][1][2],startxpoint[' numM1 '][2][2]);';];
    l{end+1} = ['fprintf(fp,"y1: %%lf \\t %%lf \\t %%lf \\t %%lf\\n",startypoint[' num '][1][1],startypoint[' num '][2][1],startypoint[' num '][1][2],startypoint[' num '][2][2]);';];
    l{end+1} = ['fprintf(fp,"y2: %%lf \\t %%lf \\t %%lf \\t %%lf\\n",startypoint[' numM1 '][1][1],startypoint[' numM1 '][2][1],startypoint[' numM1 '][1][2],startypoint[' numM1 '][2][2]);';];
    
elseif strcmp(globalinfo.modules(globalinfo.modulelist(ghost_globalinfo.real_index(index))),'Selene')%% corresponds to Selene
    
    l{1} = 'fprintf(fp,"Selene\\n");';
    l{end+1} = ['fprintf(fp,"x1: %%lf \\t %%lf \\t %%lf \\t %%lf\\n",startxpoint[' num '][1][1],startxpoint[' num '][2][1],startxpoint[' num '][1][2],startxpoint[' num '][2][2]);';];
    l{end+1} = ['fprintf(fp,"x2: %%lf \\t %%lf \\t %%lf \\t %%lf\\n",startxpoint[' numM1 '][1][1],startxpoint[' numM1 '][2][1],startxpoint[' numM1 '][1][2],startxpoint[' numM1 '][2][2]);';];
    l{end+1} = ['fprintf(fp,"y1: %%lf \\t %%lf \\t %%lf \\t %%lf\\n",startypoint[' num '][1][1],startypoint[' num '][2][1],startypoint[' num '][1][2],startypoint[' num '][2][2]);';];
    l{end+1} = ['fprintf(fp,"y2: %%lf \\t %%lf \\t %%lf \\t %%lf\\n",startypoint[' numM1 '][1][1],startypoint[' numM1 '][2][1],startypoint[' numM1 '][1][2],startypoint[' numM1 '][2][2]);';];
    l{end+1} = ['fprintf(fp,"xpos: %%lf \\t %%lf \\t %%lf \\t %%lf\\n",startXposition[' num '][1],startXposition[' num '][2],startXposition[' numM1 '][1],startXposition[' numM1 '][2]);'];
    l{end+1} = ['fprintf(fp,"ypos: %%lf \\t %%lf \\t %%lf \\t %%lf\\n",startYposition[' num '][1],startYposition[' num '][2],startYposition[' numM1 '][1],startYposition[' numM1 '][2]);'];
    l{end+1} = ['fprintf(fp,"pars: %%lf \\t %%lf \\t %%lf\\n",selene_distance' real_index_str ',selene_length' real_index_str ',selene_c' real_index_str ');'];
    l{end+1} = ['fprintf(fp,"xy: %%lf \\t %%lf\\n",selene_b_x' real_index_str ',selene_b_y' real_index_str ');'];
    
    
elseif strcmp(globalinfo.modules(globalinfo.modulelist(ghost_globalinfo.real_index(index))),'Slit')%% corresponds to Slit
    
    l{1} = 'fprintf(fp,"Slit\\n");';
    l{end+1} = ['fprintf(fp,"x1: %%lf \\t %%lf \\t %%lf \\t %%lf\\n",startxpoint[' num '][1][1],startxpoint[' num '][2][1],startxpoint[' num '][1][2],startxpoint[' num '][2][2]);';];
    l{end+1} = ['fprintf(fp,"x2: %%lf \\t %%lf \\t %%lf \\t %%lf\\n",startxpoint[' numM1 '][1][1],startxpoint[' numM1 '][2][1],startxpoint[' numM1 '][1][2],startxpoint[' numM1 '][2][2]);';];
    l{end+1} = ['fprintf(fp,"y1: %%lf \\t %%lf \\t %%lf \\t %%lf\\n",startypoint[' num '][1][1],startypoint[' num '][2][1],startypoint[' num '][1][2],startypoint[' num '][2][2]);';];
    l{end+1} = ['fprintf(fp,"y2: %%lf \\t %%lf \\t %%lf \\t %%lf\\n",startypoint[' numM1 '][1][1],startypoint[' numM1 '][2][1],startypoint[' numM1 '][1][2],startypoint[' numM1 '][2][2]);';];

%%%%%%% LELAND MODIFICATION %%%%%%%%
% Added elseif portion for M 
elseif strcmp(globalinfo.modules(globalinfo.modulelist(ghost_globalinfo.real_index(index))),'M')%% corresponds to Monochromator
    l{1} = 'fprintf(fp, "M\\n");';
    l{end+1} = ['fprintf(fp,"Dim: %%lf \\t %%lf \\t %%lf\\n",length' num ',startx' num ',starty' num ');';];
    %Mono_options = Parse_options(globalinfo.options{last-index+1});
    Mono_options = Parse_options(globalinfo.options{length(globalinfo.options)-index+1});
    if (~strcmp(Mono_options.HGeometry,'flat') && ~strcmp(Mono_options.VGeometry,'flat'))
        l{end+1} = ['fprintf(fp,"Lf: %%lf \\t %%lf\\n",HL1_' num ',VL1_' num ');';];
    elseif (strcmp(Mono_options.HGeometry,'flat') && ~strcmp(Mono_options.VGeometry,'flat'))
        l{end+1} = ['fprintf(fp,"Lf: N/A \\t %%lf\\n",VL1_' num ');';];
    elseif (~strcmp(Mono_options.HGeometry,'flat') && strcmp(Mono_options.VGeometry,'flat'))
        l{end+1} = ['fprintf(fp,"Lf: %%lf \\t N/A\\n",HL1_' num ');';];
    elseif (strcmp(Mono_options.HGeometry,'flat') && strcmp(Mono_options.VGeometry,'flat'))
        l{end+1} = 'fprintf(fp,"Lf: N/A \\t N/A\\n");';
    else
        disp('Invalid monochromator defined')
    end
    l{end+1} = ['fprintf(fp,"zp: %%lf\\n",startxpoint[' num '][1][1] + length' num '/2);';];
%%%% END OF LELAND MODIFICATION %%%%
    
else
    disp(['Visualizer_output not up to date, a new module must have beeen added' globalinfo.modules(globalinfo.modulelist(ghost_globalinfo.real_index(index)))])
end

if index == 1
   % print sample info
   l{end+1} = 'fprintf(fp,"Sample\\n");';
   l{end+1} = 'fprintf(fp,"%%lf \\t %%lf \\t %%lf\\n",sizeX,sizeY,sample_dist);';
   l{end+1} = 'fprintf(fp,"xp: %%lf \\t %%lf\\n",startXposition[0][1],startXposition[0][2]);';
   l{end+1} = 'fprintf(fp,"yp: %%lf \\t %%lf\\n",startYposition[0][1],startYposition[0][2]);';
   l{end+1} = 'fprintf(fp,"xd: %%lf \\t %%lf\\n",startXdirec[0][1],startXdirec[0][2]);';
   l{end+1} = 'fprintf(fp,"yd: %%lf \\t %%lf\\n",startYdirec[0][1],startYdirec[0][2]);';
   l{end+1} = 'fprintf(fp,"Moderator\\n");';
   l{end+1} = 'fprintf(fp,"%%lf \\t %%lf\\n",mod_x,mod_y);';
end

    

if index == 1;
    l{end+1} =['fclose(fp);'];
end


    vis_string = '';
    for i=1:length(l)
            vis_string=[vis_string l{i} '\n'];
    end
    McStasStr.visualizer_str=[McStasStr.visualizer_str '\n' vis_string];

end

% Leland Modification (included support function)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% START OF SUPPORT FUNCTION PARSE_OPTIONS() %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function StructOut = Parse_options(options)
for ops=1:length(options)
    command = options{ops};
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
StructOut = StructArg;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% END OF SUPPORT FUNCTION PARSE_OPTIONS() %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%