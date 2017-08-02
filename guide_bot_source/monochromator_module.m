function [McStasStr] = monochromator_module(McStasStr,index,last,requirements,demands,options,globalinfo,ghost_globalinfo,defaults)
%Smodule adds text to McStas bot to make a straight guide
   % Adds the component to trace 
   % Handles the start parameters
   % Handles the end parameters
   % Handles the length parameter
   % Handles the endPoint parameter
   % Handles the options given to the function
   
num=num2str(index);
numM1=num2str(index-1);

defaultnames=fieldnames(defaults);
% Input not concerning the geometry
% WARNING, the indexes are used in the part extracting data from options
linp=length(McStasStr.input);
R0Index=linp+1;
McStasStr.input{R0Index}=['R0' num];
if max(ismember(defaultnames,'R0'))
McStasStr.inputvalue(R0Index)=defaults.R0;
else    
McStasStr.inputvalue(R0Index)=0.99;
end

QcIndex=linp+2;
McStasStr.input{QcIndex}=['Qc' num];
if max(ismember(defaultnames,'Qc'))
McStasStr.inputvalue(QcIndex)=defaults.Qc;    
else
McStasStr.inputvalue(QcIndex)=0.0217;
end

alphaIndex=linp+3;
McStasStr.input{alphaIndex}=['alpha' num];
if max(ismember(defaultnames,'alpha'))
McStasStr.inputvalue(alphaIndex)=defaults.alpha;
else
McStasStr.inputvalue(alphaIndex)=6.07;
end

mIndex=linp+4;
McStasStr.input{mIndex}=['m' num];
if max(ismember(defaultnames,'m'))
McStasStr.inputvalue(mIndex)=defaults.m;
else
McStasStr.inputvalue(mIndex)=3;
end

WIndex=linp+5;
McStasStr.input{WIndex}=['W' num];
if max(ismember(defaultnames,'W'))
McStasStr.inputvalue(WIndex)=defaults.W;
else
McStasStr.inputvalue(WIndex)=0.003;
end

optimize_end_logic(1) = 0; % width
optimize_end_logic(2) = 0; % height
optimlimits=-ones(1,4);
locked=-ones(1,2); % width,height
optimlimits_end=-ones(1,4);
locked_end=-ones(1,2); % width,height
% Extracting options from the option cell
for ops=1:length(options)
    for i=1:length(options{ops})
        if strcmp(options{ops}(i),'='); 
            equals=i;
        end
    end
    
    % check for options on reflectivity curve
    if strcmp(options{ops}(1:equals-1),'R0')
       McStasStr.inputvalue(R0Index)=str2num(options{ops}(equals+1:end));
    elseif strcmp(options{ops}(1:equals-1),'Qc')
       McStasStr.inputvalue(QcIndex)=str2num(options{ops}(equals+1:end));
    elseif strcmp(options{ops}(1:equals-1),'alpha')
       McStasStr.inputvalue(alphaIndex)=str2num(options{ops}(equals+1:end));
    elseif strcmp(options{ops}(1:equals-1),'m')
       McStasStr.inputvalue(mIndex)=str2num(options{ops}(equals+1:end));
    elseif strcmp(options{ops}(1:equals-1),'W')
       McStasStr.inputvalue(WIndex)=str2num(options{ops}(equals+1:end));
    elseif strcmp(options{ops}(1:equals-1),'minStartWidth')
       if optimlimits(1)>0
           optimlimits(1)=max([optimlimits(1) str2num(options{ops}(equals+1:end))]);
       else
           optimlimits(1)=str2num(options{ops}(equals+1:end));
       end
    elseif strcmp(options{ops}(1:equals-1),'maxStartWidth')
       if optimlimits(2)>0
           optimilimts(2)=min([optimlimits(2) str2num(options{ops}(equals+1:end))]);
       else
           optimlimits(2)=str2num(options{ops}(equals+1:end));
       end
    elseif strcmp(options{ops}(1:equals-1),'minStartHeight')
       if optimlimits(3)>0
           optimilimts(3)=max([optimlimits(3) str2num(options{ops}(equals+1:end))]);
       else
           optimlimits(3)=str2num(options{ops}(equals+1:end));
       end
    elseif strcmp(options{ops}(1:equals-1),'maxStartHeight')
       if optimlimits(4)>0
           optimilimts(4)=min([optimlimits(4) str2num(options{ops}(equals+1:end))]);
       else
           optimlimits(4)=str2num(options{ops}(equals+1:end));   
       end
    elseif strcmp(options{ops}(1:equals-1),'StartWidth')
       if locked(1)>0
           locked(1)=min([locked(1) str2num(options{ops}(equals+1:end))]);
       else
           locked(1)=str2num(options{ops}(equals+1:end));   
       end
    elseif strcmp(options{ops}(1:equals-1),'StartHeight')
       if locked(2)>0
           locked(2)=min([locked(2) str2num(options{ops}(equals+1:end))]);
       else
           locked(2)=str2num(options{ops}(equals+1:end));   
       end
       
    elseif strcmp(options{ops}(1:equals-1),'minEndWidth') && index == 1
       optimize_end_logic(1) = 1;
       if optimlimits_end(1)>0
           optimlimits_end(1)=max([optimlimits_end(1) str2num(options{ops}(equals+1:end))]);
       else
           optimlimits_end(1)=str2num(options{ops}(equals+1:end));
       end
    elseif strcmp(options{ops}(1:equals-1),'maxEndWidth') && index == 1
       optimize_end_logic(1) = 1;
       if optimlimits_end(2)>0
           optimilimts_end(2)=min([optimlimits_end(2) str2num(options{ops}(equals+1:end))]);
       else
           optimlimits_end(2)=str2num(options{ops}(equals+1:end));
       end
    elseif strcmp(options{ops}(1:equals-1),'minEndHeight') && index == 1
       optimize_end_logic(2) = 1;
       if optimlimits_end(3)>0
           optimilimts_end(3)=max([optimlimits_end(3) str2num(options{ops}(equals+1:end))]);
       else
           optimlimits_end(3)=str2num(options{ops}(equals+1:end));
       end
    elseif strcmp(options{ops}(1:equals-1),'maxEndHeight') && index == 1
       optimize_end_logic(2) = 1;
       if optimlimits_end(4)>0
           optimilimts_end(4)=min([optimlimits_end(4) str2num(options{ops}(equals+1:end))]);
       else
           optimlimits_end(4)=str2num(options{ops}(equals+1:end));   
       end
    elseif strcmp(options{ops}(1:equals-1),'EndWidth') && index == 1
       optimize_end_logic(1) = 1;
       if locked_end(1)>0
           locked_end(1)=min([locked_end(1) str2num(options{ops}(equals+1:end))]);
       else
           locked_end(1)=str2num(options{ops}(equals+1:end));   
       end
    elseif strcmp(options{ops}(1:equals-1),'EndHeight') && index == 1
       optimize_end_logic(2) = 1;
       if locked_end(2)>0
           locked_end(2)=min([locked_end(2) str2num(options{ops}(equals+1:end))]);
       else
           locked_end(2)=str2num(options{ops}(equals+1:end));   
       end  
    elseif strcmp(options{ops}(1:equals-1),'Optimize_end') && index == 1
        if str2num(options{ops}(equals+1:end)) == 1
            optimize_end_logic(1) = 1;optimize_end_logic(2) = 1;
        end
    elseif strcmp(options{ops}(1:equals-1),'Optimize_end_width') && index == 1
        if str2num(options{ops}(equals+1:end)) == 1
            optimize_end_logic(1) = 1;
        end
    elseif strcmp(options{ops}(1:equals-1),'Optimize_end_height') && index == 1
        if str2num(options{ops}(equals+1:end)) == 1
            optimize_end_logic(2) = 1;
        end
       
    elseif (strcmp(options{ops}(1:equals-1),'start') || strcmp(options{ops}(1:equals-1),'length') || strcmp(options{ops}(1:equals-1),'changed') ||...
            strcmp(options{ops}(1:equals-1),'maxstart') || strcmp(options{ops}(1:equals-1),'minstart') || strcmp(options{ops}(1:equals-1),'maxlength') ||...
            strcmp(options{ops}(1:equals-1),'los_start') || strcmp(options{ops}(1:equals-1),'los_end') || strcmp(options{ops}(1:equals-1),'los_divide') ||...
            strcmp(options{ops}(1:equals-1),'minlength'))
    else
       disp(['ERROR, the input ' options{ops} ' is not valid in module S (straight guide)'])
    end
    
end

% trace string   
%l{1}=['COMPONENT straight_guide_' num '= Guide_gravity('];
%l{2}=['w1=startx' num ', h1=starty' num ', w2=endx' num ', h2=endy' num ', l=length' num ','];
%l{3}=['R0=R0' num ', Qc=Qc' num ', alpha=alpha' num ', m=m' num ', W=W' num ',G=-9.82)'];
%l{4}='AT (0,0, u) RELATIVE PREVIOUS';    
%l{5}='ROTATED (0,0,0) RELATIVE PREVIOUS';
%l{6}='';
%l{7}=['COMPONENT EndOfelement_' num '= Arm()'];
%l{8}=['AT (0,0,length' num ' + 2*u) RELATIVE PREVIOUS'];

l{1}=    ['COMPONENT mono_arm_start_' num '= Arm()'];
l{end+1}=['AT (0,0,mono_width' num ') RELATIVE PREVIOUS'];

l{end+1}=['COMPONENT mono_arm_position_' num '= Arm()'];
l{end+1}=['AT (0,0,0) RELATIVE PREVIOUS'];
l{end+1}=['ROTATED (0,mono_angle' num ',0) RELATIVE PREVIOUS'];

l{end+1}=['COMPONENT monochromator_' num '= Monochromator_flat('];
l{end+1}=['zwidth=mono_width' num ',yheight=mono_height' num ','];
l{end+1}=['mosaich=mono_mosaic_h' num ',mosaicv=mono_mosaic_v' num ','];
l{end+1}=['r0=mono_r0' num ',DM=mono_d_spacing' num ')'];
l{end+1}=['AT (0,0,0) RELATIVE PREVIOUS'];

l{end+1}=['COMPONENT mono_arm_end_' num '= Arm()'];
l{end+1}=['AT (0,0,mono_width' num ') RELATIVE PREVIOUS'];


tracestring='';
for i=1:length(l)
    tracestring=[tracestring l{i} '\n'];
end
McStasStr.trace=[McStasStr.trace '\n\n' tracestring];
clear l;

% trace string
% Et lignende stykke som overstående der laver N lige store dele med
% forskellige coating værdier. Måske skal der også lidt til i initialize.

if sum(ismember(fieldnames(McStasStr),'input_seg')) < 0.5
    McStasStr.input_seg{1}=['m' num '_' num2str(1)];
    McStasStr.inputvalue_seg(1)=0;
    McStasStr.optimize_seg(1)=0;
else
    idec = length(McStasStr.input_seg);
    McStasStr.input_seg{end+1}=['m' num '_' num2str(1)];
    McStasStr.inputvalue_seg(idec+1)=0;
    McStasStr.optimize_seg(idec+1)=0;
end

N_m_segs=5;
for seg=2:N_m_segs
    idec = length(McStasStr.input_seg);
    McStasStr.input_seg{end+1}=['m' num '_' num2str(seg)];
    McStasStr.inputvalue_seg(idec+1)=0;
    McStasStr.optimize_seg(idec+1)=0;
end

l{1} = '';
for seg = 1:N_m_segs
    % trace string   
    l{end+1}=['COMPONENT straight_guide_' num '_' num2str(seg) '= Guide_gravity('];
    l{end+1}=['w1=' num2str(seg-1) '*(endx' num '-startx' num ')*0.2+startx' num ','];
    l{end+1}=['h1=' num2str(seg-1) '*(endy' num '-starty' num ')*0.2+starty' num ','];
    l{end+1}=['w2=' num2str(seg) '*(endx' num '-startx' num ')*0.2+startx' num ','];
    l{end+1}=['h2=' num2str(seg) '*(endy' num '-starty' num ')*0.2+starty' num ','];
    l{end+1}=['l=length' num '/' num2str(N_m_segs) ','];
    l{end+1}=['R0=R0' num ', Qc=Qc' num ', alpha=alpha' num ', m=m' num '_' num2str(seg) ', W=W' num ',G=-9.82)'];
    l{end+1}='AT (0,0, u) RELATIVE PREVIOUS';    
    l{end+1}='ROTATED (0,0,0) RELATIVE PREVIOUS';
    l{end+1}='';
    l{end+1}=['COMPONENT EndOfelement_' num '_' num2str(seg) '= Arm()'];
    l{end+1}=['AT (0,0,(length' num ' + 2*u)/' num2str(N_m_segs) ') RELATIVE PREVIOUS'];
    l{end+1}='';
end

tracestring='';
for i=1:length(l)
    tracestring=[tracestring l{i} '\n'];
end
McStasStr.trace_seg=[McStasStr.trace_seg '\n\n' tracestring];
clear l;

% Limits
if optimlimits(1)>0
    xlimdata(1)=optimlimits(1);
else
    xlimdata(1)=0.005;
end
if optimlimits(2)>0
    xlimdata(2)=optimlimits(2);
else
    xlimdata(2)=0.2;
end
xlimdata(3)=0.03;
if xlimdata(3)>xlimdata(2) || xlimdata(3)<xlimdata(1)
    xlimdata(3)=mean(xlimdata(1:2));
end    

if optimlimits(3)>0
    ylimdata(1)=optimlimits(3);
else
    ylimdata(1)=0.005;
end
if optimlimits(4)>0
    ylimdata(2)=optimlimits(4);
else
    ylimdata(2)=0.2;
end
ylimdata(3)=0.03;
if ylimdata(3)>ylimdata(2) || ylimdata(3)<ylimdata(1)
    ylimdata(3)=mean(ylimdata(1:2));
end    

% Limits end (default values entered here)
if optimlimits_end(1)>0
    xlimdata_end(1)=optimlimits_end(1);
else
    xlimdata_end(1)=0.005;
end
if optimlimits_end(2)>0
    xlimdata_end(2)=optimlimits_end(2);
else
    xlimdata_end(2)=0.2;
end
xlimdata_end(3)=0.03;
if xlimdata_end(3)>xlimdata_end(2) || xlimdata_end(3)<xlimdata_end(1)
    xlimdata_end(3)=mean(xlimdata_end(1:2));
end    

if optimlimits_end(3)>0
    ylimdata_end(1)=optimlimits_end(3);
else
    ylimdata_end(1)=0.005;
end
if optimlimits_end(4)>0
    ylimdata_end(2)=optimlimits_end(4);
else
    ylimdata_end(2)=0.2;
end
ylimdata_end(3)=0.03;
if ylimdata_end(3)>ylimdata_end(2) || ylimdata_end(3)<ylimdata_end(1)
    ylimdata_end(3)=mean(ylimdata_end(1:2));
end    

% Tese next functions write the part of the instrument file for this module
% concerning start, end and length parameters.
% start
McStasStr=guide_writer_start(McStasStr,index,last,xlimdata,ylimdata,locked,globalinfo);
% The last two vectors is the [min,max,guess] of (x,y) size

% end
McStasStr=guide_writer_end(McStasStr,index,last,xlimdata_end,ylimdata_end,locked_end,globalinfo,optimize_end_logic);

% length
McStasStr=guide_writer_length(McStasStr,index,last,demands,requirements,globalinfo,[0.1 0.95]);


% Example of the component in mcstas to read of variable names
%COMPONENT name1 = Guide_gravity(
%    w1 = 1, h1 = 1, w2 = 1, h2 = 1, l = 1, R0 = 1, Qc = 1,
%    alpha = 1, m = 1, W = 1, nslit = 1, d = 1, nhslit = 1, G = 1,
%    wavy = 1, wavy_z = 1, wavy_tb = 1, wavy_lr = 1, chamfers = 0)
%  AT (0, 0, 10) RELATIVE source
%  ROTATED (0, 0, 0) RELATIVE origin
end

