function [McStasStr] = E_module(McStasStr,index,last,requirements,demands,options,globalinfo,ghost_globalinfo,defaults)
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
if ismember(defaultnames,'R0')
McStasStr.inputvalue(R0Index)=defaults.R0;
else    
McStasStr.inputvalue(R0Index)=0.99;
end

QcIndex=linp+2;
McStasStr.input{QcIndex}=['Qc' num];
if ismember(defaultnames,'Qc')
McStasStr.inputvalue(QcIndex)=defaults.Qc;    
else
McStasStr.inputvalue(QcIndex)=0.0217;
end

alphaIndex=linp+3;
McStasStr.input{alphaIndex}=['alpha' num];
if ismember(defaultnames,'alpha')
McStasStr.inputvalue(alphaIndex)=defaults.alpha;
else
McStasStr.inputvalue(alphaIndex)=6.07;
end

mIndex=linp+4;
McStasStr.input{mIndex}=['m' num];
if ismember(defaultnames,'m')
McStasStr.inputvalue(mIndex)=defaults.m;
else
McStasStr.inputvalue(mIndex)=3;
end

WIndex=linp+5;
McStasStr.input{WIndex}=['W' num];
if ismember(defaultnames,'W')
McStasStr.inputvalue(WIndex)=defaults.W;
else
McStasStr.inputvalue(WIndex)=0.003;
end

optimlimits=-ones(1,4);
locked=-ones(1,2); % width,height
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
    elseif (strcmp(options{ops}(1:equals-1),'start') || strcmp(options{ops}(1:equals-1),'length') || strcmp(options{ops}(1:equals-1),'changed') ||...
            strcmp(options{ops}(1:equals-1),'maxstart') || strcmp(options{ops}(1:equals-1),'minstart') || strcmp(options{ops}(1:equals-1),'maxlength') ||...
            strcmp(options{ops}(1:equals-1),'los_start') || strcmp(options{ops}(1:equals-1),'los_end') || strcmp(options{ops}(1:equals-1),'los_divide') ||...
            strcmp(options{ops}(1:equals-1),'minlength'))
    else
       disp(['ERROR, the input ' options{ops} ' is not valid in module E (Elliptic guide)'])
    end
    
end


% COMPONENT test = Elliptic_guide_gravity(
%     L = 1, Linh = 1, Louth = 1, Linv = 1, Loutv = 1,
%     widthOfEllipseh = 1, widthOfEllipsev = 1, R0 = 1, Qc = 1,
%     alpha = 1, m = 1, W = 1)
%   AT (0, 0, 0) RELATIVE PREVIOUS
%   ROTATED (0, 0, 0) RELATIVE PREVIOUS

% The commented lines below is for the 84beta component.
% trace string   
l{1}=['COMPONENT elliptical_guide_gravity' num '= Elliptic_guide_gravity('];
%l{1}=['COMPONENT elliptical_guide_gravity' num '= guide_elliptical_84beta('];
%l{end+1}=['l=length' num ', linxw=Linx' num ', linyh=Liny' num ', loutxw=Loutx' num ', loutyh=Louty' num ','];
%l{end+1}=['xwidth = smallaxis_x' num ', yheight = smallaxis_y' num ', dimensionsAt="mid",'];
l{end+1}=['L=length' num ', Linh=Linx' num ', Linv=Liny' num ', Louth=Loutx' num ', Loutv=Louty' num ','];
l{end+1}=['widthOfEllipseh = smallaxis_x' num ', widthOfEllipsev = smallaxis_y' num ','];
l{end+1}=['R0=R0' num ', Qc=Qc' num ', alpha=alpha' num ', m=m' num ', W=W' num ')'];
l{end+1}='AT (0,0, u) RELATIVE PREVIOUS';    
l{end+1}='ROTATED (0,0,0) RELATIVE PREVIOUS';
l{end+1}='';
l{end+1}=['COMPONENT EndOfelement_' num '= Arm()'];
l{end+1}=['AT (0,0,length' num ' + 2*u) RELATIVE PREVIOUS'];

tracestring='';
for i=1:length(l)
    tracestring=[tracestring l{i} '\n'];
end
McStasStr.trace=[McStasStr.trace '\n\n' tracestring];
clear l;

% Alternate version of the trace section for coating optimization: 
% maybe use the new version of the code with coating segments?

l{1}=['COMPONENT elliptical_guide_gravity' num '= guide_elliptical_89beta('];
l{end+1}=['l=length' num ', linxw=Linx' num ', linyh=Liny' num ', loutxw=Loutx' num ', loutyh=Louty' num ','];
l{end+1}=['xwidth = 2*smallaxis_x' num ', yheight = 2*smallaxis_y' num ', dimensionsAt="mid",'];
l{end+1}=['R0=R0' num ', Qc=Qc' num ', alpha=alpha' num ', W=W' num ','];
l{end+1}=['mvaluesright = mValues' num ', mvaluesleft = mValues' num ',mvaluestop = mValues' num ',mvaluesbottom = mValues' num ','];
l{end+1}=['seglength = elementLength' num ')'];
l{end+1}='AT (0,0, u) RELATIVE PREVIOUS';    
l{end+1}='ROTATED (0,0,0) RELATIVE PREVIOUS';
l{end+1}='';
l{end+1}=['COMPONENT EndOfelement_' num '= Arm()'];
l{end+1}=['AT (0,0,length' num ' + 2*u) RELATIVE PREVIOUS'];

tracestring='';
for i=1:length(l)
    tracestring=[tracestring l{i} '\n'];
end
McStasStr.trace_seg=[McStasStr.trace_seg '\n\n' tracestring];
clear l;

% controlls the number of segments per elliptic guide
N_m_segs=5;
if sum(ismember(fieldnames(McStasStr),'input_seg')) < 0.5
    McStasStr.input_seg{1}=['m' num '_' num2str(1)];
    McStasStr.inputvalue_seg(1)=0;
    McStasStr.optimize_seg(1)=0;
else
    idec = length(McStasStr.input_seg);
    McStasStr.input_seg{idec+1}=['m' num '_' num2str(1)];
    McStasStr.inputvalue_seg(idec+1)=0;
    McStasStr.optimize_seg(idec+1)=0;
end
for seg=2:N_m_segs
idec = length(McStasStr.input_seg);    
McStasStr.input_seg{idec+1}=['m' num '_' num2str(seg)];
McStasStr.inputvalue_seg(idec+1) = 0; % will be changed when coating is relevant
McStasStr.optimize_seg(idec+1) = 0;
end

if sum(ismember(fieldnames(McStasStr),'declare_seg')) < 0.5
McStasStr.declare_seg{1}=['mValues' num '[' num2str(N_m_segs) ']'];
McStasStr.declare_seg{2}=['elementLength' num '[' num2str(N_m_segs) ']'];
else
McStasStr.declare_seg{end+1}=['mValues' num '[' num2str(N_m_segs) ']'];
McStasStr.declare_seg{end+1}=['elementLength' num '[' num2str(N_m_segs) ']'];
end

% Tese next functions write the part of the instrument file for this module
% concerning start, end and length parameters.

% code for calculating the focus points and smallaxis based on start and
% end determined by the minimalist concept and ifit optimizer

% Problem with the smallaxis. It can not be smaller than the average of the
% openings. factor from 1 to 5 ?

    linp=length(McStasStr.input);
    McStasStr.input{linp+1}=['smallaxis_x_factor' num];
    McStasStr.optimize(linp+1)=1;    
    McStasStr.optimvals.min(linp+1)=1.001;
    McStasStr.optimvals.max(linp+1)=9;
    McStasStr.optimvals.guess(linp+1)=1.25;
    
    McStasStr.input{linp+2}=['smallaxis_y_factor' num];
    McStasStr.optimize(linp+2)=1;    
    McStasStr.optimvals.min(linp+2)=1.001;
    McStasStr.optimvals.max(linp+2)=9;
    McStasStr.optimvals.guess(linp+2)=1.25;
    
    McStasStr.declare{end+1}=['smallaxis_x' num];
    McStasStr.declare{end+1}=['smallaxis_y' num];
    McStasStr.declare{end+1}=['Linx' num];
    McStasStr.declare{end+1}=['Liny' num];
    McStasStr.declare{end+1}=['Loutx' num];
    McStasStr.declare{end+1}=['Louty' num];
    
    if sum(ismember(McStasStr.declare,'tmp_k'))<0.5
    McStasStr.declare{end+1}=['tmp_k']; % declare temp variables
    McStasStr.declare{end+1}=['tmp_L1'];
    McStasStr.declare{end+1}=['tmp_L2'];
    McStasStr.declare{end+1}=['tmp_c'];
    McStasStr.declare{end+1}=['tmp_b'];
    McStasStr.declare{end+1}=['tmp_w1'];
    McStasStr.declare{end+1}=['tmp_w2'];
    McStasStr.declare{end+1}=['tmp_L'];
    end
    
    l{1} = [''];
    l{end+1} = ['if (startx' num ' > endx' num ') {'];
    l{end+1} = ['smallaxis_x' num ' = 0.5*startx' num '*smallaxis_x_factor' num ';'];
    %l{end+1} = ['smallaxis_x' num ' = 0.5*startx' num '+ 0.5*smallaxis_x_factor' num '*(max_smallaxis_x' num '-startx' num ');'];
    l{end+1} = ['}'];
    l{end+1} = ['else {'];
    l{end+1} = ['smallaxis_x' num ' = 0.5*endx' num '*smallaxis_x_factor' num ';'];
    l{end+1} = ['}'];
    l{end+1} = ['if (starty' num ' > endy' num ') {'];
    l{end+1} = ['smallaxis_y' num ' = 0.5*starty' num '*smallaxis_y_factor' num ';'];
    l{end+1} = ['}'];
    l{end+1} = ['else {'];
    l{end+1} = ['smallaxis_y' num ' = 0.5*endy' num '*smallaxis_y_factor' num ';'];
    l{end+1} = ['}'];
    %l{end+1} = ['smallaxis_x' num ' = 0.5*( startx' num ' + endx' num ' )*smallaxis_x_factor' num ';'];
    %l{end+1} = ['smallaxis_y' num ' = 0.5*( starty' num ' + endy' num ' )*smallaxis_y_factor' num ';'];
    l{end+1} = [''];
    l{end+1} = [''];
    l{end+1} = ['// calculating focus points for elliptic guide ' num ' x direction'];
    l{end+1} = ['tmp_w1=startx' num ';'];
    l{end+1} = ['tmp_w2=endx' num ';'];
    l{end+1} = ['tmp_b=smallaxis_x' num ';'];
    l{end+1} = ['tmp_L=length' num ';'];
    l{end+1} = [''];
    l{end+1} = ['tmp_k=cos(asin(tmp_w2/(2*tmp_b)))/cos(asin(tmp_w1/(2*tmp_b)));'];
    l{end+1} = ['tmp_L1=tmp_L/(1+tmp_k);'];
    l{end+1} = ['tmp_L2=tmp_L-tmp_L1;'];
    l{end+1} = ['tmp_c=cos(asin(tmp_w1/(2*tmp_b)))*cos(asin(tmp_w1/(2*tmp_b)));'];
    l{end+1} = ['Linx' num '=-tmp_L1+sqrt(tmp_L1*tmp_L1*(1/tmp_c)-tmp_b*tmp_b);'];
    l{end+1} = ['Loutx' num '= Linx' num '+tmp_L1-tmp_L2;'];
    l{end+1} = [''];
    l{end+1} = ['// calculating focus points for elliptic guide ' num ' y direction'];
    l{end+1} = ['tmp_w1=starty' num ';'];
    l{end+1} = ['tmp_w2=endy' num ';'];
    l{end+1} = ['tmp_b=smallaxis_y' num ';'];
    l{end+1} = ['tmp_L=length' num ';'];
    l{end+1} = [''];
    l{end+1} = ['tmp_k=cos(asin(tmp_w2/(2*tmp_b)))/cos(asin(tmp_w1/(2*tmp_b)));'];
    l{end+1} = ['tmp_L1=tmp_L/(1+tmp_k);'];
    l{end+1} = ['tmp_L2=tmp_L-tmp_L1;'];
    l{end+1} = ['tmp_c=cos(asin(tmp_w1/(2*tmp_b)))*cos(asin(tmp_w1/(2*tmp_b)));'];
    l{end+1} = ['Liny' num '=-tmp_L1+sqrt(tmp_L1*tmp_L1*(1/tmp_c)-tmp_b*tmp_b);'];
    l{end+1} = ['Louty' num '= Liny' num '+tmp_L1-tmp_L2;'];
    l{end+1} = [''];
    
    initstring='';
    for i=1:length(l)
        initstring=[initstring l{i} '\n'];
    end
    clear l;
    McStasStr.initialize=[initstring '\n\n' McStasStr.initialize];
    
    l{1} = ['// Setting up the mvalue array'];
    for seg=1:N_m_segs
    l{end+1} = ['mValues' num '[' num2str(seg-1) '] = m' num '_' num2str(seg) ';'];
    l{end+1} = ['elementLength' num '[' num2str(seg-1) '] = length' num '/' num2str(N_m_segs) ';']; 
    end
    
    initstring='';
    for i=1:length(l)
        initstring=[initstring l{i} '\n'];
    end
    clear l;
    McStasStr.initialize_seg=[initstring '\n\n' McStasStr.initialize_seg];
    
    

% // guide_1_x
% w1=w_in_guide_1_x;
% w2=w_out_guide_1_x;
% b=smallaxis_guide_1_x;
% L=guide_length;
% 
% k=cos(asin(w2/(2*b)))/cos(asin(w1/(2*b)));
% L1=L/(1+k);
% L2=L-L1;
% c=cos(asin(w1/(2*b)))*cos(asin(w1/(2*b)));
% focus_in_guide_x=-L1+sqrt(L1*L1*(1/c)-b*b);
% f1=focus_in_guide_x;
% focus_out_guide_x=f1+L1-L2;
% 
% // guide_1_y
% w1=w_in_guide_1_y;
% w2=w_out_guide_1_y;
% b=smallaxis_guide_1_y;
% L=guide_length;
% 
% k=cos(asin(w2/(2*b)))/cos(asin(w1/(2*b)));
% L1=L/(1+k);
% L2=L-L1;
% c=cos(asin(w1/(2*b)))*cos(asin(w1/(2*b)));
% focus_in_guide_y=-L1+sqrt(L1*L1*(1/c)-b*b);
% f1=focus_in_guide_y;
% focus_out_guide_y=f1+L1-L2;

if optimlimits(1)>0
    xlimdata(1)=optimlimits(1);
else
    xlimdata(1)=0.005;
end
if optimlimits(2)>0
    xlimdata(2)=optimlimits(2);
else
    xlimdata(2)=0.15;
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
    ylimdata(2)=0.15;
end
ylimdata(3)=0.03;
if ylimdata(3)>ylimdata(2) || ylimdata(3)<ylimdata(1)
    ylimdata(3)=mean(ylimdata(1:2));
end    

% start
McStasStr=guide_writer_start(McStasStr,index,last,xlimdata,ylimdata,locked,globalinfo);
% The last two vectors is the [min,max,guess] of (x,y) size

% end
McStasStr=guide_writer_end(McStasStr,index,last);

% length
McStasStr=guide_writer_length(McStasStr,index,last,demands,requirements,globalinfo,[0.1 0.95]);




end

