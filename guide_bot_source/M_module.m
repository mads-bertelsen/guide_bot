function [McStasStr] = M_module(McStasStr,index,last,requirements,demands,options,globalinfo,ghost_globalinfo,defaults)

% This is the sequential number assigned to the component with num = 1 corresponding
% to the component closest to the sample and num = max the component
% closest to the source.
num=num2str(index);

% Append monochromator options as new elements in McStasStr.input and McStasStr.inputvalue
% Assign a default value from the structure defaults, if no value there,
% then use a default value embedded in the code.
defaultnames=fieldnames(defaults);
linp=length(McStasStr.input);


%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Start of Monocromator options with default filling %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Used to define correct monochromator L2 distance if BeamDir = 'transmit' (m)
if max(ismember(defaultnames,'L2_override'))
    L2_override = defaults.L2_override;
else    
    L2_override = -1; % Don't override calculated distance.
end

% Monochromator reflectivity (single value)
r0Index = linp+1;
McStasStr.input{r0Index}=['r0' num];
if max(ismember(defaultnames,'r0'))
    McStasStr.inputvalue(r0Index)=defaults.r0;
else    
    McStasStr.inputvalue(r0Index)=0.7;
end

% Monochromator reflectivty (file)
% Note: McStasStr.inputvalue is a 1D array and cannot accept data type char
if max(ismember(defaultnames,'Rfile'))
    rfile=defaults.Rfile;
else
    rfile='NULL';
end

% Monochromator d-spacing (Ang)
dMIndex=linp+2;
McStasStr.input{dMIndex}=['dM' num];
if max(ismember(defaultnames,'dM'))
    McStasStr.inputvalue(dMIndex)=defaults.dM;    
else
    McStasStr.inputvalue(dMIndex)=3.354;
end

% Monochromator vertical mosaic (minutes arc)
VmosIndex=linp+3;
McStasStr.input{VmosIndex}=['Vmos' num];
if max(ismember(defaultnames,'Vmos'))
    McStasStr.inputvalue(VmosIndex)=defaults.Vmos;
else
    McStasStr.inputvalue(VmosIndex)=30.0;
end

% Monochromator horizontal mosaic (minutes arc)
HmosIndex=linp+4;
McStasStr.input{HmosIndex}=['Hmos' num];
if max(ismember(defaultnames,'Hmos'))
    McStasStr.inputvalue(HmosIndex)=defaults.Hmos;
else
    McStasStr.inputvalue(HmosIndex)=30.0;
end

% Isotropic mosaic (minutes arc)
if max(ismember(defaultnames,'mos'))
    McStasStr.inputvalue(HmosIndex)=defaults.mos;
    McStasStr.inputvalue(VmosIndex)=defaults.mos;
end
    

% Number of vertical segments
NVIndex=linp+5;
McStasStr.input{NVIndex}=['NV' num];
if max(ismember(defaultnames,'NV'))
    McStasStr.inputvalue(NVIndex)=defaults.NV;
else
    McStasStr.inputvalue(NVIndex)=10.0;
end

% Gap between segments (m)
SegGapIndex=linp+6;
McStasStr.input{SegGapIndex}=['SegGap' num];
if max(ismember(defaultnames,'SegGap'))
    McStasStr.inputvalue(SegGapIndex)=defaults.SegGap;
else
    McStasStr.inputvalue(SegGapIndex)=0.0005;
end

% Number of blades
BladeNIndex=linp+7;
McStasStr.input{BladeNIndex}=['BladeN' num];
if max(ismember(defaultnames,'BladeN'))
    McStasStr.inputvalue(BladeNIndex)=defaults.BladeN;
else
    McStasStr.inputvalue(BladeNIndex)=11;
end

% Width of an individual blade (m)
BladeWIndex=linp+8;
McStasStr.input{BladeWIndex}=['BladeW' num];
if max(ismember(defaultnames,'BladeW'))
    McStasStr.inputvalue(BladeWIndex)=defaults.BladeW;
else
    McStasStr.inputvalue(BladeWIndex)=0.02;
end

% Height of an individual blade (m)
BladeHIndex=linp+9;
McStasStr.input{BladeHIndex}=['BladeH' num];
if max(ismember(defaultnames,'BladeH'))
    McStasStr.inputvalue(BladeHIndex)=defaults.BladeH;
else
    McStasStr.inputvalue(BladeHIndex)=0.2;
end

% Gap between blades (m)
BladeGapIndex=linp+10;
McStasStr.input{BladeGapIndex}=['BladeGap' num];
if max(ismember(defaultnames,'BladeGap'))
    McStasStr.inputvalue(BladeGapIndex)=defaults.BladeGap;
else
    McStasStr.inputvalue(BladeGapIndex)=0.0005;
end

% Monochromator Binning (0 -> No binning, -1 -> Binning, pos. val. ->
% Sample single Ei with Ei = pos. val. (meV) )
EbinIndex=linp+11;
McStasStr.input{EbinIndex}=['Ebin' num];
if max(ismember(defaultnames,'Ebin'))
    McStasStr.inputvalue(EbinIndex)=defaults.Ebin;
else
    McStasStr.inputvalue(EbinIndex)=-2;
end

% Beam direction flag ('reflect', 'transmit')
if max(ismember(defaultnames,'BeamDir'))
    BeamDir=defaults.BeamDir;
else
    BeamDir='reflect';
end

% Horizontal focusing geometry ('flat', 'lens', 'rowland')
% Note: McStasStr.inputvalue is a 1D array and cannot accept data type char
if max(ismember(defaultnames,'HGeometry'))
    HGeo=defaults.HGeometry;
else
    HGeo='rowland';
end

% Vertical focusing geometry ('flat', 'lens')
% Note: McStasStr.inputvalue is a 1D array and cannot accept data type char
if max(ismember(defaultnames,'VGeometry'))
    VGeo=defaults.VGeometry;
else
    VGeo='lens';
end



% Minimalist Principle calculation at the Monochromator
% (1 -> apply Minimalist Principle, -1 -> Do not apply
if max(ismember(defaultnames,'MinimalistMono'))
    MinimalistMono=defaults.MinimalistMono;
else
    MinimalistMono=-1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% End of Monochromator options with default filling. %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Start of Monochromator options with user filling. %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This sets optional overrides for poor limit/value choices assigned by the
% user. For now they are defined such that the overrides are always inactive.
% ie: set to -1. Note that the -1 value serves a second purpose which is to
% trigger default limits later in the code or keep an optimization
% parameter unlocked.
optimlimits_L1=-ones(1,4); % minHL1, maxHL1, minVL1, maxVL1
locked_L1=-ones(1,2); % HL1, VL1

% loop through any options defined by user
for ops=1:length(options)
    
    % parse option ops into name, value pair.
    for i=1:length(options{ops})
        if strcmp(options{ops}(i),'='); 
            equals=i;
        end
    end
    name_ops = options{ops}(1:equals-1);
    value_ops = options{ops}(equals+1:end);
    
    % Find name_ops element location in McStasStr.inputvalue and assign value_ops
    % This overrides defaults assigned above.
    switch(name_ops)
        case 'L2_override'
            L2_override = str2num(value_ops);
        case 'r0'
           McStasStr.inputvalue(r0Index)=str2num(value_ops);
        case 'Rfile'
            rfile = value_ops;
        case 'dM'
            McStasStr.inputvalue(dMIndex)=str2num(value_ops);
        case 'Vmos'
            McStasStr.inputvalue(VmosIndex)=str2num(value_ops);
        case 'Hmos'
            McStasStr.inputvalue(HmosIndex)=str2num(value_ops);
        case 'mos'
            McStasStr.inputvalue(VmosIndex)=str2num(value_ops);
            McStasStr.inputvalue(HmosIndex)=str2num(value_ops);
        case 'NV'
            McStasStr.inputvalue(NVIndex)=str2num(value_ops);
        case 'SegGap'
            McStasStr.inputvalue(SegGapIndex)=str2num(value_ops);
        case 'BladeN'
            McStasStr.inputvalue(BladeNIndex)=str2num(value_ops);
        case 'BladeW'
            McStasStr.inputvalue(BladeWIndex)=str2num(value_ops);
        case 'BladeH'
            McStasStr.inputvalue(BladeHIndex)=str2num(value_ops);
        case 'BladeGap'
            McStasStr.inputvalue(BladeGapIndex)=str2num(value_ops);
        case 'HGeometry'
            HGeo = value_ops;
        case 'VGeometry'
            VGeo = value_ops;
        case 'Ebin'
            McStasStr.inputvalue(EbinIndex) = str2num(value_ops);
        case 'BeamDir'
            BeamDir = value_ops;
        case 'MinimalistMono'
            MinimalistMono=str2num(value_ops);
        case 'minHL1'
            if optimlimits_L1(1)>0
                optimlimits_L1(1) = max([optimlimits_L1(1) str2num(value_ops)]);
            else
                optimlimits_L1(1) = str2num(value_ops);
            end
        case 'maxHL1'
            if optimlimits_L1(2)>0
                optimlimits_L1(2) = min([optimlimits_L1(2) str2num(value_ops)]);
            else
                optimlimits_L1(2) = str2num(value_ops);
            end
        case 'HL1'
            if locked_L1(1)>0
                locked_L1(1)=min([locked_L1(1) str2num(value_ops)]);
            else
                locked_L1(1)=str2num(value_ops);   
            end
        case 'minVL1'
            if optimlimits_L1(3)>0
                optimlimits_L1(3) = max([optimlimits_L1(3) str2num(value_ops)]);
            else
                optimlimits_L1(3) = str2num(value_ops);
            end
        case 'maxVL1'
            if optimlimits_L1(4)>0
                optimlimits_L1(4) = min([optimlimits_L1(4) str2num(value_ops)]);
            else
                optimlimits_L1(4) = str2num(value_ops);
            end
        case 'VL1'
            if locked_L1(2)>0
                locked_L1(2)=min([locked_L1(1) str2num(value_ops)]);
            else
                locked_L1(2)=str2num(value_ops);   
            end
        case {'minstart', 'maxstart', 'start', 'length'}
            % Do nothing, start is handled by globalinfo structure
        otherwise
            disp(['ERROR, the input ' options{ops} ' is not valid in module M (monochromator)']) 
    end
end

numBlades = McStasStr.inputvalue(BladeNIndex);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% End of monochromator options with user filling %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Add optimization parameters %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define min, max, guess of optimization parameters. If none defined, then
% set with default value embedded in code.

if optimlimits_L1(1)>0
    minHL1=optimlimits_L1(1);
else
    minHL1 = 0.1;
end

if optimlimits_L1(2)>0
    maxHL1=optimlimits_L1(2);
else
    maxHL1 = 20;
end

if ((maxHL1 > demands.Dist) && (minHL1 < demands.Dist))
    guessHL1 = demands.Dist; % Equidistant Rowland Focusing Condition
else
    guessHL1 = (maxHL1 + minHL1)/2;
end


if optimlimits_L1(3)>0
    minVL1=optimlimits_L1(3);
else
    minVL1 = 0.1;
end

if optimlimits_L1(4)>0
    maxVL1=optimlimits_L1(4);
else
    maxVL1 = 20;
end

guessVL1 = (maxVL1 + minVL1)/2;

% Flat mode has no optimization focal length (HL1 or VL1)
if ~strcmp(HGeo, 'flat')
    % Set-up optimization parameters
    linp=length(McStasStr.input);
    McStasStr.input{linp+1}=['HL1_' num];
    if locked_L1(1)<0                       % HL1 is a free parameter
        McStasStr.optimize(linp+1)=1;
        McStasStr.optimvals.min(linp+1)=minHL1;
        McStasStr.optimvals.max(linp+1)=maxHL1;
        McStasStr.optimvals.guess(linp+1)=guessHL1;
    else                                    % HL1 is a fixed parameter
        McStasStr.inputvalue(linp+1)=locked_L1(1);
    end
end

if ~strcmp(VGeo, 'flat')
    % Set-up optimization parameters
    linp=length(McStasStr.input);    
    McStasStr.input{linp+1}=['VL1_' num];
    if locked_L1(2)<0                       % VL1 is a free parameter
        McStasStr.optimize(linp+1)=1;
        McStasStr.optimvals.min(linp+1)=minVL1;
        McStasStr.optimvals.max(linp+1)=maxVL1;
        McStasStr.optimvals.guess(linp+1)=guessVL1;
    else                                    % VL1 is a fixed parameter
        McStasStr.inputvalue(linp+1)=locked_L1(2);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% End of optimization parameters %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Add DECLARE parameters %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

McStasStr.declare{end+1} = ['endPoint' num];
McStasStr.declare{end+1} = ['var_length' num];
McStasStr.declare{end+1} = ['endx' num];
McStasStr.declare{end+1} = ['endy' num];
McStasStr.declare{end+1} = ['startx' num];
McStasStr.declare{end+1} = ['starty' num];
McStasStr.declare{end+1} = ['A1_' num];
McStasStr.declare{end+1} = ['A2_' num];
McStasStr.declare{end+1} = ['Ki' num];
McStasStr.declare{end+1} = ['Vi' num];
McStasStr.declare{end+1} = ['Li' num];
McStasStr.declare{end+1} = ['Lbin' num];
McStasStr.declare{end+1} = ['hmos_rad' num];
McStasStr.declare{end+1} = ['vmos_rad' num];
McStasStr.declare{end+1} = ['mono_phi' num];
McStasStr.declare{end+1} = ['VRstar' num];
McStasStr.declare{end+1} = ['HRstar' num];
McStasStr.declare{end+1} = ['L2_' num];
McStasStr.declare{end+1} = ['var_dM' num];
McStasStr.declare{end+1} = ['var_Ebin' num];
if strcmp(HGeo, 'rowland')
    McStasStr.declare{end+1} = ['A2_max' num];
    McStasStr.declare{end+1} = ['mono_phi_max' num];
end
% Only use focus variables if not flat
if ~strcmp(HGeo, 'flat')
    McStasStr.declare{end+1} = ['var_HL1_' num];
end
if ~strcmp(VGeo, 'flat')
    McStasStr.declare{end+1} = ['var_VL1_' num];
end


McStasStr.declare{end+1} = ['var_binscale' num];
McStasStr.declareint{end+1} = ['bincheck' num];
McStasStr.declareint{end+1} = ['singlebin' num];
McStasStr.declare{end+1} = ['bmin' num];
McStasStr.declare{end+1} = ['bmax' num];
McStasStr.declare{end+1} = ['var_m' num];

% Horizontal Geometry Parameters
if ~strcmp(HGeo, 'lens') % horiz. lens does not use blades.
    % Blade positions and angles for flat, rowland
    for i = 1:numBlades
        McStasStr.declare{end+1} = ['Mpos' num '_' num2str(i)];
        McStasStr.declare{end+1} = ['PHA' num '_' num2str(i)];
    end
end

% These are McStas data types used to dynamically move the monochromator
if isfield(McStasStr, 'declaregen_name')
    McStasStr.declaregen_name{end+1} = ['mctr1_' num];
    McStasStr.declaregen_type{end+1} = 'Rotation';
    McStasStr.declaregen_name{end+1} = ['mctc1_' num];
    McStasStr.declaregen_type{end+1} = 'Coords';
    McStasStr.declaregen_name{end+1} = ['mctc2_' num];
    McStasStr.declaregen_type{end+1} = 'Coords';
else
    McStasStr.declaregen_name{1} = ['mctr1_' num];
    McStasStr.declaregen_type{1} = 'Rotation';
    McStasStr.declaregen_name{end+1} = ['mctc1_' num];
    McStasStr.declaregen_type{end+1} = 'Coords';
    McStasStr.declaregen_name{end+1} = ['mctc2_' num];
    McStasStr.declaregen_type{end+1} = 'Coords';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% End of DECLARE parameters %%%%%%    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Initialize String %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear l
l{1} = ['var_dM' num ' = dM' num ';'];

% Only need focus variable if not flat.
if ~strcmp(HGeo, 'flat')
    l{end+1} = ['var_HL1_' num ' = HL1_' num ';'];
end
if ~strcmp(VGeo, 'flat')
    l{end+1} = ['var_VL1_' num ' = VL1_' num ';'];
end

% Variable for caluclating monochromator binning
l{end+1} = ['var_Ebin' num ' = Ebin' num ';'];

% Monochromator Length
l{end+1} = ['var_length' num ' = length' num ';'];

% Note that for a monochromator endx and endy do not correspond to any
% property of the monochromator. It is defined like other components to
% match the size startx and starty of component num-1. This propogation of
% startx and starty of component num-1 is used in the minimalist 
% calculation of the gap if BeamDir == transmit & MinimalistMono == -1.
if index == 1
    l{end+1} = 'endPoint1 = Mod_sample - sample_dist;';
    l{end+1}=['endx' num ' = sizeX + 2*sample_dist*tan(divreq_x*DEG2RAD);'];
    l{end+1}=['endy' num ' = sizeY + 2*sample_dist*tan(divreq_y*DEG2RAD);'];
else
    l{end+1} = ['endPoint' num ' = endPoint' num2str(index-1) ' - length' num2str(index-1) ' - u;'];
    l{end+1}=['endx' num ' = startx' num2str(index-1) ';'];
    l{end+1}=['endy' num ' = starty' num2str(index-1) ';'];
end

% Define mono to sample distance
if L2_override == -1 % If BeamDir = reflect then this is what you want
    l{end+1} = ['L2_' num ' = Mod_sample - endPoint' num ' + (length' num ')/2;'];
else % If BeamDir = transmit then the monochromator sample distance must be put in by hand.
    l{end+1} = ['L2_' num ' = ' num2str(L2_override) ';'];
end


% startx, starty of the monochromator correspond to its height (starty) and
% its maximal projected width (startx).
if strcmp(HGeo, 'rowland')
    l{end+1} = '// Positive Ebin is a single energy so calculate Rowland tangent angle at that energy.';
    l{end+1} = ['if (Ebin' num ' > 0){'];
    l{end+1} = ['     A2_max' num ' = 2*asin(sqrt(81.81/Ebin' num ')/(2*dM' num '));'];
    l{end+1} = ['     mono_phi_max' num ' = atan2(var_HL1_' num '*sin(A2_max' num '),(var_HL1_' num '*cos(A2_max' num ')+L2_' num '));'];
    l{end+1} = '}';
    l{end+1} = '// Negative or zero Ebin is tougher. Need to find the largest projection over the energy range probed (WaveMin to WaveMax)';
    l{end+1} = 'else {';
    l{end+1} = ['     A2_max' num ' = 2*asin(WaveMax/(2*dM' num '));'];
    % Special exception where Rowland blade group angle (phi) is not
    % monotonic. See Validation powerpoint for explanation of
    % calculation.
    l{end+1} = ['     if (HL1_' num '< L2_' num '){'];
    % Maxima of function phi(A2) is at A2 = acos(-L1/L2)
    l{end+1} = ['          if (A2_max' num ' > acos(-HL1_' num '/L2_' num ')){'];
    l{end+1} = ['               A2_max' num ' = acos(-HL1_' num '/L2_' num ');'];
    % Make sure that if you adjust A2 down to the phi maxima that it is not smaller than the minimum A2 angle.
    l{end+1} = ['               if (A2_max' num '< 2*asin(WaveMin/(2*dM' num '))){'];
    l{end+1} = ['                    A2_max' num ' = 2*asin(WaveMin/(2*dM' num '));'];
    l{end+1} = '     }}}';       
    l{end+1} = '';
    l{end+1} = ['     mono_phi_max' num ' = atan2(var_HL1_' num '*sin(A2_max' num '),(var_HL1_' num '*cos(A2_max' num ')+L2_' num '));'];
    % Projected monochromator width is proportional to sin(phi), thus if
    % phi > pi/2 then the you are actually shrinking the projected width.
    l{end+1} = ['     if (mono_phi_max' num '> PI/2){'];
    l{end+1} = ['          mono_phi_max' num ' = PI/2;'];
    % Make sure that if you adjust A2 down to phi = pi/2 that it is not smaller than the minimum A2 angle.
    l{end+1} = ['               if (2*asin(WaveMin/(2*dM' num ')) > acos(-L2_' num '/HL1_' num ')){'];
    l{end+1} = ['                    A2_max' num ' = acos(-L2_' num '/HL1_' num ');'];
    l{end+1} = ['                    mono_phi_max' num ' = atan2(var_HL1_' num '*sin(A2_max' num '),(var_HL1_' num '*cos(A2_max' num ')+L2_' num '));'];
    l{end+1} = '     }}';
    l{end+1} = '}';
    l{end+1} = '';
    l{end+1} = ['// Use calcuated Rowland angle to find maximum mono projection (startx' num ')'];
    l{end+1} = ['startx' num ' = (BladeW' num '*BladeN' num ' + BladeGap' num '*(BladeN' num '-1))*sin(mono_phi_max' num ');'];
else
    l{end+1} = ['if (Ebin' num ' > 0){'];
    l{end+1} = ['     startx' num ' = (BladeW' num '*BladeN' num ' + BladeGap' num '*(BladeN' num '-1))*sqrt(81.81/Ebin' num ')/(2*dM' num ');'];
    l{end+1} = '}';
    l{end+1} = 'else {';
    l{end+1} = ['     startx' num ' = (BladeW' num '*BladeN' num ' + BladeGap' num '*(BladeN' num '-1))*WaveMax/(2*dM' num ');'];
    l{end+1} = '}';
end
l{end+1} = ['starty' num ' = BladeH' num ';'];


l{end+1} = ['hmos_rad' num ' = Hmos' num '/60.0*PI/180.0;'];
l{end+1} = ['vmos_rad' num ' = Vmos' num '/60.0*PI/180.0;'];

l{end+1} = '';
l{end+1} = '// Variables for calculating monochromator binning';
l{end+1} = '// Only used for dynamic monochromator (Ebin < 0)';
l{end+1} = ['var_binscale' num ' = -1*Ebin' num ';'];
if isfield(defaults,'m')
    l{end+1} = ['var_m' num ' = ' num2str(defaults.m) ';'];
else
    disp('WARNING: Energy binning uses defaults.m value which is currently not defined, using m=2 for determining bin sizes')
    l{end+1} = ['     var_m' num ' = 2;'];
end

l{end+1} = '';
l{end+1} = '// Initialize standard definitions of A1 and A2';
l{end+1} = '// A1 is average wavelength, this is what is used if monochromator is static';
l{end+1} = ['if (Ebin' num ' > 0){  // Use central E value set by Ebin'];
l{end+1} = ['     A1_' num ' = asin(sqrt(81.81/Ebin' num ')/(2*var_dM' num '))*RAD2DEG;'];
l{end+1} = '}';
l{end+1} = 'else {  // Use average wavelength if mono is dynamic (this value is not used since it updates dynamically.) ';
l{end+1} = ['     A1_' num ' = asin(Lambda0/(2*var_dM' num '))*RAD2DEG;'];
l{end+1} = '}';
l{end+1} = ['A2_' num ' = 2*A1_' num ';'];


% Initialize blade positions/angles (lens does not use blades)
if ~strcmp(HGeo, 'lens')
    % Positions
    l{end+1} = '';
    l{end+1} = '/****** MONOCHROMATOR BLADE PLACEMENT ******/';
    for i = 1:numBlades
        l{end+1} = ['Mpos' num '_' num2str(i) ' = ' num2str(i-(numBlades+1)/2) '*(BladeW' num ' + BladeGap' num ');'];
    end
    l{end+1} = '';
    % Angles for rowland (for flat, they are set to zero directly in the COMPONENT)
    if strcmp(HGeo, 'rowland')
        l{end+1} = '/****** MONOCHROMATOR ANGLE PLACEMENT ******/';
        l{end+1} = ['mono_phi' num ' = atan2(var_HL1_' num '*sin(A2_' num '*DEG2RAD),(var_HL1_' num '*cos(A2_' num '*DEG2RAD)+L2_' num '));'];
        for i = 1:numBlades
            l{end+1} = ['PHA' num '_' num2str(i) ' = -(atan2(var_HL1_' num ',(var_HL1_' num '/tan(mono_phi' num ') - Mpos' num '_' num2str(numBlades+1-i) '/sin(mono_phi' num ')))*RAD2DEG - A1_' num ');'];
        end
        l{end+1}='';
    end
end
        
% Initialize lens focusing geometry (lens, flat)
if  strcmp(HGeo, 'lens')
    l{end+1} = ['HRstar' num ' = 2/(1/var_HL1_' num ' + 1/L2_' num ');'];
else
    l{end+1} = ['HRstar' num ' = 0;'];
end

if strcmp(VGeo, 'lens') 
    l{end+1} = ['VRstar' num ' = 2/(1/var_VL1_' num ' + 1/L2_' num ');'];
else
    l{end+1} = ['VRstar' num ' = 0;'];
end
    

% Add initstring to McStasStr.initialize
initstring='';
for i=1:length(l)
    initstring=[initstring l{i} '\n'];
end
clear l;
McStasStr.initialize=[initstring '\n\n' McStasStr.initialize];
McStasStr.initialize_seg=[initstring '\n\n' McStasStr.initialize_seg];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% End of Initialize String %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Trace String %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

% LastBladeIndex is used to properly link the Mono_out component to the last blade
% when dynamic positioning is used. It is initialized to the value
% corresponding to HGeo = 'rowland' or 'flat'.
LastBladeIndex = numBlades;


clear l;

l{1} = ['COMPONENT calculation_arm' num ' = Arm()'];
l{end+1} = 'AT (0,0,u) RELATIVE PREVIOUS';
l{end+1} = 'ROTATED (0,0,0) RELATIVE PREVIOUS';

% Set up dynamic monochromator using an EXTEND block

l{end+1} = 'EXTEND %%{';
l{end+1} = ['if (var_Ebin' num ' <= 0){   // dynamic mono'];      
l{end+1} = ['     Vi' num ' = sqrt(vx*vx + vy*vy + vz*vz);'];
l{end+1} = ['     Ki' num ' = V2K*Vi' num ';'];
l{end+1} = ['     Li' num ' = 2*PI/Ki' num ';'];
l{end+1} = '';
    
% Set up energy binning
l{end+1} = ['     if (var_Ebin' num ' < 0){  // Binning'];    
l{end+1} = '          // For a given Li, set mono to corresponding bin wavelength Lbin';
l{end+1} = ['          bincheck' num ' = 0;'];
l{end+1} = ['          singlebin' num ' = 0;'];
l{end+1} = ['          bmin' num ' = Lambda0 - dLambda - 0.0001;'];
l{end+1} = ['          bmax' num ' = bmin' num ' + var_binscale' num '*2*var_dM' num '*(hmos_rad' num ' + var_m' num '*bmin' num '*.0024)/2*cos(asin(bmin' num '/2/var_dM' num '));'];
l{end+1} = ['          if((bmax' num ' - bmin' num ') > (2*dLambda)){'];
l{end+1} = ['               Lbin' num ' = Lambda0;'];
l{end+1} = ['               singlebin' num ' = 1;'];
l{end+1} = ['               bincheck' num ' = 1;'];
l{end+1} = '          }';
l{end+1} = ['          while ((bmin' num ' <= (Lambda0 + dLambda)) && ((Li' num '/2/var_dM' num ') <= 1) && (singlebin' num ' == 0)) {'];
l{end+1} = ['               if ((bmin' num ' <= Li' num ') && (bmax' num ' >= Li' num ')) {'];
l{end+1} = ['                    Lbin' num ' = (bmin' num ' + bmax' num ')/2;'];
l{end+1} = ['                    bincheck' num ' = 1;'];
l{end+1} =  '                    break;';
l{end+1} =  '               }';
l{end+1} =  '               else {';
l{end+1} = ['                    bmin' num ' = bmax' num ';'];
l{end+1} = ['                    bmax' num ' = bmax' num ' +  var_binscale' num '*2*var_dM' num '*(hmos_rad' num ' + var_m' num '*bmin' num '*.0024)/2*cos(asin(bmin' num '/2/var_dM' num '));'];
l{end+1} =  '               }';
l{end+1} = '          }';            
l{end+1} = '';
l{end+1} = ['          if (bincheck' num ' == 0) {'];
l{end+1} = ['               if(Li' num '/2/var_dM' num ' > 1) {'];
l{end+1} = ['                    Lbin' num ' = var_dM' num ';'];
l{end+1} = '               }';
l{end+1} = '               else {';
l{end+1} = ['                    printf("No Bin Found for Wavelength = %%0.4g\\n", Li' num ');'];
l{end+1} = '               }';
l{end+1} = '          }';
l{end+1} = '';
l{end+1} = '     }';
l{end+1} =  '     else {  // No binning';
l{end+1} = ['          Lbin' num ' = Li' num ';'];
l{end+1} = '     }';
    
% Standard definitions of A1 and A2
l{end+1} = ['     A1_' num ' = asin(Lbin' num '/(2*var_dM' num '))*RAD2DEG;'];
l{end+1} = ['     A2_' num ' = 2*A1_' num ';'];
    
% Select horizontal focusing geometry (rowland, lens, flat)
switch (HGeo)
    case 'rowland'        
        % Set up Monochromator table angles            
        l{end+1} = ['     mono_phi' num ' = atan2(var_HL1_' num '*sin(A2_' num '*DEG2RAD),(var_HL1_' num '*cos(A2_' num '*DEG2RAD)+L2_' num '));'];
        l{end+1} = '';

        % Set up blade angles
        l{end+1} = '     /* Define blade angles */';
        for i = 1:numBlades
            l{end+1} = ['     PHA' num '_' num2str(i) ' = -(atan2(var_HL1_' num ',(var_HL1_' num '/tan(mono_phi' num ') - Mpos' num '_' num2str(numBlades+1-i) '/sin(mono_phi' num ')))*RAD2DEG - A1_' num ');'];
        end
        l{end+1} = '';

        % Dynamically position the Monocromator 
        % First  position the table
        l{end+1} = ['     /* Component mono_arm' num ' (Table Rotation) */'];
        l{end+1} = ['     rot_set_rotation(mctr1_' num ', (0)*DEG2RAD, mono_phi' num ', (0)*DEG2RAD);'];
        l{end+1} = ['     rot_mul(mctr1_' num ', mcrotacalculation_arm' num ', mcrotamono_arm' num ');'];
        l{end+1} = ['     rot_transpose(mcrotacalculation_arm' num ', mctr1_' num ');'];
        l{end+1} = ['     rot_mul(mcrotamono_arm' num ', mctr1_' num ', mcrotrmono_arm' num ');'];
        l{end+1} = ['     mctc1_' num ' = coords_set(0, 0, var_length' num '/2);'];
        l{end+1} = ['     rot_transpose(mcrotacalculation_arm' num ', mctr1_' num ');'];
        l{end+1} = ['     mctc2_' num ' = rot_apply(mctr1_' num ', mctc1_' num ');'];
        l{end+1} = ['     mcposamono_arm' num ' = coords_add(mcposacalculation_arm' num ', mctc2_' num ');'];
        l{end+1} = ['     mctc1_' num ' = coords_sub(mcposacalculation_arm' num ', mcposamono_arm' num ');'];
        l{end+1} = ['     mcposrmono_arm' num ' = rot_apply(mcrotamono_arm' num ', mctc1_' num ');'];
        l{end+1} = '';
        l{end+1} = '     /* Not sure what McStas uses the array version of the transformations for. */';
        l{end+1} = '     /* It is hard to get the correct index values in guide_bot and so the array update is commented out */';
        l{end+1} = '     /* Only the non-array transformation appear to be used but beware. If so desired you can count up to each component */';
        l{end+1} = '     /* and put in the index number by hand. You can also check the .c file to see what number McStas was writing in. */';
        l{end+1} = '';
        l{end+1} = ['     /* mccomp_posa[?] = mcposamono_arm' num '; */'];
        l{end+1} = ['     /* mccomp_posr[?] = mcposrmono_arm' num '; */'];
        l{end+1} = '';
        
        % Need to place the first blade without using a loop
        l{end+1} = ['     /* Component mono' num '_1 (Blade 1) */'];
        l{end+1} = ['     rot_set_rotation(mctr1_' num ', (0)*DEG2RAD, (PHA' num '_1)*DEG2RAD, (0)*DEG2RAD);'];
        l{end+1} = ['     rot_mul(mctr1_' num ', mcrotamono_arm' num ', mcrotamono' num '_1);'];
        l{end+1} = ['     rot_transpose(mcrotamono_arm' num ', mctr1_' num ');'];
        l{end+1} = ['     rot_mul(mcrotamono' num '_1, mctr1_' num ', mcrotrmono' num '_1);'];
        l{end+1} = ['     mctc1_' num ' = coords_set(0, 0, Mpos' num '_1);'];
        l{end+1} = ['     rot_transpose(mcrotamono_arm' num ', mctr1_' num ');'];
        l{end+1} = ['     mctc2_' num ' = rot_apply(mctr1_' num ', mctc1_' num ');'];
        l{end+1} = ['     mcposamono' num '_1 = coords_add(mcposamono_arm' num ', mctc2_' num ');'];
        l{end+1} = ['     mctc1_' num ' = coords_sub(mcposamono_arm' num ', mcposamono' num '_1);'];
        l{end+1} = ['     mcposrmono' num '_1 = rot_apply(mcrotamono' num '_1, mctc1_' num ');'];
        l{end+1} = ['     /* mccomp_posa[?+1] = mcposamono' num '_1; */'];
        l{end+1} = ['     /* mccomp_posr[?+1] = mcposrmono' num '_1; */'];
        l{end+1} = '';

        % Now position all other blades with a loop
        for i = 2:numBlades
            l{end+1} = ['     /* Component mono' num '_' num2str(i) ' (Blade ' num2str(i) ') */'];
            l{end+1} = ['     rot_set_rotation(mctr1_' num ', (0)*DEG2RAD, (PHA' num '_' num2str(i) ')*DEG2RAD, (0)*DEG2RAD);'];
            l{end+1} = ['     rot_mul(mctr1_' num ', mcrotamono_arm' num ', mcrotamono' num '_' num2str(i) ');'];
            l{end+1} = ['     rot_transpose(mcrotamono' num '_' num2str(i-1) ', mctr1_' num ');'];
            l{end+1} = ['     rot_mul(mcrotamono' num '_' num2str(i) ', mctr1_' num ', mcrotrmono' num '_' num2str(i) ');'];
            l{end+1} = ['     mctc1_' num ' = coords_set(0, 0, Mpos' num '_' num2str(i) ');'];
            l{end+1} = ['     rot_transpose(mcrotamono_arm' num ', mctr1_' num ');'];
            l{end+1} = ['     mctc2_' num ' = rot_apply(mctr1_' num ', mctc1_' num ');'];
            l{end+1} = ['     mcposamono' num '_' num2str(i) ' = coords_add(mcposamono_arm' num ', mctc2_' num ');'];
            l{end+1} = ['     mctc1_' num ' = coords_sub(mcposamono' num '_' num2str(i-1) ', mcposamono' num '_' num2str(i) ');'];
            l{end+1} = ['     mcposrmono' num '_' num2str(i) ' = rot_apply(mcrotamono' num '_' num2str(i) ', mctc1_' num ');'];
            l{end+1} = ['     /* mccomp_posa[?+' num2str(i) '] = mcposamono' num '_' num2str(i) '; */'];
            l{end+1} = ['     /* mccomp_posr[?+' num2str(i) '] = mcposrmono' num '_' num2str(i) '; */'];
            l{end+1} = '';
        end
        l{end+1} = '';
        
    case 'flat'
        % First  position the table
        l{end+1} = ['     /* Component mono_arm' num ' (Table Rotation) */'];
        l{end+1} = ['     rot_set_rotation(mctr1_' num ', (0)*DEG2RAD, A1_' num '*DEG2RAD, (0)*DEG2RAD);'];
        l{end+1} = ['     rot_mul(mctr1_' num ', mcrotacalculation_arm' num ', mcrotamono_arm' num ');'];
        l{end+1} = ['     rot_transpose(mcrotacalculation_arm' num ', mctr1_' num ');'];
        l{end+1} = ['     rot_mul(mcrotamono_arm' num ', mctr1_' num ', mcrotrmono_arm' num ');'];
        l{end+1} = ['     mctc1_' num ' = coords_set(0, 0, var_length' num '/2);'];
        l{end+1} = ['     rot_transpose(mcrotacalculation_arm' num ', mctr1_' num ');'];
        l{end+1} = ['     mctc2_' num ' = rot_apply(mctr1_' num ', mctc1_' num ');'];
        l{end+1} = ['     mcposamono_arm' num ' = coords_add(mcposacalculation_arm' num ', mctc2_' num ');'];
        l{end+1} = ['     mctc1_' num ' = coords_sub(mcposacalculation_arm' num ', mcposamono_arm' num ');'];
        l{end+1} = ['     mcposrmono_arm' num ' = rot_apply(mcrotamono_arm' num ', mctc1_' num ');'];
        l{end+1} = '';
        l{end+1} = '     /* Not sure what McStas uses the array version of the transformations for. */';
        l{end+1} = '     /* It is hard to get the correct index values in guide_bot and so the array update is commented out */';
        l{end+1} = '     /* Only the non-array transformation appear to be used but beware. If so desired you can count up to each component */';
        l{end+1} = '     /* and put in the index number by hand. You can also check the .c file to see what number McStas was writing in. */';
        l{end+1} = '';
        l{end+1} = ['     /* mccomp_posa[?] = mcposamono_arm' num '; */'];
        l{end+1} = ['     /* mccomp_posr[?] = mcposrmono_arm' num '; */'];
        l{end+1} = '';
        
        % Need to place the first blade without using a loop
        l{end+1} = ['     /* Component mono' num '_1 (Blade 1) */'];
        l{end+1} = ['     rot_set_rotation(mctr1_' num ', (0)*DEG2RAD, 0*DEG2RAD, (0)*DEG2RAD);'];
        l{end+1} = ['     rot_mul(mctr1_' num ', mcrotamono_arm' num ', mcrotamono' num '_1);'];
        l{end+1} = ['     rot_transpose(mcrotamono_arm' num ', mctr1_' num ');'];
        l{end+1} = ['     rot_mul(mcrotamono' num '_1, mctr1_' num ', mcrotrmono' num '_1);'];
        l{end+1} = ['     mctc1_' num ' = coords_set(0, 0, Mpos' num '_1);'];
        l{end+1} = ['     rot_transpose(mcrotamono_arm' num ', mctr1_' num ');'];
        l{end+1} = ['     mctc2_' num ' = rot_apply(mctr1_' num ', mctc1_' num ');'];
        l{end+1} = ['     mcposamono' num '_1 = coords_add(mcposamono_arm' num ', mctc2_' num ');'];
        l{end+1} = ['     mctc1_' num ' = coords_sub(mcposamono_arm' num ', mcposamono' num '_1);'];
        l{end+1} = ['     mcposrmono' num '_1 = rot_apply(mcrotamono' num '_1, mctc1_' num ');'];
        l{end+1} = ['     /* mccomp_posa[?+1] = mcposamono' num '_1; */'];
        l{end+1} = ['     /* mccomp_posr[?+1] = mcposrmono' num '_1; */'];
        l{end+1} = '';

        % Now position all other blades with a loop
        for i = 2:numBlades
            l{end+1} = ['     /* Component mono' num '_' num2str(i) ' (Blade ' num2str(i) ') */'];
            l{end+1} = ['     rot_set_rotation(mctr1_' num ', (0)*DEG2RAD, (0)*DEG2RAD, (0)*DEG2RAD);'];
            l{end+1} = ['     rot_mul(mctr1_' num ', mcrotamono_arm' num ', mcrotamono' num '_' num2str(i) ');'];
            l{end+1} = ['     rot_transpose(mcrotamono' num '_' num2str(i-1) ', mctr1_' num ');'];
            l{end+1} = ['     rot_mul(mcrotamono' num '_' num2str(i) ', mctr1_' num ', mcrotrmono' num '_' num2str(i) ');'];
            l{end+1} = ['     mctc1_' num ' = coords_set(0, 0, Mpos' num '_' num2str(i) ');'];
            l{end+1} = ['     rot_transpose(mcrotamono_arm' num ', mctr1_' num ');'];
            l{end+1} = ['     mctc2_' num ' = rot_apply(mctr1_' num ', mctc1_' num ');'];
            l{end+1} = ['     mcposamono' num '_' num2str(i) ' = coords_add(mcposamono_arm' num ', mctc2_' num ');'];
            l{end+1} = ['     mctc1_' num ' = coords_sub(mcposamono' num '_' num2str(i-1) ', mcposamono' num '_' num2str(i) ');'];
            l{end+1} = ['     mcposrmono' num '_' num2str(i) ' = rot_apply(mcrotamono' num '_' num2str(i) ', mctc1_' num ');'];
            l{end+1} = ['     /* mccomp_posa[?+' num2str(i) '] = mcposamono' num '_' num2str(i) '; */'];
            l{end+1} = ['     /* mccomp_posr[?+' num2str(i) '] = mcposrmono' num '_' num2str(i) '; */'];
            l{end+1} = '';
        end
        l{end+1} = '';        
          
    case 'lens'        
        % Now  position the table
        l{end+1} = ['     /* Component mono_arm' num ' (Table Rotation) */'];
        l{end+1} = ['     rot_set_rotation(mctr1_' num ', (0)*DEG2RAD, A1_' num '*DEG2RAD, (0)*DEG2RAD);'];
        l{end+1} = ['     rot_mul(mctr1_' num ', mcrotacalculation_arm' num ', mcrotamono_arm' num ');'];
        l{end+1} = ['     rot_transpose(mcrotacalculation_arm' num ', mctr1_' num ');'];
        l{end+1} = ['     rot_mul(mcrotamono_arm' num ', mctr1_' num ', mcrotrmono_arm' num ');'];
        l{end+1} = ['     mctc1_' num ' = coords_set(0, 0, var_length' num '/2);'];
        l{end+1} = ['     rot_transpose(mcrotacalculation_arm' num ', mctr1_' num ');'];
        l{end+1} = ['     mctc2_' num ' = rot_apply(mctr1_' num ', mctc1_' num ');'];
        l{end+1} = ['     mcposamono_arm' num ' = coords_add(mcposacalculation_arm' num ', mctc2_' num ');'];
        l{end+1} = ['     mctc1_' num ' = coords_sub(mcposacalculation_arm' num ', mcposamono_arm' num ');'];
        l{end+1} = ['     mcposrmono_arm' num ' = rot_apply(mcrotamono_arm' num ', mctc1_' num ');'];
        l{end+1} = '';
        l{end+1} = '     /* Not sure what McStas uses the array version of the transformations for. */';
        l{end+1} = '     /* It is hard to get the correct index values in guide_bot and so the array update is commented out */';
        l{end+1} = '     /* Only the non-array transformation appear to be used but beware. If so desired you can count up to each component */';
        l{end+1} = '     /* and put in the index number by hand. You can also check the .c file to see what number McStas was writing in. */';
        l{end+1} = '';
        l{end+1} = ['     /* mccomp_posa[?] = mcposamono_arm' num '; */'];
        l{end+1} = ['     /* mccomp_posr[?] = mcposrmono_arm' num '; */'];
        l{end+1} = '';
          
        % Place Monochromator
        l{end+1} = ['     /* Component mono' num '_1 (Blade 1) */'];
        l{end+1} = ['     rot_set_rotation(mctr1_' num ', (0)*DEG2RAD, (0)*DEG2RAD, (0)*DEG2RAD);'];
        l{end+1} = ['     rot_mul(mctr1_' num ', mcrotamono_arm' num ', mcrotamono' num '_1);'];
        l{end+1} = ['     rot_transpose(mcrotamono_arm' num ', mctr1_' num ');'];
        l{end+1} = ['     rot_mul(mcrotamono' num '_1, mctr1_' num ', mcrotrmono' num '_1);'];
        l{end+1} = ['     mctc1_' num ' = coords_set(0, 0, 0);'];
        l{end+1} = ['     rot_transpose(mcrotamono_arm' num ', mctr1_' num ');'];
        l{end+1} = ['     mctc2_' num ' = rot_apply(mctr1_' num ', mctc1_' num ');'];
        l{end+1} = ['     mcposamono' num '_1 = coords_add(mcposamono_arm' num ', mctc2_' num ');'];
        l{end+1} = ['     mctc1_' num ' = coords_sub(mcposamono_arm' num ', mcposamono' num '_1);'];
        l{end+1} = ['     mcposrmono' num '_1 = rot_apply(mcrotamono' num '_1, mctc1_' num ');'];
        l{end+1} = ['     /* mccomp_posa[?+1] = mcposamono' num '_1; */'];
        l{end+1} = ['     /* mccomp_posr[?+1] = mcposrmono' num '_1; */'];
        l{end+1} = '';
          
        % LastBladeIndex is used to properly link the next component,
        % Mono_out, to the last monochromator blade, for lens there is
        % only one blade.
        LastBladeIndex = 1;
  
    otherwise
        disp('Invalid Horizontal Monochromator Geometry')
end
    
% Finally position the beam direction after the monochromator.    
l{end+1} = ['     /* Component Mono_out' num ' (Beam Direction)*/'];
if strcmp(BeamDir, 'reflect') % A2 changes for each Ei        
    l{end+1} = ['     rot_set_rotation(mctr1_' num ', (0)*DEG2RAD, (A2_' num ')*DEG2RAD,(0)*DEG2RAD);'];
    l{end+1} = ['     rot_mul(mctr1_' num ', mcrotacalculation_arm' num ', mcrotaMono_out' num ');'];
    l{end+1} = ['     rot_transpose(mcrotamono' num '_' num2str(LastBladeIndex) ', mctr1_' num ');'];
    l{end+1} = ['     rot_mul(mcrotaMono_out' num ', mctr1_' num ', mcrotrMono_out' num ');'];
    l{end+1} = ['     mctc1_' num ' = coords_set(0, 0, var_length' num '/2);'];
    l{end+1} = ['     rot_transpose(mcrotacalculation_arm' num ', mctr1_' num ');'];
    l{end+1} = ['     mctc2_' num ' = rot_apply(mctr1_' num ', mctc1_' num ');'];
    l{end+1} = ['     mcposaMono_out' num ' = coords_add(mcposacalculation_arm' num ', mctc2_' num ');'];
    l{end+1} = ['     mctc1_' num ' = coords_sub(mcposamono' num '_' num2str(LastBladeIndex) ', mcposaMono_out' num ');'];
    l{end+1} = ['     mcposrMono_out' num ' = rot_apply(mcrotaMono_out' num ', mctc1_' num ');'];
    l{end+1} = ['     /* mccomp_posa[?+' num2str(i+1) '] = mcposaMono_out' num '; */'];
    l{end+1} = ['     /* mccomp_posr[?+' num2str(i+1) '] = mcposrMono_out' num '; */'];
elseif strcmp(BeamDir, 'transmit')  % A2 = 0 for all Ei, but still must dynamically update (not sure why...)
    l{end+1} = ['     rot_set_rotation(mctr1_' num ', (0)*DEG2RAD, (0)*DEG2RAD,(0)*DEG2RAD);'];
    l{end+1} = ['     rot_mul(mctr1_' num ', mcrotacalculation_arm' num ', mcrotaMono_out' num ');'];
    l{end+1} = ['     rot_transpose(mcrotamono' num '_' num2str(LastBladeIndex) ', mctr1_' num ');'];
    l{end+1} = ['     rot_mul(mcrotaMono_out' num ', mctr1_' num ', mcrotrMono_out' num ');'];
    l{end+1} = ['     mctc1_' num ' = coords_set(0, 0, var_length' num '/2);'];
    l{end+1} = ['     rot_transpose(mcrotacalculation_arm' num ', mctr1_' num ');'];
    l{end+1} = ['     mctc2_' num ' = rot_apply(mctr1_' num ', mctc1_' num ');'];
    l{end+1} = ['     mcposaMono_out' num ' = coords_add(mcposacalculation_arm' num ', mctc2_' num ');'];
    l{end+1} = ['     mctc1_' num ' = coords_sub(mcposamono' num '_' num2str(LastBladeIndex) ', mcposaMono_out' num ');'];
    l{end+1} = ['     mcposrMono_out' num ' = rot_apply(mcrotaMono_out' num ', mctc1_' num ');'];
    l{end+1} = ['     /* mccomp_posa[?+' num2str(i+1) '] = mcposaMono_out' num '; */'];
    l{end+1} = ['     /* mccomp_posr[?+' num2str(i+1) '] = mcposrMono_out' num '; */'];
else
    disp('WARNING: BeamDir must be set to transmit or reflect!')
    return
end
l{end+1} = '';                         
l{end+1} = '}'; 
l{end+1} = '%%}'; % End of EXTEND module for dynamic positioning 
l{end+1} = '';   



% Setup COMPONENT definitions

% COMPONONT definition is written different for each horizontal geometry.
switch (HGeo)
    case 'rowland'
        % Table
        l{end+1} = ['COMPONENT mono_arm' num ' = Arm()'];
        l{end+1} = ['  AT (0, 0, (length' num ')/2) RELATIVE calculation_arm' num];
        l{end+1} = ['  ROTATED (0, mono_phi' num '*RAD2DEG, 0) RELATIVE calculation_arm' num];
        l{end+1} = '';
        
        % Blades
        for i = 1:numBlades
            l{end+1} = ['COMPONENT mono' num '_' num2str(i) ' = Monochromator_curved_dynamic(DM = dM' num ', mosaich = Hmos' num ', mosaicv = Vmos' num ','];
            l{end+1} = ['width = BladeW' num ', height = BladeH' num ', gap = SegGap' num ', binning = Ebin' num ', WaveMin = WaveMin, WaveMax = WaveMax, RstarV=VRstar' num ', RstarH=HRstar' num ', NV = NV' num ', NH = 1, r0 = r0' num ', reflect = "' rfile '", binscale = var_binscale' num ', mvalue = var_m' num ')'];
            l{end+1} = ['AT (0, 0, Mpos' num '_' num2str(i) ') RELATIVE mono_arm' num];
            l{end+1} = ['ROTATED (0, PHA' num '_' num2str(i) ', 0) RELATIVE mono_arm' num];
            l{end+1} = '';
        end

        McStasStr.declare{end+1} = ['Mpos' num '_' num2str(i)];
        McStasStr.declare{end+1} = ['PHA' num '_' num2str(i)];
        % Beam direction after monochromator
    case 'flat'
        % Table
        l{end+1} = ['COMPONENT mono_arm' num ' = Arm()'];
        l{end+1} = ['  AT(0, 0, (length' num ')/2) RELATIVE calculation_arm' num];
        l{end+1} = ['  ROTATED (0, A1_' num ', 0) RELATIVE calculation_arm' num];
        l{end+1} = '';
        
        % Blades
        for i = 1:numBlades
            l{end+1} = ['COMPONENT mono' num '_' num2str(i) ' = Monochromator_curved_dynamic(DM = dM' num ', mosaich = Hmos' num ', mosaicv = Vmos' num ','];
            l{end+1} = ['width = BladeW' num ', height = BladeH' num ', gap = SegGap' num ', binning = Ebin' num ', WaveMin = WaveMin, WaveMax = WaveMax, RstarV=VRstar' num ', RstarH=HRstar' num ', NV = NV' num ', NH = 1, r0 = r0' num ', reflect = "' rfile '", binscale = var_binscale' num ', mvalue = var_m' num ')'];
            l{end+1} = ['AT (0, 0, Mpos' num '_' num2str(i) ') RELATIVE mono_arm' num];
            l{end+1} = ['ROTATED (0, 0, 0) RELATIVE mono_arm' num];
            l{end+1} = 'GROUP MONO';
            l{end+1} = '';
        end
    case 'lens'
        % Table
        l{end+1} = ['COMPONENT mono_arm' num ' = Arm()'];
        l{end+1} = ['  AT(0, 0, (length' num ')/2) RELATIVE calculation_arm' num];
        l{end+1} = ['  ROTATED (0, A1_' num ', 0) RELATIVE calculation_arm' num];
        l{end+1} = '';
        
        % Monochromator
        l{end+1} = ['COMPONENT mono' num '_1 = Monochromator_curved_dynamic(DM = dM' num ', mosaich = Hmos' num ', mosaicv = Vmos' num ','];
        l{end+1} = ['width = length' num ', height = BladeH' num ', gap = SegGap' num ', binning = Ebin' num ', WaveMin = WaveMin, WaveMax = WaveMax, RstarV=VRstar' num ', RstarH=HRstar' num ', NH = BladeN' num ', NV = NV' num ', r0 = r0' num ', reflect = "' rfile '", binscale = var_binscale' num ', mvalue = var_m' num ')'];
        l{end+1} = ['AT (0, 0, 0) RELATIVE mono_arm' num];
        l{end+1} = ['ROTATED (0, 0, 0) RELATIVE mono_arm' num];
        l{end+1} = '';
    otherwise
        disp('Invalid Horizontal Monochromator Geometry')
end


% Define beam direction after monochromator
if strcmp(BeamDir, 'reflect')
    l{end+1} = ['COMPONENT Mono_out' num ' = Arm()'];
    l{end+1} = ['  AT (0, 0, (length' num ')/2) RELATIVE calculation_arm' num];
    l{end+1} = ['  ROTATED (0, A2_' num ', 0) RELATIVE calculation_arm' num];
else
    l{end+1} = ['COMPONENT Mono_out' num ' = Arm()'];
    l{end+1} = ['  AT (0, 0, (length' num ')/2) RELATIVE calculation_arm' num];
    l{end+1} = ['  ROTATED (0, 0, 0) RELATIVE calculation_arm' num];
end

% End of element COMPONENT
l{end+1}='';
l{end+1}=['COMPONENT EndOfelement_' num '= Arm()'];
l{end+1}=['AT (0,0,length' num '/2) RELATIVE PREVIOUS'];
l{end+1}='';


% Add trace string to McStasStr.trace
tracestring='';
for i=1:length(l)
    tracestring=[tracestring l{i} '\n'];
end
McStasStr.trace=[McStasStr.trace '\n\n' tracestring];
McStasStr.trace_seg=[McStasStr.trace_seg '\n\n' tracestring];
clear l;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%%%%%% End of Trace String %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




