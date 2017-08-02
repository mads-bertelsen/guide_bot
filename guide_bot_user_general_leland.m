clear; clc;

% NAME OF PROJECT:
name='Ebinneg2_mos45';

addpath('..')
% --- Input / demands / figure of merit ---
addpath(genpath('/home/lwh1/iFit/ifit-1.8-4'));
% Horizontal divergence [deg (+/-)]
demands.Hdiv=3;
% Vertical divergence [deg (+/-)]
demands.Vdiv=3;

%options.locked.Hdiv='Vdiv';

% Sample size horizontal [cm]
demands.Hsize=2;
% Sample size vertical [cm]
demands.Vsize=2;
% Lowest wavelength needed [Å]
demands.WaveLmin=2.335;
%demands.WaveLmin=3.96645;
% Highest wavelength needed [Å]
demands.WaveLmax=5.72;
%demands.WaveLmax=4.12841;

% Distance between guide end and sample [m]
% Distance between moderator and sample [m]
demands.Mod_sample=43.2097; demands.Dist=1.6; % Standard
% demands.Mod_sample=45.2097; demands.Dist=0.3;

% --- requirements made by facility / external factors ---

% The minimum distance between moderator and guide allowed
requirements.closest_element=1.605;
% The maximum distance between moderator and guide allowed
requirements.latest_start=1.605;
% Horizontal moderator size (rectangular)
requirements.moderator_size_x=0.15;
% Vertical moderator size
requirements.moderator_size_y=0.2;
% Source spectrum (ESS), cold, thermal or bispectral
% Notice this does not have an impact on the optimization
% The bispectral option is the maximum of the cold and thermal at all wavelengths
requirements.source='Source_gen';
options.source_gen_file = 'components/NCNR_coldsource.dat';

options.optimizer_mode = 'realistic_source';

% --- options which affects the optimization ---
% Controll how much should be calculated by the minimalist concept
%  1: Extraction, gaps and exit
%  0: Gaps and exit
% -1: Exit
% The lower the number, the more degrees of freedom for the optimizer
options.minimalist=-1;
% Option to include figures on absolute flux performance on ESS
% 1: Enable. 0:Disable.
options.absolute_intensity_run = 1;
% Beamport, 0 is the center, goes from -30 to 30 deg.
options.beamport=-13.5;
% Orthogonal distance from moderator face to LOS origin of guide. [m]
% Only used for 'Source_gen' with a non-zero beamport angle. 
options.Lc = 0.02895;

% When a monochromator in reflection mode is present, then a figure of
% merit table is generated during analysis, along with some associated
% plots. It can be useful to compare this to some other instument (ie: a
% baseline). To do this, just copy the baselines FOM table into the 
% guide_bot main directory (for instance: /home/username/guide_bot_v1.02)
% and set it to options.baseline
% Remember to use a baseline with same brilliance window!
options.baseline = 'Current_SPINS.mat';

%%%%%%%%%%%%% Intermediate Brilliance with Custom Window %%%%%%%%%%%%%%

% Generate scripts for analyzing brilliance at guide locations other than
% just the sample. Guide locations are identified by component index. For
% instance SEGMS has indicies S=5 E=4 G=3 M=2 S=1. Thus:
% options.Intermediate_Brilliance = [2 3];
% will produce .inst and ifit scripts at the start position of M and G.
% Note that the value 0 can be used to produce analysis scripts at the
% sample position (which can be useful for analyzing different
% brilliance windows other than the FOM)
options.Intermediate_Brilliance = [0 1 2 3 4];

% The Billiance_Window option only applies to Intermedieate_Brilliance
% analysis scripts. This allows the user to define a brilliance window other 
% than the FOM used to run the optimization at the sample position. There
% are three options:
% 'FOM' :       Set to the sample FOM used for optimization
% 'component' : Set to the optimized start value of the component where the 
%               intermediate brilliance is to be calculated. Only valid for
%               Hsize and Vsize.
% value :       Set to user defined value; Hsize and Vsize [cm],
%               Hdiv and Vdiv [deg].
options.Brilliance_Window.Hsize = 'component';
options.Brilliance_Window.Vsize = 'component';
options.Brilliance_Window.Hdiv = 1;
options.Brilliance_Window.Vdiv = 1;

% HOW TO USE: 
% After the optimization is done, navigate to the instrument
% folder of interest, in the above example this would be:
% ~/guide_bot_v1.02/<project_name>/SEGMS
% The intermediate brilliance ifit files are written here and follow the
% naming convention <instrument>_<module><component index>_ifit.m
% Thus in the above example, the names would be:
% SEGMS_G_module3_ifit.m
% SEGMS_M_module2_ifit.m
% Simply run the one of interest. (either directly or with the .bat file) 
% It will automatically write a brilliance ifit file at the moderator for 
% normalization and run this in tandem.
% When done, all results are sent to the analysis folder and a
% similar naming convention is used to write an analyze_all_ifit.m file:
% G_module_SEGMS1_analyze_all_ifit.m
% E_module_SEGMS1_analyze_all_ifit.m
% Run these to produce analysis plots

% TIP 1:
% If an intermediate brilliance analysis is needed for an optimization, but
% it was not set up to produce the ifit files at the time, then this is not
% a problem. Just move the guide_bot record copy from the project folder to
% the guide_bot main directory, modify the brilliance options you want and
% then rename and rerun. You can then copy the data and parameter files over from the
% optimization (eg: SEGMS1_all.mat, and SEGMS1_geometry.par). This is the 
% only external input needed by the intermediate ifit files to perform the analysis.

% TIP 2:
% If you have intermediate brilliance files and the corresponding optimization
% in place but you would like to look at different brilliance windows, then
% you can do this by opening the intermediate brilliance ifit file and
% changing the values of Hsize, Vsize, Hdiv, and Vdiv under the % primary
% parameters section. Then rerun. Be careful that you don't copy over other
% intermediate brilliance runs (you should rename them first).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% --- computing options ---
options.cluster = 'NCNR';

% Allow project folders to be run back to back. 
% Namely, all instruments in a project folder are run simultaneously. Once
% they have all completed then the next project in the project list is run.
% Note, this requires valid launch_all.sh scripts to work.
% Valid values are as follows:
% 'new':     delete old project list and add this project as the first new one.
% 'add':     add this project to the current project list
% 'single':  make this project but do not add it to the project list
options.projectlist = 'new'; % ('new', 'single', or 'add')

%options.queue = 'long'; 
%options.enable_nested_compile = 1;
options.mpi=100;
options.machines='/home/lwh1/machines';

% --- defaults ---
% mainly used for reflectivity options R0,Qc,m,alpha,W.
defaults.m=4;

% setting defaults.alpha = 0 and defaults.W = 0 will have mcstas pick the
% correct values associated with the defaults.m value chosen. (These are
% taken from Swiss Neutrons and paper by Henrik Jacobson)
defaults.alpha = 0;
defaults.W = 0;

clear MStruct

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Monochromator Parameters %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MStruct.start = 41.3052;        % Fixed start position (m)
% MStruct.minstart = 10;        % Lower limit if start is fit (m)
% MStruct.maxstart = 25;        % Upper limit if start is fit (m)
% MStruct.r0 = 0.7;             % Reflectivity constanst [0-1]
MStruct.Rfile = '../../data/HOPG.rfl';     % Reflectivity file
MStruct.dM = 3.354;             % Monochromator d-spacing (Ang)
% MStruct.Vmos = 60;            % Monochromator vertical mosaic (minutes arc)
% MStruct.Hmos = 60;            % Monochromator horizontal mosaic (minutes arc)
MStruct.mos = 45;               % Isotropic mosaic (minutes arc)
MStruct.NV = 10;                % Number of vertical monochromator segments
%MStruct.NV = 15;                % Number of vertical monochromator segments
MStruct.SegGap = 0.0005;        % Gap between monochromator segments (m)
MStruct.BladeN = 10;            % Number of monochromator blades
%MStruct.BladeN = 15;            % Number of monochromator blades
MStruct.BladeW = 0.02;          % Width of an individual blade (m)
MStruct.BladeH = 0.2;           % Height of an individual blade (m)
%MStruct.BladeH = 0.3;           % Height of an individual blade (m)
MStruct.BladeGap = 0.0005;      % Gap between blades (m)
MStruct.BeamDir = 'reflect';    % Beam direction flag ('reflect', 'transmit')
MStruct.HGeometry = 'rowland';  % Horizontal focusing geometry ('rowland', 'lens', 'flat')
MStruct.VGeometry = 'lens';     % Vertical focusing geometry ('lens', 'flat')
MStruct.Ebin = -2;            % Dynamic Monochromator Binning (neg. val. -> Binning with bin sizes scaled in units of energy resolution FWHM (meV);  0 -> No Binning;  pos. val. -> Sample single Ei with Ei = pos. val. (meV))
MStruct.MinimalistMono = -1;    % Minimalist Calculation for Monochromator (if Monochromaor is preceded by a gap). (1 -> apply Minimalist Principle;  -1 -> Do not apply)
MStruct.L2_override = -1;       % Overrides mono to sample calculated distance, use when BeamDir = 'transmit'. Only needed to get focusing correct. (set to -1 to disable) (m).

% This option makes sure that guide_bot only analyzes the performance for
% wavelengths less than 2.1*dM when a monochromator in reflection mode is
% present.
if strcmp(MStruct.BeamDir, 'reflect')
    options.max_wavelength_investigated_multiplier = 2.1*MStruct.dM/demands.WaveLmax;
end

% This option allows a user to investigate different monochromator energy binning widths before running a full optimization. 
% Note that it can be used with any instrument parameter if the name is known. For instance SEGM has monochromator
% with index 1 and the monochromator parameter Ebin. Thus, guidebot defines this as the instrument parameter Ebin1. 
% An associated array of values for this parameter must also be defined. All other parameter values are set to their initialization value. 
% To run, navigate to the instrument definition folder and run the [NAME]_show.inst file.
% General definition: options.show = {parameter_name, [parameter_values]}
options.show = {'Ebin1', [0, -0.5, -1, -2, -3]};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% End of Monochromator Parameters %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Straight Guide Parameters %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SStruct.start = 1.605;
SStruct.StartWidth = 0.05;
SStruct.StartHeight = 0.12;
SStruct.EndWidth = 0.05;
SStruct.EndHeight = 0.12;
SStruct.length = 0.8854;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% End of Straight Guide Parameters %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% ---- input strings ---
% each line here is a seperate guide to be optimized by guide_bot!


% input{1} = 'S(SStruct) E E M(MStruct)';
% input{1} = 'S(SStruct) E G M(MStruct)';
% input{end+1} = 'S(SStruct) E P G M(MStruct)';
% input{end+1} = 'S(SStruct) P E S G M(MStruct)';
% input{end+1} = 'S(SStruct) S G M(MStruct)';
%input{end+1} = 'S(SStruct) E E E G M(MStruct)';
%input{end+1} = 'S(SStruct) S S S G M(MStruct)';

%%%%%%%% CURRENT SPINS CONFIGURATION %%%%%%%
%load('SPINS_Params.mat')
%demands.Mod_sample=33.894; 
%demands.Dist=1.60;
%defaults.m = 1;
%defaults.alpha = 0;
%defaults.W = 0.000251;
%defaults.Qc = 0.0257;
%input{1} = 'S(SSPINS) G(GSPINS) M(MSPINS)';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Guide input section ---
% - Guides:
% S: Straight E: Elliptic P: Parabolic
% - Line of sight breakers:
% K: Kink C: Curved
% - Misc:
% G: Gap
% Needs to start and end with an element from the guide list.

% - General options (optional)
% length,StartWidth,StartHeight,start: will fix the value for the module
% start sets the absolute distance between the moderator and module start
% max or min in front of these will change the optimization interval
% e.g. maxStartWidth=0.04 will limit the width to maximum 4cm

% - Specific options (optional)
% Guide elements can have reflectivity options, R0,m,Qc,alpha.
% Line of sight breakers can be controlled diretly by:
% rot: the total angle the beam will be turned in degrees
% rotd: 'h' horizontal (standard), 'v' vertical turn
% rots: 1: rotates to the right, -1: rotates to the left


input{1} = 'S(SStruct) E G M(MStruct)';
 input{end+1} = 'S(SStruct) E G M(MStruct)';
 input{end+1} = 'S(SStruct) E G M(MStruct)';
 input{end+1} = 'S(SStruct) E G M(MStruct)';
 input{end+1} = 'S(SStruct) E G M(MStruct)';
%input{end+1} = 'S S S';
% input{1} = 'E(start=1.605,StartWidth=0.05,StartHeight=0.12) G M(MStruct)';
%input{end+1} = 'P K(maxlength=0.05) C E(los_end=2m_e)';
%input{end+1} = 'P K(maxlength=0.05) C E(Optimize_end=1,los_end=2m_e)';
%input{end+1} = 'P C K(maxlength=0.05) E(los_end=2m_e)';
%input{end+1} = 'P C K(maxlength=0.05) E(Optimize_end=1,los_end=2m_e)';
%input{end+1} = 'E K(maxlength=0.05) E(los_end=2m_e)';
%input{end+1} = 'E K(maxlength=0.05) E(Optimize_end=1,los_end=2m_e)';

%input{1}='P G(start=6.5,length=0.1) E(maxStartWidth=0.035) K(maxlength=2) E';
%input{1+end}='P G(start=6.5,length=0.1) E(maxStartWidth=0.035) K(maxlength=2) E P';
%input{1+end}='P G(start=6.5,length=0.1) E(maxStartWidth=0.035) C(minlength=1) E P';
%input{1+end}='P G(start=6.5,length=0.1) E(maxStartWidth=0.035) P C(minlength=1) E';
%input{1+end}='P G(start=6.5,length=0.1) E(maxStartWidth=0.035) S C(minlength=1) S E';
%input{1+end}='P G(start=6.5,length=0.1) E(maxStartWidth=0.035) E';

guide_bot(name,demands,requirements,defaults,input,options);

% Send a copy of this script to the run folder for record keeping
sourcefile = mfilename('fullpath');
[source_dir, source_name, source_ext] = fileparts(sourcefile);
copyfile([sourcefile, '.m'], ['./' name '/' source_name '_RECORDCOPY.m']);

