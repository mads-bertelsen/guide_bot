clear; clc;

% NAME OF PROJECT:
name='SPINS';

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
% demands.WaveLmin=3.96645;
% Highest wavelength needed [Å]
demands.WaveLmax=5.72;
% demands.WaveLmax=4.12841;

% Distance between guide end and sample [m]
% Distance between moderator and sample [m]
demands.Mod_sample=33.894; demands.Dist=1.6; % Standard
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

% --- computing options ---
options.cluster = 'NCNR';
%options.queue = 'long'; 
%options.enable_nested_compile = 1;
options.mpi=60;
options.machines='/home/lwh1/machines';

options.projectlist = 'add'; % ('new', 'single', or 'add')

% --- defaults ---
% mainly used for reflectivity options R0,Qc,m,alpha,W.
defaults.m = 1;
defaults.alpha = 0;
defaults.W = 0.000251;
defaults.Qc = 0.0257;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Monochromator Parameters %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mos = 30;
AveReflec = 0.75;
rcalc = AdjustReflectivityTable(mos,2,30,2,AveReflec);

MSPINS.r0 = 0.7;             % Reflectivity constanst [0-1]
MSPINS.Rfile = 'HOPG.rfl';% Reflectivity file
MSPINS.dM = 3.354;             % Monochromator d-spacing (Ang)
MSPINS.mos = mos;              % Isotropic mosaic (minutes arc)
MSPINS.NV = 5;                 % Number of vertical monochromator segments
MSPINS.SegGap = 0.0005;        % Gap between monochromator segments (m)
MSPINS.BladeN = 3;             % Number of monochromator blades
MSPINS.BladeW = 0.0492;        % Width of an individual blade (m)
MSPINS.BladeH = 0.1206;        % Height of an individual blade (m)
MSPINS.BladeGap = 0.0005;      % Gap between blades (m)
MSPINS.BeamDir = 'reflect';    % Beam direction flag ('reflect', 'transmit')
MSPINS.HGeometry = 'flat';     % Horizontal focusing geometry ('rowland', 'lens', 'flat')
MSPINS.VGeometry = 'lens';     % Vertical focusing geometry ('lens', 'flat')
MSPINS.Ebin = -2;              % Dynamic Monochromator Binning (neg. val. -> Binning with bin sizes scaled in units of energy resolution FWHM (meV);  0 -> No Binning;  pos. val. -> Sample single Ei with Ei = pos. val. (meV))
MSPINS.MinimalistMono = -1;    % Minimalist Calculation for Monochromator (if Monochromaor is preceded by a gap). (1 -> apply Minimalist Principle;  -1 -> Do not apply)
MSPINS.L2_override = -1;       % Overrides mono to sample calculated distance, use when BeamDir = 'transmit'. Only needed to get focusing correct. (set to -1 to disable) (m).

% This option makes sure that guide_bot only analyzes the performance for
% wavelengths less than 2.1*dM when a monochromator in reflection mode is
% present.
if strcmp(MSPINS.BeamDir, 'reflect')
    options.max_wavelength_investigated_multiplier = 2.1*MSPINS.dM/demands.WaveLmax;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% End of Monochromator Parameters %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Straight Guide Parameters (SSPINS) %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SSPINS.start = 1.605;
SSPINS.StartWidth = 0.05;
SSPINS.StartHeight = 0.12;
SSPINS.EndWidth = 0.05;
SSPINS.EndHeight = 0.12;
SSPINS.length = 30.318;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% End of Straight Guide Parameters %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Gap Parameters (GSPINS) %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

GSPINS.length = 0.0786;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% End of Gap Guide Parameters %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% ---- input strings ---
% each line here is a seperate guide to be optimized by guide_bot!



%%%%%%%% CURRENT SPINS CONFIGURATION %%%%%%%
defaults.m = 1;
defaults.alpha = 0;
defaults.W = 0.000251;
defaults.Qc = 0.0257;
input{1} = 'S(SSPINS) G(GSPINS) M(MSPINS)';
input{end+1} = 'S(SSPINS) G(GSPINS) M(MSPINS)';
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


%input{end+1} = 'S C E(Optimize_end=1,los_end=2m_e)';
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

