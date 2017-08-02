clear all; clc;

% NAME OF PROJECT:
name='First_guide_bot_run';

% Need to add the path of the folder in which guide_bot.m is located, but
% not necessary if the work directory is in the same folder.
%addpath('.') 

% --- Input / demands / figure of merit ---

% Horizontal divergence [deg (+/-)]
demands.Hdiv=1.5;
% Vertical divergence [deg (+/-)]
demands.Vdiv=1.0;
% Sample size horizontal [cm]
demands.Hsize=[2 2 2]; % repeats the optimization three times
% Sample size vertical [cm]
demands.Vsize=3;
% Lowest wavelength needed [Å]
demands.WaveLmin=3.0;
% Highest wavelength needed [Å]
demands.WaveLmax=4.0;
% Distance between guide end and sample [m]
demands.Dist=0.22;
% Distance between moderator and sample [m]
demands.Mod_sample=35;


% --- requirements made by facility / external factors ---

% The minimum distance between moderator and guide allowed
requirements.closest_element=2.0;
% The maximum distance between moderator and guide allowed
requirements.latest_start=4.0;
% Horizontal moderator size (rectangular)
requirements.moderator_size_x=0.10;
% Vertical moderator size
requirements.moderator_size_y=0.03;
% Source spectrum (ESS), cold, thermal or bispectral
% Notice this does not have an impact on the optimization
% The bispectral option is the maximum of the cold and thermal at all wavelengths
requirements.source='butterfly';
options.beamline_number = 5;
options.butterfly_face = 'E';

% uses ideal source for optimization (flat wavelength distribution)
options.optimizer_mode='ideal_source';
% uses ESS source for optimization (realistic wavelength distribtution)
%options.optimizer_mode='realistic_source'; 
% optimization on ideal source follow by optimization of moderator view on
% realistic source
%options.optimizer_mode='combined';

% --- options which affects the optimization ---
% Controll how much should be calculated by the minimalist concept
%  1: Extraction, gaps and exit
%  0: Gaps and exit
% -1: Exit
% The lower the number, the more degrees of freedom for the optimizer
options.minimalist=0;

% Option to include figures on absolute flux performance on ESS
% 1: Enable. 0:Disable. % Currently bugged, always keep enabled.
options.absolute_intensity_run = 1;

% --- computing options ---
options.cluster = 'ESSS';
options.queue = 'long'; 
options.mpi=2; % number of cores used when running locally
%options.enable_nested_compile = 1;

% --- defaults ---
% mainly used for reflectivity options R0,Qc,m,alpha,W.
defaults.m=3.0;

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

% input strings, each string describes a guide that is optimized individually 
input{1}='S C S';
input{end+1}='S(maxlength=10) C S';
input{end+1}='S C(start=12) S';


guide_bot(name,demands,requirements,defaults,input,options);
