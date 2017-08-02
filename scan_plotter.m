% script for reading and plotting a general guide_bot scan.
clear all;close all;clc;
% EXPLANATION
clear all

% This script can be used to plot large amounts of data made with guide_bot
% in order to get an overview. Because it is made for large data sets, it
% is "uncrashable"^TM. Every time something is done, it is saved to disk,
% and if the process crash or is killed it will start from where it was
% killed the next time it is started. I usually run this script over night,
% as it can take 5 hours or more for reasonable data sets. 
% If the process somehow fails and you want to start over, delete the
% progress.mat file in this folder.

% HOW TO USE IT:
% Simply navigate to the output/analysis folder and run this script. It
% will automaticly detect any files made by guide_bot and run the analyzse
% script while recording the overall performance from each script. Remember
% to set meaningful wavelength ranges, if Wmax is above what is recoreded
% it will fail.

% START INPUT SECTION

% Option to plot the scans at the end. 1 = Yes, 0 = No.
plot_options.plot_results = 1;
% Using this script with plot_results = 0 is like analyze_all with a system to save progress.

%plot_options.nickname = 'complex';

% Which of the above wavelength ranges to plot performance for?
plot_options.wavelength_range = 7;

% Calculate performance in these intervals:
%calculate.Wmax = [8 0.5 1.0 2.0 3.0 4];
%calculate.Wmin = [4 0.1 0.5 1.0 2.0 1];


% Pascale
%calculate.Wmax = [9.0 5.0 5.0 2.0 3.00 4.00 5.00 6.00];% 7.00 8.00 9.00];
%calculate.Wmin = [1.0 1.0 1.5 1.0 2.00 3.00 4.00 5.00];% 6.00 7.00 8.00];
% CAMEA
%calculate.Wmax = [6.40 3.00 2.0 3.00 4.00 5.00 6.00 7.00];
%calculate.Wmin = [1.65 1.65 1.0 2.00 3.00 4.00 5.00 6.00];
% ODIN
%calculate.Wmax = [7 4 7 2.0 3.00 4.00 5.00 6.00 7.00 7.99];
%calculate.Wmin = [1 1 4 1.0 2.00 3.00 4.00 5.00 6.00 6.99];
% Werner Thermal
%calculate.Wmax = [2.4 1.5 1.0 1.5 2.0 2.5 3.0 3.5 4.0 4.5]
%calculate.Wmin = [0.8 0.8 0.5 1.0 1.5 2.0 2.5 3.0 3.5 4.0]
% Werner Cold
%calculate.Wmax = [10.0 4.0 2.0 3.00 4.00 5.00 6.00 7.00 8.00 9.00 10.00 11.00 12.00 13.00 14.00 15.00 16.00 17.00];
%calculate.Wmin = [2.40 2.4 1.0 2.00 3.00 4.00 5.00 6.00 7.00 8.00 9.000 10.00 11.00 12.00 13.00 14.00 15.00 16.00];
% Heimdal
%calculate.Wmax = [2.27 2.40 1.00 2.0 2.4 3.00 4.00 5.00 6.00];
%calculate.Wmin = [0.60 0.60 0.60 1.0 2.0 2.00 3.00 4.00 5.00];
% Heimdal cold
%calculate.Wmax = [10.0 10.0 5.00 6.00 7.00 8.00 9.00 10.00];
%calculate.Wmin = [4.00 6.00 4.00 5.00 6.00 7.00 8.00 9.000];
% Fundamental Old
%calculate.Wmax = [16.0 4.00 2.0 3.00 4.00 5.00 6.00 7.00 8.00 9.00 10.00 11.00 12.00 13.00 14.00 15.00 16.00];
%calculate.Wmin = [2.00 2.00 1.0 2.00 3.00 4.00 5.00 6.00 7.00 8.00 9.000 10.00 11.00 12.00 13.00 14.00 15.00];
% Fundamental
%calculate.Wmax = [8.00 4.00 2.0 3.00 4.00 5.00 6.00 7.00 8.00 9.00 10.00 11.00 12.00 13.00 14.00 15.00 16.00];
%calculate.Wmin = [3.00 2.00 1.0 2.00 3.00 4.00 5.00 6.00 7.00 8.00 9.000 10.00 11.00 12.00 13.00 14.00 15.00];
% Minimalist
%calculate.Wmax = [4 8 8 8 8 2.0 3.00 4.00 5.00 6.00 7.00 8];
%calculate.Wmin = [1 1 3 4 5 1.0 2.00 3.00 4.00 5.00 6.00 7];
% Thesis complicated 1
%calculate.Wmax = [4.00 8.00 2.0 3.00 4.00 5.00 6.00 7.00 8.00];
%calculate.Wmin = [1.50 4.00 1.0 2.00 3.00 4.00 5.00 6.00 7.00];
% Thesis complicated feeder
%calculate.Wmax = [10.00 10.00 2.0 3.00 4.00 5.00 6.00 7.00 8.00];
%calculate.Wmin = [2.00  4.00  1.0 2.00 3.00 4.00 5.00 6.00 7.00];
% Thesis complicated 1
%calculate.Wmax = [4.00 8.00 2.0 3.00 4.00 5.00 6.00 7.00 8.00];
%calculate.Wmin = [1.25 4.00 1.0 2.00 3.00 4.00 5.00 6.00 7.00];
% Wiebke CSPEC
%calculate.Wmax = [10.0 10.0 5.0 5.0 2.0 3.00 4.00 5.00 6.00 7.00 8.00 9.00 10.00];
%calculate.Wmin = [1.00 2.00 1.0 2.0 1.0 2.00 3.00 4.00 5.00 6.00 7.00 8.00 9.000];
% ESSENSE
%calculate.Wmax = [25.0 25.0 6 8 10.0 12 16];
%calculate.Wmin = [4.00 6.00 4 6 8.00 10 12];
% Werner singleXmagnetism
%calculate.Wmax = [8.00 3.00 1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0]
%calculate.Wmin = [0.70 0.70 0.7 1.0 2.0 3.0 4.0 5.0 6.0 7.0]
% WANSE
%calculate.Wmax = [10.0 10.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0]
%calculate.Wmin = [2.00 4.00 1.0 2.0 3.0 4.0 5.0 6.0 7.0]
% MIRACLES
%calculate.Wmax = [8.00 16.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0 11.0]
%calculate.Wmin = [2.00 8.00 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0]
% LOKI
%calculate.Wmax = [13.00 13.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0 11.0]
%calculate.Wmin = [2.000 4.00 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0]
% BiSpecChopper TREX
%calculate.Wmax = [7.20 2.00 7.2 1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0];
%calculate.Wmin = [0.80 0.80 2.0 0.8 1.0 2.0 3.0 4.0 5.0 6.0 7.0];
% NMX
calculate.Wmax = [3.3 2.0 3.3 2.00 3.00 4.00 5.00 6.00 7.00];
calculate.Wmin = [1.5 1.5 2.0 1.00 2.00 3.00 4.00 5.00 6.00];
% DREAM
%calculate.Wmax = [4.6 2.0 4.6 2.00 3.00 4.00 5.00 6.00 7.00];
%calculate.Wmin = [0.8 0.8 2.0 1.00 2.00 3.00 4.00 5.00 6.00];
% HOD
%calculate.Wmax = [2.0 3.0 4.0 5.0 6.0 7.0 8.0]
%calculate.Wmin = [1.2 2.0 3.0 4.0 5.0 6.0 7.0]
% VERITAS
%calculate.Wmax = [10.0 10.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0]
%calculate.Wmin = [2.00 4.00 1.0 2.0 3.0 4.0 5.0 6.0 7.0]
% Thermal Chopper Spectrometer
%calculate.Wmax = [3.0 1.5 1.0 1.5 2.0 2.5 3.0]
%calculate.Wmin = [0.6 0.6 0.5 1.0 1.5 2.0 2.5]
% SLEIPNIR
%calculate.Wmax = [19.00 16.5 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0 12.0 14 16 18]
%calculate.Wmin = [3.000 3.00 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0  10.0 12 14 16]
% SKADI
%calculate.Wmax = [10.00 9.5 9.5 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0]
%calculate.Wmin = [2.000 2.0 4.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.00]
% Vibrational
%calculate.Wmax = [6.0 3.0 6.0 1.0 1.5 2.0 2.5 3.0 4 5 6]
%calculate.Wmin = [0.4 0.4 3.0 0.5 1.0 1.5 2.0 2.5 3 4 5]
% VESPA
%calculate.Wmax = [4.7 2.0 4.7 1.0 1.5 2.0 2.5 3.0 4 5]
%calculate.Wmin = [0.6 0.6 2.0 0.5 1.0 1.5 2.0 2.5 3 4]
% Micron
%calculate.Wmax = [8.5 8.5 1.0 1.5 2.0 2.5 3.0 4 5 6 7 8 9];
%calculate.Wmin = [0.5 1.0 0.5 1.0 1.5 2.0 2.5 3 4 5 6 7 8];
% FREIA
%calculate.Wmax = [10.0 5.0 3.5 4.5 5.5 6.5 7.5 8.5 10 ];
%calculate.Wmin = [2.5  2.5 2.5 3.5 4.5 5.5 6.5 7.5 8.5];
% RESPECT
%calculate.Wmax = [6.0 10.0 3 4 5 6 7 8 9 10];
%calculate.Wmin = [2.0 2.00 2 3 4 5 6 7 8 9];
% D7
%calculate.Wmax = [6.0 3.0 3 4 5 6]
%calculate.Wmin = [2.0 2.0 2 3 4 5]
% Mark
%calculate.Wmax = [3.5 6.8 4 5 6 7 8];
%calculate.Wmin = [3.1 6.2 3 4 5 6 7];
% Undervisning
%calculate.Wmax = [4.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0];
%calculate.Wmin = [2.0 1.0 2.0 3.0 4.0 5.0 6.0 7.0];
% Katsuaki
%calculate.Wmax = [8.0 2.0 8.0 1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0]
%calculate.Wmin = [0.5 0.5 2.0 0.5 1.0 2.0 3.0 4.0 5.0 6.0 7.0]


%calculate.Wmax = [2.0 1.0 1.5 2.0 3.0 4.0 5.0 6.0];
%calculate.Wmin = [1.0 0.5 1.0 1.5 2.0 3.0 4.0 5.0];

% What figures do you want? [3 5] recomended.
% 1: Performance from simple source
% 2: 1 with degraded guide
% 3: ESS 
% 4: 3 with degraded guide
% 5: Brilliance transfer from uniform source (equal I on all wavelengths)  
% 6: 5 with degraded guide
plot_options.fig_choice = [5 3 7];

% Colormap to distinguish different guides, hot and jet are ok.
plot_options.colormap_choice = jet;    
    
% If nxm 2D plot, plot n 1D graphs with m points instead if force_1d = 1.
plot_options.force_1D = 1;
% If reverse_order = 1 and force_1D = 1, it will plot m 1D graphs with n points. 
plot_options.reverse_order = 1;

% Text size, FSL (large) and FS (small)
plot_options.FSL = 14;
plot_options.FS = 12;

% Algorithm for cleaning up in data from optimizer which can fail.
% cleanup_low_points should be zero when looking at the data for the first time
% when one, it will remove points which is below tollerance*data(i-1) and
% tollerance*data(i+1), meaning failled optimizations.here 0
plot_options.cleanup_low_points = 1;
plot_options.tollerance = 0.70;

% There is hardcoded data for the unperterped moderator gains for pancake
% moderators in this script. It will only be applied if moderator_size_y is
% scanned. Enable this by setting the following command to one.
plot_options.apply_flux = 0;


% the keyword field can be used to add comments to some of the tested
% guides. These will apear in the legend.
%plot_options.keyword.EGSlitGE = '50mm';
%plot_options.keyword.PGSlitGE = '15mm';
%plot_options.keyword.PGSlitGE_alt1 = '20mm';
%plot_options.keyword.PGSlitGE_alt2 = '25mm';
%plot_options.keyword.PGSlitGE_alt3 = '50mm';

% Werner cold2
%plot_options.keyword.EGSelene = ' no m limit';
%plot_options.keyword.EGSelene_alt1 = 'm 2.5 limit';
%plot_options.keyword.EGSelene_alt2 = 'm 2.0 limit';
%plot_options.keyword.SCGSelene = ' no m limit';
%plot_options.keyword.SCGSelene_alt1 = 'm 2.5 limit';
%plot_options.keyword.SCGSelene_alt2 = 'm 2.0 limit';
%plot_options.keyword.Selene = ' no m limit';
%plot_options.keyword.Selene = '2 m slit m=4.3';
%plot_options.keyword.Selene_alt1 = '2 m slit m=3.0 limit';
%plot_options.keyword.Selene_alt1 = 'moderator focus m=4.3';
%plot_options.keyword.Selene_alt3 = 'moderator focus m=3.0 limit';


%plot_options.keyword.PGSCSCGP = '2 m start focusing end';
%plot_options.keyword.PGSCSCGP_alt1 = '3.5 m start focusing end';
%plot_options.keyword.PGSCSCGS = '2 m start straight end';
%plot_options.keyword.PGSCSCGS_alt1 = '3.5 m start straight end';


plot_options.keyword.PGESKSE = 'Proposal guide';
%plot_options.keyword.PGPSP = '(STAP) 0 los';
%plot_options.keyword.PGPCP = '(STAP) 1 los';
%plot_options.keyword.PGPCP_alt1 = '(STAP) 2 los';

%plot_options.keyword.ECE = '6.1 m start - 1 los';
%plot_options.keyword.ECP = '6.1 m start - 1 los';
%plot_options.keyword.PCP = '6.1 m start - 1 los';
%plot_options.keyword.PCE = '6.1 m start - 1 los';

%plot_options.keyword.PGECE = 'free start - 1 los';
%plot_options.keyword.PGECP = 'free start - 1 los';
%plot_options.keyword.PGPCP = 'free start - 1 los';
%plot_options.keyword.PGPCE = 'free start - 1 los';

%plot_options.keyword.PGPCP = 'no pinhole - 1 los';
%plot_options.keyword.PGPCP_alt1 = 'no pinhole - 2 los';
%plot_options.keyword.PGSlitGPCP = '15mm pinhole - 1 los';
%plot_options.keyword.PGSlitGPCP_alt1 = '20mm pinhole - 1 los';
%plot_options.keyword.PGSlitGPCP_alt2 = '25mm pinhole - 1 los';
%plot_options.keyword.PGSlitGPSP = '15mm pinhole - 0 los';
%plot_options.keyword.PGSlitGPSP_alt1 = '20mm pinhole - 0 los';
%plot_options.keyword.PGSlitGPSP_alt2 = '25mm pinhole - 0 los';


%plot_options.keyword.PGECE = '3.0 cm pinhole - 1 los';
%plot_options.keyword.PGECE_alt1 = '3.0 cm pinhole - 2 los';
%plot_options.keyword.PGPCP = '3.0 cm pinhole - 1 los';
%plot_options.keyword.PGECE_alt2 = '4.5 cm pinhole - 1 los';
%plot_options.keyword.PGECE_alt3 = '4.5 cm pinhole - 2 los';
%plot_options.keyword.PGPCP_alt1 = '4.5 cm pinhole - 1 los';


% ODIN6 input names
% plot_options.keyword.PGSlitGECE = ' 1.5 cm slit 1 channel';
% plot_options.keyword.PGSlitGECE_alt1 = ' 1.5 cm slit 3 channels';
% plot_options.keyword.PGSlitGECE_alt2 = ' 1.5 cm slit 5 channels';
% plot_options.keyword.PGSlitGECE_alt3 = ' 2.0 cm slit 1 channel';
% plot_options.keyword.PGSlitGECE_alt4 = ' 2.0 cm slit 3 channels';
% plot_options.keyword.PGSlitGECE_alt5 = ' 2.0 cm slit 5 channels';
% plot_options.keyword.PGSlitGECE_alt6 = ' 2.5 cm slit 1 channel';
% plot_options.keyword.PGSlitGECE_alt7 = ' 2.5 cm slit 3 channels';
% plot_options.keyword.PGSlitGECE_alt8 = ' 2.5 cm slit 5 channels';
% plot_options.keyword.PGSlitGPCP = ' 1.5 cm slit 1 channel';
% plot_options.keyword.PGSlitGPCP_alt1 = ' 1.5 cm slit 3 channels';
% plot_options.keyword.PGSlitGPCP_alt2 = ' 1.5 cm slit 5 channels';
% plot_options.keyword.PGSlitGPCP_alt3 = ' 2.0 cm slit 1 channel';
% plot_options.keyword.PGSlitGPCP_alt4 = ' 2.0 cm slit 3 channels';
% plot_options.keyword.PGSlitGPCP_alt5 = ' 2.0 cm slit 5 channels';
% plot_options.keyword.PGSlitGPCP_alt6 = ' 2.5 cm slit 1 channel';
% plot_options.keyword.PGSlitGPCP_alt7 = ' 2.5 cm slit 3 channels';
% plot_options.keyword.PGSlitGPCP_alt8 = ' 2.5 cm slit 5 channels';


%plot_options.keyword.EGECE = '        2.0 m from moderator';
%plot_options.keyword.EGECE_alt1 = '1.0 m from moderator';
%plot_options.keyword.EGE_alt2 = '1.0 m from moderator';

plot_options.keyword.EGESE = '        2.0 m from moderator';
plot_options.keyword.EGESE_alt1 = '1.5 m from moderator';
plot_options.keyword.EGESE_alt2 = '1.0 m from moderator';

plot_options.keyword.EGPCEGP='        2.0 m from moderator';
plot_options.keyword.EGPCEGP_alt1= '1.0 m from moderator';

%plot_options.keyword.SCSCS = '        2.0 m from moderator';
%plot_options.keyword.SCSCS_alt1 = '1.5 m from moderator';
%plot_options.keyword.SCSCS_alt2 = '1.0 m from moderator';

%plot_options.keyword.EGSlitGESE = '        2.0 m from moderator';
%plot_options.keyword.EGSlitGESE_alt1 = '1.0 m from moderator';
%plot_options.keyword.EGSlitGE = '        2.0 m from moderator';
%plot_options.keyword.EGSlitGE_alt1 = '1.0 m from moderator';

%plot_options.plot_logical.ESE = 0;
%plot_options.plot_logical.ESE_alt1 = 0;
%plot_options.plot_logical.ES_alt3 = 0;


%plot_options.plot_logical.PSP = 0;
%plot_options.plot_logical.PSP_alt1 = 0;
%plot_options.plot_logical.ESE = 0;
%plot_options.plot_logical.ESE_alt1 = 0;

plot_options.keyword.SCCS = '        1.0 m from moderator vertical bend';
plot_options.keyword.SCCS_alt1 = '2.0 m from moderator vertical bend';
plot_options.keyword.SCCS_alt2 = '1.0 m from moderator horizontal bend';
plot_options.keyword.SCCS_alt3 = '2.0 m from moderator horizontal bend';

% END INPUT SECTION


plot_options.keyword.ES = '        2.0 m from moderator';
plot_options.keyword.ES_alt1 = '1.0 m from moderator';
plot_options.keyword.S = '        2.0 m from moderator';
plot_options.keyword.S_alt1 = '1.0 m from moderator';
plot_options.keyword.SSS = '        2.0 m from moderator';
plot_options.keyword.SSS_alt1 = '1.0 m from moderator';


plot_options.keyword.ESESelene = '4 m Selene';
plot_options.keyword.ESKSESelene = '4 m Selene';
plot_options.keyword.PGECESelene = '4 m Selene';

plot_options.keyword.ESESelene_alt1 = '8 m Selene';
plot_options.keyword.ESKSESelene_alt1 = '8 m Selene';
plot_options.keyword.PGECESelene_alt1 = '8 m Selene';


plot_options.keyword.SCSCSP = '       4 channels - 2 los at 11.8 m';
plot_options.keyword.SCSCSP_alt1 = '4 channels - 2 los at 15.0 m';


plot_options.keyword.SCSCSP_alt2 = '- 2 channels - 2 los at 11.8';
plot_options.keyword.SCSCSP_alt3 = '- 1 channels - 2 los at 11.8';
plot_options.keyword.SCSCSP_alt4 = '- 4 channels - 2 los at 15.0';
plot_options.keyword.SCSCSP_alt5 = '- 3 channels - 2 los at 15.0';
plot_options.keyword.SCSCSP_alt6 = '- 2 channels - 2 los at 15.0';
plot_options.keyword.SCSCSP_alt7 = '- 1 channels - 2 los at 15.0';
plot_options.keyword.SCCSP = '- No middle S - 3 channels - 2 los at 11.8';
plot_options.keyword.SCCSP_alt1 = '- No middle S - 3 channels - 2 los at 15.0';

plot_options.keyword.S = 'Straight guide';
plot_options.keyword.E = 'Elliptic guide';
plot_options.keyword.SCS = 'Curved guide';
plot_options.keyword.PCP = 'Balistic guide: parabolic';
plot_options.keyword.ECE = 'Balistic guide: elliptic';

% Code for identifying the different instruments
    % find all .mat files in the directory
    % look for unique names before the first underscore

% Code for identifying the different scanned variables
    % count the number of underscores in the names
    % check against a list of possible scan names (Hdiv, moderator_size_y)
    % conclude which values were scanned


    % could add wavemax / wavemin to the list.
    possible_scan_names={'Hdiv' 'Vdiv' 'Hsize' 'Vsize' 'WaveLmin' 'WaveLmax' 'moderator_size_x' 'moderator_size_y' 'minimalist_factor'};
    possible_scan_names_mcstas={'divreq_x' 'divreq_y' 'sizeX' 'sizeY' 'WaveMin' 'WaveMax' 'mod_x' 'mod_y' 'minimalist_factor'};
    possible_scan_units={'[Deg]' '[Deg]' '[m]' '[m]' '[AA]' '[AA]' '[m]' '[m]' ''};
    
    dir_output = dir;
    verbose = 0;
    
    no_scan = 0;
    % This code does not work for _alt versions!
    for ii = 3:length(dir_output) % . and .. are ellements.
        if length(dir_output(ii).name) > 4
            if strcmp(dir_output(ii).name(end-3:end),'.mat')
                disp(ii)
                disp(dir_output(ii).name)
                %if strcmp(dir_output(ii).name,'progress.mat'); verbose = 1; else; verbose = 0; end;
                
                for j = 1:length(dir_output(ii).name)
                    if strcmp(dir_output(ii).name(j),'_')
                        underscore{ii}(j) = 1;
                    else
                        underscore{ii}(j) = 0;
                    end
                end
                
                if verbose; underscore{ii}
                end;

                if sum(underscore{ii}) ~= 0
                    if verbose; disp('in underscore search'); end;

                    underscore_index{ii} = find(underscore{ii});
                    
                    if strcmp(dir_output(ii).name(underscore_index{ii}(1)+1:underscore_index{ii}(1)+3),'alt')
                        pre{ii} = dir_output(ii).name(1:underscore_index{ii}(2)-1);
                    else
                        pre{ii} = dir_output(ii).name(1:underscore_index{ii}(1)-1);
                    end

                    not_found = 1;
                    interesting = 1;
                    limit = 0;
                    while not_found
                        limit = limit + 1; 
                        for k = 1:length(possible_scan_names)
                            lower = underscore_index{ii}(end)-limit-1;
                            if lower > 2
                                if strcmp(dir_output(ii).name(lower:underscore_index{ii}(end)-1),possible_scan_names{k})
                                    not_found = 0;
                                    scan_dim(1).name = possible_scan_names{k};
                                    scan_dim(1).name_mcstas = possible_scan_names_mcstas{k};
                                    scan_dim(1).index = k;
                                end
                            else
                                not_found = 0;
                                interesting = 0;
                                % not interesting
                            end
                        end
                    end
                    
                    if strcmp(dir_output(ii).name(length(pre{ii})+1:end),'_all.mat')
                       % In here if the file is made by guide_bot and not
                       % in a scan.
                       
                       % Safe to assume there is no scan in this folder
                       no_scan = 1;
                       relevant(ii) = true;
                       
                       if exist('scan_pre')
                           scan_pre{end+1} = pre{ii};
                       else
                           scan_pre{1} = pre{ii};
                       end
                       
                    end


                    if interesting

                        relevant(ii) = true;
                        %disp(ii)
                        % need to extract the scan number, what number is t
                        lower_underscore = underscore_index{ii} < lower;
                        previous_underscore = max(underscore_index{ii}(lower_underscore));

                        value1 = str2num(dir_output(ii).name(previous_underscore+1:lower-1));
                        if sum(ismember(fieldnames(scan_dim(1)),'values')) > 0.5
                            scan_dim(1).values(end+1) = value1;
                            scan_pre{end+1}=pre{ii};
                        else
                            scan_dim(1).values(1) = value1;
                            scan_pre{1}=pre{ii};
                        end


                        % Need to check if there is another dimension or not.
                        % is there underscore between length(pre) and first letter in scan name? 
                        logic1 = underscore_index{ii} > length(pre{ii}) + 1; % check
                        logic2 = underscore_index{ii} < underscore_index{ii}(end) - length(scan_dim(1).name); % check

                        logic_combine = logic1 & logic2;
                        sum(logic_combine);

                        if sum(logic_combine) > 0.5 
                            % this is a two dimensional scan!
                            % store name in scan_dim(2).name and index in scan_dim(2).index
                            % store somewhere that this is a two dimensional scan
                            underscore_after = max(underscore_index{ii}(logic2));

                            not_found = 1;
                            limit = 0;
                            interesting = 0;
                            while not_found
                                limit = limit + 1; 
                                for k = 1:length(possible_scan_names)
                                    lower = underscore_after-limit-1;
                                    if lower > 2
                                        if strcmp(dir_output(ii).name(lower:underscore_after-1),possible_scan_names{k})
                                            not_found = 0;interesting = 1;
                                            scan_dim(2).name = possible_scan_names{k};
                                            scan_dim(2).name_mcstas = possible_scan_names_mcstas{k};
                                            scan_dim(2).index = k;
                                        end
                                    else
                                       interesting = 0;
                                       not_found = 0;
                                    end 
                                end
                            end

                            if interesting

                                lower_underscore = underscore_index{ii} < lower;
                                previous_underscore = max(underscore_index{ii}(lower_underscore));

                                value2 = str2num(dir_output(ii).name(previous_underscore+1:lower-1));
                                if sum(ismember(fieldnames(scan_dim(2)),'values')) > 0.5
                                    scan_dim(2).values(end+1) = value2;
                                else
                                    scan_dim(2).values(1) = value2;
                                end


                            end
                        end
                    end
                end
            end
        end
    end
    
% As the names are taken from the back of the string, they are swithced if
% there were two.



if no_scan == 1
    size_select = 0;
    % need to build up nessecary variables
    list{1}=1;
    list{2}=1;
    name{1} = scan_pre{1};
    for ii = 1:length(scan_pre)
        new = 1;
        for jj = 1:length(name)
            if strcmp(name{jj},scan_pre{ii})
               new = 0; 
            end
        end
        if new
        name{end+1} = scan_pre{ii};
        end
    end
else
    if length(scan_dim) == 2
        temp = scan_dim(1);
        scan_dim(1) = scan_dim(2);
        scan_dim(2) = temp;
        size_select = 2;
    else
        size_select = 1;
    end


    % Need to tie pre names together with values
    name{1} = scan_pre{1};
    for ii = 1:length(scan_pre)
        new = 1;
        for jj = 1:length(name)
            if strcmp(name{jj},scan_pre{ii})
               new = 0; 
            end
        end
        if new
        name{end+1} = scan_pre{ii};
        end
    end
    
    list{1}(1) = scan_dim(1).values(1);
    for ii = 1:length(scan_dim(1).values)
        new = 1;
        for jj = 1:length(list{1})
            if list{1}(jj) == scan_dim(1).values(ii)
                new = 0;
            end
        end
        if new
        list{1}(end+1) = scan_dim(1).values(ii);
        end
    end
    
    % wrong
    if length(scan_dim) == 2
        list{2}(1) = scan_dim(2).values(1);
        for ii = 1:length(scan_dim(2).values)
            new = 1;
            for jj = 1:length(list{2})
                if list{2}(jj) == scan_dim(2).values(ii)
                    new = 0;
                end
            end
            if new
            list{2}(end+1) = scan_dim(2).values(ii);
            end
        end
    else
        list{2} = 1;
    end
end



% Code for looping over the analyse scripts (can read from analyse_all)
    % check for program state file.
    % for i = all instruments
    %  for j = all 1st_scan_variable
    %   for k = all 2nd_scan_variable
    %    check if the data, script and brill_ref exists
    %    run the analyse script
    %    save relevant data (B,B_deg,ESS,ESS_deg)
    %    clear irrelevant data
    %    save new program state
    %   end
    %  end
    % end
    
    % Update the following code to use brill transfer and user defined
    % wavelength interval. 
    
    % Length of Wmax should equal Wmin.
    %Wmax = [6.4 10 2.27 2.27 8 3];
    %Wmin = [1.65 1.0 0.6 0.4 2 2];
    
    % CAMEA
    %Wmax = [6.4 6.40 2.00 1.65];
    %Wmin = [1.65 1.00 1.65 1.00];
    
    % ODIN
    %Wmax = [8  2 3 4 10];
    %Wmin = [1.0 1 2 3 2];
    
    % Selene
    %Wmax = [6.00 2.00 3.00 4.00 5.00 6.00 7.00 8.00 8.00];
    %Wmin = [1.00 1.00 2.00 3.00 4.00 5.00 6.00 7.00 4.00];
    
    % Applying the inputs
    Wmax = calculate.Wmax;
    Wmin = calculate.Wmin;
    
% Before we start, check if any progress is saved
DONE = 0;
if exist('progress.mat') > 0.5;
    load('progress.mat')
    
    instr = 1;
    while DONE == 0 && sum(sum(progress_logic(instr,:,:))) == length(progress_logic(1,:,1))*length(progress_logic(1,1,:));
        instr = instr + 1;
        if instr > length(name); disp('Allready done!'); DONE = 1; end;
    end
    
    
    
    if DONE == 0;
    current_instr = instr;
    instr_range = current_instr:length(name);
    
    % make standard to do list first, then remove what have been done.
    for i = 1:length(name); iindex_range{i} = 1:length(list{1}); end;
    for i = 1:length(name); for j = 1:length(list{1}); jindex_range{i,j} = 1:length(list{2}); end; end;
    
    % current_instr needs to be updated
    
    current_iindex = 1;
    while sum(progress_logic(current_instr,current_iindex,:)) == length(list{2})
        current_iindex = current_iindex + 1;
    end
    iindex_range{current_instr} = current_iindex:length(list{1});
    
    current_jindex = sum(progress_logic(current_instr,current_iindex,:))+1;
    jindex_range{current_instr,current_iindex} = current_jindex:length(list{2});
    end
    
    disp('Using saved progress, to disable, delete progress.mat')
    disp(['Starting from instrument number: ' num2str(instr_range(1))]);
    disp(['Starting from iindex     number: ' num2str(iindex_range{instr_range(1)}(1))]);
    disp(['Starting from jindex     number: ' num2str( jindex_range{instr_range(1),iindex_range{instr_range(1)}(1)})]);
    disp(['Using Wmax/Wmin'])
    disp(Wmax)
    disp(Wmin)
    
else
    instr_range = 1:length(name);
    for i = 1:length(name); iindex_range{i} = 1:length(list{1}); end;
    for i = 1:length(name); for j = 1:length(list{1}); jindex_range{i,j} = 1:length(list{2}); end; end;
    progress_logic = zeros(length(name),length(list{1}),length(list{2}));
end

% Init variables
ALLW_B_wave(1,1).startup=0;
ALLW_d_B_wave(1,1).startup = 0;

save('before_main.mat')

% REMOVE WHEN DONE
% TEMP CODE
%iindex_range{1} = 1:5;
%jindex_range{1} = 1:10;

if DONE == 0;
for instrument=instr_range
    for iindex=iindex_range{instrument}
        for jindex=jindex_range{instrument,iindex}
            % size_select = 0 => no scan. = 1 => one dimensional scan. = 2 => two dimensional scan
            if size_select == 2
              logic_req1 = exist([name{instrument} '_' num2str(iindex) scan_dim(1).name '_' num2str(jindex) scan_dim(2).name '_ifit_analyse']) > 0.5;
              logic_req2 = exist([name{instrument} '_' num2str(iindex) scan_dim(1).name '_' num2str(jindex) scan_dim(2).name '_all.mat']) > 0.5;
            elseif size_select == 1  
              logic_req1 = exist([name{instrument} '_' num2str(iindex) scan_dim(1).name '_ifit_analyse']) > 0.5;
              logic_req2 = exist([name{instrument} '_' num2str(iindex) scan_dim(1).name '_all.mat']) > 0.5;
            elseif size_select == 0
              logic_req1 = exist([name{instrument} '_ifit_analyse']) > 0.5;
              logic_req2 = exist([name{instrument} '_all.mat']) > 0.5;
            end
            if (logic_req1 + logic_req2 == 2)
            close all;
            clearvars -except ALLW_wave ALLW_d_wave ESS_wave ESS_d_wave ESS_lim_wave ESS_d_lim_wave ALLW_B_wave ALLW_d_B_wave scan_dim name list iindex jindex instrument Wmax Wmin wavband progress_logic jindex_range iindex_range instr_range size_select calculate plot_options
            if size_select == 2
              eval([name{instrument} '_' num2str(iindex) scan_dim(1).name '_' num2str(jindex) scan_dim(2).name '_ifit_analyse']);
              temp_name = [name{instrument} '_' num2str(iindex) scan_dim(1).name '_' num2str(jindex) scan_dim(2).name];
            elseif size_select == 1
              eval([name{instrument} '_' num2str(iindex) scan_dim(1).name '_ifit_analyse']);
              temp_name = [name{instrument} '_' num2str(iindex) scan_dim(1).name ];
            elseif size_select == 0
              eval([name{instrument} '_ifit_analyse']);
              temp_name = [name{instrument}];
            end
            % This will be done when all local variables of the script are
            % accesable. Extract usefull info and move on to the next.
                %I_ALLW(iindex,jindex,instrument) = monitor_ALLW(10).Data.values(1);
                %I_ALLW_d(iindex,jindex,instrument) = monitor_ALLW_degraded(10).Data.values(1);
                %I_ESS(iindex,jindex,instrument) = monitor_ESSW(10).Data.values(1);
                %I_ESS_d(iindex,jindex,instrument) = monitor_ESSW_degraded(10).Data.values(1);
                
                
                if 1==2 
                    % hardcoded plotting extend for ODIN

                    monitor_num_image = 20;
                    close all;
                    figure_handle = figure(1);
                    set(figure_handle, 'Position', [0 0 1800 1800])
                    axes1 = axes('Parent',figure_handle,'YDir','reverse','Position',[0 0 1 1],'Layer','top');
                    object = monitor_ESSW(monitor_num_image);
                    signal_data = object.signal;
                    %x_data = object.x;
                    %y_data = object.y;
                    x_data = object{2};
                    y_data = object{1};

                    imagesc(x_data,y_data,signal_data)
                    colormap('gray')

                    axis off

                    %surf(peaks)
                    % control the image pixel size by manipulating the paper size and number of dots per inch
                    output_size = [1250 1250];%Size in pixels
                    resolution = 300;%Resolution in DPI
                    set(gcf,'paperunits','inches','paperposition',[0 0 output_size/resolution]);
                    % use 300 DPI
                    print([filename '_image_allW.tif'],'-dtiff',['-r' num2str(resolution)]);


                    for aa = 1:length(fieldnames(monitor_W_ess))
                    close all;
                    figure_handle = figure(1);
                    set(figure_handle, 'Position', [0 0 1800 1800])
                    axes1 = axes('Parent',figure_handle,'YDir','reverse','Position',[0 0 1 1],'Layer','top');
                    names = fieldnames(monitor_W_ess);
                    object = monitor_W_ess.(names{aa})(monitor_num_image);
                    signal_data = object.signal;
                    %x_data = object.x;
                    %y_data = object.y;
                    x_data = object{2};
                    y_data = object{1};

                    imagesc(x_data,y_data,signal_data)
                    colormap('gray')

                    axis off

                    %surf(peaks)
                    % control the image pixel size by manipulating the paper size and number of dots per inch
                    output_size = [1250 1250];%Size in pixels
                    resolution = 300;%Resolution in DPI
                    set(gcf,'paperunits','inches','paperposition',[0 0 output_size/resolution]);
                    % use 300 DPI
                    print([filename '_image_W' num2str(aa) '.tif'],'-dtiff',['-r' num2str(resolution)]);    
                    end
                end
                
                
                
                
                disp(['reflogic = ' num2str(reflogic)])
                %monitor_num = 11;
                %monitor_num_ess = 20;
                % Write txt files with wavelength and intensity /
                % brilliance transfer
                
                % Lmon naming
                ALLW_ifit = assign_by_title('Lmon_sample_B.',monitor_ALLW);
                ALLW_d_ifit = assign_by_title('Lmon_sample_B.',monitor_ALLW_degraded);
                ESS_ifit = assign_by_title('Lmon_sample.',monitor_ESSW);
                ESS_d_ifit = assign_by_title('Lmon_sample.',monitor_ESSW_degraded);
                ESS_lim_ifit = assign_by_title('Lmon_sample_B.',monitor_ESSW);
                ESS_d_lim_ifit = assign_by_title('Lmon_sample_B.',monitor_ESSW_degraded);
                if reflogic
                %brilliance_ifit = monitor_ALLW(monitor_num)/monitor_ALLW_ref(monitor_num);
                
                LAMBDA_B_ref=assign_by_title('Lmon_sample_B.',monitor_ALLW_ref);
                ALLW_B_ifit = ALLW_ifit/LAMBDA_B_ref.Data.Mean;
                ALLW_d_B_ifit = ALLW_d_ifit/LAMBDA_B_ref.Data.Mean;
                
                %LAMBDA_B_ref=assign_by_title('Lmon_sample_B.',monitor_ALLW_ref);
                %brilliance_ifit=assign_by_title('Lmon_sample_B.',monitor_ALLW)/LAMBDA_B_ref.Data.Mean;
                
                wavelength = ALLW_B_ifit{1};
                brilliance = ALLW_B_ifit.signal;
                
                %ESS_fom_ifit = monitor_ESSW(monitor_num);
                ESS_fom_ifit= assign_by_title('Lmon_sample_B.',monitor_ESSW);
                ESS_fom = ESS_fom_ifit.signal;
                
                %ESS_all_ifit = monitor_ESSW(monitor_num);
                ESS_all_ifit = assign_by_title('Lmon_sample.',monitor_ESSW);
                ESS_all = ESS_all_ifit.signal;
                
                
                % open file
                fid = fopen([temp_name '.txt'], 'w');
                
                % header line
                % Works for all monitor_num
                %mod_x = monitor_ALLW(monitor_num).data.mod_x;
                %mod_y = monitor_ALLW(monitor_num).data.mod_y;
                
                mod_x = ALLW_ifit.data.mod_x;
                mod_y = ALLW_ifit.data.mod_y;
                %sizeX = monitor_ALLW(monitor_num).data.sizeX;
                %sizeY = monitor_ALLW(monitor_num).data.sizeY;
                sizeX = ALLW_ifit.data.sizeX;
                sizeY = ALLW_ifit.data.sizeY;
                %divreq_x = monitor_ALLW(monitor_num).data.divreq_x;
                %divreq_y = monitor_ALLW(monitor_num).data.divreq_y;
                divreq_x = ALLW_ifit.data.divreq_x;
                divreq_y = ALLW_ifit.data.divreq_y;

                h_line = ['guide_bot inputstring = ' inputstring '\n'];
                fprintf(fid,h_line);
                h_line = ['moderator height [cm] = ' num2str(mod_y*100) '\n'];
                fprintf(fid,h_line);
                h_line = ['moderator width  [cm] = ' num2str(mod_x*100) '\n'];
                fprintf(fid,h_line);
                h_line = ['sample width  [cm] = ' num2str(sizeX*100) '\n'];
                fprintf(fid,h_line);
                h_line = ['sample height  [cm] = ' num2str(sizeY*100) '\n'];
                fprintf(fid,h_line);
                h_line = ['horizontal divergence (plus/minus) [deg] = ' num2str(divreq_x) '\n'];
                fprintf(fid,h_line);
                h_line = ['vertical divergence (plus/minus) [deg] = ' num2str(divreq_y) '\n'];
                fprintf(fid,h_line);
                
                h_line = ['-------- DATA ------------ \n'];
                fprintf(fid,h_line);
                h_line = ['wavelength [AA], brilliance transfer [unitless], intensity within fom [n/s] (in wavelength bin), intensity on sample [n/s] (in wavelength bin)\n'];
                fprintf(fid,h_line);
                
                % write data
                for ii = 1:length(wavelength)
                   fprintf(fid,[num2str(wavelength(ii),'%10.6f') '\t' num2str(brilliance(ii),'%6.5f') '\t' num2str(ESS_fom(ii),'%10e') '\t' num2str(ESS_all(ii),'%10e') '\n'])
                end
                
                fclose(fid)
                
                end
                
                for wavband = 1:length(Wmax)
                    %ALLW=xlim(monitor_ALLW(monitor_num),[Wmin(wavband) Wmax(wavband)]);
                    ALLW=xlim(ALLW_ifit,[Wmin(wavband) Wmax(wavband)]);
                    
                    ALLW_wave(instrument,wavband).signal(iindex,jindex)=sum(ALLW);
                    ALLW_wave(instrument,wavband).error(iindex,jindex)=sqrt(sum(ALLW.error.^2));
                    ALLW_wave(instrument,wavband).waveinfo = [Wmin(wavband) Wmax(wavband)];
                    if size_select > 0.5;
                     ALLW_wave(instrument,wavband).scan(1).name = scan_dim(1).name;
                     ALLW_wave(instrument,wavband).scan(1).value(iindex) = str2num(p.(scan_dim(1).name_mcstas));
                     if size_select == 2
                        ALLW_wave(instrument,wavband).scan(2).name = scan_dim(2).name;
                        ALLW_wave(instrument,wavband).scan(2).value(jindex) = str2num(p.(scan_dim(2).name_mcstas));
                     end
                    end
                    

                    if reflogic
                    ALLW_B = xlim(ALLW_B_ifit,[Wmin(wavband) Wmax(wavband)]);
                    ALLW_B_wave(instrument,wavband).signal(iindex,jindex)=mean(ALLW_B);
                    ALLW_B_wave(instrument,wavband).error(iindex,jindex)=sqrt(sum(ALLW_B.error.^2))/length(ALLW_B.error);
                    ALLW_B_wave(instrument,wavband).waveinfo = [Wmin(wavband) Wmax(wavband)];
                    if size_select > 0.5;
                     ALLW_B_wave(instrument,wavband).scan(1).name = scan_dim(1).name;
                     ALLW_B_wave(instrument,wavband).scan(1).value(iindex) = str2num(p.(scan_dim(1).name_mcstas));
                     if size_select == 2
                        ALLW_B_wave(instrument,wavband).scan(2).name = scan_dim(2).name;
                        ALLW_B_wave(instrument,wavband).scan(2).value(jindex) = str2num(p.(scan_dim(2).name_mcstas));
                     end
                    end
                    end
                    %ALLW_d=xlim(monitor_ALLW_degraded(monitor_num),[Wmin(wavband) Wmax(wavband)]);
                    ALLW_d=xlim(ALLW_d_ifit,[Wmin(wavband) Wmax(wavband)]);
                    ALLW_d_wave(instrument,wavband).signal(iindex,jindex)=sum(ALLW_d);
                    ALLW_d_wave(instrument,wavband).error(iindex,jindex)=sqrt(sum(ALLW_d.error.^2));
                    ALLW_d_wave(instrument,wavband).waveinfo = [Wmin(wavband) Wmax(wavband)];
                    if size_select > 0.5;
                     ALLW_d_wave(instrument,wavband).scan(1).name = scan_dim(1).name;
                     ALLW_d_wave(instrument,wavband).scan(1).value(iindex) = str2num(p.(scan_dim(1).name_mcstas));
                     if size_select == 2
                        ALLW_d_wave(instrument,wavband).scan(2).name = scan_dim(2).name;
                        ALLW_d_wave(instrument,wavband).scan(2).value(jindex) = str2num(p.(scan_dim(2).name_mcstas));
                     end
                    end

                    if reflogic
                    ALLW_d_B = xlim(ALLW_d_B_ifit,[Wmin(wavband) Wmax(wavband)]);
                    ALLW_d_B_wave(instrument,wavband).signal(iindex,jindex)=mean(ALLW_d_B);
                    ALLW_d_B_wave(instrument,wavband).error(iindex,jindex)=sqrt(sum(ALLW_d_B.error.^2))/length(ALLW_d_B.error);
                    ALLW_d_B_wave(instrument,wavband).waveinfo = [Wmin(wavband) Wmax(wavband)];
                    if size_select < 0.5;
                     ALLW_d_B_wave(instrument,wavband).name = name{instrument};
                        
                    else
                     ALLW_d_B_wave(instrument,wavband).scan(1).name = scan_dim(1).name;
                     ALLW_d_B_wave(instrument,wavband).scan(1).value(iindex) = str2num(p.(scan_dim(1).name_mcstas));
                     if size_select == 2
                        ALLW_d_B_wave(instrument,wavband).scan(2).name = scan_dim(2).name;
                        ALLW_d_B_wave(instrument,wavband).scan(2).value(jindex) = str2num(p.(scan_dim(2).name_mcstas));
                     end
                    end
                    end
                    ESS=xlim(ESS_ifit,[Wmin(wavband) Wmax(wavband)]);
                    ESS_wave(instrument,wavband).signal(iindex,jindex)=sum(ESS);
                    ESS_wave(instrument,wavband).error(iindex,jindex)=sqrt(sum(ESS.error.^2));
                    ESS_wave(instrument,wavband).waveinfo = [Wmin(wavband) Wmax(wavband)];
                    if size_select > 0.5;
                     ESS_wave(instrument,wavband).scan(1).name = scan_dim(1).name;
                     ESS_wave(instrument,wavband).scan(1).value(iindex) = str2num(p.(scan_dim(1).name_mcstas));
                     if size_select == 2
                        ESS_wave(instrument,wavband).scan(2).name = scan_dim(2).name;
                        ESS_wave(instrument,wavband).scan(2).value(jindex) = str2num(p.(scan_dim(2).name_mcstas));
                     end
                    end
                    
                    ESS_d=xlim(ESS_d_ifit,[Wmin(wavband) Wmax(wavband)]);
                    ESS_d_wave(instrument,wavband).signal(iindex,jindex)=sum(ESS_d);
                    ESS_d_wave(instrument,wavband).error(iindex,jindex)=sqrt(sum(ESS_d.error.^2));
                    ESS_d_wave(instrument,wavband).waveinfo = [Wmin(wavband) Wmax(wavband)];
                    if size_select > 0.5;
                     ESS_d_wave(instrument,wavband).scan(1).name = scan_dim(1).name;
                     ESS_d_wave(instrument,wavband).scan(1).value(iindex) = str2num(p.(scan_dim(1).name_mcstas));
                     if size_select == 2
                        ESS_d_wave(instrument,wavband).scan(2).name = scan_dim(2).name;
                        ESS_d_wave(instrument,wavband).scan(2).value(jindex) = str2num(p.(scan_dim(2).name_mcstas));
                     end
                    end
                    
                    ESS_lim=xlim(ESS_lim_ifit,[Wmin(wavband) Wmax(wavband)]);
                    ESS_lim_wave(instrument,wavband).signal(iindex,jindex)=sum(ESS_lim);
                    ESS_lim_wave(instrument,wavband).error(iindex,jindex)=sqrt(sum(ESS_lim.error.^2));
                    ESS_lim_wave(instrument,wavband).waveinfo = [Wmin(wavband) Wmax(wavband)];
                    if size_select > 0.5;
                     ESS_lim_wave(instrument,wavband).scan(1).name = scan_dim(1).name;
                     ESS_lim_wave(instrument,wavband).scan(1).value(iindex) = str2num(p.(scan_dim(1).name_mcstas));
                     if size_select == 2
                        ESS_lim_wave(instrument,wavband).scan(2).name = scan_dim(2).name;
                        ESS_lim_wave(instrument,wavband).scan(2).value(jindex) = str2num(p.(scan_dim(2).name_mcstas));
                     end
                    end
                    
                    ESS_d_lim=xlim(ESS_d_lim_ifit,[Wmin(wavband) Wmax(wavband)]);
                    ESS_d_lim_wave(instrument,wavband).signal(iindex,jindex)=sum(ESS_d_lim);
                    ESS_d_lim_wave(instrument,wavband).error(iindex,jindex)=sqrt(sum(ESS_d_lim.error.^2));
                    ESS_d_lim_wave(instrument,wavband).waveinfo = [Wmin(wavband) Wmax(wavband)];
                    if size_select > 0.5;
                     ESS_d_lim_wave(instrument,wavband).scan(1).name = scan_dim(1).name;
                     ESS_d_lim_wave(instrument,wavband).scan(1).value(iindex) = str2num(p.(scan_dim(1).name_mcstas));
                     if size_select == 2
                        ESS_d_lim_wave(instrument,wavband).scan(2).name = scan_dim(2).name;
                        ESS_d_lim_wave(instrument,wavband).scan(2).value(jindex) = str2num(p.(scan_dim(2).name_mcstas));
                     end
                    end
                end
                
            else
                % Do I need to do anything specific to ensure zeros in the
                % scan?
                
            end 
            progress_logic(instrument,iindex,jindex) = 1;
            if exist('ALLW_wave') > 0.5
            if exist('scan_dim') > 0.5 
            save('progress.mat','ALLW_wave','ALLW_d_wave','ESS_wave','ESS_d_wave','ESS_lim_wave','ESS_d_lim_wave','ALLW_B_wave','ALLW_d_B_wave','scan_dim','name','list','iindex','jindex','instrument','Wmax','Wmin','wavband','progress_logic','jindex_range','iindex_range','instr_range','size_select');
            else
            save('progress.mat','ALLW_wave','ALLW_d_wave','ESS_wave','ESS_d_wave','ESS_lim_wave','ESS_d_lim_wave','ALLW_B_wave','ALLW_d_B_wave','name','list','iindex','jindex','instrument','Wmax','Wmin','wavband','progress_logic','jindex_range','iindex_range','instr_range','size_select');    
            end
            end
        end
    end
end
end
% Includes all for debugging, can be reduced heavily.
clearvars -except ALLW_wave ALLW_d_wave ESS_wave ESS_d_wave ESS_lim_wave ESS_d_lim_wave ALLW_B_wave ALLW_d_B_wave scan_dim name list iindex jindex instrument Wmax Wmin wavband progress_logic jindex_range iindex_range instr_range size_select calculate plot_options


% Checking the .signal matrices are all the same size. If an entire row
% fails, it can lead to nasty results.

for instr = 1:length(ALLW_wave(:,1))
   size_test_array(instr,:) = size(ALLW_wave(instr,1).signal); 
end

max_first = max(size_test_array(:,1));
max_second = max(size_test_array(:,2));

for instr = 1:length(ALLW_wave(:,1))
   if size_test_array(instr,1) < max_first || size_test_array(instr,2) < max_second 
       for w_band = 1:length(ALLW_wave(1,:))
        ALLW_wave(instr,w_band).signal(max_first,max_second)=0;
        ALLW_B_wave(instr,w_band).signal(max_first,max_second)=0;
        ALLW_d_wave(instr,w_band).signal(max_first,max_second)=0;
        ALLW_d_B_wave(instr,w_band).signal(max_first,max_second)=0;
        ESS_wave(instr,w_band).signal(max_first,max_second)=0;
        ESS_d_wave(instr,w_band).signal(max_first,max_second)=0;
        ESS_lim_wave(instr,w_band).signal(max_first,max_second)=0;
        ESS_d_lim_wave(instr,w_band).signal(max_first,max_second)=0;
        
        ALLW_wave(instr,w_band).error(max_first,max_second)=0;
        ALLW_B_wave(instr,w_band).error(max_first,max_second)=0;
        ALLW_d_wave(instr,w_band).error(max_first,max_second)=0;
        ALLW_d_B_wave(instr,w_band).error(max_first,max_second)=0;
        ESS_wave(instr,w_band).error(max_first,max_second)=0;
        ESS_d_wave(instr,w_band).error(max_first,max_second)=0;
        ESS_lim_wave(instr,w_band).error(max_first,max_second)=0;
        ESS_d_lim_wave(instr,w_band).error(max_first,max_second)=0;
        
        ALLW_wave(instr,w_band).scan(2).value=[0.01 0.02 0.03 0.03 0.04];
        ALLW_B_wave(instr,w_band).scan(2).value=[0.01 0.02 0.03 0.03 0.04];
        ALLW_d_wave(instr,w_band).scan(2).value=[0.01 0.02 0.03 0.03 0.04];
        ALLW_d_B_wave(instr,w_band).scan(2).value=[0.01 0.02 0.03 0.03 0.04];
        ESS_wave(instr,w_band).scan(2).value=[0.01 0.02 0.03 0.03 0.04];
        ESS_d_wave(instr,w_band).scan(2).value=[0.01 0.02 0.03 0.03 0.04];
        ESS_lim_wave(instr,w_band).scan(2).value=[0.01 0.02 0.03 0.03 0.04];
        ESS_d_lim_wave(instr,w_band).scan(2).value=[0.01 0.02 0.03 0.03 0.04];
       end
   end
end



save('debug_pre_redundant.mat')

% Check for a scan dim which does not change.

if size_select == 2
    
   not_done = 1;
   instr_check = 1;
   instr_length = length(ALLW_wave(:,1));
   for instr = 1:instr_length
       scan_length(instr) = length(ALLW_wave(1,1).scan(1).value);
   end
   min_length = max(scan_length);
   while not_done
   values1 = ALLW_wave(instr_check,1).scan(1).value;
   check = values1 == 0;
   if (sum(check) == 0 && min_length == length(values1)) || instr_check == instr_length
        not_done = 0;
   else
       instr_check = instr_check + 1;
   end
   end
   
   if length(values1) > 1
       duplicates1 = values1(1)*ones(1,length(values1)) == values1;
       if sum(duplicates1) == length(values1)
           redundant(1)=1;
       else
           redundant(1)=0;
       end
   else
       redundant(1)=0;
   end
   
   %values2 = ALLW_wave(1,1).scan(2).value;
   not_done = 1;
   instr_check = 1;
   instr_length = length(ALLW_wave(:,1));
   for instr = 1:instr_length
       scan_length(instr) = length(ALLW_wave(1,1).scan(2).value);
   end
   min_length = max(scan_length);
   while not_done
   values2 = ALLW_wave(instr_check,1).scan(2).value;
   check = values2 == 0;
   if (sum(check) == 0 && min_length == length(values1)) || instr_check == instr_length
        not_done = 0;
   else
       instr_check = instr_check + 1;
   end
   end
   if length(values2) > 1
       duplicates2 = values2(1)*ones(1,length(values2)) == values2;
       if sum(duplicates2) == length(values2)
           redundant(2)=1;
       else
           redundant(2)=0;
       end
   else
       redundant(2)=0;
   end
   
elseif size_select == 1
   values1 = ALLW_wave(1,1).scan(1).value; 
   % BUG, if first is not complete this will not work
   
   not_done = 1;
   instr_check = 1;
   instr_length = length(ALLW_wave(:,1));
   for instr = 1:instr_length
       scan_length(instr) = length(ALLW_wave(1,1).scan.value);
   end
   min_length = max(scan_length);
   while not_done
   values1 = ALLW_wave(instr_check,1).scan(1).value;
   check = values1 == 0;
   if (sum(check) == 0 && min_length == length(values1)) || instr_check == instr_length
        not_done = 0;
   else
       instr_check = instr_check + 1;
   end
   end
   
   if length(values1) > 1
       duplicates1 = values1(1)*ones(1,length(values1)) == values1;
       if sum(duplicates1) == length(values1)
           redundant(1)=1;
       else
           redundant(1)=0;
       end
   else
       redundant(1)=0;
   end
   redundant(2)=0;
elseif size_select == 0
    redundant(1) = 0;
    redundant(2) = 0;
end 



if redundant(1)
    % for now only sort by absolute intensities
    sort_band = 1;
    if exist('ALLW_B_wave')>0.5 % sort by briliance transfer 
        for instr = 1:length(ALLW_wave(:,1))
            for jindex = 1:length(ALLW_wave(instr,1).signal(1,:))
                zero_test = ALLW_B_wave(instr,sort_band).signal(:,jindex) == 0;
                if sum(zero_test) == length(zero_test)
                [val(instr,jindex) ind(instr,jindex)] = max(ALLW_wave(instr,sort_band).signal(:,jindex));
                else
                [val(instr,jindex) ind(instr,jindex)] = max(ALLW_B_wave(instr,sort_band).signal(:,jindex));
                end
                if ind(instr,jindex) == 0
                    disp('problem index')
                    disp(instr)
                    disp(jindex)
                end 
            end
        end
    else % sort by absolute intensities
        for instr = 1:length(ALLW_wave(:,1))
            for jindex = 1:length(ALLW_wave(instr,1).signal(1,:))
                [val(instr,jindex) ind(instr,jindex)] = max(ALLW_wave(instr,sort_band).signal(:,jindex));
            end
        end
    end
    for w_band = 1:length(ALLW_wave(1,:))
        for instr=1:length(ALLW_wave(:,1))
            for jindex = 1:length(ALLW_wave(1,1).signal(1,:))
                ALLW_wave_reduced(instr,w_band).signal(jindex) = ALLW_wave(instr,w_band).signal(ind(instr,jindex),jindex);
                ALLW_wave_reduced(instr,w_band).error(jindex) = ALLW_wave(instr,w_band).error(ind(instr,jindex),jindex);
                ALLW_wave_reduced(instr,w_band).waveinfo = ALLW_wave(instr,w_band).waveinfo;
                if size_select == 2
                   ALLW_wave_reduced(instr,w_band).scan(1).name = ALLW_wave(instr,w_band).scan(2).name;
                   ALLW_wave_reduced(instr,w_band).scan(1).value(jindex) = ALLW_wave(instr,w_band).scan(2).value(jindex);
                end
                ALLW_d_wave_reduced(instr,w_band).signal(jindex) = ALLW_d_wave(instr,w_band).signal(ind(instr,jindex),jindex);
                ALLW_d_wave_reduced(instr,w_band).error(jindex) = ALLW_d_wave(instr,w_band).error(ind(instr,jindex),jindex);
                ALLW_d_wave_reduced(instr,w_band).waveinfo = ALLW_d_wave(instr,w_band).waveinfo;
                if size_select == 2
                   ALLW_d_wave_reduced(instr,w_band).scan(1).name = ALLW_d_wave(instr,w_band).scan(2).name;
                   ALLW_d_wave_reduced(instr,w_band).scan(1).value(jindex) = ALLW_d_wave(instr,w_band).scan(2).value(jindex);
                end
                ESS_wave_reduced(instr,w_band).signal(jindex) = ESS_wave(instr,w_band).signal(ind(instr,jindex),jindex);
                ESS_wave_reduced(instr,w_band).error(jindex) = ESS_wave(instr,w_band).error(ind(instr,jindex),jindex);
                ESS_wave_reduced(instr,w_band).waveinfo = ESS_wave(instr,w_band).waveinfo;
                if size_select == 2
                   ESS_wave_reduced(instr,w_band).scan(1).name = ESS_wave(instr,w_band).scan(2).name;
                   ESS_wave_reduced(instr,w_band).scan(1).value(jindex) = ESS_wave(instr,w_band).scan(2).value(jindex);
                end
                ESS_d_wave_reduced(instr,w_band).signal(jindex) = ESS_d_wave(instr,w_band).signal(ind(instr,jindex),jindex);
                ESS_d_wave_reduced(instr,w_band).error(jindex) = ESS_d_wave(instr,w_band).error(ind(instr,jindex),jindex);
                ESS_d_wave_reduced(instr,w_band).waveinfo = ESS_d_wave(instr,w_band).waveinfo;
                if size_select == 2
                   ESS_d_wave_reduced(instr,w_band).scan(1).name = ESS_d_wave(instr,w_band).scan(2).name;
                   ESS_d_wave_reduced(instr,w_band).scan(1).value(jindex) = ESS_d_wave(instr,w_band).scan(2).value(jindex);
                end
                ESS_lim_wave_reduced(instr,w_band).signal(jindex) = ESS_lim_wave(instr,w_band).signal(ind(instr,jindex),jindex);
                ESS_lim_wave_reduced(instr,w_band).error(jindex) = ESS_lim_wave(instr,w_band).error(ind(instr,jindex),jindex);
                ESS_lim_wave_reduced(instr,w_band).waveinfo = ESS_lim_wave(instr,w_band).waveinfo;
                if size_select == 2
                   ESS_lim_wave_reduced(instr,w_band).scan(1).name = ESS_lim_wave(instr,w_band).scan(2).name;
                   ESS_lim_wave_reduced(instr,w_band).scan(1).value(jindex) = ESS_lim_wave(instr,w_band).scan(2).value(jindex);
                end
                ESS_d_lim_wave_reduced(instr,w_band).signal(jindex) = ESS_d_lim_wave(instr,w_band).signal(ind(instr,jindex),jindex);
                ESS_d_lim_wave_reduced(instr,w_band).error(jindex) = ESS_d_lim_wave(instr,w_band).error(ind(instr,jindex),jindex);
                ESS_d_lim_wave_reduced(instr,w_band).waveinfo = ESS_d_lim_wave(instr,w_band).waveinfo;
                if size_select == 2
                   ESS_d_lim_wave_reduced(instr,w_band).scan(1).name = ESS_d_lim_wave(instr,w_band).scan(2).name;
                   ESS_d_lim_wave_reduced(instr,w_band).scan(1).value(jindex) = ESS_d_lim_wave(instr,w_band).scan(2).value(jindex);
                end
                if exist('ALLW_B_wave')>0.5
                   ALLW_B_wave_reduced(instr,w_band).signal(jindex) = ALLW_B_wave(instr,w_band).signal(ind(instr,jindex),jindex);
                   ALLW_B_wave_reduced(instr,w_band).error(jindex) = ALLW_B_wave(instr,w_band).error(ind(instr,jindex),jindex);
                   ALLW_B_wave_reduced(instr,w_band).waveinfo = ALLW_B_wave(instr,w_band).waveinfo;
                   if size_select == 2
                      ALLW_B_wave_reduced(instr,w_band).scan(1).name = ALLW_B_wave(instr,w_band).scan(2).name;
                      ALLW_B_wave_reduced(instr,w_band).scan(1).value(jindex) = ALLW_B_wave(instr,w_band).scan(2).value(jindex);
                   end
                   ALLW_d_B_wave_reduced(instr,w_band).signal(jindex) = ALLW_d_B_wave(instr,w_band).signal(ind(instr,jindex),jindex);
                   ALLW_d_B_wave_reduced(instr,w_band).error(jindex) = ALLW_d_B_wave(instr,w_band).error(ind(instr,jindex),jindex);
                   ALLW_d_B_wave_reduced(instr,w_band).waveinfo = ALLW_d_B_wave(instr,w_band).waveinfo;
                   if size_select == 2
                      ALLW_d_B_wave_reduced(instr,w_band).scan(1).name = ALLW_d_B_wave(instr,w_band).scan(2).name;
                      ALLW_d_B_wave_reduced(instr,w_band).scan(1).value(jindex) = ALLW_d_B_wave(instr,w_band).scan(2).value(jindex);
                   end
                end 
            end
        end
    end
elseif redundant(2)
    
    sort_band = 1;
    if exist('ALLW_B') % sort by briliance transfer 
        for instr = 1:length(ALLW_wave(:,1))
            for iindex = 1:length(ALLW_wave(instr,1).signal(:,1))
                % Needs a check to see if all the values are zero!
                zero_test = ALLW_B_wave(instr,sort_band).signal(iindex,:) == 0;
                if sum(zero_test) == length(zero_test)
                [val(instr,iindex) ind(instr,iindex)] = max(ALLW_wave(instr,sort_band).signal(iindex,:));
                else
                [val(instr,iindex) ind(instr,iindex)] = max(ALLW_B_wave(instr,sort_band).signal(iindex,:));
                end
            end
        end
    else % sort by absolute intensities
        for instr = 1:length(ALLW_wave(:,1))
            for iindex = 1:length(ALLW_wave(instr,1).signal(:,1))
                [val(instr,iindex) ind(instr,iindex)] = max(ALLW_wave(instr,sort_band).signal(iindex,:));
            end
        end
    end
    for w_band = 1:length(ALLW_wave(1,:))
        for instr=1:length(ALLW_wave(:,1))
            for iindex = 1:length(ALLW_wave(1,1).signal(:,1))
                ALLW_wave_reduced(instr,w_band).signal(iindex) = ALLW_wave(instr,w_band).signal(iindex,ind(instr,iindex));
                ALLW_wave_reduced(instr,w_band).error(iindex) = ALLW_wave(instr,w_band).error(iindex,ind(instr,iindex));
                ALLW_wave_reduced(instr,w_band).waveinfo = ALLW_wave(instr,w_band).waveinfo;
                if size_select == 2
                   ALLW_wave_reduced(instr,w_band).scan(1).name = ALLW_wave(instr,w_band).scan(1).name;
                   ALLW_wave_reduced(instr,w_band).scan(1).value(iindex) = ALLW_wave(instr,w_band).scan(1).value(iindex);
                end
                ALLW_d_wave_reduced(instr,w_band).signal(iindex) = ALLW_d_wave(instr,w_band).signal(iindex,ind(instr,iindex));
                ALLW_d_wave_reduced(instr,w_band).error(iindex) = ALLW_d_wave(instr,w_band).error(iindex,ind(instr,iindex));
                ALLW_d_wave_reduced(instr,w_band).waveinfo = ALLW_d_wave(instr,w_band).waveinfo;
                if size_select == 2
                   ALLW_d_wave_reduced(instr,w_band).scan(1).name = ALLW_d_wave(instr,w_band).scan(1).name;
                   ALLW_d_wave_reduced(instr,w_band).scan(1).value(iindex) = ALLW_d_wave(instr,w_band).scan(1).value(iindex);
                end
                ESS_wave_reduced(instr,w_band).signal(iindex) = ESS_wave(instr,w_band).signal(iindex,ind(instr,iindex));
                ESS_wave_reduced(instr,w_band).error(iindex) = ESS_wave(instr,w_band).error(iindex,ind(instr,iindex));
                ESS_wave_reduced(instr,w_band).waveinfo = ESS_wave(instr,w_band).waveinfo;
                if size_select == 2
                   ESS_wave_reduced(instr,w_band).scan(1).name = ESS_wave(instr,w_band).scan(1).name;
                   ESS_wave_reduced(instr,w_band).scan(1).value(iindex) = ESS_wave(instr,w_band).scan(1).value(iindex);
                end
                ESS_d_wave_reduced(instr,w_band).signal(iindex) = ESS_d_wave(instr,w_band).signal(iindex,ind(instr,iindex));
                ESS_d_wave_reduced(instr,w_band).error(iindex) = ESS_d_wave(instr,w_band).error(iindex,ind(instr,iindex));
                ESS_d_wave_reduced(instr,w_band).waveinfo = ESS_d_wave(instr,w_band).waveinfo;
                if size_select == 2
                   ESS_d_wave_reduced(instr,w_band).scan(1).name = ESS_d_wave(instr,w_band).scan(1).name;
                   ESS_d_wave_reduced(instr,w_band).scan(1).value(iindex) = ESS_d_wave(instr,w_band).scan(1).value(iindex);
                end
                ESS_lim_wave_reduced(instr,w_band).signal(iindex) = ESS_lim_wave(instr,w_band).signal(iindex,ind(instr,iindex));
                ESS_lim_wave_reduced(instr,w_band).error(iindex) = ESS_lim_wave(instr,w_band).error(iindex,ind(instr,iindex));
                ESS_lim_wave_reduced(instr,w_band).waveinfo = ESS_lim_wave(instr,w_band).waveinfo;
                if size_select == 2
                   ESS_lim_wave_reduced(instr,w_band).scan(1).name = ESS_lim_wave(instr,w_band).scan(1).name;
                   ESS_lim_wave_reduced(instr,w_band).scan(1).value(iindex) = ESS_lim_wave(instr,w_band).scan(1).value(iindex);
                end
                ESS_d_lim_wave_reduced(instr,w_band).signal(iindex) = ESS_d_lim_wave(instr,w_band).signal(iindex,ind(instr,iindex));
                ESS_d_lim_wave_reduced(instr,w_band).error(iindex) = ESS_d_lim_wave(instr,w_band).error(iindex,ind(instr,iindex));
                ESS_d_lim_wave_reduced(instr,w_band).waveinfo = ESS_d_lim_wave(instr,w_band).waveinfo;
                if size_select == 2
                   ESS_d_lim_wave_reduced(instr,w_band).scan(1).name = ESS_d_lim_wave(instr,w_band).scan(1).name;
                   ESS_d_lim_wave_reduced(instr,w_band).scan(1).value(iindex) = ESS_d_lim_wave(instr,w_band).scan(1).value(iindex);
                end
                if exist('ALLW_B_wave')>0.5
                   ALLW_B_wave_reduced(instr,w_band).signal(iindex) = ALLW_B_wave(instr,w_band).signal(iindex,ind(instr,iindex));
                   ALLW_B_wave_reduced(instr,w_band).error(iindex) = ALLW_B_wave(instr,w_band).error(iindex,ind(instr,iindex));
                   ALLW_B_wave_reduced(instr,w_band).waveinfo = ALLW_B_wave(instr,w_band).waveinfo;
                   if size_select == 2
                      ALLW_B_wave_reduced(instr,w_band).scan(1).name = ALLW_B_wave(instr,w_band).scan(1).name;
                      ALLW_B_wave_reduced(instr,w_band).scan(1).value(iindex) = ALLW_B_wave(instr,w_band).scan(1).value(iindex);
                   end
                   ALLW_d_B_wave_reduced(instr,w_band).signal(iindex) = ALLW_d_B_wave(instr,w_band).signal(iindex,ind(instr,iindex));
                   ALLW_d_B_wave_reduced(instr,w_band).error(iindex) = ALLW_d_B_wave(instr,w_band).error(iindex,ind(instr,iindex));
                   ALLW_d_B_wave_reduced(instr,w_band).waveinfo = ALLW_d_B_wave(instr,w_band).waveinfo;
                   if size_select == 2
                      ALLW_d_B_wave_reduced(instr,w_band).scan(1).name = ALLW_d_B_wave(instr,w_band).scan(1).name;
                      ALLW_d_B_wave_reduced(instr,w_band).scan(1).value(iindex) = ALLW_d_B_wave(instr,w_band).scan(1).value(iindex);
                   end
                end 
            end
        end
    end
    
end



% code for writing which files were selected when reducing the redundant
% dimension.

if exist('ind') > 0
    % open file
    fid = fopen(['selected_by_scan_plotter.txt'], 'w');

    
    for ii = 1:length(ind(:,1))
        h_line = [ name{ii} ' index ' num2str(ind(ii,:))  '\n'];
        fprintf(fid,h_line);
    end

    fclose(fid)
    
end



reduce_redundant = 1;
if reduce_redundant && (redundant(1) || redundant(2))
    ALLW_wave = ALLW_wave_reduced;
    ALLW_d_wave = ALLW_d_wave_reduced;
    ESS_wave = ESS_wave_reduced;
    ESS_d_wave = ESS_d_wave_reduced;
    ESS_lim_wave = ESS_lim_wave_reduced;
    ESS_d_lim_wave = ESS_d_lim_wave_reduced;
    
    for w_band = 1:length(ALLW_wave(1,:))
      for instr=1:length(ALLW_wave(:,1))
        ALLW_wave(instr,w_band).signal = ALLW_wave(instr,w_band).signal';
        ALLW_wave(instr,w_band).error = ALLW_wave(instr,w_band).error';
        ALLW_d_wave(instr,w_band).signal = ALLW_d_wave(instr,w_band).signal';
        ALLW_d_wave(instr,w_band).error = ALLW_d_wave(instr,w_band).error';
        ESS_wave(instr,w_band).signal = ESS_wave(instr,w_band).signal';
        ESS_wave(instr,w_band).error = ESS_wave(instr,w_band).error';
        ESS_d_wave(instr,w_band).signal = ESS_d_wave(instr,w_band).signal';
        ESS_d_wave(instr,w_band).error = ESS_d_wave(instr,w_band).error';
        ESS_lim_wave(instr,w_band).signal = ESS_lim_wave(instr,w_band).signal';
        ESS_lim_wave(instr,w_band).error = ESS_lim_wave(instr,w_band).error';
        ESS_d_lim_wave(instr,w_band).signal = ESS_d_lim_wave(instr,w_band).signal';
        ESS_d_lim_wave(instr,w_band).error = ESS_d_lim_wave(instr,w_band).error';
      end
    end
    
    if exist('ALLW_B_wave')>0.5
        ALLW_B_wave = ALLW_B_wave_reduced;
        ALLW_d_B_wave = ALLW_d_B_wave_reduced;
        
      for w_band = 1:length(ALLW_wave(1,:))
        for instr=1:length(ALLW_wave(:,1))
        ALLW_B_wave(instr,w_band).signal = ALLW_B_wave(instr,w_band).signal';
        ALLW_B_wave(instr,w_band).error = ALLW_B_wave(instr,w_band).error';
        end
      end
        disp('overwrote')
    end
    if redundant(1)
        if size_select == 2
        scan_dim_tmp = scan_dim(2);
        clear scan_dim
        scan_dim = scan_dim_tmp;
        end
    elseif redundant(2)
        scan_dim_tmp = scan_dim(1);
        clear scan_dim
        scan_dim = scan_dim_tmp;
    end
    size_select = size_select - 1;
end

% flux hardcoded for CAMEA
%    for w_band = 1:length(ALLW_wave(1,:))
%      for instr=1:length(ALLW_wave(:,1))
%        ALLW_wave(instr,w_band).signal(1,:) = ALLW_wave(instr,w_band).signal(1,:)./(1.5*0.5)';
%        ALLW_wave(instr,w_band).signal(2,:) = ALLW_wave(instr,w_band).signal(2,:)./(1.5*1.0)';
%        ALLW_wave(instr,w_band).signal(3,:) = ALLW_wave(instr,w_band).signal(3,:)./(1.5*1.5)';
%        ALLW_wave(instr,w_band).signal(4,:) = ALLW_wave(instr,w_band).signal(4,:)./(1.5*2.0)';%
%        ESS_wave(instr,w_band).signal(1,:) = ESS_wave(instr,w_band).signal(1,:)./(1.5*0.5)';
%        ESS_wave(instr,w_band).signal(2,:) = ESS_wave(instr,w_band).signal(2,:)./(1.5*1.0)';
%        ESS_wave(instr,w_band).signal(3,:) = ESS_wave(instr,w_band).signal(3,:)./(1.5*1.5)';
%        ESS_wave(instr,w_band).signal(4,:) = ESS_wave(instr,w_band).signal(4,:)./(1.5*2.0)';
%      end
%    end
    


save('debug_post_redundant.mat')


clc;close all;
% Code for plotting
    % How many paramaters was scanned?
    %   if 1
    %       show data as lines as a function of scanned variable
    %       line color corresponds to instrument (legend)
    %       different plots for the different sources / degraded.
    %   if 2
    %       can be shown as n*4 2D plots where n is the number of
    %       instruments.
    %       Can also be shown as above with n*4 plots with line color
    %       showing the 2nd scanned variable instead of the instrument
    %       Can also have different 1d graphs for on of the two scanned
    %       variables.

    % For a certain case:

wavelength_range = plot_options.wavelength_range;    
    
colormap_choice = plot_options.colormap_choice;    
    
plot_results = plot_options.plot_results;

if size_select == 2
    force_1D = plot_options.force_1D;
else
    force_1D = 0;
end
    

reverse_order = plot_options.reverse_order;

if size_select == 2;
    reverse_order = plot_options.reverse_order;
else
    reverse_order = 0;
end

FSL = plot_options.FSL;
FS = plot_options.FS;

if isfield(plot_options,'keyword')
keyword = plot_options.keyword;
end

if isfield(plot_options,'plot_logical')
plot_logical = plot_options.plot_logical;
end

if isfield(plot_options,'nickname')
    nickname = plot_options.nickname;
end


cleanup_low_points = plot_options.cleanup_low_points;
tollerance = plot_options.tollerance;

apply_flux = plot_options.apply_flux;


num_wavelength_intervals = length(ALLW_wave(1,:));

if wavelength_range > num_wavelength_intervals
   if ispc == 0
       disp('collecting results to pdf')
       if exist('nickname') > 0.5
        line{1}='#!/bin/bash';
        line{2} = ['convert -page A4 scan_plotter_ESS_*_' nickname '.png scan_plotter_collected_ESS_' nickname '.pdf'];
        line{3} = ['convert -page A4 scan_plotter_ESS_div_lim_*_' nickname '.png scan_plotter_collected_ESS_div_lim_' nickname '.pdf'];
        line{4} = ['convert -page A4 scan_plotter_Brill_*_' nickname '.png scan_plotter_collected_Brill_' nickname '.pdf'];
       else
        line{1}='#!/bin/bash';
        line{2} = ['convert -page A4 scan_plotter_ESS_*.png scan_plotter_collected_ESS.pdf'];
        line{3} = ['convert -page A4 scan_plotter_ESS_div_lim_*.png scan_plotter_collected_ESS_div_lim.pdf'];
        line{4} = ['convert -page A4 scan_plotter_Brill_*.png scan_plotter_collected_Brill.pdf'];
       end
       
        file_to_write='';
        
        for i=1:length(line)
          file_to_write = [file_to_write line{i} '\n'];
        end
        

        fid = fopen(['./generate_pdf.sh'], 'w');
        fprintf(fid,file_to_write);
        fclose(fid);
        unix(['chmod 744 ./generate_pdf.sh']);
   end
end


% Code for plotting, this can be changed according to taste.
moderator_flux = [2.4899 2.6280 2.5536 2.4474 2.3677 2.2747 2.2030 2.0968 2.0198 1.9852 1.9029 1.8392 1.7993 1.7223 1.6958 1.6320 1.5683 1.5311 1.5205 1.4860 1.4355 1.3957 1.3664 1.3107 1.2841 1.2629 1.2230];
moderator_flux = moderator_flux./moderator_flux(23);
moderator_heights = linspace(0.01,0.14,27);
plot_names = {'arb_unit' 'arb_unit_deg' 'ESS' 'ESS_deg' 'Brill' 'Brill_deg' 'ESS_div_lim' 'ESS_deg_div_lim'};

set(0,'DefaultAxesFontName', 'Times New Roman')
set(0,'DefaultAxesFontSize',FS)


if plot_results    
if (size_select == 1 || force_1D == 1) && size_select ~= 0
    num_instruments = length(name);
    temp_size = size(ALLW_wave(instrument,1).signal);
    if reverse_order == 1
    parts = temp_size(1);%1
    else
    parts = temp_size(2);%2ne
    end
    
    if num_instruments == 1
        collect = 1;
        colormatrix = colormap(colormap_choice);
        colorlength = length(colormatrix);
        skip = floor(colorlength/parts);
        colorpoint = skip*(1:parts);
        colorpoint = [1 colorpoint];
    else
        collect = 0;
        colormatrix = colormap(colormap_choice);
        colorlength = length(colormatrix);
        skip = floor(colorlength/num_instruments);
        colorpoint = skip*(1:num_instruments);
        colorpoint = [1 colorpoint];
    end

    for plots = plot_options.fig_choice% 5:6
        figure_handle = figure;
        max_val_ever = 0;
        switch plots
            case 1
                monitor = ALLW_wave;
                title_string = 'ALLW';
                yaxis = 'Intensity on sample, ARB unit';
            case 2
                monitor = ALLW_d_wave;
                title_string = 'ALLW_d';
                yaxis = 'Intensity on sample, ARB unit';
            case 3
                monitor = ESS_wave;
                title_string = 'ESS';
                yaxis = 'Intensity on sample n/s';
                %yaxis = 'Flux on sample n/s/cm^2';
            case 4
                monitor = ESS_d_wave;
                title_string = 'ESS_d';
                yaxis = 'Intensity on sample n/s';
            case 5
                monitor = ALLW_B_wave;
                title_string = 'Brilliance.';
                yaxis = 'Brilliance transfer';
                for i = 1:num_instruments
                    for j = 1:1 %wavelength ranges
                        %monitor(i,j).scan(2) = ALLW_wave(i,j).scan(2);
                    end
                end
                hold on
                %plot([0 1],[0 1],'--','linewidth',2,'color',[0.45 0.45 0.45])
                %plot([1 1],[0 1],'-.','linewidth',2,'color',[0.45 0.45 0.45])
            case 6
                monitor = ALLW_d_B_wave;
                title_string = 'Brilliance d';
                yaxis = 'Brilliance transfer';
                for i = 1:num_instruments
                    for j = 1:1 %wavelength ranges
                        %monitor(i,j).scan(2) = ALLW_wave(i,j).scan(2);
                    end
                end
            case 7
                monitor = ESS_lim_wave;
                title_string = 'ESS (div limit)';
                yaxis = 'Intensity on sample n/s (within div bounds)';
                %yaxis = 'Flux on sample n/s/cm^2';
            case 8
                monitor = ESS_d_lim_wave;
                title_string = 'ESS_d (div limit)';
                yaxis = 'Intensity on sample n/s (within div bounds)';
        end
        
    for part = 1:parts
        %if collect == 0 && part > 1; figure; else; hold on; end;
        if collect == 0;  subplot(parts,1,part); box on;else; hold on; end;
        
        for instrument = 1:num_instruments
            hold on
            if collect == 1
                if reverse_order
                cut_zero = monitor(instrument,wavelength_range).signal(part,:) ~= 0;
                %cut_zero = monitor(instrument,1).signal(:,part) ~= -1;
                scan_number = 2;
                if strcmp(monitor(instrument,wavelength_range).scan(scan_number).name,'moderator_size_y') && apply_flux
                    for index = 1:length(monitor(instrument,wavelength_range).scan(scan_number).value)
                        if cut_zero(index)
                        smaller_part = monitor(instrument,wavelength_range).scan(scan_number).value(index) >= moderator_heights; % effectivly rounds off
                        flux_index = sum(smaller_part);
                        %monitor(instrument,wavelength_range).signal(part,index) = monitor(instrument,wavelength_range).signal(part,index)*moderator_flux(flux_index);
                        end
                    end
                end
                
                X_plot = monitor(instrument,wavelength_range).scan(2).value(cut_zero);
                Y_plot = monitor(instrument,wavelength_range).signal(part,cut_zero);
                Y_error_plot = monitor(instrument,wavelength_range).error(part,cut_zero).*0;
                color_input = colormatrix(colorpoint(part),:);
                
                %handle = errorbar(monitor(instrument,1).scan(2).value(cut_zero),monitor(instrument,1).signal(part,cut_zero),monitor(instrument,1).error(part,cut_zero).*0,'-','color',colormatrix(colorpoint(part),:));    
                %handle = errorbar(X_plot,Y_plot,Y_error_plot,'-','color',color_input);
                else
                cut_zero = monitor(instrument,1).signal(:,part) ~= 0;
                %cut_zero = monitor(instrument,1).signal(:,part) ~= -1;
                
                scan_number = 1;
                if strcmp(monitor(instrument,wavelength_range).scan(scan_number).name,'moderator_size_y') && apply_flux
                    for index = 1:length(monitor(instrument,wavelength_range).scan(scan_number).value)
                        if cut_zero(index)
                        smaller_part = monitor(instrument,wavelength_range).scan(scan_number).value(index) >= moderator_heights; % effectivly rounds off
                        flux_index = sum(smaller_part);
                        %monitor(instrument,wavelength_range).signal(index,part) = monitor(instrument,wavelength_range).signal(index,part)*moderator_flux(flux_index);
                        end
                    end
                end
                
                X_plot = monitor(instrument,wavelength_range).scan(1).value(cut_zero);
                Y_plot = monitor(instrument,wavelength_range).signal(cut_zero,part);
                Y_error_plot = monitor(instrument,wavelength_range).error(cut_zero,part).*0;
                color_input = colormatrix(colorpoint(part),:);
                
                %handle = errorbar(monitor(instrument,1).scan(1).value(cut_zero),monitor(instrument,1).signal(cut_zero,part),monitor(instrument,1).error(cut_zero,part).*0,'-','color',colormatrix(colorpoint(part),:));
                %handle = errorbar(X_plot,Y_plot,Y_error_plot,'-','color',color_input);
                end
                %set(handle,'Marker','o');
                %set(handle,'MarkerSize',5)
                %set(handle,'LineWidth',2);
            else
                
                if reverse_order
                
                cut_zero = monitor(instrument,wavelength_range).signal(part,:) ~= 0;
                %cut_zero = monitor(instrument,1).signal(:,part) ~= -1;
                
                scan_number = 2;
                if strcmp(monitor(instrument,wavelength_range).scan(scan_number).name,'moderator_size_y') && apply_flux
                    for index = 1:length(monitor(instrument,wavelength_range).scan(scan_number).value)
                        if cut_zero(index)
                        smaller_part = monitor(instrument,wavelength_range).scan(scan_number).value(index) >= moderator_heights; % effectivly rounds off
                        flux_index = sum(smaller_part);
                        %monitor(instrument,wavelength_range).signal(part,index) = monitor(instrument,wavelength_range).signal(part,index)*moderator_flux(flux_index);
                        end
                    end
                end
                
                X_plot = monitor(instrument,wavelength_range).scan(2).value(cut_zero);
                Y_plot = monitor(instrument,wavelength_range).signal(part,cut_zero);
                Y_error_plot = monitor(instrument,wavelength_range).error(part,cut_zero);%.*0;
                color_input = colormatrix(colorpoint(instrument),:);
                
                %handle = errorbar(monitor(instrument,1).scan(2).value(cut_zero),monitor(instrument,1).signal(part,cut_zero),monitor(instrument,1).error(part,cut_zero).*0,'-','color',colormatrix(colorpoint(instrument),:));    
                %handle = errorbar(X_plot,Y_plot,Y_error_plot,'-','color',color_input);    
                else
                cut_zero = monitor(instrument,wavelength_range).signal(:,part) ~= 0;
                %cut_zero = monitor(instrument,1).signal(:,part) ~= -1;
                
                
                scan_number = 1;
                if strcmp(monitor(instrument,wavelength_range).scan(scan_number).name,'moderator_size_y') && apply_flux
                    for index = 1:length(monitor(instrument,wavelength_range).scan(scan_number).value)
                        if cut_zero(index)
                        smaller_part = monitor(instrument,wavelength_range).scan(scan_number).value(index) >= moderator_heights-0.0001; % effectivly rounds off
                        flux_index = sum(smaller_part);
                        %monitor(instrument,wavelength_range).signal(index,part) = monitor(instrument,wavelength_range).signal(index,part)*moderator_flux(flux_index);
                        end
                    end
                end
                
                X_plot = monitor(instrument,wavelength_range).scan(1).value(cut_zero);
                Y_plot = monitor(instrument,wavelength_range).signal(cut_zero,part);
                Y_error_plot = monitor(instrument,wavelength_range).error(cut_zero,part);%.*0;
                color_input = colormatrix(colorpoint(instrument),:);
                
                %handle = errorbar(monitor(instrument,1).scan(1).value(cut_zero),monitor(instrument,1).signal(cut_zero,part),monitor(instrument,1).error(cut_zero,part).*0,'-','color',colormatrix(colorpoint(instrument),:));
                %handle = errorbar(X_plot,Y_plot,Y_error_plot,'-','color',color_input);
                end
                %set(handle,'Marker','o')
                %set(handle,'MarkerSize',5)
                %set(handle,'LineWidth',2);
            end
            
           
            
            if cleanup_low_points
            clear good_point
            plot_length = length(X_plot);
            good_point(1:plot_length) = true;
            
            for index = 2:plot_length-1
                good_point(index) = Y_plot(index) > Y_plot(index-1)*tollerance || Y_plot(index) > Y_plot(index+1)*tollerance;
            end
            
            else    
               good_point = true(size(X_plot));
            end
            
            global_Y_plot{plots,instrument} = Y_plot;    
            
            X_plot = X_plot(good_point);
            Y_plot = Y_plot(good_point);
            Y_error_plot = Y_error_plot(good_point);
            
            
            %X_plot_save{part} = X_plot;
            %Y_plot_save{part} = Y_plot;
            
            ERROR_C_X = size(X_plot);
            ERROR_C_Y = size(Y_plot);
            
            Dont_plot = 0;
            if (ERROR_C_X(1) == 0 || ERROR_C_X(2) == 0)
                Dont_plot = 1;
            end
            if (ERROR_C_Y(1) == 0 || ERROR_C_Y(2) == 0)
                Dont_plot = 1;
            end
            
            if Dont_plot == 0
            %handle = errorbar(X_plot(good_point),Y_plot(good_point),Y_error_plot(good_point),'-','color',color_input);
            
            %if instrument == 2
            %    hold off
            %else
                hold on
            %end
            
            hide_this_plot = 0;
            if exist('plot_logical')>0.5
                if isfield(plot_logical,name{instrument})
                    if plot_logical.(name{instrument}) == 0
                        hide_this_plot = 1;
                    end
                end
            end
            
            
            if hide_this_plot == 0
                % surpress my wrong error calculations
            handle = errorbar(X_plot,Y_plot,Y_error_plot.*0,'-','color',color_input);
            box on
            %handle = errorbar(X_plot{part},Y_plot{part},Y_error_plot{part},'-','color',color_input);
            set(handle,'Marker','o')
            set(handle,'MarkerSize',5)
            set(handle,'LineWidth',2);
            %set(gca,'FontSize',13)
            end
            
            max_val_this = max(Y_plot);
            
            if max_val_this > max_val_ever
                max_val_ever = max_val_this;
            end
            if max_val_this ~= 0
            ylim([0 max_val_ever*1.05]);
            %ylim([0 6*10^10]);
            end
            end
        end
        
        
        
        hold off
    
        set(figure_handle, 'Position', [0 0 600 700])
        set(figure_handle, 'paperpositionmode', 'auto');
        
        if collect == 1
            if force_1D == 1
            disp('here0 sure')
            if reverse_order
                base_name = scan_dim(1).name;
                %base_name = scan_dim(2).name;
            else
                base_name = scan_dim(2).name;
                %base_name = scan_dim(1).name;
            end
            disp('here1')
            for i = 1:parts
                if reverse_order
                    legend_cell{i} = [base_name ' = ' num2str(ALLW_wave(instrument,wavelength_range).scan(1).value(i))];
                    %legend_cell{i} = [base_name ' = ' num2str(ALLW_wave(instrument,wavelength_range).scan(2).value(i))];
                else
                    legend_cell{i} = [base_name ' = ' num2str(ALLW_wave(instrument,wavelength_range).scan(2).value(i))];
                    %legend_cell{i} = [base_name ' = ' num2str(ALLW_wave(instrument,wavelength_range).scan(1).value(i))];
                end
            end
            if part == parts
            
            
            % temp code
            %new_legend_cell{1} = 'k=B';
            %for i = 1:parts
            %   new_legend_cell{i+1} = legend_cell{i}
            %end
            
            legend_handle = legend(legend_cell,'location','Best','interpreter','none','GridLineStyle','-','fontsize',14);
            end
            end
            title([title_string ' Integrated over ' num2str(Wmin(wavelength_range)) '-' num2str(Wmax(wavelength_range)) ' AA'] ,'interpreter','none','FontSize',FSL)
        else
            
            
            % NOT DONE
            plot_logical_local = true(1,num_instruments);
            if exist('plot_logical') > 0.5
                for ii = 1:num_instruments
                    if isfield(plot_logical,name{ii})
                        if plot_logical.(name{ii}) == 0 
                            plot_logical_local(ii) = 0;
                        end
                    end
                end
            end
            
            clear legelnd_list
            for ii = 1:num_instruments
                if plot_logical_local(ii) == 1
                    corres_index = sum(plot_logical_local(1:ii));
                    legend_list{corres_index} = name{ii};
                end
            end
            
            if exist('keyword')>0.5
            for ii = 1:length(name)
                if isfield(keyword,name{ii}) && plot_logical_local(ii) == 1
                    corres_index = sum(plot_logical_local(1:ii));
                    legend_list{corres_index} = [legend_list{corres_index} ' ' keyword.(name{ii})];
                end
            end    
            end
            
             % temp code
            %new_legend_cell{1} = 'k=B';
            %new_legend_cell{2} = 'k=1';
            %for i = 1:3
            %   new_legend_cell{i+2} = legend_list{i}
            %end
            
            if part == parts
            legend_handle = legend(legend_list,'location','Best','interpreter','none','fontsize',FS);
            %legend_handle = legend(legend_cell,'location','Best','interpreter','none','fontsize',22);
            end
            
            if reverse_order
                if size_select == 2
                title([title_string ' ' monitor(instrument,wavelength_range).scan(1).name ' = ' num2str(monitor(instrument,wavelength_range).scan(1).value(part)) ' Integrated over ' num2str(Wmin(wavelength_range)) '-' num2str(Wmax(wavelength_range)) ' Å'],'interpreter','none','FontSize',FSL)
                %title([monitor(instrument,wavelength_range).scan(1).name ' = ' num2str(monitor(instrument,wavelength_range).scan(1).value(part)) ' Wavelengths = ' num2str(Wmin(wavelength_range)) '-' num2str(Wmax(wavelength_range)) ' Å'],'interpreter','none','FontSize',FSL)
                %title(['Wavelengths = ' num2str(Wmin(wavelength_range)) '-' num2str(Wmax(wavelength_range)) ' Å'],'interpreter','none','FontSize',FSL)
                else
                title([title_string ' Wavelengths = ' num2str(Wmin(wavelength_range)) '-' num2str(Wmax(wavelength_range)) ' AA'],'interpreter','none','FontSize',FSL)    
                end
            else
                if size_select == 2
                title([title_string ' ' monitor(instrument,wavelength_range).scan(2).name ' = ' num2str(monitor(instrument,wavelength_range).scan(2).value(part)) ' Integrated over = ' num2str(Wmin(wavelength_range)) '-' num2str(Wmax(wavelength_range)) ' Å'],'interpreter','none','FontSize',FSL)
                %title([monitor(instrument,wavelength_range).scan(2).name ' = ' num2str(monitor(instrument,wavelength_range).scan(2).value(part)) ' Wavelengths = ' num2str(Wmin(wavelength_range)) '-' num2str(Wmax(wavelength_range)) ' Å'],'interpreter','none','FontSize',FSL)
                %title(['Wavelengths = ' num2str(Wmin(wavelength_range)) '-' num2str(Wmax(wavelength_range)) ' Å'],'interpreter','none','FontSize',FSL)
                else
                title([title_string ' Wavelengths = ' num2str(Wmin(wavelength_range)) '-' num2str(Wmax(wavelength_range)) ' AA'],'interpreter','none','FontSize',FSL)    
                end
            end
            
        end
        if reverse_order
            xlabel(monitor(instrument,wavelength_range).scan(2).name,'interpreter','none','FontSize',FS)
            used_xlabel = monitor(instrument,wavelength_range).scan(2).name;
        else
            xlabel([monitor(instrument,wavelength_range).scan(1).name],'interpreter','none','FontSize',FS)
            used_xlabel = monitor(instrument,wavelength_range).scan(1).name;
        end
        
        %xlabel('Moderator height [m]','FontSize',FS)
        ylabel(yaxis,'FontSize',FS)
        %set(gca,'fontsize',18)
        hold off
        
        if plots == plot_options.fig_choice(end)
            for instrument = 1:num_instruments
               % write a file for each instrument with all y_plot selected
               indexer=floor(wavelength_range*0.1);
               if collect == 1
               fid = fopen(['scan_plotter_' num2str(indexer) '_' num2str(wavelength_range) '_' name{instrument} '.txt'], 'w');    
               else
                   % name was legend list, I do not know why that was
                   % needed here. I will keep it, so that if I find a
                   % reason it can be changed back
               fid = fopen(['scan_plotter_' num2str(indexer) '_' num2str(wavelength_range) '_' name{instrument} '.txt'], 'w');
               %fid = fopen(['scan_plotter_' num2str(indexer) '_' num2str(wavelength_range) '_' legend_list{instrument} '.txt'], 'w');
               end
                
                % header lines
                fprintf(fid,[num2str(Wmin(wavelength_range)) '-' num2str(Wmax(wavelength_range)) 'AA\n'])
                
                header_line = used_xlabel;
                for ii = plot_options.fig_choice
                   switch ii
                       case 1
                          header_line = [header_line ' - homogen source'];
                       case 2
                          header_line = [header_line ' - homogen source (degraded)'];
                       case 3
                          header_line = [header_line ' - ESS source'];
                       case 4
                          header_line = [header_line ' - ESS source (degraded)'];
                       case 5
                          header_line = [header_line ' - brilliance t.'];
                       case 6
                          header_line = [header_line ' - brilliance t. (degraded)'];
                       case 7
                          header_line = [header_line ' - ESS source (div limited)'];
                       case 8
                          header_line = [header_line ' - ESS source (div limited, degraded)'];
                   end
                end
                fprintf(fid,[header_line '\n'])
                
                % write data
                for ii = 1:length(X_plot)
                   line_string = [ num2str(X_plot(ii),'%10.6f') ];
                   for jj = plot_options.fig_choice
                       if jj == 5 || jj == 6
                            if length(global_Y_plot{jj,instrument}) >= ii
                                line_string = [ line_string '\t' num2str(global_Y_plot{jj,instrument}(ii),'%6.5f')];
                            else
                                line_string = [ line_string '\t' num2str(0)];
                            end
                       else
                           if length(global_Y_plot{jj,instrument}) >= ii
                            line_string = [ line_string '\t' num2str(global_Y_plot{jj,instrument}(ii),'%10e')];
                           else
                            line_string = [ line_string '\t' num2str(0)];   
                           end
                       end
                   end
                   line_string = [ line_string '\n'];
                   fprintf(fid,line_string)
                end
                
                fclose(fid)
            end
        end
         
    end

        hold off
        
        %xlim([0 1.75])
        %set(gca,'Xtick',[0 0.25 0.5 0.75 1 1.25 1.5])
        if plots==5
        %ylim([0 1])
        end
        %print(figure_handle,'-dpng','-r150',['scan_plotter_' plot_names(plots) '.png'])
        
        W_min = ALLW_wave(1,wavelength_range).waveinfo(1);
        W_max = ALLW_wave(1,wavelength_range).waveinfo(2);
        %title(['Integrated over wavelength band:' num2str(W_min) '-' num2str(W_max) ' Å'],'fontsize',FSL)
        box on
        indexer=floor(wavelength_range*0.1);
        if exist('nickname') > 0.5 
            print(figure_handle,'-dpng','-r300',['scan_plotter_' plot_names{plots} '_' num2str(indexer) '_' num2str(wavelength_range) '_' nickname '.png'])
        else
            print(figure_handle,'-dpng','-r300',['scan_plotter_' plot_names{plots} '_' num2str(indexer) '_' num2str(wavelength_range) '.png'])
        end
        %print(figure_handle,'-depsc',['scan_plotter_' plot_names{plots} '_' num2str(indexer) '_' num2str(wavelength_range) '.eps'])
    
    end
    
elseif size_select == 0
    
    num_instruments = length(name);
    temp_size = size(ALLW_wave(instrument,1).signal);
    
    %if num_instruments == 1
    %    collect = 1;
    %    colormatrix = colormap(colormap_choice);
    %    colorlength = length(colormatrix);
    %    skip = floor(colorlength/parts);
    %    colorpoint = skip*(1:parts);
    %    colorpoint = [1 colorpoint];
%     else
%         collect = 0;
%         colormatrix = colormap(colormap_choice);
%         colorlength = length(colormatrix);
%         skip = floor(colorlength/num_instruments);
%         colorpoint = skip*(1:num_instruments);
%         colorpoint = [1 colorpoint];
%     end

    for plots = plot_options.fig_choice% 5:6
        figure_handle = figure;
        max_val_ever = 0;
        switch plots
            case 1
                monitor = ALLW_wave;
                title_string = 'ALLW';
                yaxis = 'Intensity on sample, ARB unit';
            case 2
                monitor = ALLW_d_wave;
                title_string = 'ALLW_d';
                yaxis = 'Intensity on sample, ARB unit';
            case 3
                monitor = ESS_wave;
                title_string = 'ESS';
                yaxis = 'Intensity on sample n/s';
            case 4
                monitor = ESS_d_wave;
                title_string = 'ESS_d';
                yaxis = 'Intensity on sample n/s';
            case 5
                monitor = ALLW_B_wave;
                title_string = 'ALLW_B';
                yaxis = 'Brilliance transfer';
            case 6
                monitor = ALLW_d_B_wave;
                title_string = 'ALLW_d_B';
                yaxis = 'Brilliance transfer';
            case 7
                monitor = ESS_lim_wave;
                title_string = 'ESS (div limit)';
                yaxis = 'Intensity on sample n/s (within div bounds)';
                %yaxis = 'Flux on sample n/s/cm^2';
            case 8
                monitor = ESS_d_lim_wave;
                title_string = 'ESS_d (div limit)';
                yaxis = 'Intensity on sample n/s (within div bounds)';
        end
                
        for instrument = 1:num_instruments
                
                %X_plot(instrument) = instrument;
                %Y_plot(instrument) = monitor(instrument,wavelength_range).signal;
                %Y_error_plot(instrument) = monitor(instrument,wavelength_range).error*0;
                Y_plot(instrument) = instrument;
                X_plot(instrument) = monitor(instrument,wavelength_range).signal;
                
                
                
                %handle = errorbar(monitor(instrument,1).scan(2).value(cut_zero),monitor(instrument,1).signal(part,cut_zero),monitor(instrument,1).error(part,cut_zero).*0,'-','color',colormatrix(colorpoint(part),:));    
                %handle = errorbar(X_plot,Y_plot,Y_error_plot,'-','color',color_input);
                        
        end
        
%             if cleanup_low_points
%             clear good_point
%             plot_length = length(X_plot);
%             good_point(1:plot_length) = true;
%             for index = 2:plot_length-1
%                 good_point(index) = Y_plot(index) > Y_plot(index-1)*tollerance && Y_plot(index) > Y_plot(index+1)*tollerance;
%             end
%             
%             else    
%                good_point = ones(size(X_plot));
%             end
            
%            X_plot = X_plot(good_point);
%            Y_plot = Y_plot(good_point);
%            Y_error_plot = Y_error_plot(good_point);
            
            %X_plot_save{part} = X_plot;
            %Y_plot_save{part} = Y_plot;
            
            %handle = errorbar(X_plot(good_point),Y_plot(good_point),Y_error_plot(good_point),'-','color',color_input);
            %handle = errorbar(X_plot,Y_plot,Y_error_plot,'k+');
            
            handle = plot(X_plot,Y_plot,'k+');
            %handle = errorbar(X_plot{part},Y_plot{part},Y_error_plot{part},'-','color',color_input);
            set(handle,'Marker','+')
            set(handle,'MarkerSize',10)
            %set(handle,'LineWidth',2);
            %set(gca,'FontSize',13)
            
            
            
            %set(gca,'XTick',X_plot)
            %set(gca,'XTickLabel',name)
            
            keyword_name_list = name;
            if exist('keyword')>0.5
            for ii = 1:length(name)
                if isfield(keyword,name{ii})
                    keyword_name_list{ii} = [keyword_name_list{ii} ' ' keyword.(name{ii})];
                end
            end    
            end
            
            set(gca,'YTick',Y_plot)
            set(gca,'YTickLabel',keyword_name_list)
            
            %max_val_this = max(Y_plot);
            %if max_val_this > max_val_ever
            %    max_val_ever = max_val_this;
            %end
            %ylim([0 max_val_ever*1.05]);

        
            max_val_this = max(X_plot);
            if max_val_this > max_val_ever
                max_val_ever = max_val_this;
            end
            xlim([0 max_val_ever*1.05]);
            ylim([0.5 length(name)+0.5]);

                
%         
%         hold off
%     
         set(figure_handle, 'Position', [0 0 900 600])
         set(figure_handle, 'paperpositionmode', 'auto');
%         
%         if collect == 1
%             if reverse_order
%                 base_name = scan_dim(1).name;
%             else
%                 base_name = scan_dim(2).name;
%             end
%             for i = 1:parts
%                 if reverse_order
%                     legend_cell{i} = [base_name ' = ' num2str(ALLW_wave(instrument,wavelength_range).scan(1).value(i))];
%                 else
%                     legend_cell{i} = [base_name ' = ' num2str(ALLW_wave(instrument,wavelength_range).scan(2).value(i))];
%                 end
%             end
%             if part == parts
%             legend_handle = legend(legend_cell,'location','Best','interpreter','none','GridLineStyle','-');
%             end
%             title([title_string ' Wavelengths = ' num2str(Wmin(wavelength_range)) '-' num2str(Wmax(wavelength_range)) ' Å'] ,'interpreter','none','FontSize',FSL)
%         else
%             if part == parts
%             legend_handle = legend(name,'location','Best','interpreter','none');
%             end
%             if reverse_order
%                 %title([title_string ' ' monitor(instrument,wavelength_range).scan(1).name ' = ' num2str(monitor(instrument,wavelength_range).scan(1).value(part)) ' Wavelengths = ' num2str(Wmin(wavelength_range)) '-' num2str(Wmax(wavelength_range)) ' Å'],'interpreter','none','FontSize',FSL)
%                 %title([monitor(instrument,wavelength_range).scan(1).name ' = ' num2str(monitor(instrument,wavelength_range).scan(1).value(part)) ' Wavelengths = ' num2str(Wmin(wavelength_range)) '-' num2str(Wmax(wavelength_range)) ' Å'],'interpreter','none','FontSize',FSL)
%                 title(['Wavelengths = ' num2str(Wmin(wavelength_range)) '-' num2str(Wmax(wavelength_range)) ' Å'],'interpreter','none','FontSize',FSL)
%             else
%                 %title([title_string ' ' monitor(instrument,wavelength_range).scan(2).name ' = ' num2str(monitor(instrument,wavelength_range).scan(2).value(part)) ' Wavelengths = ' num2str(Wmin(wavelength_range)) '-' num2str(Wmax(wavelength_range)) ' Å'],'interpreter','none','FontSize',FSL)
%                 %title([monitor(instrument,wavelength_range).scan(2).name ' = ' num2str(monitor(instrument,wavelength_range).scan(2).value(part)) ' Wavelengths = ' num2str(Wmin(wavelength_range)) '-' num2str(Wmax(wavelength_range)) ' Å'],'interpreter','none','FontSize',FSL)
%                 title(['Wavelengths = ' num2str(Wmin(wavelength_range)) '-' num2str(Wmax(wavelength_range)) ' Å'],'interpreter','none','FontSize',FSL)
%             end
%             
%         end
%         if reverse_order
%             xlabel(monitor(instrument,wavelength_range).scan(2).name,'interpreter','none','FontSize',FS)
%         else
%             xlabel(monitor(instrument,wavelength_range).scan(1).name,'interpreter','none','FontSize',FS)
%         end
%         
%         xlabel('Moderator height [m]','FontSize',FS)
%         ylabel(yaxis,'FontSize',FS)
%         hold off
        
        W_min = ALLW_wave(1,wavelength_range).waveinfo(1);
        W_max = ALLW_wave(1,wavelength_range).waveinfo(2);
        title(['Integrated over wavelength band:' num2str(W_min) '-' num2str(W_max) ' AA'],'fontsize',FSL)
        xlabel(yaxis,'FontSize',FS)
        
        %print(figure_handle,'-dpng','-r150',['scan_plotter_' plot_names(plots) '.png'])
        box on
        indexer=floor(wavelength_range);
        print(figure_handle,'-dpng','-r300',['scan_plotter_' plot_names{plots} '_' num2str(indexer) '_' num2str(wavelength_range) '.png'])
        %print(figure_handle,'-depsc',['scan_plotter_' plot_names{plots} '_' num2str(indexer) '_' num2str(wavelength_range) '.eps'])
        
    
    end
    
elseif size_select == 2    
    
    disp('Sadly the 2D histogram have not been written yet, use force 1D and reverse_order instead!')
    something = 1;
    
     for plots = plot_options.fig_choice% 5:6
        figure_handle = figure;
        max_val_ever = 0;
        switch plots
            case 1
                monitor = ALLW_wave;
                title_string = 'ALLW ARB';
                yaxis = 'Intensity on sample, ARB unit';
            case 2
                monitor = ALLW_d_wave;
                title_string = 'ALLW_d ARB';
                yaxis = 'Intensity on sample, ARB unit';
            case 3
                monitor = ESS_wave;
                title_string = 'ESS flux';
                yaxis = 'Intensity on sample n/s';
                %yaxis = 'Flux on sample n/s/cm^2';
            case 4
                monitor = ESS_d_wave;
                title_string = 'ESS_d flux';
                yaxis = 'Intensity on sample n/s';
            case 5
                monitor = ALLW_B_wave;
                title_string = 'Brilliance transfer.';
                yaxis = 'Brilliance transfer';
            case 6
                monitor = ALLW_d_B_wave;
                title_string = 'Brilliance d transfer';
                yaxis = 'Brilliance transfer';
            case 7
                monitor = ESS_lim_wave;
                title_string = 'ESS (div limit)';
                yaxis = 'Intensity on sample n/s (within div bounds)';
            case 8
                monitor = ESS_d_lim_wave;
                title_string = 'ESS_d (div limit)';
                yaxis = 'Intensity on sample n/s (within div bounds)';
        end
        
        num_instruments = length(name);
        fig_handle = figure
        subplot(num_instruments,1,1)
        set(fig_handle, 'Position', [0 0 600 700])
        %set(fig_handle, 'paperpositionmode', 'auto');
        for instr = 1:num_instruments
            switch instr
                case 1
                    subplot(num_instruments,1,2)
                case 2
                    subplot(num_instruments,1,3)
                case 3
                    subplot(num_instruments,1,1)
            end
             
%             if plots < 3.5
%             % flux special case
%             T(1,:) = monitor(instr,wavelength_range).signal(1,:)./((monitor(instr,wavelength_range).scan(2).value*100).^2)
%             T(2,:) = monitor(instr,wavelength_range).signal(2,:)./((monitor(instr,wavelength_range).scan(2).value*100).^2)
%             T(3,:) = monitor(instr,wavelength_range).signal(3,:)./((monitor(instr,wavelength_range).scan(2).value*100).^2)
%             T(4,:) = monitor(instr,wavelength_range).signal(4,:)./((monitor(instr,wavelength_range).scan(2).value*100).^2)
%             T(5,:) = monitor(instr,wavelength_range).signal(5,:)./((monitor(instr,wavelength_range).scan(2).value*100).^2)
%             else
%                 T = monitor(instr,wavelength_range).signal;
%             end
            
            imagesc(monitor(instr,wavelength_range).scan(2).value,monitor(instr,wavelength_range).scan(1).value,monitor(instr,wavelength_range).signal)
            %imagesc(monitor(instr,wavelength_range).scan(2).value,monitor(instr,wavelength_range).scan(1).value,T)
            set(gca,'Ydir','normal')
            
            max_val = max(max(T));
            min_val = min(min(T));
            
            set(gca,'CLim',[0 max_val])
            
            set(gca,'Xtick',monitor(instr,wavelength_range).scan(2).value)
            set(gca,'Ytick',monitor(instr,wavelength_range).scan(1).value)
            xlabel(monitor(instr,wavelength_range).scan(2).name,'fontsize',FS)
            ylabel(monitor(instr,wavelength_range).scan(1).name,'fontsize',FS)
            title([name{instr} ' ' title_string ' over ' num2str(Wmin(wavelength_range)) '-' num2str(Wmax(wavelength_range)) ' AA'] ,'interpreter','none','FontSize',FSL)
            colormap jet
            colorbar
        end
     end
    
    
    
        
else
    disp('ERROR, corrpupted data?')
end

%if plots == 5
%    hold on
%plot([0 1],[0 1],'--','linewidth',2,'color',[0.45 0.45 0.45])
%plot([1 1],[0 1],'-.','linewidth',2,'color',[0.45 0.45 0.45])
%print(figure_handle,'-dpng','-r300',['scan_plotter_' plot_names{plots} '_' num2str(indexer) '_' num2str(wavelength_range) '.png'])
%hold off
%end


end

%%
if 1==2
load('thermal_heimdal_source_flux.mat')
%correction = 6e9/flux_in_band(10);
hold on
%handle = plot(index2_thermal.*0.01,flux_in_band.*correction,'k','linewidth',2)
correction = 7e9/flux_in_band(10);
handle = plot(index2_thermal.*0.01,flux_in_band.*correction,'k','linewidth',2)
correction = 5.8e9/flux_in_band(10);
handle = plot(index2_thermal.*0.01,flux_in_band.*correction,'k','linewidth',2)
correction = 7.6e9/flux_in_band(10);
handle = plot(index2_thermal.*0.01,flux_in_band.*correction,'k','linewidth',2)
end

%%
if 1==2
load('cold_camea_source_flux.mat')
%correction = 6e9/flux_in_band(10);
hold on
%handle = plot(index2_thermal.*0.01,flux_in_band.*correction,'k','linewidth',2)
correction = 8.4e10/flux_in_band(9);
handle = plot(index2.*0.01,flux_in_band.*correction,'--b','linewidth',2)
correction = 9.8e10/flux_in_band(9);
handle = plot(index2.*0.01,flux_in_band.*correction,'--b','linewidth',2)
end
































