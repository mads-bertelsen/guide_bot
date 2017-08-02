function IntermediateBrilliance_Analysis(Hdiv,Vdiv, Hsize, Vsize, div1d, div2d, psd1d, psd2d, apsd, adiv, matfile, brillfile, savename)
% This function is equivelent to the analysis script that Mads tailor
% generates for each optimization within a project. All values hardwired
% into the tailored scripts have been recast as inputs in this function, thereby
% generalizing this single analysis script to support any optimization run.

% This function is written to be run inconjuction with the function
% IntermediateBrilliance_ifit. They should of course both use the same brilliance
% window. Taken together, they allow any <name>_analyze.inst and <name>_analyze_ess.inst
% file pair written with the Mads detectors ensemble to be analyzed and
% plotted with guide_bot analysis.


%%%% INPUTS %%%%

% Brilliance Window that analysis is to be performed over:
% Hdiv:       Horizontal Divergence
% Vdiv:       Vertical Divergence
% Hsize:      Horizontal width
% Vsize:      Vertical height

% Scale factors for plotting:
% div1d:      1d plots of divergence
% div2d:      2d plots of divergence (Vdiv x Hdiv)
% psd1d:      1d plots of beam width/height
% psd2d:      2d plots of area (Vsize x Hsize)
% apsd/adiv:  2d plots of (Hdiv x Hsize) or (Vdiv x Vsize)

% matfile:    .mat file that contains all data to be analyzed
% brillfile:  .mat file that contains brilliance reference data, just set
% to 'none' if there is none.

% savename:   name prefix for saving plots, if set to 'default' then it
%             will plot using the name in the .mat matfile

% Any input (exept matfile) can be set to 'default, in which case the input
% will be set to the .mat matfile value
%%%%%%%%%%%%%%%%


% Note, may need to add detailed_absolute_run to the inputs, waiting for
% Mads to select a solution to a current mcstas_bot bug. This input will
% require an additional line, which I have added but commented out below
load(matfile);


if ~strcmp(savename,'default')
    filename = savename;
end

cpath=pwd;

MaxIndex=length(wavecenters)+1;
Foverall=figure
set(Foverall, 'Position', [0 0 700 1000])
subplot(4,2,1)
if select==1; set(gcf, 'Renderer', 'painters'); end;
Fposdiv=figure
set(Fposdiv, 'Position', [0 0 700 1000])
subplot(MaxIndex,2,1)
if select==1; set(gcf, 'Renderer', 'painters'); end;
Facceptance=figure
set(Facceptance, 'Position', [0 0 700 1000])
subplot(MaxIndex,2,1)
if select==1; set(gcf, 'Renderer', 'painters'); end;
background=figure
set(background, 'Position', [0 0 700 1000])
subplot(4,2,1)
if select==1; set(gcf, 'Renderer', 'painters'); end;

colors={'r' 'g' 'b' 'k' 'm'};
linecolor='k';
linethick=2;
FS=14;
FS_medium=12;
FS_small=10;

reflogic=exist(['../brill_ref/' brillfile ])>0.5;

if reflogic
load(['../brill_ref/' brillfile]);
end

for i=1:MaxIndex
   if i==MaxIndex
    monitors=monitor_fom;
    if reflogic; monitors_ref=monitor_fom_ref; end;
   else
    Wnames=fieldnames(monitor_W);
    monitors=monitor_W.(Wnames{i});
    if reflogic; monitors_ref=monitor_W_ref.(Wnames{i}); end;
   end
   
   if reflogic
   guide_end_lambda = assign_by_title('Lmon_guide_end.',monitors);
   DIV2D_ref=assign_by_title('Div2d_sample.',monitors_ref)
   DIV2D=assign_by_title('Div2d_sample.',monitors)/DIV2D_ref.Data.Mean;
   PSD2D_ref=assign_by_title('PSD_sample.',monitors_ref);
   PSD2D=assign_by_title('PSD_sample.',monitors)/PSD2D_ref.Data.Mean;
   HPSD_ref=assign_by_title('HPSD_sample.',monitors_ref);
   HPSD=assign_by_title('HPSD_sample.',monitors)/HPSD_ref.Data.Mean;
   HPSD{1} = HPSD{1}.*100; %Conversion to cm
   VPSD_ref=assign_by_title('VPSD_sample.',monitors_ref);
   VPSD=assign_by_title('VPSD_sample.',monitors)/VPSD_ref.Data.Mean;
   VPSD{1} = VPSD{1}.*100; %Conversion to cm
   HDIV_ref=assign_by_title('Hdiv_sample.',monitors_ref);
   HDIV=assign_by_title('Hdiv_sample.',monitors)/HDIV_ref.Data.Mean;
   VDIV_ref=assign_by_title('Vdiv_sample.',monitors_ref);
   VDIV=assign_by_title('Vdiv_sample.',monitors)/VDIV_ref.Data.Mean;
   HACCP_ref=assign_by_title('acceptance_x_divx.',monitors_ref);
   HACCP=assign_by_title('acceptance_x_divx.',monitors)/HACCP_ref.Data.Mean;
   HACCP{2} = HACCP{2}.*100; %Conversion to cm
   VACCP_ref=assign_by_title('acceptance_y_divy.',monitors_ref);
   VACCP=assign_by_title('acceptance_y_divy.',monitors)/HACCP_ref.Data.Mean;
   VACCP{2} = VACCP{2}.*100; %Conversion to cm
   LAMBDA_B_RAW_monitor=assign_by_title('Lmon_sample_B.',monitors);
   LAMBDA_B_RAW=LAMBDA_B_RAW_monitor.Data.values(1);
   LAMBDA_B_ref=assign_by_title('Lmon_sample_B.',monitors_ref);
   LAMBDA_B=assign_by_title('Lmon_sample_B.',monitors)/LAMBDA_B_ref.Data.Mean;
   DIV2D_B_ref=assign_by_title('Div2d_sample_maxdiv.',monitors_ref);
   DIV2D_B=assign_by_title('Div2d_sample_maxdiv.',monitors)/DIV2D_B_ref.Data.Mean;
   PSD2D_B_ref=assign_by_title('PSD_sample_maxdiv.',monitors_ref);
   PSD2D_B=assign_by_title('PSD_sample_maxdiv.',monitors)/PSD2D_B_ref.Data.Mean;
   HPSD_B_ref=assign_by_title('HPSD_sample_maxdiv.',monitors_ref);
   HPSD_B=assign_by_title('HPSD_sample_maxdiv.',monitors)/HPSD_B_ref.Data.Mean;
   HPSD_B{1} = HPSD_B{1}.*100; %Conversion to cm
   VPSD_B_ref=assign_by_title('VPSD_sample_maxdiv.',monitors_ref);
   VPSD_B=assign_by_title('VPSD_sample_maxdiv.',monitors)/VPSD_B_ref.Data.Mean;
   VPSD_B{1} = VPSD_B{1}.*100; %Conversion to cm
   HDIV_B_ref=assign_by_title('Hdiv_sample_maxdiv.',monitors_ref);
   HDIV_B=assign_by_title('Hdiv_sample_maxdiv.',monitors)/HDIV_B_ref.Data.Mean;
   VDIV_B_ref=assign_by_title('Vdiv_sample_maxdiv.',monitors_ref);
   VDIV_B=assign_by_title('Vdiv_sample_maxdiv.',monitors)/VDIV_B_ref.Data.Mean;
   HACCP_B_ref=assign_by_title('acceptance_x_divx_maxdiv.',monitors_ref);
   HACCP_B=assign_by_title('acceptance_x_divx_maxdiv.',monitors)/HACCP_B_ref.Data.Mean;
   HACCP_B{2} = HACCP_B{2}.*100; %Conversion to cm
   VACCP_B_ref=assign_by_title('acceptance_y_divy_maxdiv.',monitors_ref);
   VACCP_B=assign_by_title('acceptance_y_divy_maxdiv.',monitors)/VACCP_B_ref.Data.Mean;
   VACCP_B{2} = VACCP_B{2}.*100; %Conversion to cm
   LAMBDA_ref=assign_by_title('Lmon_sample.',monitors_ref);
   LAMBDA=assign_by_title('Lmon_sample.',monitors)/LAMBDA_ref.Data.Mean;
   LAMBDA_RAW_monitor=assign_by_title('Lmon_sample.',monitors);
   AROUND_SAMPLE_BT_ref=assign_by_title('Lmon_sample.',monitors_ref);
   AROUND_SAMPLE_BT_val=assign_by_title('Lmon_sample.',monitors)/AROUND_SAMPLE_BT_ref.Data.Mean;
   AROUND_SAMPLE_BT=AROUND_SAMPLE_BT_val.Data.values(1)/AROUND_SAMPLE_BT_ref.Data.values(1);
   else
   guide_end_lambda = assign_by_title('Lmon_guide_end.',monitors);
   DIV2D=assign_by_title('Div2d_sample.',monitors);
   PSD2D=assign_by_title('PSD_sample.',monitors);
   HPSD=assign_by_title('HPSD_sample.',monitors);
   HPSD{1} = HPSD{1}.*100; %Conversion to cm
   VPSD=assign_by_title('VPSD_sample.',monitors);
   VPSD{1} = VPSD{1}.*100; %Conversion to cm
   HDIV=assign_by_title('Hdiv_sample.',monitors);
   VDIV=assign_by_title('Vdiv_sample.',monitors);
   HACCP=assign_by_title('acceptance_x_divx.',monitors);
   HACCP{2} = HACCP{2}.*100; %Conversion to cm
   VACCP=assign_by_title('acceptance_y_divy.',monitors);
   VACCP{2} = VACCP{2}*100; %Conversion to cm
   LAMBDA_B=assign_by_title('Lmon_sample_B.',monitors);
   LAMBDA_B_RAW_monitor=assign_by_title('Lmon_sample_B.',monitors);
   LAMBDA_B_RAW=LAMBDA_B_RAW_monitor.Data.values(1);
   DIV2D_B=assign_by_title('Div2d_sample_maxdiv.',monitors);
   PSD2D_B=assign_by_title('PSD_sample_maxdiv.',monitors);
   HPSD_B=assign_by_title('HPSD_sample_maxdiv.',monitors);
   HPSD_B{1} = HPSD_B{1}*100; %Conversion to cm
   VPSD_B=assign_by_title('VPSD_sample_maxdiv.',monitors);
   VPSD_B{1} = VPSD_B{1}*100; %Conversion to cm
   HDIV_B=assign_by_title('Hdiv_sample_maxdiv.',monitors);
   VDIV_B=assign_by_title('Vdiv_sample_maxdiv.',monitors);
   HACCP_B=assign_by_title('acceptance_x_divx_maxdiv.',monitors);
   HACCP_B{2} = HACCP_B{2}*100; %Conversion to cm
   VACCP_B=assign_by_title('acceptance_y_divy_maxdiv.',monitors);
   VACCP_B{2} = VACCP_B{2}*100; %Conversion to cm
   LAMBDA=assign_by_title('Lmon_sample.',monitors);
   LAMBDA_RAW=assign_by_title('Lmon_sample.',monitors);
   end
   
   if i==MaxIndex
      figure(Foverall);
      subplot(4,2,7:8)
      axis off;
      text(0.5,1,Project_name,'interpreter','none')
      text(0,0.83,[instrument_name ' - ' filename ' - '  inputstring],'interpreter','none')
      text(0,0.66,['Sample size: Horizontal=' num2str(Hsize) 'cm, Vertical=' num2str(Vsize) 'cm'],'interpreter','none');
      text(0,0.5,['Divergence requirement: Horizontal=' num2str(Hdiv) 'deg, Vertical= ' num2str(Vdiv) 'deg'],'interpreter','none');
      text(0,0.33,['Intensity on sample of 100 emitted = ' num2str(LAMBDA.Data.values(1)) ' (no divergence limit) ' num2str(LAMBDA_B_RAW) ' (width divergence limits)' ],'interpreter','none');
      text(0,0.16,['Intensity at guide end of 100 emitted = ' num2str(guide_end_lambda.Data.values(1)) ' (no divergence limit)' ],'interpreter','none');
   if reflogic
      text(0,0,['BT on sample = ' num2str(LAMBDA_B.Data.values(1)) ' (width divergence limits)' ],'interpreter','none');
      text(0,-0.16,['BT near sample = ' num2str(AROUND_SAMPLE_BT) ' (width divergence limits)' ],'interpreter','none');
   end
   
      figure(Foverall)
      subplot(4,2,1:2) 
      
      %plot(LAMBDA_B,'k')
    if reflogic
      axis([0 MaxWB 0 1]) % may need to comment out for monochromator
    else
      axis([0 MaxWB 0 LAMBDA_B.Data.Max*1.1])
    end
      %set(gca,'XTick',[0 0.5 1.0 1.5 2.0 2.5 3.0])
      title(['Wavelength dependence, + are wavelengths for 1d graphs. I=' num2str(LAMBDA_B_RAW)])
      ylabel('Brilliance transfer')
      set(gca,'Fontsize',FS)
    if reflogic
      markerheight=0.1;
    else
      markerheight=LAMBDA_B.Data.Max*0.2;
    end
        for j=1:length(wavecenters)
                hold on
                plot(wavecenters(j),markerheight,['+' colors{j}],'MarkerSize',12)
                hold off
        end
   else
   
    figure(Foverall)
    subplot(4,2,3)
    hold on
    plot(HDIV_B,colors{i})
    box
    if reflogic
    plot([-Hdiv -Hdiv],[0 1],'--k');
    plot([Hdiv Hdiv],[0 1],'--k');
    axis([-div1d*Hdiv div1d*Hdiv 0 1]) % may need to comment out for monochromator
    else
    if i==1; max_HDIV_B = HDIV_B.Data.Max;
    else max_HDIV_B = max([max_HDIV_B HDIV_B.Data.Max]); end;
    if i==MaxIndex-1;
      plot([-Hdiv -Hdiv],[0 max_HDIV_B*1.1],'--k');
      plot([Hdiv Hdiv],[0 max_HDIV_B*1.1],'--k');
      axis([-div1d*Hdiv div1d*Hdiv 0 max_HDIV_B*1.1]);
     end
    end
    %set(gca,'XTick',[-1 -0.5 -0.3 0 0.3 0.5 1])
    title('')
    xlabel('Horizontal divergence [deg]')
    if reflogic
      ylabel('Brilliance transfer')
    else
      ylabel('Intensity in FOM [arb]')
    end
    set(gca,'Fontsize',FS)

    subplot(4,2,4)
    hold on
    plot(VDIV_B,colors{i})
    box
    if reflogic
    plot([-Vdiv -Vdiv],[0 1],'--k');
    plot([Vdiv Vdiv],[0 1],'--k');
    axis([-div1d*Vdiv div1d*Vdiv 0 1]) % may need to comment out for monochromator
    else
    if i==1; max_VDIV_B = VDIV_B.Data.Max;
    else max_VDIV_B = max([max_VDIV_B VDIV_B.Data.Max]); end;
    if i==MaxIndex-1;
      plot([-Vdiv -Vdiv],[0 max_VDIV_B*1.1],'--k');
      plot([Vdiv Vdiv],[0 max_VDIV_B*1.1],'--k');
      axis([-div1d*Vdiv div1d*Vdiv 0 max_VDIV_B*1.1]);
     end
    end
    %set(gca,'XTick',[-1 -0.5 -0.3 0 0.3 0.5 1])
    title('')
    xlabel('Vertical divergence [deg]')
    if reflogic
      ylabel('Brilliance transfer')
    else
      ylabel('Intensity in FOM [arb]')
    end
    set(gca,'Fontsize',FS)

    subplot(4,2,5)
    hold on
    plot(HPSD_B,colors{i})
    box
    %axis([-0.75 0.75 0 1])
    if reflogic
    plot([-0.5*Hsize -0.5*Hsize],[0 1],'--k');
    plot([0.5*Hsize 0.5*Hsize],[0 1],'--k');
    axis([-0.5*psd1d*Hsize 0.5*psd1d*Hsize 0 1]) % may need to comment out for monochromator
    else
    if i==1; max_HPSD_B= HPSD_B.Data.Max;
    else max_HPSD_B= max([max_HPSD_B HPSD_B.Data.Max]); end;
    if i==MaxIndex-1;
      plot([-0.5*Hsize -0.5*Hsize],[0 max_HPSD_B*1.1],'--k');
      plot([0.5*Hsize 0.5*Hsize],[0 max_HPSD_B*1.1],'--k');
      axis([-0.5*psd1d*Hsize 0.5*psd1d*Hsize 0 max_HPSD_B*1.1]); 
     end
    end
    %set(gca,'XTick',[-0.75 -0.5 -0.25 0 0.25 0.5 0.75])
    title('')
    xlabel('Horizontal position [cm]')
    if reflogic
      ylabel('Brilliance transfer')
    else
      ylabel('Intensity in FOM [arb]')
    end
    set(gca,'Fontsize',FS)

    subplot(4,2,6)
    hold on
    plot(VPSD_B,colors{i})
    box
    %axis([-1 1 0 1])
    if reflogic
    plot([-0.5*Vsize -0.5*Vsize],[0 1],'--k');
    plot([0.5*Vsize 0.5*Vsize],[0 1],'--k');
    axis([-0.5*psd1d*Vsize 0.5*psd1d*Vsize 0 1]) % may need to comment out for monochromator
    else
    if i==1; max_VPSD_B= VPSD_B.Data.Max;
    else max_VPSD_B= max([max_VPSD_B VPSD_B.Data.Max]); end;
    if i==MaxIndex-1;
      plot([-0.5*Vsize -0.5*Vsize],[0 max_VPSD_B*1.1],'--k');
      plot([0.5*Vsize 0.5*Vsize],[0 max_VPSD_B*1.1],'--k');
      axis([-0.5*psd1d*Vsize 0.5*psd1d*Vsize 0 max_VPSD_B*1.1]);
     end
    end
    %set(gca,'XTick',[-1 -0.5 -0.25 0 0.25 0.5 1])
    title('')
    xlabel('Vertical position [cm]')
    if reflogic
      ylabel('Brilliance transfer')
    else
      ylabel('Intensity in FOM [arb]')
    end
    set(gca,'Fontsize',FS)
   
   end
   
   
    figure(Fposdiv)
    subplot(MaxIndex,2,(i-1)*2+1)
    plot(DIV2D_B)
    box
    view([-0.5 90]);
    %axis([-1 1 -1 1])
   axis([-div2d*Hdiv div2d*Hdiv -div2d*Vdiv div2d*Vdiv])
    maxi=max(DIV2D_B)+100;
    hold on
    x=Hdiv;
    y=Vdiv;
    plot3([-x x], [-y -y], [maxi maxi],'color',linecolor,'LineWidth',linethick)
    plot3([-x x], [y y], [maxi maxi],'color',linecolor,'LineWidth',linethick)
    plot3([-x -x], [-y y], [maxi maxi],'color',linecolor,'LineWidth',linethick)
    plot3([x x], [-y y], [maxi maxi],'color',linecolor,'LineWidth',linethick)
    hold off
    %set(gca,'XTick',[-1 -0.5 -0.3 0 0.3 0.5 1])
    %set(gca,'YTick',[-1 -0.5 -0.3 0 0.3 0.5 1])
    xlabel('Horizontal divergence [deg]');
    set(gca,'Fontsize',FS_medium)
    ylabel('Vertical divergence [deg]','fontsize',FS_medium);
    if i==MaxIndex
        w_min = DIV2D_B.Data.WaveMin;w_max=DIV2D_B.Data.WaveMax;
        title(['2d div ' num2str(w_min) '-' num2str(w_max) ' Å'],'fontsize',FS_small)
    else
        title(['2d div Lambda=' num2str(wavecenters(i)) 'A'],'fontsize',FS_small)
    end


    subplot(MaxIndex,2,(i-1)*2+2)
    plot(PSD2D_B)
    box
    view([-0.5 90]);
    %axis([-1 1 -1 1])
   axis([-0.5*psd2d*Hsize 0.5*psd2d*Hsize -0.5*psd2d*Vsize 0.5*psd2d*Vsize])
    maxi=max(PSD2D_B)+100;
    hold on
    x=Hsize*0.5;
    y=Vsize*0.5;
    plot3([-x x], [-y -y], [maxi maxi],'color',linecolor,'LineWidth',linethick)
    plot3([-x x], [y y], [maxi maxi],'color',linecolor,'LineWidth',linethick)
    plot3([-x -x], [-y y], [maxi maxi],'color',linecolor,'LineWidth',linethick)
    plot3([x x], [-y y], [maxi maxi],'color',linecolor,'LineWidth',linethick)
    hold off
    %set(gca,'XTick',[-1 -0.5 -0.3 0 0.3 0.5 1])
    %set(gca,'YTick',[-1 -0.5 -0.3 0 0.3 0.5 1])
    xlabel('Horizontal position [cm]');
    ylabel('Vertical position [cm]','fontsize',FS_medium);
    set(gca,'Fontsize',FS_medium)
    if i==MaxIndex
        w_min = PSD2D_B.Data.WaveMin;w_max=PSD2D_B.Data.WaveMax;
        title(['2d psd ' num2str(w_min) '-' num2str(w_max) ' Å'],'fontsize',FS_small)
    else
        title(['2d psd Lambda=' num2str(wavecenters(i)) 'A'],'fontsize',FS_small)
    end
   
    
    figure(Facceptance)
    subplot(MaxIndex,2,(i-1)*2+1)
    box
    plot(HACCP_B)
    view([-0.5 90]);
    %axis([-1 1 -1 1])
   axis([-0.5*apsd*Hsize 0.5*apsd*Hsize -adiv*Hdiv adiv*Hdiv])
    maxi=max(HACCP_B)+100;
    hold on
    x=Hsize*0.5;
    y=Hdiv;
    plot3([-x x], [-y -y], [maxi maxi],'color',linecolor,'LineWidth',linethick)
    plot3([-x x], [y y], [maxi maxi],'color',linecolor,'LineWidth',linethick)
    plot3([-x -x], [-y y], [maxi maxi],'color',linecolor,'LineWidth',linethick)
    plot3([x x], [-y y], [maxi maxi],'color',linecolor,'LineWidth',linethick)
    hold off
    %set(gca,'XTick',[-1 -0.5 -0.3 0 0.3 0.5 1])
    %set(gca,'YTick',[-1 -0.5 -0.3 0 0.3 0.5 1])
    xlabel('Horizontal position [cm]')
    set(gca,'Fontsize',FS_medium)
    ylabel('Hor. divergence [deg]','fontsize',FS_medium)
    if i==MaxIndex
    w_min = HACCP_B.Data.WaveMin;w_max=HACCP_B.Data.WaveMax;
    title(['Horizontal accep. dia. ' num2str(w_min) '-' num2str(w_max) ' Å'],'fontsize',FS_small)
    else
    title(['Horizontal accep. dia. Wavelength=' num2str(wavecenters(i)) 'A'],'fontsize',FS_small);
    end    

    subplot(MaxIndex,2,(i-1)*2+2)
    plot(VACCP_B)
    box
    view([-0.5 90]);
    %axis([-1 1 -1 1])
   axis([-0.5*apsd*Vsize 0.5*apsd*Vsize -adiv*Vdiv adiv*Vdiv])
    maxi=max(VACCP_B)+100;
    hold on
    x=Vsize*0.5;
    y=Vdiv;
    plot3([-x x], [-y -y], [maxi maxi],'color',linecolor,'LineWidth',linethick)
    plot3([-x x], [y y], [maxi maxi],'color',linecolor,'LineWidth',linethick)
    plot3([-x -x], [-y y], [maxi maxi],'color',linecolor,'LineWidth',linethick)
    plot3([x x], [-y y], [maxi maxi],'color',linecolor,'LineWidth',linethick)
    hold off
    %set(gca,'XTick',[-1 -0.5 -0.3 0 0.3 0.5 1])
    %set(gca,'YTick',[-1 -0.5 -0.3 0 0.3 0.5 1])
    xlabel('Vertical position [cm]')
    set(gca,'Fontsize',FS_medium)
    ylabel('Ver. divergence [deg]','fontsize',FS_medium)
    if i==MaxIndex
    w_min = VACCP_B.Data.WaveMin;w_max=VACCP_B.Data.WaveMax;
    title(['Vertical accep. dia. ' num2str(w_min) '-' num2str(w_max) ' Å'],'fontsize',FS_small)
    else
    title(['Vertical accep. dia. Wavelength=' num2str(wavecenters(i)) 'A'],'fontsize',FS_small);
    end
    
end

if reflogic
flux=monitor_fom(11).Data.values(1)/monitor_fom_ref(11).Data.values(1);
fluxtext=' BT ';
else
flux=monitor_fom(11).Data.values(1);
fluxtext=' Absolute ';
end

fid = fopen([cpath '/master_record-analyzed' scanname '.txt'],'a');
fprintf(fid,[num2str(flux) fluxtext '= ' filename ' - ' inputstring '\n'])
fclose(fid);

set(Foverall, 'paperpositionmode', 'auto');
set(Fposdiv, 'paperpositionmode', 'auto');
set(Facceptance, 'paperpositionmode', 'auto');

if reflogic
hold on
figure(Foverall)
subplot(4,2,1:2)
hold on
LAMBDA_B_ALLW_deg=assign_by_title('Lmon_sample_B.',monitor_ALLW_degraded);
LAMBDA_B_ALLW_deg_ref=assign_by_title('Lmon_sample_B.',monitor_ALLW_ref);
LAMBDA_B_ALLW=assign_by_title('Lmon_sample_B.',monitor_ALLW);
LAMBDA_B_ALLW_ref=assign_by_title('Lmon_sample_B.',monitor_ALLW_ref);
plot(LAMBDA_B_ALLW_deg/LAMBDA_B_ALLW_deg_ref,'r')
plot(LAMBDA_B_ALLW/LAMBDA_B_ALLW_ref,'k')
ylabel('Brilliance transfer')
xlabel('Wavelength [AA]')
set(gca,'Fontsize',FS)
title(['Wavelength dependence, + are wavelengths for 1d graphs. I=' num2str(LAMBDA_B_RAW)],'fontsize',FS_medium)
hold off
box
else
figure(Foverall)
subplot(4,2,1:2)
hold on
LAMBDA_B_ALLW_deg=assign_by_title('Lmon_sample_B.',monitor_ALLW_degraded);
LAMBDA_B_ALLW=assign_by_title('Lmon_sample_B.',monitor_ALLW);
plot(LAMBDA_B_ALLW_deg,'r')
plot(LAMBDA_B_ALLW,'k')
ylim([0 LAMBDA_B_ALLW.Data.Max*1.1+1e-10])
set(gca,'Fontsize',FS)
ylabel('Intensity in FOM [Arb]')
xlabel('Wavelength [AA]')
title(['Wavelength dependence, + are wavelengths for 1d graphs. I=' num2str(LAMBDA_B_RAW)],'fontsize',FS_medium)
box on
hold off
end

figure(background)
subplot(4,2,1:4)
LAMBDA_ALLW=assign_by_title('Lmon_sample.',monitor_ALLW);
guide_end_lambda_ALLW = assign_by_title('Lmon_guide_end.',monitor_ALLW);
plot(LAMBDA_ALLW/guide_end_lambda_ALLW,'k');
hold on
plot(LAMBDA_B_ALLW/guide_end_lambda_ALLW,'b');
hold off
box on
title('Signal compared to total guide output') 
xlabel('Wavelength [AA]')
ylabel('Signal fraction [Unitless]') 
set(gca,'Fontsize',FS)

print(Foverall,'-dpng','-r300',[cpath '/' filename '_overall_pure.png'])
print(Fposdiv,'-dpng','-r300',[cpath '/' filename '_posdiv_pure.png'])
print(Facceptance,'-dpng','-r300',[cpath '/' filename '_acceptance_pure.png'])


if exist('OnlyPureFigs') < 0.5
close(Foverall);close(Fposdiv);close(Facceptance);
MaxIndex=length(wavecenters)+1;
Foverall_ess=figure
set(Foverall_ess, 'Position', [0 0 700 1000])
subplot(4,2,1)
if select==1; set(gcf, 'Renderer', 'painters'); end;
Fposdiv_ess=figure
set(Fposdiv_ess, 'Position', [0 0 700 1000])
subplot(MaxIndex,2,1)
if select==1; set(gcf, 'Renderer', 'painters'); end;
Facceptance_ess=figure
set(Facceptance_ess, 'Position', [0 0 700 1000])
subplot(MaxIndex,2,1)
if select==1; set(gcf, 'Renderer', 'painters'); end;

colors={'r' 'g' 'b' 'k' 'm'};
linecolor='k';
linethick=2;

Wnames=fieldnames(monitor_W_ess);
max_HPSD = 0;
max_VPSD = 0;
max_HDIV = 0;
max_VDIV = 0;

% if detailed_absolute_run  % Commented out for now, end statement needs to
                            % be added as well though. Not sure where to
                            % put it until Mads patches mcstas_bot bug

for i=1:MaxIndex
   if i==MaxIndex
    monitors=monitor_ESSW;
   else
    monitors=monitor_W_ess.(Wnames{i});
   end
   
   guide_end_lambda = assign_by_title('Lmon_guide_end.',monitors);
   DIV2D=assign_by_title('Div2d_sample.',monitors)
   PSD2D=assign_by_title('PSD_sample.',monitors);
   HPSD=assign_by_title('HPSD_sample.',monitors);
   HPSD{1} = HPSD{1}.*100; %Conversion to cm
   VPSD=assign_by_title('VPSD_sample.',monitors);
   VPSD{1} = VPSD{1}.*100; %Conversion to cm
   HDIV=assign_by_title('Hdiv_sample.',monitors);
   VDIV=assign_by_title('Vdiv_sample.',monitors);
   HACCP=assign_by_title('acceptance_x_divx.',monitors);
   HACCP{2} = HACCP{2}.*100; %Conversion to cm
   VACCP=assign_by_title('acceptance_y_divy.',monitors);
   VACCP{2} = VACCP{2}.*100; %Conversion to cm
   LAMBDA_B=assign_by_title('Lmon_sample_B.',monitors);
   LAMBDA_B_RAW_monitor=assign_by_title('Lmon_sample_B.',monitors);
   LAMBDA_B_RAW=LAMBDA_B_RAW_monitor.Data.values(1);
   DIV2D_B=assign_by_title('Div2d_sample_maxdiv.',monitors);
   PSD2D_B=assign_by_title('PSD_sample_maxdiv.',monitors);
   HPSD_B=assign_by_title('HPSD_sample_maxdiv.',monitors);
   HPSD_B{1} = HPSD_B{1}.*100; %Conversion to cm
   VPSD_B=assign_by_title('VPSD_sample_maxdiv.',monitors);
   VPSD_B{1} = VPSD_B{1}.*100; %Conversion to cm
   HDIV_B=assign_by_title('Hdiv_sample_maxdiv.',monitors);
   VDIV_B=assign_by_title('Vdiv_sample_maxdiv.',monitors);
   HACCP_B=assign_by_title('acceptance_x_divx_maxdiv.',monitors);
   HACCP_B{2} = HACCP_B{2}.*100; %Conversion to cm
   VACCP_B=assign_by_title('acceptance_y_divy_maxdiv.',monitors);
   VACCP_B{2} = VACCP_B{2}.*100; %Conversion to cm
   LAMBDA=assign_by_title('Lmon_sample.',monitors);
   LAMBDA_RAW_monitor=assign_by_title('Lmon_sample.',monitors);
   
   if i==MaxIndex
      figure(Foverall_ess);
      subplot(4,2,7:8)
      axis off;
      text(0.5,1,Project_name,'interpreter','none')
      text(0,0.83,[instrument_name ' - ' filename ' - '  inputstring],'interpreter','none')
      text(0,0.66,['Sample size: Horizontal=' num2str(Hsize) 'cm, Vertical=' num2str(Vsize) 'cm'],'interpreter','none');
      text(0,0.5,['Divergence requirement: Horizontal=' num2str(Hdiv) 'deg, Vertical= ' num2str(Vdiv) 'deg'],'interpreter','none');
      text(0,0.33,['Intensity on sample of 100 emitted = ' num2str(LAMBDA.Data.values(1)) ' (no divergence limit) ' num2str(LAMBDA_B_RAW) ' (width divergence limits)' ],'interpreter','none');
      text(0,0.16,['Intensity near sample of 100 emitted = ' num2str(PSD2D.Data.values(1)) ' (no divergence limit)' ],'interpreter','none');
   
      figure(Foverall_ess)
      subplot(4,2,1:2) 
     figure(Foverall_ess)
     subplot(4,2,1:2)
     hold on
     sample_area_cm2 = 10000*str2num_safe(LAMBDA.Data.Parameters.sizeX)*str2num_safe(LAMBDA.Data.Parameters.sizeY);
     norm_factor = LAMBDA.Data.array_1d/(LAMBDA.Data.xlimits(2)-LAMBDA.Data.xlimits(1))/sample_area_cm2;

     plot(LAMBDA*norm_factor,'k')
     ylim([0 LAMBDA.Data.Max*norm_factor*1.1+1e-10]);
     hold off
      %set(gca,'XTick',[0 0.5 1.0 1.5 2.0 2.5 3.0])
      ylabel('Flux n/s/cm^2/AA')
      xlabel('Wavelength [AA]')
      set(gca,'Fontsize',FS)
      title(['Wavelength dependence, + are wavelengths for 1d graphs. I=' num2str(LAMBDA_B_RAW)],'fontsize',FS_medium)
      box
      markerheight=LAMBDA.Data.Max*norm_factor*0.1;
        for j=1:length(wavecenters)
                hold on
                plot(wavecenters(j),markerheight,['+' colors{j}],'MarkerSize',12)
                hold off
        end
   else
   
    figure(Foverall_ess)
    subplot(4,2,3)
    hold on
    plot(HDIV,colors{i})
    box
     max_HDIV = max([max_HDIV HDIV.Data.Max]);
    if i==MaxIndex-1;
      plot([-Hdiv -Hdiv],[0 max_HDIV*1.1],'--k');
      plot([Hdiv Hdiv],[0 max_HDIV*1.1],'--k');
      axis([-div1d*Hdiv div1d*Hdiv 0 max_HDIV*1.1]);
     end
    %set(gca,'XTick',[-1 -0.5 -0.3 0 0.3 0.5 1])
    set(gca,'Fontsize',FS)
    title('')
    xlabel('Horizontal divergence [deg]')
    ylabel('Arb intensity unit')
    set(gca,'Fontsize',FS)

    subplot(4,2,4)
    hold on
    plot(VDIV,colors{i})
    box
     max_VDIV = max([max_VDIV VDIV.Data.Max]);
    if i==MaxIndex-1;
      plot([-Vdiv -Vdiv],[0 max_VDIV*1.1],'--k');
      plot([Vdiv Vdiv],[0 max_VDIV*1.1],'--k');
      axis([-div1d*Vdiv div1d*Vdiv 0 max_VDIV*1.1]);
     end
    %set(gca,'XTick',[-1 -0.5 -0.3 0 0.3 0.5 1])
    title('')
    xlabel('Vertical divergence [deg]')
    ylabel('Arb intensity unit')
    set(gca,'Fontsize',FS)

    subplot(4,2,5)
    hold on
    monitor_area_cm2 = 10000*str2num_safe(HPSD.Data.Parameters.sizeY)*(HPSD.Data.xlimits(2)-HPSD.Data.xlimits(1));
    number_of_bins = HPSD.Data.array_1d;
    wavelength_band = str2num_safe(HPSD.Data.Parameters.WaveMax)-str2num_safe(HPSD.Data.Parameters.WaveMin);
    norm_factor = number_of_bins/monitor_area_cm2/wavelength_band;
    plot(HPSD*norm_factor,colors{i})
    box
    %axis([-0.75 0.75 0 1])
     max_HPSD = max([max_HPSD HPSD.Data.Max]);
    if i==MaxIndex-1;
      plot([-0.5*Hsize -0.5*Hsize],[0 max_HPSD*norm_factor*1.1],'--k');
      plot([0.5*Hsize 0.5*Hsize],[0 max_HPSD*norm_factor*1.1],'--k');
      axis([-0.5*psd1d*Hsize 0.5*psd1d*Hsize 0 max_HPSD*norm_factor*1.1]);
     end

    %set(gca,'XTick',[-0.75 -0.5 -0.25 0 0.25 0.5 0.75])
    title('')
    xlabel('Horizontal position [cm]')
    ylabel('Flux n/s/cm^2/AA')
    set(gca,'Fontsize',FS)

    subplot(4,2,6)
    hold on
    monitor_area_cm2 = 10000*str2num_safe(VPSD.Data.Parameters.sizeX)*(VPSD.Data.xlimits(2)-VPSD.Data.xlimits(1));
    number_of_bins = VPSD.Data.array_1d;
    wavelength_band = str2num_safe(VPSD.Data.Parameters.WaveMax)-str2num_safe(VPSD.Data.Parameters.WaveMin);
    norm_factor = number_of_bins/monitor_area_cm2/wavelength_band;
    plot(VPSD*norm_factor,colors{i})
    box
    %axis([-1 1 0 1])
     max_VPSD = max([max_VPSD VPSD.Data.Max]);
    if i==MaxIndex-1;
      plot([-0.5*Vsize -0.5*Vsize],[0 max_VPSD*norm_factor*1.1],'--k');
      plot([0.5*Vsize 0.5*Vsize],[0 max_VPSD*norm_factor*1.1],'--k');
      axis([-0.5*psd1d*Vsize 0.5*psd1d*Vsize 0 max_VPSD*norm_factor*1.1]);
     end
    %set(gca,'XTick',[-1 -0.5 -0.25 0 0.25 0.5 1])
    title('')
    xlabel('Vertical position [cm]')
    ylabel('Flux n/s/cm^2/AA') 
    set(gca,'Fontsize',FS)
   
   end
   
   
    figure(Fposdiv_ess)
    subplot(MaxIndex,2,(i-1)*2+1)
    plot(DIV2D)
    box
    view([-0.5 90]);
    %axis([-1 1 -1 1])
   axis([-div2d*Hdiv div2d*Hdiv -div2d*Vdiv div2d*Vdiv])
    maxi=max(DIV2D)+100;
    hold on
    x=Hdiv;
    y=Vdiv;
    plot3([-x x], [-y -y], [maxi maxi],'color',linecolor,'LineWidth',linethick)
    plot3([-x x], [y y], [maxi maxi],'color',linecolor,'LineWidth',linethick)
    plot3([-x -x], [-y y], [maxi maxi],'color',linecolor,'LineWidth',linethick)
    plot3([x x], [-y y], [maxi maxi],'color',linecolor,'LineWidth',linethick)
    hold off
    %set(gca,'XTick',[-1 -0.5 -0.3 0 0.3 0.5 1])
    %set(gca,'YTick',[-1 -0.5 -0.3 0 0.3 0.5 1])
    xlabel('Horizontal divergence [deg]');
    set(gca,'Fontsize',FS_medium)
    ylabel('Vertical divergence [deg]','fontsize',FS_medium);
    if i==MaxIndex
        w_min = DIV2D.Data.WaveMin;w_max=DIV2D.Data.WaveMax;
        title(['2d div ' num2str(w_min) '-' num2str(w_max) ' Å'],'fontsize',FS_small)
    else
        title(['2d div Lambda=' num2str(wavecenters(i)) 'A'],'fontsize',FS_small)
    end


    subplot(MaxIndex,2,(i-1)*2+2)
    plot(PSD2D)
    box
    view([-0.5 90]);
    %axis([-1 1 -1 1])
   axis([-0.5*psd2d*Hsize 0.5*psd2d*Hsize -0.5*psd2d*Vsize 0.5*psd2d*Vsize])
    maxi=max(PSD2D)+100;
    hold on
    x=0.5*Hsize;
    y=0.5*Vsize;
    plot3([-x x], [-y -y], [maxi maxi],'color',linecolor,'LineWidth',linethick)
    plot3([-x x], [y y], [maxi maxi],'color',linecolor,'LineWidth',linethick)
    plot3([-x -x], [-y y], [maxi maxi],'color',linecolor,'LineWidth',linethick)
    plot3([x x], [-y y], [maxi maxi],'color',linecolor,'LineWidth',linethick)
    hold off
    %set(gca,'XTick',[-1 -0.5 -0.3 0 0.3 0.5 1])
    %set(gca,'YTick',[-1 -0.5 -0.3 0 0.3 0.5 1])
    xlabel('Horizontal position [cm]');
    ylabel('Vertical position [cm]','fontsize',FS_medium);
    set(gca,'Fontsize',FS_medium)
    if i==MaxIndex
        w_min = PSD2D.Data.WaveMin;w_max=PSD2D.Data.WaveMax;
        title(['2d psd ' num2str(w_min) '-' num2str(w_max) ' Å'],'fontsize',FS_small)
    else
        title(['2d psd Lambda=' num2str(wavecenters(i)) 'A'],'fontsize',FS_small)
    end
   
    
    figure(Facceptance_ess)
    subplot(MaxIndex,2,(i-1)*2+1)
    plot(HACCP)
    box
    view([-0.5 90]);
    %axis([-1 1 -1 1])
   axis([-0.5*apsd*Hsize 0.5*apsd*Hsize -adiv*Hdiv adiv*Hdiv])
    maxi=max(HACCP)+100;
    hold on
    x=0.5*Hsize;
    y=Hdiv;
    plot3([-x x], [-y -y], [maxi maxi],'color',linecolor,'LineWidth',linethick)
    plot3([-x x], [y y], [maxi maxi],'color',linecolor,'LineWidth',linethick)
    plot3([-x -x], [-y y], [maxi maxi],'color',linecolor,'LineWidth',linethick)
    plot3([x x], [-y y], [maxi maxi],'color',linecolor,'LineWidth',linethick)
    hold off
    %set(gca,'XTick',[-1 -0.5 -0.3 0 0.3 0.5 1])
    %set(gca,'YTick',[-1 -0.5 -0.3 0 0.3 0.5 1])
    xlabel('Horizontal position [cm]')
    set(gca,'Fontsize',FS_medium)
    ylabel('Hor. divergence [deg]','fontsize',FS_medium)
    if i==MaxIndex
    w_min = HACCP.Data.WaveMin;w_max=HACCP.Data.WaveMax;
    title(['Horizontal accep. dia. ' num2str(w_min) '-' num2str(w_max) ' Å'],'fontsize',FS_small)
    else
    title(['Horizontal accep. dia. Wavelength=' num2str(wavecenters(i)) 'A'],'fontsize',FS_small)
    end    

    subplot(MaxIndex,2,(i-1)*2+2)
    plot(VACCP)
    box
    view([-0.5 90]);
    %axis([-1 1 -1 1])
   axis([-0.5*apsd*Vsize 0.5*apsd*Vsize -adiv*Vdiv adiv*Vdiv])
    maxi=max(VACCP)+100;
    hold on
    x=0.5*Vsize;
    y=Vdiv;
    plot3([-x x], [-y -y], [maxi maxi],'color',linecolor,'LineWidth',linethick)
    plot3([-x x], [y y], [maxi maxi],'color',linecolor,'LineWidth',linethick)
    plot3([-x -x], [-y y], [maxi maxi],'color',linecolor,'LineWidth',linethick)
    plot3([x x], [-y y], [maxi maxi],'color',linecolor,'LineWidth',linethick)
    hold off
    %set(gca,'XTick',[-1 -0.5 -0.3 0 0.3 0.5 1])
    %set(gca,'YTick',[-1 -0.5 -0.3 0 0.3 0.5 1])
    xlabel('Vertical position [cm]')
    set(gca,'Fontsize',FS_medium)
    ylabel('Ver. divergence [deg]','fontsize',FS_medium)
    if i==MaxIndex
    w_min = VACCP.Data.WaveMin;w_max=VACCP.Data.WaveMax;
    title(['Vertical accep. dia. ' num2str(w_min) '-' num2str(w_max) ' Å'],'fontsize',FS_small)
    else
    title(['Vertical accep. dia. Wavelength=' num2str(wavecenters(i)) 'A'],'fontsize',FS_small)
    end
    
    
   if i==MaxIndex
      figure(background)
      subplot(4,2,5:8)
      p1 = plot(LAMBDA_RAW_monitor/guide_end_lambda,'k');
      hold on
      p2 = plot(LAMBDA_B_RAW_monitor/guide_end_lambda,'b');
      hold off
      box on
      title('Signal compared to total guide output')
      legend([p1 p2], {'No Divergence Limits', 'Divergence Limits'})
      xlabel('Wavelength [AA]')
      ylabel('Signal fraction [Unitless]') 
      set(gca,'Fontsize',FS)
    end
    
end

flux=monitor_fom(11).Data.values(1);
fluxtext=' Absolute ';


set(Foverall_ess, 'paperpositionmode', 'auto');
set(Fposdiv_ess, 'paperpositionmode', 'auto');
set(Facceptance_ess, 'paperpositionmode', 'auto');
set(background, 'paperpositionmode', 'auto');


print(Foverall_ess,'-dpng','-r300',[cpath '/' filename '_overall_ess.png'])
print(Fposdiv_ess,'-dpng','-r300',[cpath '/' filename '_posdiv_ess.png'])
print(Facceptance_ess,'-dpng','-r300',[cpath '/' filename '_acceptance_ess.png'])
print(background,'-dpng','-r300',[cpath '/' filename '_background.png'])
end


if exist('FOM_table.m')
  close all;
  FOM_table;
end

