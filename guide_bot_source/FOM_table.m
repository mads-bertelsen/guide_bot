try
    baseStruct = load('Baseline_Table.mat');
    baseline = baseStruct.t_data;
end

f_temp = figure; sp6 = subplot(3,3,8); pos6 = get(sp6, 'Position'); close(f_temp);
f_table = figure('units','normalized','outerposition',[0.2 0.2 0.75 0.75]);
%f_table = figure('Position', [0 0 1000 1000]);
t = uitable(f_table, 'Units', 'normalized', 'FontSize', 11, 'ColumnName', {'Ei (meV)', 'FOM (n/s/cm2)', '3x Div Limit (n/s/cm2)', 'Resolution (meV)', 'FOM^2/Resolution'});
set(t, 'Position', [0.5 - pos6(3) pos6(2) 2*pos6(3) pos6(4)]);

allwaves = fields(monitor_W_ess);
for i = 1:numel(allwaves)
    wavei = monitor_W_ess.(allwaves{i});
    wmon = assign_by_title('Lmon_sample.dat',wavei);
    wavemon = wmon;
    wavemon.Error = 1;
    center_e = 81.81/((wavemon.WaveMax + wavemon.WaveMin)/2)^2;
    amplitude = 0.5*max(wavemon.Signal);
    wp = fits(wavemon, gauss, [amplitude,wavemon.X0,wavemon.dX,0], 'fminsimplex', [0,0,0,1]);
    final_fit = wavemon(gauss,wp);
    fwhm_e = 81.81*4*wp(2)*wp(3)/(wp(2)^2-wp(3)^2)^2;
    
    wavemon_area_cm2 = 10000*str2num_safe(wavemon.Data.Parameters.sizeX)*str2num_safe(wavemon.Data.Parameters.sizeY);
    wavemon_norm_factor = wavemon.Data.array_1d/(81.81/wavemon.Data.xlimits(1)^2-81.81/wavemon.Data.xlimits(2)^2)/wavemon_area_cm2;
    
    wx = 81.81./wavemon.x.^2;
    wy = wavemon.Signal*wavemon_norm_factor;
    err = wmon.Error*wavemon_norm_factor;
    
    fx = 81.81./final_fit.x.^2;
    fy = final_fit.Signal*wavemon_norm_factor;
    
    subplot(3,3,i)
    errorbar(wx,wy,err,'b')
    hold on
    p=plot(fx,fy,'r');
    hold off
    set(p,'LineWidth',2)
    xlabel('Energy (meV)')
    ylabel('Flux at Sample (n/s/cm^{2}/meV)')
    title(['E_{i} = ' num2str(center_e) 'meV'])
    
    fomdivmon = assign_by_title('Div2d_sample_B.dat', wavei);
    fulldivmon = assign_by_title('Div2d_sample.dat', wavei);
    
    
    t_data{i,1} = center_e;
    p_data(i,1) = center_e;
    t_data{i,2} = fomdivmon.values(1)/10000/fomdivmon.sizeX/fomdivmon.sizeY;
    p_data(i,2) = fomdivmon.values(1)/10000/fomdivmon.sizeX/fomdivmon.sizeY;
    t_data{i,3} = fulldivmon.values(1)/10000/fomdivmon.sizeX/fomdivmon.sizeY;
    p_data(i,3) = fulldivmon.values(1)/10000/fomdivmon.sizeX/fomdivmon.sizeY;
    t_data{i,4} = fwhm_e;
    p_data(i,4) = fwhm_e;
    t_data{i,5} = t_data{i,2}^2/fwhm_e;
    
end
t_data{6,1} = 'Average';
t_data{6,2} = (t_data{1,2} + t_data{2,2} + t_data{3,2} + t_data{4,2} + t_data{5,2})/5;
t_data{6,3} = (t_data{1,3} + t_data{2,3} + t_data{3,3} + t_data{4,3} + t_data{5,3})/5;
t_data{6,4} = (t_data{1,4} + t_data{2,4} + t_data{3,4} + t_data{4,4} + t_data{5,4})/5;
t_data{6,5} = (t_data{1,5} + t_data{2,5} + t_data{3,5} + t_data{4,5} + t_data{5,5})/5;


subplot(3,3,6)
%[fax,l1,l2] = plotyy(p_data(:,1), p_data(:,2),p_data(:,1),p_data(:,4));
[fax,l1,l2] = plotyy([p_data(:,1), p_data(:,1)], [p_data(:,2), p_data(:,3)], p_data(:,1), p_data(:,4));
set(l1, 'Marker', 'o')
set(l2, 'Marker', 'o')
f_leg = legend('FOM', '3x Div', 'Location', 'SouthEast');
set(get(fax(1),'YLabel'),'String', 'Flux at Sample (n/s/cm^{2})')
set(get(fax(2),'YLabel'),'String', 'Resolution (meV)')
xlabel('E_{i} (meV)');
set(t, 'Data', t_data);
set(f_table, 'paperpositionmode', 'auto');
print(f_table,'-dpng','-r300',[filename '_FOM_table.png'])
tablename = [filename '_FOM_table_data.mat'];
save([filename '_FOM_table_data'], 't_data', 'tablename')

if exist('baseStruct', 'var')
    f_temp = figure; 
    sp12 = subplot(3,3,1:2); pos12 = get(sp12, 'Position');
    sp45 = subplot(3,3,4:5); pos45 = get(sp45, 'Position');
    sp78 = subplot(3,3,7:8); pos78 = get(sp78, 'Position');
    close(f_temp);
   
    norm_table = figure('units','normalized','outerposition',[0.2 0.2 0.75 0.75]);
    bk_col = get(norm_table, 'Color');
    uih = uicontrol('Style', 'text', 'Units', 'normalized', 'Position', [pos12(1) 0.92 0.3 0.05], 'String', ['Baseline: ' baseStruct.tablename], 'FontSize', 12, 'BackgroundColor', bk_col);
    t_base = uitable('Units', 'normalized', 'FontSize', 11, 'Position', pos12, 'ColumnName', {'Ei (meV)', 'FOM (n/s/cm2)', '3x Div Limit (n/s/cm2)', 'Resolution (meV)', 'FOM^2/Resolution'});
    t_orig = uitable('Units', 'normalized', 'FontSize', 11, 'Position', pos45, 'ColumnName', {'Ei (meV)', 'FOM (n/s/cm2)', '3x Div Limit (n/s/cm2)', 'Resolution (meV)', 'FOM^2/Resolution'});
    t_norm = uitable('Units', 'normalized', 'FontSize', 11, 'Position', pos78, 'ColumnName', {'Ei (meV)', 'FOM (n/s/cm2)', '3x Div Limit (n/s/cm2)', 'Resolution (meV)', 'FOM^2/Resolution'});
    
    for i = 2:5
        for j = 1:5
            t_data_norm{j,i} = t_data{j,i}/baseline{j,i};
            p_data_norm(j,i) = t_data{j,i}/baseline{j,i};
            p_baseline(j,i) = baseline{j,i};
        end
    end
    
    for i = 1:5
        t_data_norm{i,1} = t_data{i,1};
        p_data_norm(i,1) = t_data{i,1};
        p_baseline(i,1) = baseline{i,1};
    end
    
    t_data_norm{6,1} = 'Average';
    t_data_norm{6,2} = (t_data_norm{1,2} + t_data_norm{2,2} + t_data_norm{3,2} + t_data_norm{4,2} + t_data_norm{5,2})/5;
    t_data_norm{6,3} = (t_data_norm{1,3} + t_data_norm{2,3} + t_data_norm{3,3} + t_data_norm{4,3} + t_data_norm{5,3})/5;
    t_data_norm{6,4} = (t_data_norm{1,4} + t_data_norm{2,4} + t_data_norm{3,4} + t_data_norm{4,4} + t_data_norm{5,4})/5;
    t_data_norm{6,5} = (t_data_norm{1,5} + t_data_norm{2,5} + t_data_norm{3,5} + t_data_norm{4,5} + t_data_norm{5,5})/5;
    
    
    t_data{7,1} = 'OPT';
    baseline{7,1} = 'BASELINE';
    t_data_norm{7,1} = 'RATIO';
    
    set(t_orig, 'Data', t_data);
    set(t_base, 'Data', baseline);
    set(t_norm, 'Data', t_data_norm);
    
    subplot(3,3,3)
    [fnax,ln1,ln2] = plotyy([p_data(:,1), p_baseline(:,1)], [p_data(:,2), p_baseline(:,2)], p_data_norm(:,1), p_data_norm(:,2));
    set(ln1, 'Marker', 'o')
    set(ln2, 'Marker', 'o')
    fn_leg = legend('Optimization', 'Baseline', 'Ratio', 'Location', 'SouthEast');
    set(get(fnax(1),'YLabel'),'String', 'Flux at Sample (n/s/cm^{2})')
    set(get(fnax(2),'YLabel'),'String', 'Ratio')
    xlabel('E_{i} (meV)');
    title('FOM')
    
    subplot(3,3,6)
    [fnax,ln1,ln2] = plotyy([p_data(:,1), p_baseline(:,1)], [p_data(:,3), p_baseline(:,3)], p_data_norm(:,1), p_data_norm(:,3));
    set(ln1, 'Marker', 'o')
    set(ln2, 'Marker', 'o')
    fn_leg = legend('Optimization', 'Baseline', 'Ratio', 'Location', 'SouthEast');
    set(get(fnax(1),'YLabel'),'String', 'Flux at Sample (n/s/cm^{2})')
    set(get(fnax(2),'YLabel'),'String', 'Ratio')
    xlabel('E_{i} (meV)');
    title('3x Divergence Limit')
    
    subplot(3,3,9)
    [fnax,ln1,ln2] = plotyy([p_data(:,1), p_baseline(:,1)], [p_data(:,4), p_baseline(:,4)], p_data_norm(:,1), p_data_norm(:,4));
    set(ln1, 'Marker', 'o')
    set(ln2, 'Marker', 'o')
    fn_leg = legend('Optimization', 'Baseline', 'Ratio', 'Location', 'SouthEast');
    set(get(fnax(1),'YLabel'),'String', 'Resolution (meV)')
    set(get(fnax(2),'YLabel'),'String', 'Ratio')
    xlabel('E_{i} (meV)');
    title('Resolution')
    
    set(norm_table, 'paperpositionmode', 'auto');
    print(norm_table,'-dpng','-r300',[filename '_FOM_table_NORM.png'])
end