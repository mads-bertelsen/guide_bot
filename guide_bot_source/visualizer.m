function visualizer(file_name)


data_name = [file_name '_geometry.dat'];
if exist(file_name) > 0.5
   data_name = file_name;
end

if exist(data_name) > 0.5
% ------------------------------------------------------------------------
% First part reads the data from file 
fid = fopen(data_name);
tline = fgets(fid); % gets the first line
tline = tline(1:end-1);
module_num = 0;

%%%%%%%%%%%% LELAND MODIFICATION %%%%%%%%%%%%
% Getting constuctor lines with more than 1 extra character
% tline = tline(1:end-1); % Original line commented out by Leland
tline = strtrim(tline); % New line added by Leland
%%%%%%%%%%%% LELAND MODIFICATION %%%%%%%%%%%%

while ischar(tline);
    module_num = module_num + 1;
    disp(tline)
    if strcmp(tline,'S') % Working
    
    tline = fgets(fid);tline = tline(1:end-1);    
    Xdata{1} = cell2mat(textscan(tline,'x1: %n \t %n \t %n \t %n'));
    tline = fgets(fid);tline = tline(1:end-1);
    Xdata{2} = cell2mat(textscan(tline,'x2: %n \t %n \t %n \t %n'));
    tline = fgets(fid);tline = tline(1:end-1);
    Ydata{1} = cell2mat(textscan(tline,'y1: %n \t %n \t %n \t %n'));
    tline = fgets(fid);tline = tline(1:end-1);
    Ydata{2} = cell2mat(textscan(tline,'y2: %n \t %n \t %n \t %n'));
        
    geometry{module_num}.name = 'S';
    geometry{module_num}.Xdata = Xdata;
    geometry{module_num}.Ydata = Ydata;
    clear Xdata;clear Ydata;

    
    %%%%%%%%%%%%%% LELAND MODIFICATION %%%%%%%%%%%%%%%
    % Added Monochromator option
    elseif strcmp(tline, 'M')
        % Read in params (z0, monLength, startx, starty, Lx, Ly)
        tline = fgets(fid);tline = tline(1:end-1);
        linevals1 = textscan(tline, 'Dim: %f \t %f \t %f');
        tline = fgets(fid);tline = tline(1:end-1);
        % linevals2 = textscan(tline, 'Lf: %f \t %f');
        linevals2 = textscan(tline, 'Lf: %s \t %s');
        linevals2 = {str2num(linevals2{1}{1}) str2num(linevals2{2}{1})};
        tline = fgets(fid);tline = tline(1:end-1);
        linevals3 = textscan(tline, 'zp: %f');

        % Assign parameters to fields
        geometry{module_num}.name = 'M';
        geometry{module_num}.monLength = linevals1{1};
        geometry{module_num}.startx = linevals1{2};
        geometry{module_num}.starty = linevals1{3};
        geometry{module_num}.Lx = linevals2{1};
        geometry{module_num}.Ly = linevals2{2};
        geometry{module_num}.z0 = linevals3{1};

        % Remove temporary variables
        clear linevals1 linevals2 linevals3
    %%%%%%%%%%% END OF LELAND MODIFICATION %%%%%%%%%%%
        
    
    elseif strcmp(tline,'P')
    
    tline = fgets(fid);tline = tline(1:end-1);
    Xdata{1} = cell2mat(textscan(tline,'x1: %n \t %n \t %n \t %n'));
    tline = fgets(fid);tline = tline(1:end-1);
    Xdata{2} = cell2mat(textscan(tline,'x2: %n \t %n \t %n \t %n'));
    tline = fgets(fid);tline = tline(1:end-1);
    Ydata{1} = cell2mat(textscan(tline,'y1: %n \t %n \t %n \t %n'));
    tline = fgets(fid);tline = tline(1:end-1);
    Ydata{2} = cell2mat(textscan(tline,'y2: %n \t %n \t %n \t %n'));
    tline = fgets(fid);tline = tline(1:end-1);
    Xpos = cell2mat(textscan(tline,'xpos: %n \t %n'));
    tline = fgets(fid);tline = tline(1:end-1);
    Ypos = cell2mat(textscan(tline,'ypos: %n \t %n'));
    tline = fgets(fid);tline = tline(1:end-1);
    Xdir = cell2mat(textscan(tline,'xdir: %n \t %n'));
    tline = fgets(fid);tline = tline(1:end-1);
    Ydir = cell2mat(textscan(tline,'ydir: %n \t %n'));
    tline = fgets(fid);tline = tline(1:end-1);
    Xpars = cell2mat(textscan(tline,'xpars: %n \t %n \t %n \t %n'));
    tline = fgets(fid);tline = tline(1:end-1);
    Ypars = cell2mat(textscan(tline,'ypars: %n \t %n \t %n \t %n'));
    
    geometry{module_num}.name = 'P';
    geometry{module_num}.Xdata = Xdata;
    geometry{module_num}.Ydata = Ydata;    
    geometry{module_num}.Xpos = Xpos;
    geometry{module_num}.Ypos = Ypos;
    geometry{module_num}.Xdir = Xdir;
    geometry{module_num}.Ydir = Ydir;
    geometry{module_num}.Xpars = Xpars;
    geometry{module_num}.Ypars = Ypars;    
    
    elseif strcmp(tline,'E')

    tline = fgets(fid);tline = tline(1:end-1);    
    Xdata{1} = cell2mat(textscan(tline,'x1: %n \t %n \t %n \t %n'));
    tline = fgets(fid);tline = tline(1:end-1);
    Xdata{2} = cell2mat(textscan(tline,'x2: %n \t %n \t %n \t %n'));
    tline = fgets(fid);tline = tline(1:end-1);
    Ydata{1} = cell2mat(textscan(tline,'y1: %n \t %n \t %n \t %n'));
    tline = fgets(fid);tline = tline(1:end-1);
    Ydata{2} = cell2mat(textscan(tline,'y2: %n \t %n \t %n \t %n'));
    tline = fgets(fid);tline = tline(1:end-1);
    Xpos = cell2mat(textscan(tline,'xpos: %n \t %n'));
    tline = fgets(fid);tline = tline(1:end-1);
    Ypos = cell2mat(textscan(tline,'ypos: %n \t %n'));
    tline = fgets(fid);tline = tline(1:end-1);
    Xdir = cell2mat(textscan(tline,'xdir: %n \t %n'));
    tline = fgets(fid);tline = tline(1:end-1);
    Ydir = cell2mat(textscan(tline,'ydir: %n \t %n'));
    tline = fgets(fid);tline = tline(1:end-1);
    Xpars = cell2mat(textscan(tline,'xpars: %n \t %n \t %n \t %n'));
    tline = fgets(fid);tline = tline(1:end-1);
    Ypars = cell2mat(textscan(tline,'ypars: %n \t %n \t %n \t %n'));
        
    geometry{module_num}.name = 'E';
    geometry{module_num}.Xdata = Xdata;
    geometry{module_num}.Ydata = Ydata;
    geometry{module_num}.Xpos = Xpos;
    geometry{module_num}.Ypos = Ypos;
    geometry{module_num}.Xdir = Xdir;
    geometry{module_num}.Ydir = Ydir;
    geometry{module_num}.Xpars = Xpars;
    geometry{module_num}.Ypars = Ypars;
    
    elseif strcmp(tline,'K') % Working

    tline = fgets(fid);tline = tline(1:end-1);    
    Xdata{1} = cell2mat(textscan(tline,'x1: %n \t %n \t %n \t %n'));
    tline = fgets(fid);tline = tline(1:end-1);
    Xdata{2} = cell2mat(textscan(tline,'x2: %n \t %n \t %n \t %n'));
    tline = fgets(fid);tline = tline(1:end-1);
    Ydata{1} = cell2mat(textscan(tline,'y1: %n \t %n \t %n \t %n'));
    tline = fgets(fid);tline = tline(1:end-1);
    Ydata{2} = cell2mat(textscan(tline,'y2: %n \t %n \t %n \t %n'));
        
    geometry{module_num}.name = 'K';
    geometry{module_num}.Xdata = Xdata;
    geometry{module_num}.Ydata = Ydata;    
    
    
    elseif strcmp(tline,'C') % Working

    tline = fgets(fid);tline = tline(1:end-1);   
    Xdata{1} = cell2mat(textscan(tline,'x1: %n \t %n \t %n \t %n'));
    tline = fgets(fid);tline = tline(1:end-1);
    Xdata{2} = cell2mat(textscan(tline,'x2: %n \t %n \t %n \t %n'));
    tline = fgets(fid);tline = tline(1:end-1);
    Ydata{1} = cell2mat(textscan(tline,'y1: %n \t %n \t %n \t %n'));
    tline = fgets(fid);tline = tline(1:end-1);
    Ydata{2} = cell2mat(textscan(tline,'y2: %n \t %n \t %n \t %n'));
        
    tline = fgets(fid);tline = tline(1:end-1);
    circle_data = textscan(tline,'%s %n \t %n \t %n \t %n \t %n');

    if strcmp(circle_data{1}(1),'y:')
     geometry{module_num}.circle_orientation = 'y';
    else
     geometry{module_num}.circle_orientation = 'x';
    end
   
    geometry{module_num}.name = 'C';
    geometry{module_num}.Xdata = Xdata;
    geometry{module_num}.Ydata = Ydata;
    geometry{module_num}.circle_data = cell2mat(circle_data(2:end));
    clear Xdata;clear Ydata;   
  
    
    elseif strcmp(tline,'Selene')

    tline = fgets(fid);tline = tline(1:end-1);    
    Xdata{1} = cell2mat(textscan(tline,'x1: %n \t %n \t %n \t %n'));
    tline = fgets(fid);tline = tline(1:end-1);
    Xdata{2} = cell2mat(textscan(tline,'x2: %n \t %n \t %n \t %n'));
    tline = fgets(fid);tline = tline(1:end-1);
    Ydata{1} = cell2mat(textscan(tline,'y1: %n \t %n \t %n \t %n'));
    tline = fgets(fid);tline = tline(1:end-1);
    Ydata{2} = cell2mat(textscan(tline,'y2: %n \t %n \t %n \t %n'));
    tline = fgets(fid);tline = tline(1:end-1);
    Xpos = cell2mat(textscan(tline,'xpos: %n \t %n \t %n \t %n'));
    tline = fgets(fid);tline = tline(1:end-1);
    Ypos = cell2mat(textscan(tline,'ypos: %n \t %n \t %n \t %n'));
    tline = fgets(fid);tline = tline(1:end-1);
    Pars = cell2mat(textscan(tline,'pars: %n \t %n \t %n'));
    tline = fgets(fid);tline = tline(1:end-1);
    Bxy = cell2mat(textscan(tline,'xy: %n \t %n'));
    
    geometry{module_num}.name = 'Selene';
    geometry{module_num}.Xdata = Xdata;
    geometry{module_num}.Ydata = Ydata;
    geometry{module_num}.Xpos = Xpos;
    geometry{module_num}.Ypos = Ypos;
    geometry{module_num}.Pars = Pars;
    geometry{module_num}.Bxy = Bxy;
    
    elseif strcmp(tline,'G')

    tline = fgets(fid);tline = tline(1:end-1);
    Xdata{1} = cell2mat(textscan(tline,'x1: %n \t %n \t %n \t %n'));
    tline = fgets(fid);tline = tline(1:end-1);
    Xdata{2} = cell2mat(textscan(tline,'x2: %n \t %n \t %n \t %n'));
    tline = fgets(fid);tline = tline(1:end-1);
    Ydata{1} = cell2mat(textscan(tline,'y1: %n \t %n \t %n \t %n'));
    tline = fgets(fid);tline = tline(1:end-1);
    Ydata{2} = cell2mat(textscan(tline,'y2: %n \t %n \t %n \t %n'));
        
    geometry{module_num}.name = 'G';
    geometry{module_num}.Xdata = Xdata;
    geometry{module_num}.Ydata = Ydata;    
        
    elseif strcmp(tline,'Slit')

    tline = fgets(fid);tline = tline(1:end-1);
    Xdata{1} = cell2mat(textscan(tline,'x1: %n \t %n \t %n \t %n'));
    tline = fgets(fid);tline = tline(1:end-1);
    Xdata{2} = cell2mat(textscan(tline,'x2: %n \t %n \t %n \t %n'));
    tline = fgets(fid);tline = tline(1:end-1);
    Ydata{1} = cell2mat(textscan(tline,'y1: %n \t %n \t %n \t %n'));
    tline = fgets(fid);tline = tline(1:end-1);
    Ydata{2} = cell2mat(textscan(tline,'y2: %n \t %n \t %n \t %n'));
        
    geometry{module_num}.name = 'Slit';
    geometry{module_num}.Xdata = Xdata;
    geometry{module_num}.Ydata = Ydata;    
    
    elseif strcmp(tline,'Sample')
        
    tline = fgets(fid);tline = tline(1:end-1);
    data = cell2mat(textscan(tline,'%n \t %n \t %n'));
    tline = fgets(fid);tline = tline(1:end-1);
    data_pos_x = cell2mat(textscan(tline,'xp: %n \t %n'));
    tline = fgets(fid);tline = tline(1:end-1);
    data_pos_y = cell2mat(textscan(tline,'yp: %n \t %n'));
    tline = fgets(fid);tline = tline(1:end-1);
    data_div_x = cell2mat(textscan(tline,'xd: %n \t %n'));
    tline = fgets(fid);tline = tline(1:end-1);
    data_div_y = cell2mat(textscan(tline,'yd: %n \t %n'));
    
    geometry{module_num}.name = 'Sample';
    geometry{module_num}.data = data;
    geometry{module_num}.start_pos_x = data_pos_x;
    geometry{module_num}.start_pos_y = data_pos_y;
    geometry{module_num}.start_div_x = data_div_x;
    geometry{module_num}.start_div_y = data_div_y;
    
    
    elseif strcmp(tline,'Moderator')
        
    tline = fgets(fid);tline = tline(1:end-1);
    data = cell2mat(textscan(tline,'%n \t %n'));
    
    geometry{module_num}.name = 'Moderator';
    geometry{module_num}.data = data;
    
    else
        disp('visualizer error')
    end
    
    %%%%%%%%%%%% LELAND MODIFICATION %%%%%%%%%%%%
    % Was getting constuctor lines with more than 1 extra character
    % tline = fgets(fid);tline = tline(1:end-1); % Original line commented out by Leland (only removes a single extra character)
    tline = fgets(fid);
    if ischar(tline) % strtrim only works on type char, so EOF value -1 will throw an error.
        tline = strtrim(tline);
    else % EOF
        tline = tline(1:end-1);
    end
    %%%%%%%%%%%% LELAND MODIFICATION %%%%%%%%%%%%
    if ischar(tline)
     if strcmp(tline(1),'L')
        if strcmp(tline(end),'s')
          geometry{module_num}.type = 'S';
        elseif strcmp(tline(end),'e')
          geometry{module_num}.type = 'E';
        elseif strcmp(tline(end),'t')
          geometry{module_num}.type = 'T';
        else
            disp('bug, los keyword misused?')
        end
        tline = fgets(fid);tline = tline(1:end-1);
     else
        geometry{module_num}.type = 'N';
     end
    else
        geometry{module_num}.type = 'N';
    end

end

fclose(fid);


% ------------------------------------------------------------------------
% Second part plots the data in matlab

disp('drawing part')
colormatrix = colormap(hot);
colorlength = length(colormatrix);
skip = floor(colorlength/(module_num*2));
colorpoint = skip*(1:module_num);
colorpoint = [1 colorpoint];

lw = 1.5;
fsize = 26;

overview = figure(1);
axes('fontsize',20)
set(overview, 'Position', [200 200 1000 700])
set(overview, 'paperpositionmode', 'auto');

maximum_z = 0;
        
for index = 1:module_num
    
    %color_now = colormatrix(colorpoint(index),:);
    if mod(index,2)
        color_now = 'b';
    else
        color_now = 'k';
    end
    los_start_color = 'g';
    los_end_color = 'm';
    
   switch geometry{index}.name
       case 'S' 
           line1x{1}=[geometry{index}.Xdata{1}(1) geometry{index}.Xdata{2}(1)];
           line1x{2}=[geometry{index}.Xdata{1}(2) geometry{index}.Xdata{2}(2)];
           line2x{1}=[geometry{index}.Xdata{1}(3) geometry{index}.Xdata{2}(3)];
           line2x{2}=[geometry{index}.Xdata{1}(4) geometry{index}.Xdata{2}(4)];
           subplot(2,2,1:2) %X
           hold on
           plot(line1x{1},line1x{2},'color',color_now,'linewidth',lw)
           plot(line2x{1},line2x{2},'color',color_now,'linewidth',lw)
           hold off
           
           line1y{1}=[geometry{index}.Ydata{1}(1) geometry{index}.Ydata{2}(1)];
           line1y{2}=[geometry{index}.Ydata{1}(2) geometry{index}.Ydata{2}(2)];
           line2y{1}=[geometry{index}.Ydata{1}(3) geometry{index}.Ydata{2}(3)];
           line2y{2}=[geometry{index}.Ydata{1}(4) geometry{index}.Ydata{2}(4)];
           subplot(2,2,3:4) %Y
           hold on
           plot(line1y{1},line1y{2},'color',color_now,'linewidth',lw)
           plot(line2y{1},line2y{2},'color',color_now,'linewidth',lw)
           hold off
           
           subplot(2,2,1:2)
           hold on
           if strcmp(geometry{index}.type,'S')
              loswall{1} = [geometry{index}.Xdata{1}(1) geometry{index}.Xdata{1}(3)];
              loswall{2} = [geometry{index}.Xdata{1}(2) geometry{index}.Xdata{1}(4)];
              plot(loswall{1},loswall{2},los_start_color) 
           elseif strcmp(geometry{index}.type,'E')
              loswall{1} = [geometry{index}.Xdata{1}(1) geometry{index}.Xdata{1}(3)];
              loswall{2} = [geometry{index}.Xdata{1}(2) geometry{index}.Xdata{1}(4)];       
              plot(loswall{1},loswall{2},los_end_color) 
           elseif strcmp(geometry{index}.type,'T')
              loswall{1} = [geometry{index}.Xdata{2}(1) geometry{index}.Xdata{2}(3)];
              loswall{2} = [geometry{index}.Xdata{2}(2) geometry{index}.Xdata{2}(4)];     
              plot(loswall{1},loswall{2},los_end_color) 
           end
           hold off
           
           subplot(2,2,3:4)
           hold on
           if strcmp(geometry{index}.type,'S')
              loswall{1} = [geometry{index}.Ydata{1}(1) geometry{index}.Ydata{1}(3)];
              loswall{2} = [geometry{index}.Ydata{1}(2) geometry{index}.Ydata{1}(4)];
              plot(loswall{1},loswall{2},los_start_color) 
           elseif strcmp(geometry{index}.type,'E')
              loswall{1} = [geometry{index}.Ydata{1}(1) geometry{index}.Ydata{1}(3)];
              loswall{2} = [geometry{index}.Ydata{1}(2) geometry{index}.Ydata{1}(4)];       
              plot(loswall{1},loswall{2},los_end_color) 
           elseif strcmp(geometry{index}.type,'T')
              loswall{1} = [geometry{index}.Ydata{2}(1) geometry{index}.Ydata{2}(3)];
              loswall{2} = [geometry{index}.Ydata{2}(2) geometry{index}.Ydata{2}(4)];     
              plot(loswall{1},loswall{2},los_end_color) 
           end
           hold off
           
       %%%%%%%%%%% LELAND MODIFICATION %%%%%%%%%%%%
       % Added Monochromator
       case 'M'
           ang = 0:0.01:2*pi;

           % x-projection plot
           subplot(2,2,1:2)
           hold on
           % Plot Circle
           zp=geometry{index}.monLength/2*sin(ang);
           xp=geometry{index}.monLength/2*cos(ang);
           plot(geometry{index}.z0+zp,xp,':','LineWidth', lw/2)
           % Plot Projection
           plot([geometry{index}.z0 geometry{index}.z0], [-geometry{index}.startx/2 geometry{index}.startx/2], 'r', 'LineWidth', lw)
           % Plot Focal Length
           if ~isempty(geometry{index}.Lx)
               plot([geometry{index}.z0-geometry{index}.Lx geometry{index}.z0],[0 0], 'r', 'LineWidth', lw/2)
           end
           hold off

           % y=projection plot
           subplot(2,2,3:4)
           hold on
           % Plot Projection
           plot([geometry{index}.z0 geometry{index}.z0], [-geometry{index}.starty/2 geometry{index}.starty/2], 'r', 'LineWidth', lw)
           % Plot Focal Length
           if ~isempty(geometry{index}.Ly)
               plot([geometry{index}.z0-geometry{index}.Ly geometry{index}.z0],[0 0], 'r', 'LineWidth', lw/2)
           end
           hold off
           %%%%%%%%%%%%% END OF LELAND MODIFICATION %%%%%%%%%%%%%%
                      
       case 'C'
           % curve in x or in y?
           mean_radius = geometry{index}.circle_data(1);
           inner_radius = geometry{index}.circle_data(2);
           outer_radius = (mean_radius-inner_radius)+mean_radius;
           circle_center_z = geometry{index}.circle_data(3);
           channels = geometry{index}.circle_data(5);
           walls = channels - 1;
           
           bender_mode = 0;
           if walls > 0
               bender_mode = 1;
               for jj = 1:walls
                  wall_radius(jj) = inner_radius + jj/channels*(outer_radius-inner_radius);
               end
           end
                      
           
           if strcmp(geometry{index}.circle_orientation,'x')
               % X  
               circle_center_x = geometry{index}.circle_data(4);
               circle(:,1) = linspace(geometry{index}.Xdata{1}(1),geometry{index}.Xdata{2}(1),100);
               
                if circle_center_x > 0
                     circle(:,2)=circle_center_x-sqrt(inner_radius^2-(circle(:,1)-circle_center_z).^2);
                     circle(:,3)=circle_center_x-sqrt(outer_radius^2-(circle(:,1)-circle_center_z).^2);
                     if bender_mode
                        for jj=1:walls 
                            circle(:,3+jj) =circle_center_x-sqrt(wall_radius(jj)^2-(circle(:,1)-circle_center_z).^2);
                        end
                     end
                else
                     circle(:,2)=circle_center_x+sqrt(inner_radius^2-(circle(:,1)-circle_center_z).^2);
                     circle(:,3)=circle_center_x+sqrt(outer_radius^2-(circle(:,1)-circle_center_z).^2);    
                     if bender_mode
                        for jj=1:walls 
                            circle(:,3+jj) =circle_center_x+sqrt(wall_radius(jj)^2-(circle(:,1)-circle_center_z).^2);    
                        end
                     end
                end

               subplot(2,2,1:2) %X
               hold on
               plot(circle(:,1),circle(:,2),'color',color_now,'linewidth',lw)
               plot(circle(:,1),circle(:,3),'color',color_now,'linewidth',lw)
               if bender_mode
                  for jj=1:walls 
                     eval(['plot(circle(:,1),circle(:,' num2str(3+jj) '),''color'',color_now,''linewidth'',lw)']);
                  end
               end
               hold off
                              
               % Y part normal
               line1y{1}=[geometry{index}.Ydata{1}(1) geometry{index}.Ydata{2}(1)];
               line1y{2}=[geometry{index}.Ydata{1}(2) geometry{index}.Ydata{2}(2)];
               line2y{1}=[geometry{index}.Ydata{1}(3) geometry{index}.Ydata{2}(3)];
               line2y{2}=[geometry{index}.Ydata{1}(4) geometry{index}.Ydata{2}(4)];
               subplot(2,2,3:4) %Y
               hold on
               plot(line1y{1},line1y{2},'color',color_now,'linewidth',lw)
               plot(line2y{1},line2y{2},'color',color_now,'linewidth',lw)
               hold off
               
               subplot(2,2,1:2)
               hold on
               if strcmp(geometry{index}.type,'S')
                  loswall{1} = [geometry{index}.Xdata{1}(1) geometry{index}.Xdata{1}(3)];
                  loswall{2} = [geometry{index}.Xdata{1}(2) geometry{index}.Xdata{1}(4)];
                  plot(loswall{1},loswall{2},los_start_color) 
               elseif strcmp(geometry{index}.type,'E')
                  loswall{1} = [geometry{index}.Xdata{1}(1) geometry{index}.Xdata{1}(3)];
                  loswall{2} = [geometry{index}.Xdata{1}(2) geometry{index}.Xdata{1}(4)];       
                  plot(loswall{1},loswall{2},los_end_color) 
               elseif strcmp(geometry{index}.type,'T')
                  loswall{1} = [geometry{index}.Xdata{2}(1) geometry{index}.Xdata{2}(3)];
                  loswall{2} = [geometry{index}.Xdata{2}(2) geometry{index}.Xdata{2}(4)];     
                  plot(loswall{1},loswall{2},los_end_color) 
               end
               hold off
               subplot(2,2,3:4)
               hold on
               if strcmp(geometry{index}.type,'S')
                  loswall{1} = [geometry{index}.Ydata{1}(1) geometry{index}.Ydata{1}(3)];
                  loswall{2} = [geometry{index}.Ydata{1}(2) geometry{index}.Ydata{1}(4)];
                  plot(loswall{1},loswall{2},los_start_color) 
               elseif strcmp(geometry{index}.type,'E')
                  loswall{1} = [geometry{index}.Ydata{1}(1) geometry{index}.Ydata{1}(3)];
                  loswall{2} = [geometry{index}.Ydata{1}(2) geometry{index}.Ydata{1}(4)];       
                  plot(loswall{1},loswall{2},los_end_color) 
               elseif strcmp(geometry{index}.type,'T')
                  loswall{1} = [geometry{index}.Ydata{2}(1) geometry{index}.Ydata{2}(3)];
                  loswall{2} = [geometry{index}.Ydata{2}(2) geometry{index}.Ydata{2}(4)];     
                  plot(loswall{1},loswall{2},los_end_color) 
               end
               hold off
               
           else
               % Y
               circle_center_y = geometry{index}.circle_data(4);
               circle(:,1) = linspace(geometry{index}.Xdata{1}(1),geometry{index}.Xdata{2}(1),100);
               
                if circle_center_y > 0
                     circle(:,2)=circle_center_y-sqrt(inner_radius^2-(circle(:,1)-circle_center_z).^2);
                     circle(:,3)=circle_center_y-sqrt(outer_radius^2-(circle(:,1)-circle_center_z).^2);
                     if bender_mode
                        for jj=1:walls 
                            circle(:,3+jj) =circle_center_y-sqrt(wall_radius(jj)^2-(circle(:,1)-circle_center_z).^2);
                        end
                     end
                else
                     circle(:,2)=circle_center_y+sqrt(inner_radius^2-(circle(:,1)-circle_center_z).^2);
                     circle(:,3)=circle_center_y+sqrt(outer_radius^2-(circle(:,1)-circle_center_z).^2);    
                     if bender_mode
                        for jj=1:walls 
                            circle(:,3+jj) =circle_center_y+sqrt(wall_radius(jj)^2-(circle(:,1)-circle_center_z).^2);    
                        end
                     end
                end

               subplot(2,2,3:4) %Y
               hold on
               plot(circle(:,1),circle(:,2),'color',color_now,'linewidth',lw)
               plot(circle(:,1),circle(:,3),'color',color_now,'linewidth',lw)
               if bender_mode
                  for jj=1:walls 
                     eval(['plot(circle(:,1),circle(:,' num2str(3+jj) '),''color'',color_now,''linewidth'',lw)']);
                  end
               end
               hold off
                              
               % X part normal
               line1x{1}=[geometry{index}.Xdata{1}(1) geometry{index}.Xdata{2}(1)];
               line1x{2}=[geometry{index}.Xdata{1}(2) geometry{index}.Xdata{2}(2)];
               line2x{1}=[geometry{index}.Xdata{1}(3) geometry{index}.Xdata{2}(3)];
               line2x{2}=[geometry{index}.Xdata{1}(4) geometry{index}.Xdata{2}(4)];
               subplot(2,2,1:2) %X
               hold on
               plot(line1x{1},line1x{2},'color',color_now,'linewidth',lw)
               plot(line2x{1},line2x{2},'color',color_now,'linewidth',lw)
               hold off
               
               subplot(2,2,1:2)
               hold on
               if strcmp(geometry{index}.type,'S')
                  loswall{1} = [geometry{index}.Xdata{1}(1) geometry{index}.Xdata{1}(3)];
                  loswall{2} = [geometry{index}.Xdata{1}(2) geometry{index}.Xdata{1}(4)];
                  plot(loswall{1},loswall{2},los_start_color) 
               elseif strcmp(geometry{index}.type,'E')
                  loswall{1} = [geometry{index}.Xdata{1}(1) geometry{index}.Xdata{1}(3)];
                  loswall{2} = [geometry{index}.Xdata{1}(2) geometry{index}.Xdata{1}(4)];       
                  plot(loswall{1},loswall{2},los_end_color) 
               elseif strcmp(geometry{index}.type,'T')
                  loswall{1} = [geometry{index}.Xdata{2}(1) geometry{index}.Xdata{2}(3)];
                  loswall{2} = [geometry{index}.Xdata{2}(2) geometry{index}.Xdata{2}(4)];     
                  plot(loswall{1},loswall{2},los_end_color) 
               end
               hold off
               subplot(2,2,3:4)
               hold on
               if strcmp(geometry{index}.type,'S')
                  loswall{1} = [geometry{index}.Ydata{1}(1) geometry{index}.Ydata{1}(3)];
                  loswall{2} = [geometry{index}.Ydata{1}(2) geometry{index}.Ydata{1}(4)];
                  plot(loswall{1},loswall{2},los_start_color) 
               elseif strcmp(geometry{index}.type,'E')
                  loswall{1} = [geometry{index}.Ydata{1}(1) geometry{index}.Ydata{1}(3)];
                  loswall{2} = [geometry{index}.Ydata{1}(2) geometry{index}.Ydata{1}(4)];       
                  plot(loswall{1},loswall{2},los_end_color) 
               elseif strcmp(geometry{index}.type,'T')
                  loswall{1} = [geometry{index}.Ydata{2}(1) geometry{index}.Ydata{2}(3)];
                  loswall{2} = [geometry{index}.Ydata{2}(2) geometry{index}.Ydata{2}(4)];     
                  plot(loswall{1},loswall{2},los_end_color) 
               end
               hold off
               
           end
           
       case 'G'
           % Do nothing!
       case 'K'
           % Do nothing!
       case 'E'    
           % find original E
           test_index = index;
           test_logic = test_index > 0 && strcmp(geometry{test_index}.name,'E');
           while test_logic
               if test_index > 0 && strcmp(geometry{test_index}.name,'E')
                   if geometry{test_index}.Xpars(3) == geometry{index}.Xpars(3)
                       test_index = test_index - 1;
                       test_logic = 1;
                   else
                       test_logic = 0;
                   end
               else
                   test_logic = 0;
               end
           end
           orig_index = test_index + 1;
           % find last E
           test_index = index + 1;
           %while test_index <= module_num && strcmp(geometry{test_index}.name,'E') % && ~strcmp(geometry{test_index}.type,'N')
           %     test_index = test_index + 1;
           %end
           test_logic = test_index > 0 && strcmp(geometry{test_index}.name,'E');
           while test_logic
               if test_index > 0 && strcmp(geometry{test_index}.name,'E')
                   if geometry{test_index}.Xpars(3) == geometry{index}.Xpars(3)
                       test_index = test_index + 1;
                       test_logic = 1;
                   else
                       test_logic = 0;
                   end
               else
                   test_logic = 0;
               end
           end
           last_index = test_index - 1;
           
           Elength = geometry{orig_index}.Xpars(1);
           focus_s_v = - geometry{orig_index}.Ypars(3);
           focus_e_v = geometry{orig_index}.Ypars(1) + geometry{orig_index}.Ypars(4);
           smallaxis_v = geometry{orig_index}.Ypars(2);
           elength_v = focus_e_v - focus_s_v;
           focus_s_h = - geometry{orig_index}.Xpars(3);
           focus_e_h = geometry{orig_index}.Xpars(1) + geometry{orig_index}.Xpars(4);
           smallaxis_h = geometry{orig_index}.Xpars(2);
           elength_h = focus_e_h - focus_s_h;
                
            if index ~= last_index
                end_length = sqrt((geometry{orig_index}.Xpos(1) - geometry{index+1}.Xpos(1))^2 + (geometry{orig_index}.Xpos(2) - geometry{index+1}.Xpos(2))^2);
            else
                end_length = Elength;
            end
            if index ~= orig_index
                start_length = sqrt((geometry{index}.Xpos(1) - geometry{orig_index}.Xpos(1))^2 + (geometry{index}.Xpos(2) - geometry{orig_index}.Xpos(2))^2);
            else
                start_length = 0;
            end
                    
           xx = linspace(start_length,end_length,300);
           for ii = 1:length(xx);
                 y_height(ii) = smallaxis_v*sqrt(1-(((xx(ii)-focus_s_v)-elength_v/2)/(elength_v/2))*(((xx(ii)-focus_s_v)-elength_v/2)/(elength_v/2)));
                 x_width(ii)  = smallaxis_h*sqrt(1-(((xx(ii)-focus_s_h)-elength_h/2)/(elength_h/2))*(((xx(ii)-focus_s_h)-elength_h/2)/(elength_h/2)));
           end
           
           start_dir_x = geometry{orig_index}.Xdir;
           start_dir_y = geometry{orig_index}.Ydir;
           
           hat_dir_x(1) =  - start_dir_x(2);
           hat_dir_x(2) =  start_dir_x(1);
           hat_dir_y(1) =  - start_dir_y(2);
           hat_dir_y(2) =  start_dir_y(1);
           
           start_pos_x = geometry{orig_index}.Xpos;
           start_pos_y = geometry{orig_index}.Ypos;
           
           for ii=1:length(xx)
             upper_x_pos(ii,:) = start_pos_x + start_dir_x.*xx(ii) + hat_dir_x.*x_width(ii);
             lower_x_pos(ii,:) = start_pos_x + start_dir_x.*xx(ii) - hat_dir_x.*x_width(ii);
             upper_y_pos(ii,:) = start_pos_y + start_dir_y.*xx(ii) + hat_dir_y.*y_height(ii);
             lower_y_pos(ii,:) = start_pos_y + start_dir_y.*xx(ii) - hat_dir_y.*y_height(ii);
           end
           
           subplot(2,2,1:2)
           hold on
           plot(upper_x_pos(:,1),upper_x_pos(:,2),'-','color',color_now,'linewidth',lw)
           plot(lower_x_pos(:,1),lower_x_pos(:,2),'-','color',color_now,'linewidth',lw)
           hold off
           
           subplot(2,2,3:4)
           hold on
           plot(upper_y_pos(:,1),upper_y_pos(:,2),'-','color',color_now,'linewidth',lw)
           plot(lower_y_pos(:,1),lower_y_pos(:,2),'-','color',color_now,'linewidth',lw)
           hold off
           
               subplot(2,2,1:2)
               hold on
               if strcmp(geometry{index}.type,'S')
                  loswall{1} = [geometry{index}.Xdata{1}(1) geometry{index}.Xdata{1}(3)];
                  loswall{2} = [geometry{index}.Xdata{1}(2) geometry{index}.Xdata{1}(4)];
                  plot(loswall{1},loswall{2},los_start_color) 
               elseif strcmp(geometry{index}.type,'E')
                  loswall{1} = [geometry{index}.Xdata{1}(1) geometry{index}.Xdata{1}(3)];
                  loswall{2} = [geometry{index}.Xdata{1}(2) geometry{index}.Xdata{1}(4)];       
                  plot(loswall{1},loswall{2},los_end_color) 
               elseif strcmp(geometry{index}.type,'T')
                  loswall{1} = [geometry{index}.Xdata{2}(1) geometry{index}.Xdata{2}(3)];
                  loswall{2} = [geometry{index}.Xdata{2}(2) geometry{index}.Xdata{2}(4)];     
                  plot(loswall{1},loswall{2},los_end_color) 
               end
               hold off
               
               subplot(2,2,3:4)
               hold on
               if strcmp(geometry{index}.type,'S')
                  loswall{1} = [geometry{index}.Ydata{1}(1) geometry{index}.Ydata{1}(3)];
                  loswall{2} = [geometry{index}.Ydata{1}(2) geometry{index}.Ydata{1}(4)];
                  plot(loswall{1},loswall{2},los_start_color) 
               elseif strcmp(geometry{index}.type,'E')
                  loswall{1} = [geometry{index}.Ydata{1}(1) geometry{index}.Ydata{1}(3)];
                  loswall{2} = [geometry{index}.Ydata{1}(2) geometry{index}.Ydata{1}(4)];       
                  plot(loswall{1},loswall{2},los_end_color) 
               elseif strcmp(geometry{index}.type,'T')
                  loswall{1} = [geometry{index}.Ydata{2}(1) geometry{index}.Ydata{2}(3)];
                  loswall{2} = [geometry{index}.Ydata{2}(2) geometry{index}.Ydata{2}(4)];     
                  plot(loswall{1},loswall{2},los_end_color) 
               end
               hold off
           
           
       case 'P'
           
           % find original P
           %test_index = index;
           %while test_index > 0 && strcmp(geometry{test_index}.name,'P') %&& ~strcmp(geometry{test_index}.type,'N')
           %     test_index = test_index - 1;
           %end
           test_index = index;
           test_logic = test_index > 0 && strcmp(geometry{test_index}.name,'P');
           while test_logic
               if test_index > 0 && strcmp(geometry{test_index}.name,'P')
                   if geometry{test_index}.Xpars(3) == geometry{index}.Xpars(3)
                       test_index = test_index - 1;
                       test_logic = 1;
                   else
                       test_logic = 0;
                   end
               else
                  test_logic = 0;
               end
           end
           orig_index = test_index+1;
           % find last E
           test_index = index + 1;
           %while test_index <= module_num && strcmp(geometry{test_index}.name,'P') %&& ~strcmp(geometry{test_index}.type,'N')
           %     test_index = test_index + 1;
           %end
           test_logic = test_index > 0 && strcmp(geometry{test_index}.name,'P');
           while test_logic
               if test_index > 0 && strcmp(geometry{test_index}.name,'P')
                   if geometry{test_index}.Xpars(3) == geometry{index}.Xpars(3)
                       test_index = test_index + 1;
                       test_logic = 1;
                   else
                       test_logic = 0;
                   end
               else
                   test_logic = 0;
               end
           end
           last_index = test_index - 1;
           
           Elength = geometry{orig_index}.Xpars(1);
           focus_s_v = - geometry{orig_index}.Ypars(3);
           focus_e_v = geometry{orig_index}.Ypars(1) + geometry{orig_index}.Ypars(4);
           smallaxis_v = geometry{orig_index}.Ypars(2);
           elength_v = focus_e_v - focus_s_v;
           focus_s_h = - geometry{orig_index}.Xpars(3);
           focus_e_h = geometry{orig_index}.Xpars(1) + geometry{orig_index}.Xpars(4);
           smallaxis_h = geometry{orig_index}.Xpars(2);
           elength_h = focus_e_h - focus_s_h;
                
            if index ~= last_index
                end_length = sqrt((geometry{orig_index}.Xpos(1) - geometry{index+1}.Xpos(1))^2 + (geometry{orig_index}.Xpos(2) - geometry{index+1}.Xpos(2))^2);
            else
                end_length = Elength;
            end
            if index ~= orig_index
                start_length = sqrt((geometry{index}.Xpos(1) - geometry{orig_index}.Xpos(1))^2 + (geometry{index}.Xpos(2) - geometry{orig_index}.Xpos(2))^2);
            else
                start_length = 0;
            end
                    
           xx = linspace(start_length,end_length,300);
           for ii = 1:length(xx);
                 y_height(ii) = smallaxis_v*sqrt(1-(((xx(ii)-focus_s_v)-elength_v/2)/(elength_v/2))*(((xx(ii)-focus_s_v)-elength_v/2)/(elength_v/2)));
                 x_width(ii)  = smallaxis_h*sqrt(1-(((xx(ii)-focus_s_h)-elength_h/2)/(elength_h/2))*(((xx(ii)-focus_s_h)-elength_h/2)/(elength_h/2)));
           end
           
           start_dir_x = geometry{orig_index}.Xdir;
           start_dir_y = geometry{orig_index}.Ydir;
           
           hat_dir_x(1) =  - start_dir_x(2);
           hat_dir_x(2) =  start_dir_x(1);
           hat_dir_y(1) =  - start_dir_y(2);
           hat_dir_y(2) =  start_dir_y(1);
           
           start_pos_x = geometry{orig_index}.Xpos;
           start_pos_y = geometry{orig_index}.Ypos;
           
           for ii=1:length(xx)
             upper_x_pos(ii,:) = start_pos_x + start_dir_x.*xx(ii) + hat_dir_x.*x_width(ii);
             lower_x_pos(ii,:) = start_pos_x + start_dir_x.*xx(ii) - hat_dir_x.*x_width(ii);
             upper_y_pos(ii,:) = start_pos_y + start_dir_y.*xx(ii) + hat_dir_y.*y_height(ii);
             lower_y_pos(ii,:) = start_pos_y + start_dir_y.*xx(ii) - hat_dir_y.*y_height(ii);
           end
           
           subplot(2,2,1:2)
           hold on
           plot(upper_x_pos(:,1),upper_x_pos(:,2),'-','color',color_now,'linewidth',lw)
           plot(lower_x_pos(:,1),lower_x_pos(:,2),'-','color',color_now,'linewidth',lw)
           hold off
           
           subplot(2,2,3:4)
           hold on
           plot(upper_y_pos(:,1),upper_y_pos(:,2),'-','color',color_now,'linewidth',lw)
           plot(lower_y_pos(:,1),lower_y_pos(:,2),'-','color',color_now,'linewidth',lw)
           hold off
           
               subplot(2,2,1:2)
               hold on
               if strcmp(geometry{index}.type,'S')
                  loswall{1} = [geometry{index}.Xdata{1}(1) geometry{index}.Xdata{1}(3)];
                  loswall{2} = [geometry{index}.Xdata{1}(2) geometry{index}.Xdata{1}(4)];
                  plot(loswall{1},loswall{2},los_start_color) 
               elseif strcmp(geometry{index}.type,'E')
                  loswall{1} = [geometry{index}.Xdata{1}(1) geometry{index}.Xdata{1}(3)];
                  loswall{2} = [geometry{index}.Xdata{1}(2) geometry{index}.Xdata{1}(4)];       
                  plot(loswall{1},loswall{2},los_end_color) 
               elseif strcmp(geometry{index}.type,'T')
                  loswall{1} = [geometry{index}.Xdata{2}(1) geometry{index}.Xdata{2}(3)];
                  loswall{2} = [geometry{index}.Xdata{2}(2) geometry{index}.Xdata{2}(4)];     
                  plot(loswall{1},loswall{2},los_end_color) 
               end
               hold off
               
               subplot(2,2,3:4)
               hold on
               if strcmp(geometry{index}.type,'S')
                  loswall{1} = [geometry{index}.Ydata{1}(1) geometry{index}.Ydata{1}(3)];
                  loswall{2} = [geometry{index}.Ydata{1}(2) geometry{index}.Ydata{1}(4)];
                  plot(loswall{1},loswall{2},los_start_color) 
               elseif strcmp(geometry{index}.type,'E')
                  loswall{1} = [geometry{index}.Ydata{1}(1) geometry{index}.Ydata{1}(3)];
                  loswall{2} = [geometry{index}.Ydata{1}(2) geometry{index}.Ydata{1}(4)];       
                  plot(loswall{1},loswall{2},los_end_color) 
               elseif strcmp(geometry{index}.type,'T')
                  loswall{1} = [geometry{index}.Ydata{2}(1) geometry{index}.Ydata{2}(3)];
                  loswall{2} = [geometry{index}.Ydata{2}(2) geometry{index}.Ydata{2}(4)];     
                  plot(loswall{1},loswall{2},los_end_color) 
               end
               hold off
               
               
       case 'Selene'
          
           slit1_x_down=[geometry{index}.Xdata{1}(1) geometry{index}.Xdata{1}(2)];
           slit1_x_up=[geometry{index}.Xdata{1}(3) geometry{index}.Xdata{1}(4)];
           delta_slit_x = slit1_x_up - slit1_x_down;
           delta_slit_x = delta_slit_x./sqrt(delta_slit_x(1)^2+delta_slit_x(2)^2);
           slit1_x_up_outside = slit1_x_up + 0.1 * delta_slit_x;
           slit1_x_down_outside = slit1_x_up - 0.1 * delta_slit_x;
           subplot(2,2,1:2) %X
           hold on
           plot([slit1_x_up(1) slit1_x_up_outside(1)],[slit1_x_up(2) slit1_x_up_outside(2)],'color',color_now,'linewidth',1)
           plot([slit1_x_down(1) slit1_x_down_outside(1)],[slit1_x_down(2) slit1_x_down_outside(2)],'color',color_now,'linewidth',1)
           hold off

           slit1_y_down=[geometry{index}.Ydata{1}(1) geometry{index}.Ydata{1}(2)];
           slit1_y_up=[geometry{index}.Ydata{1}(3) geometry{index}.Ydata{1}(4)];
           delta_slit_y = slit1_y_up - slit1_y_down;
           delta_slit_y = delta_slit_y./sqrt(delta_slit_y(1)^2+delta_slit_y(2)^2);
           slit1_y_up_outside = slit1_y_up + 0.1 * delta_slit_y;
           slit1_y_down_outside = slit1_y_up - 0.1 * delta_slit_y;
           subplot(2,2,3:4) %X
           hold on
           plot([slit1_y_up(1) slit1_y_up_outside(1)],[slit1_y_up(2) slit1_y_up_outside(2)],'color',color_now,'linewidth',1)
           plot([slit1_y_down(1) slit1_y_down_outside(1)],[slit1_y_down(2) slit1_y_down_outside(2)],'color',color_now,'linewidth',1)
           hold off

           
           selene_distance = geometry{index}.Pars(1); % distance between focal point and actual guide
           selene_length = geometry{index}.Pars(2); % length of one physical guide element
           selene_c = geometry{index}.Pars(3); % quater of the total distance between start and end focal points
           
           focus_s_h = - selene_distance;
           focus_e_h = selene_distance + selene_length;
           smallaxis_h = geometry{index}.Bxy(1);
           elength_h = focus_e_h - focus_s_h;
           
           Xpos1 = geometry{index}.Xpos(1:2); %z x mcstas, x y here.
           Xpos2 = geometry{index}.Xpos(3:4);
           Xdir = (Xpos2 - Xpos1);
           Xdir = Xdir./sqrt(Xdir(1)^2 + Xdir(2)^2);
           
           hat_Xdir(1) =  - Xdir(2);
           hat_Xdir(2) =  Xdir(1);
           
           focus_s_v = - selene_distance;
           focus_e_v = selene_distance + selene_length;
           smallaxis_v = geometry{index}.Bxy(2);
           elength_v = focus_e_v - focus_s_v;
           
           Ypos1 = geometry{index}.Ypos(1:2); %z x mcstas, x y here.
           Ypos2 = geometry{index}.Ypos(3:4);
           Ydir = (Ypos2 - Ypos1);
           Ydir = Ydir./sqrt(Ydir(1)^2 + Ydir(2)^2);
           
           hat_Ydir(1) =  - Ydir(2);
           hat_Ydir(2) =  Ydir(1);
           
           xx = linspace(0,selene_length,200);
           x_width=zeros(200,1);y_height=zeros(200,1);
           for ii = 1:length(xx)
             x_width(ii)  = smallaxis_h*sqrt(1-(((xx(ii)-focus_s_h)-elength_h/2)/(elength_h/2))*(((xx(ii)-focus_s_h)-elength_h/2)/(elength_h/2)));
             y_height(ii)  = smallaxis_v*sqrt(1-(((xx(ii)-focus_s_v)-elength_v/2)/(elength_v/2))*(((xx(ii)-focus_s_v)-elength_v/2)/(elength_v/2)));
           end
             
             ellipse1_x_pos=zeros(200,2);ellipse2_x_pos=zeros(200,2);ellipse1_y_pos=zeros(200,2);ellipse2_y_pos=zeros(200,2);
           for ii=1:length(xx)
             ellipse1_x_pos(ii,:) = Xpos1 + Xdir.*(xx(ii)+selene_distance) - hat_Xdir.*x_width(ii);
             ellipse2_x_pos(ii,:) = Xpos1 + Xdir.*(xx(ii)+selene_distance+2*selene_c) + hat_Xdir.*x_width(ii);
             ellipse1_y_pos(ii,:) = Ypos1 + Ydir.*(xx(ii)+selene_distance) + hat_Xdir.*y_height(ii);
             ellipse2_y_pos(ii,:) = Ypos1 + Ydir.*(xx(ii)+selene_distance+2*selene_c) - hat_Xdir.*y_height(ii);
           end
           
           subplot(2,2,1:2)
           hold on
           plot(ellipse1_x_pos(:,1),ellipse1_x_pos(:,2),'-','color',color_now,'linewidth',lw)
           plot(ellipse2_x_pos(:,1),ellipse2_x_pos(:,2),'-','color',color_now,'linewidth',lw)
           hold off
           subplot(2,2,3:4)
           hold on
           plot(ellipse1_y_pos(:,1),ellipse1_y_pos(:,2),'-','color',color_now,'linewidth',lw)
           plot(ellipse2_y_pos(:,1),ellipse2_y_pos(:,2),'-','color',color_now,'linewidth',lw)
           hold off
           
       case 'Slit'

           slit1_x_down=[geometry{index}.Xdata{1}(1) geometry{index}.Xdata{1}(2)];
           slit1_x_up=[geometry{index}.Xdata{1}(3) geometry{index}.Xdata{1}(4)];
           delta_slit_x = slit1_x_up - slit1_x_down;
           delta_slit_x = delta_slit_x./sqrt(delta_slit_x(1)^2+delta_slit_x(2)^2);
           slit1_x_up_outside = slit1_x_up + 0.1 * delta_slit_x;
           slit1_x_down_outside = slit1_x_up - 0.1 * delta_slit_x;
           subplot(2,2,1:2) %X
           hold on
           plot([slit1_x_up(1) slit1_x_up_outside(1)],[slit1_x_up(2) slit1_x_up_outside(2)],'color',color_now,'linewidth',1)
           plot([slit1_x_down(1) slit1_x_down_outside(1)],[slit1_x_down(2) slit1_x_down_outside(2)],'color',color_now,'linewidth',1)
           hold off

           slit1_y_down=[geometry{index}.Ydata{1}(1) geometry{index}.Ydata{1}(2)];
           slit1_y_up=[geometry{index}.Ydata{1}(3) geometry{index}.Ydata{1}(4)];
           delta_slit_y = slit1_y_up - slit1_y_down;
           delta_slit_y = delta_slit_y./sqrt(delta_slit_y(1)^2+delta_slit_y(2)^2);
           slit1_y_up_outside = slit1_y_up + 0.1 * delta_slit_y;
           slit1_y_down_outside = slit1_y_up - 0.1 * delta_slit_y;
           subplot(2,2,3:4) %X
           hold on
           plot([slit1_y_up(1) slit1_y_up_outside(1)],[slit1_y_up(2) slit1_y_up_outside(2)],'color',color_now,'linewidth',1)
           plot([slit1_y_down(1) slit1_y_down_outside(1)],[slit1_y_down(2) slit1_y_down_outside(2)],'color',color_now,'linewidth',1)
           hold off
           
       case 'Sample'
           % Draw the sample
           sample_pos_x = geometry{index}.start_pos_x + geometry{index}.data(3).*geometry{index}.start_div_x;
           sample_pos_y = geometry{index}.start_pos_y + geometry{index}.data(3).*geometry{index}.start_div_y;
           
           hat_dir_x(1) =  - geometry{index}.start_div_x(2);
           hat_dir_x(2) =  geometry{index}.start_div_x(1);
           hat_dir_y(1) =  - geometry{index}.start_div_y(2);
           hat_dir_y(2) =  geometry{index}.start_div_y(1);
           
           sample_up_x   = sample_pos_x + 0.5*geometry{index}.data(1).*hat_dir_x;
           sample_down_x = sample_pos_x - 0.5*geometry{index}.data(1).*hat_dir_x;
           
           sample_up_y   = sample_pos_y + 0.5*geometry{index}.data(2).*hat_dir_y;
           sample_down_y = sample_pos_y - 0.5*geometry{index}.data(2).*hat_dir_y;
           
           subplot(2,2,1:2)
           hold on
           plot([sample_up_x(1) sample_down_x(1)],[sample_up_x(2) sample_down_x(2)],'linewidth',lw)
           hold off
           subplot(2,2,3:4)
           hold on
           plot([sample_up_y(1) sample_down_y(1)],[sample_up_y(2) sample_down_y(2)],'linewidth',lw)
           hold off
           
       case 'Moderator'
           
           Z = [0 0];
           X = 0.5.*[geometry{index}.data(1) -geometry{index}.data(1)];
           Y = 0.5.*[geometry{index}.data(2) -geometry{index}.data(2)];
           
           moderator_x = X(1);
           moderator_y = Y(1);
           
           subplot(2,2,1:2)
           hold on
           plot(Z,X,'k','linewidth',lw)
           hold off
           subplot(2,2,3:4)
           hold on
           plot(Z,Y,'k','linewidth',lw)
           hold off
           
       otherwise
           disp('ERROR, no case for this case! Visualizer error')
               
   end
    
end

    figure(1)
    subplot(2,2,1:2)
    old_xlim = get(gca,'XLim');
    new_xlim = old_xlim;
    new_xlim(1) = -1.5;
    new_xlim(2) = old_xlim(2)*1.025;
    set(gca,'XLim',new_xlim)
      
    old_ylim = get(gca,'YLim');
    new_ylim = old_ylim;
    if old_ylim(1) > -1.15*moderator_x
        new_ylim(1) = -1.15*moderator_x;
    end
    if old_ylim(2) < 1.15*moderator_x
        new_ylim(2) = 1.15*moderator_x;
    end
    set(gca,'YLim',new_ylim)
    subplot(2,2,3:4)
    old_xlim = get(gca,'XLim');
    new_xlim = old_xlim;
    new_xlim(1) = -1.5;
    new_xlim(2) = old_xlim(2)*1.025;
    set(gca,'XLim',new_xlim)
    old_ylim = get(gca,'YLim');
    new_ylim = old_ylim;
    if old_ylim(1) > -1.15*moderator_y
        new_ylim(1) = -1.15*moderator_y;
    end
    if old_ylim(2) < 1.15*moderator_y
        new_ylim(2) = 1.15*moderator_y;
    end
    set(gca,'YLim',new_ylim)
    
    subplot(2,2,1:2)
    set(gca,'fontsize',20)
    box
    xlabel('[m]','FontSize',fsize)
    ylabel('[m]','FontSize',fsize)
    title('Horizontal plane','FontSize',fsize+4)
    box on
    
    subplot(2,2,3:4)
    set(gca,'fontsize',20)
    box
    xlabel('[m]','FontSize',fsize)
    ylabel('[m]','FontSize',fsize)
    title('Vertical plane','FontSize',fsize+4)
    box on
    
    print(overview,'-dpng','-r300',[file_name '_geometry.png'])
    
    disp('fsize');disp(fsize);
    
else
    disp('ERROR, visualizer could not find the specified file')
end
