function [McStasStr] = point_calc_writer(McStasStr,index,last,globalinfo,ghost_globalinfo,defaults,enable_raytracer);
% Function which writes any needed raytracing algorithms to the mcstas
% instrument file.

num = num2str(index);
numM1 = num2str(index-1);
numP1 = num2str(index+1);


if globalinfo.old_point_calc
    
if (isfield(globalinfo,'kinks') && globalinfo.kinkoverwrite(globalinfo.kinks(1))==0 && globalinfo.raytracereq==1);   
    %allocate variables needed
    if sum(ismember(McStasStr.declare,['startxpoint[' num2str(last+1) '][3][3]']))<0.5 % not always needed
                McStasStr.declare{end+1}=['startxpoint[' num2str(last+1) '][3][3]'];
    end
    if sum(ismember(McStasStr.declare,['startypoint[' num2str(last+1) '][3][3]']))<0.5 % not always needed
                McStasStr.declare{end+1}=['startypoint[' num2str(last+1) '][3][3]'];
    end
    if globalinfo.rotlogic(globalinfo.kinks(1))>0; % horizontal kink
        if sum(ismember(McStasStr.declare,'kinkendpointx[3]'))<0.5 % not always needed
                McStasStr.declare{end+1}='kinkendpointx[3]';
        end
    else
        if sum(ismember(McStasStr.declare,'kinkendpointy[3]'))<0.5 % not always needed
                McStasStr.declare{end+1}='kinkendpointy[3]';
        end
    end

    if index > globalinfo.kinks
    % before kink (easier!)
    l{1}='';
    l{end+1} =['startxpoint[' num '][1][1] = endPoint' num ' - length' num ';'];
    l{end+1} =['startxpoint[' num '][2][1] = -0.5*startx' num ';'];
    l{end+1} =['startxpoint[' num '][1][2] = endPoint' num ' - length' num ';'];
    l{end+1} =['startxpoint[' num '][2][2] = 0.5*startx' num ';'];

    l{end+1} =['startypoint[' num '][1][1] = endPoint' num ' - length' num ';'];
    l{end+1} =['startypoint[' num '][2][1] = -0.5*starty' num ';'];
    l{end+1} =['startypoint[' num '][1][2] = endPoint' num ' - length' num ';'];
    l{end+1} =['startypoint[' num '][2][2] = 0.5*starty' num ';'];

    elseif index==globalinfo.kinks(1) % need to be generalized to have more kinks
    % at kink (may need to be done in the kink module instead)
    l{1}='';
    l{end+1} =['startxpoint[' num '][1][1] = endPoint' num ' - length' num ';'];
    l{end+1} =['startxpoint[' num '][2][1] = -0.5*startx' num ';'];
    l{end+1} =['startxpoint[' num '][1][2] = endPoint' num ' - length' num ';'];
    l{end+1} =['startxpoint[' num '][2][2] = 0.5*startx' num ';'];

    l{end+1} =['startypoint[' num '][1][1] = endPoint' num ' - length' num ';'];
    l{end+1} =['startypoint[' num '][2][1] = -0.5*starty' num ';'];
    l{end+1} =['startypoint[' num '][1][2] = endPoint' num ' - length' num ';'];
    l{end+1} =['startypoint[' num '][2][2] = 0.5*starty' num ';'];
    l{end+1}='';
    
    
    if globalinfo.rotlogic(globalinfo.kinks(1))>0;
    l{end+1} =['    kinkendpointx[1] = endPoint' num ';'];
    if globalinfo.rotsign(index)>0;
    l{end+1} =['        kinkendpointx[2] = translation' num ';'];
    else
    l{end+1} =['        kinkendpointx[2] = -translation' num ';'];
    end
    %l{end+1} =['    kinkendpointy[1] = endPoint' num ';'];
    %l{end+1} =['    kinkendpointy[2] = 0;'];
    else
    l{end+1} =['    kinkendpointy[1] = endPoint' num ';'];
    if globalinfo.rotsign(index)>0;
    l{end+1} =['        kinkendpointy[2] = translation' num ';'];
    else
    l{end+1} =['        kinkendpointy[2] = -translation' num ';'];
    end
    %l{end+1} =['    kinkendpointx[1] = endPoint' num ';'];
    %l{end+1} =['    kinkendpointx[2] = 0;'];
    end
    
    else
    % after kink  
    % doh, I know before hand if the needed code is for horizontal or
    % vertical, there is no need to check for that in the McStas code!
    kinknumstr=num2str(globalinfo.kinks(1));
    l{1}='';
    if globalinfo.rotlogic(globalinfo.kinks(1))>0
    l{end+1} =['startxpoint[' num '][1][1] = kinkendpointx[1] + (endPoint' num ' - length' num ' - endPoint' kinknumstr ')*cos(rot' kinknumstr '*DEG2RAD)+0.5*sin(rot' kinknumstr '*DEG2RAD)*startx' num ';'];
    l{end+1} =['startxpoint[' num '][2][1] = kinkendpointx[2] + (endPoint' num ' - length' num ' - endPoint' kinknumstr ')*sin(rot' kinknumstr '*DEG2RAD)-0.5*cos(rot' kinknumstr '*DEG2RAD)*startx' num ';'];
    l{end+1} =['startxpoint[' num '][1][2] = kinkendpointx[1] + (endPoint' num ' - length' num ' - endPoint' kinknumstr ')*cos(rot' kinknumstr '*DEG2RAD)-0.5*sin(rot' kinknumstr '*DEG2RAD)*startx' num ';'];
    l{end+1} =['startxpoint[' num '][2][2] = kinkendpointx[2] + (endPoint' num ' - length' num ' - endPoint' kinknumstr ')*sin(rot' kinknumstr '*DEG2RAD)+0.5*cos(rot' kinknumstr '*DEG2RAD)*startx' num ';'];

    l{end+1} =['startypoint[' num '][1][1] = endPoint' num ' - length' num ';'];
    l{end+1} =['startypoint[' num '][2][1] = -0.5*starty' num ';'];
    l{end+1} =['startypoint[' num '][1][2] = endPoint' num ' - length' num ';'];
    l{end+1} =['startypoint[' num '][2][2] = 0.5*starty' num ';'];
    else
    l{end+1} =['startxpoint[' num '][1][1] = endPoint' num ' - length' num ';'];
    l{end+1} =['startxpoint[' num '][2][1] = -0.5*startx' num ';'];
    l{end+1} =['startxpoint[' num '][1][2] = endPoint' num ' - length' num ';'];
    l{end+1} =['startxpoint[' num '][2][2] = 0.5*startx' num ';'];

    l{end+1} =['startypoint[' num '][1][1] = kinkendpointy[1] + (endPoint' num ' - length' num ' - endPoint' kinknumstr ')*cos(rot' kinknumstr '*DEG2RAD)+0.5*sin(rot' kinknumstr '*DEG2RAD)*starty' num ';'];
    l{end+1} =['startypoint[' num '][2][1] = kinkendpointy[2] + (endPoint' num ' - length' num ' - endPoint' kinknumstr ')*sin(rot' kinknumstr '*DEG2RAD)-0.5*cos(rot' kinknumstr '*DEG2RAD)*starty' num ';'];
    l{end+1} =['startypoint[' num '][1][1] = kinkendpointy[1] + (endPoint' num ' - length' num ' - endPoint' kinknumstr ')*cos(rot' kinknumstr '*DEG2RAD)-0.5*sin(rot' kinknumstr '*DEG2RAD)*starty' num ';'];
    l{end+1} =['startypoint[' num '][2][2] = kinkendpointy[2] + (endPoint' num ' - length' num ' - endPoint' kinknumstr ')*sin(rot' kinknumstr '*DEG2RAD)+0.5*cos(rot' kinknumstr '*DEG2RAD)*starty' num ';'];
    end
    end
    
    if index==1 % end of the guide
%         if sum(ismember(McStasStr.declare,'endpointx[3][3]'))<0.5 % not always needed
%                     McStasStr.declare{end+1}='endpointx[3][3]';
%         end
%         if sum(ismember(McStasStr.declare,'endpointy[3][3]'))<0.5 % not always needed
%                     McStasStr.declare{end+1}='endpointy[3][3]';
%         end        
        % Always after the kink
        kinknumstr=num2str(globalinfo.kinks(1));
        
        if globalinfo.rotlogic(globalinfo.kinks(1))>0
        l{end+1} =['startxpoint[0][1][1] = kinkendpointx[1] + (endPoint' num ' - endPoint' kinknumstr ')*cos(rot' kinknumstr '*DEG2RAD)+0.5*sin(rot' kinknumstr '*DEG2RAD)*endx' num ';'];
        l{end+1} =['startxpoint[0][2][1] = kinkendpointx[2] + (endPoint' num ' - endPoint' kinknumstr ')*sin(rot' kinknumstr '*DEG2RAD)-0.5*cos(rot' kinknumstr '*DEG2RAD)*endx' num ';'];
        l{end+1} =['startxpoint[0][1][2] = kinkendpointx[1] + (endPoint' num ' - endPoint' kinknumstr ')*cos(rot' kinknumstr '*DEG2RAD)-0.5*sin(rot' kinknumstr '*DEG2RAD)*endx' num ';'];
        l{end+1} =['startxpoint[0][2][2] = kinkendpointx[2] + (endPoint' num ' - endPoint' kinknumstr ')*sin(rot' kinknumstr '*DEG2RAD)+0.5*cos(rot' kinknumstr '*DEG2RAD)*endx' num ';'];

        l{end+1} =['startypoint[0][1][1] = endPoint' num ';'];
        l{end+1} =['startypoint[0][2][1] = -0.5*endy' num ';'];
        l{end+1} =['startypoint[0][1][2] = endPoint' num ';'];
        l{end+1} =['startypoint[0][2][2] = 0.5*endy' num ';'];
        else
        l{end+1} =['startxpoint[0][1][1] = endPoint' num ';'];
        l{end+1} =['startxpoint[0][2][1] = -0.5*endx' num ';'];
        l{end+1} =['startxpoint[0][1][2] = endPoint' num ';'];
        l{end+1} =['startxpoint[0][2][2] = 0.5*endx' num ';'];

        l{end+1} =['startypoint[0][1][1] = kinkendpointy[1] + (endPoint' num ' - endPoint' kinknumstr ')*cos(rot' kinknumstr '*DEG2RAD)+0.5*sin(rot' kinknumstr '*DEG2RAD)*endy' num ';'];
        l{end+1} =['startypoint[0][2][1] = kinkendpointy[2] + (endPoint' num ' - endPoint' kinknumstr ')*sin(rot' kinknumstr '*DEG2RAD)-0.5*cos(rot' kinknumstr '*DEG2RAD)*endy' num ';'];
        l{end+1} =['startypoint[0][1][2] = kinkendpointy[1] + (endPoint' num ' - endPoint' kinknumstr ')*cos(rot' kinknumstr '*DEG2RAD)-0.5*sin(rot' kinknumstr '*DEG2RAD)*endy' num ';'];
        l{end+1} =['startypoint[0][2][2] = kinkendpointy[2] + (endPoint' num ' - endPoint' kinknumstr ')*sin(rot' kinknumstr '*DEG2RAD)+0.5*cos(rot' kinknumstr '*DEG2RAD)*endy' num ';'];
        end
    end

    if index==1 && isfield(globalinfo,'kinks') && globalinfo.kinkoverwrite(globalinfo.kinks(1))==0;   
    % write if los here, completely in McStas
    if sum(ismember(McStasStr.declareint,'n_check'))<0.5 % not always needed
             McStasStr.declareint{end+1}='n_check';
             McStasStr.declareint{end+1}='los_logic';
             McStasStr.declareint{end+1}='n1';
             McStasStr.declareint{end+1}='n2';
             McStasStr.declareint{end+1}='line';
             McStasStr.declareint{end+1}='los_tmp[5]';
             McStasStr.declareint{end+1}='n_check';
             McStasStr.declareint{end+1}='los_logic';
             McStasStr.declareint{end+1}='los_check';
             McStasStr.declareint{end+1}=['los_logic_single[' kinknumstr '][' num2str(last+1) ']'];
    end
    if globalinfo.rotlogic(globalinfo.kinks(1))>0
        if sum(ismember(McStasStr.declare,'X1[5]'))<0.5 % not always needed
                 McStasStr.declare{end+1}='X1[5]';
                 McStasStr.declare{end+1}='X2[5]';
                 McStasStr.declare{end+1}='Z1[5]';
                 McStasStr.declare{end+1}='Z2[5]';
                 McStasStr.declare{end+1}='a[5]';
                 McStasStr.declare{end+1}='b[5]';
                 McStasStr.declare{end+1}='tmp_double';
        end
    else
        if sum(ismember(McStasStr.declare,'Y1[5]'))<0.5 % not always needed
                 McStasStr.declare{end+1}='Y1[5]';
                 McStasStr.declare{end+1}='Y2[5]';
                 McStasStr.declare{end+1}='Z1[5]';
                 McStasStr.declare{end+1}='Z2[5]';
                 McStasStr.declare{end+1}='a[5]';
                 McStasStr.declare{end+1}='b[5]';
        end
    end
    
    % Temp code, if this works, the endpointx and y are not needed
%     l{end+1}='';
%     l{end+1} =['startxpoint[0][1][1] = endpointx[1][1];'];
%     l{end+1} =['startxpoint[0][2][1] = endpointx[2][1];'];
%     l{end+1} =['startxpoint[0][1][2] = endpointx[1][2];'];
%     l{end+1} =['startxpoint[0][2][2] = endpointx[2][2];'];
%     l{end+1}='';
%     l{end+1} =['startypoint[0][1][1] = endpointy[1][1];'];
%     l{end+1} =['startypoint[0][2][1] = endpointy[2][1];'];
%     l{end+1} =['startypoint[0][1][1] = endpointy[1][2];'];
%     l{end+1} =['startypoint[0][2][2] = endpointy[2][2];'];
%     l{end+1}='';
             
    if globalinfo.rotlogic(globalinfo.kinks(1))>0         
    kinknumstr=num2str(globalinfo.kinks(1));
    l{end+1} =['los_logic = 1; // assume line of sight'];
    l{end+1} =['for (n1=0;n1<' kinknumstr ';++n1) {'];
    l{end+1} =['    for (n2=' kinknumstr ';n2<' num2str(last+1) ';++n2) {'];
    l{end+1} =['        Z1[1]=startxpoint[n1][1][1];'];
    l{end+1} =['        X1[1]=startxpoint[n1][2][1];'];
    l{end+1} =['        Z2[1]=startxpoint[n2][1][1];'];
    l{end+1} =['        X2[1]=startxpoint[n2][2][1];'];
    l{end+1} =['        Z1[2]=startxpoint[n1][1][2];'];
    l{end+1} =['        X1[2]=startxpoint[n1][2][2];'];
    l{end+1} =['        Z2[2]=startxpoint[n2][1][2];'];
    l{end+1} =['        X2[2]=startxpoint[n2][2][2];'];
    l{end+1} =['        Z1[3]=startxpoint[n1][1][1];'];
    l{end+1} =['        X1[3]=startxpoint[n1][2][1];'];
    l{end+1} =['        Z2[3]=startxpoint[n2][1][2];'];
    l{end+1} =['        X2[3]=startxpoint[n2][2][2];'];
    l{end+1} =['        Z1[4]=startxpoint[n1][1][2];'];
    l{end+1} =['        X1[4]=startxpoint[n1][2][2];'];
    l{end+1} =['        Z2[4]=startxpoint[n2][1][1];'];
    l{end+1} =['        X2[4]=startxpoint[n2][2][1];'];
    l{end+1} =['        for (line=1;line<5;++line) { // 4 lines need to be checked'];
    l{end+1} =['            tmp_double=(Z1[line]-Z2[line])/100;'];
    %l{end+1} =['            a[l]=(X1[l]-X2[l])/(Z1[l]-Z2[l]);'];
    l{end+1} =['            a[line]=(X1[line]-X2[line])/tmp_double;'];
    l{end+1} =['            tmp_double=a[line]*Z1[line];'];
    %l{end+1} =['            b[l]=X1[l]-a[l]*Z1[l];'];
    l{end+1} =['            b[line]=X1[line]-tmp_double/100;'];
    l{end+1} =['            los_tmp[line]=1;'];
    l{end+1} =['            for (n_check=0;n_check<' num2str(last+1) ';++n_check) {'];
    l{end+1} =['                // check if the line goes between the points'];
    l{end+1} =['                if (n_check != n1 && n_check != n2 && ((a[line]*startxpoint[n_check][1][1])/100+b[line]<startxpoint[n_check][2][1] || (a[line]*startxpoint[n_check][1][2])/100+b[line]>startxpoint[n_check][2][2])) {'];
    l{end+1} =['                    los_tmp[line]=0; // line of sight blocked for line l!'];
    l{end+1} =['                }'];
    l{end+1} =['            }'];
    l{end+1} =['        }'];
    l{end+1} =['        if (los_tmp[1]+los_tmp[2]+los_tmp[3]+los_tmp[4]==0) {'];
    l{end+1} =['            los_logic_single[n1][n2]=0;'];
    %l{end+1} =['            printf("los_logic_single[%%d][%%d]=0.\\n",n1,n2);'];
    l{end+1} =['        }'];
    l{end+1} =['        else {'];
    l{end+1} =['            los_logic_single[n1][n2]=1;'];
    %l{end+1} =['            printf("los_logic_single[%%d][%%d]=1.\\n",n1,n2);'];
    %l{end+1} =['            printf("kink rot = %%lf.\\n",rot' kinknumstr ');'];
    %l{end+1} =['            if (los_tmp[1]==1) printf("line=1 had los.\\n");'];
    %l{end+1} =['            if (los_tmp[2]==1) printf("line=2 had los.\\n");'];
    %l{end+1} =['            if (los_tmp[3]==1) printf("line=3 had los.\\n");'];
    %l{end+1} =['            if (los_tmp[4]==1) printf("line=4 had los.\\n");'];
    %l{end+1} =['            printf("\\n");'];
    l{end+1} =['        }'];
    l{end+1} =['    }'];
    l{end+1} =['}'];
    l{end+1} =['los_check=0;'];
    l{end+1} =['for (n1=0;n1<' kinknumstr ';++n1) {'];
    l{end+1} =['    for (n2=' kinknumstr ';n2<' num2str(last+1) ';++n2) {'];
    l{end+1} =['        los_check=los_check+los_logic_single[n1][n2];'];
    l{end+1} =['    }'];
    l{end+1} =['}'];
    l{end+1} =['if (los_check==0) {'];
    l{end+1} =['    los_logic = 0;'];
    l{end+1} =['    printf("kink rot = %%lf.\\n\\n",rot' kinknumstr ');'];
    l{end+1} =['}'];
    l{end+1} =['}']; % to end the while loop
    
    %l{end+1}=['printf("los_logic=%%d \\n",los_logic);'];
    else
    kinknumstr=num2str(globalinfo.kinks(1));
    l{end+1} =['los_logic = 1; // assume line of sight'];
    l{end+1} =['for (n1=0;n1<' kinknumstr ';++n1) {'];
    l{end+1} =['    for (n2=' kinknumstr ';n2<' num2str(last+1) ';++n2) {'];
    l{end+1} =['        Z1[1]=startypoint[n1][1][1];'];
    l{end+1} =['        Y1[1]=startypoint[n1][2][1];'];
    l{end+1} =['        Z2[1]=startypoint[n2][1][1];'];
    l{end+1} =['        Y2[1]=startypoint[n2][2][1];'];
    l{end+1} =['        Z1[2]=startypoint[n1][1][2];'];
    l{end+1} =['        Y1[2]=startypoint[n1][2][2];'];
    l{end+1} =['        Z2[2]=startypoint[n2][1][2];'];
    l{end+1} =['        Y2[2]=startypoint[n2][2][2];'];
    l{end+1} =['        Z1[3]=startypoint[n1][1][1];'];
    l{end+1} =['        Y1[3]=startypoint[n1][2][1];'];
    l{end+1} =['        Z2[3]=startypoint[n2][1][2];'];
    l{end+1} =['        Y2[3]=startypoint[n2][2][2];'];
    l{end+1} =['        Z1[4]=startypoint[n1][1][2];'];
    l{end+1} =['        Y1[4]=startypoint[n1][2][2];'];
    l{end+1} =['        Z2[4]=startypoint[n2][1][1];'];
    l{end+1} =['        Y2[4]=startypoint[n2][2][1];'];
    l{end+1} =['        for (line=1;line<5;++line) { // 4 lines need to be checked'];
    l{end+1} =['            a[line]=(Y1[line]-Y2[line])/(Z1[line]-Z2[line]);'];
    l{end+1} =['            b[line]=Y1[line]-a[line]*Z1[line];'];
    l{end+1} =['            los_tmp[line]=1;'];
    l{end+1} =['            for (n_check=0;n_check<' num2str(last+1) ';++n_check) {'];
    l{end+1} =['                // check if the line goes between the points'];
    l{end+1} =['                if (n_check != n1 && n_check != n2 && (a[line]*startxpoint[n_check][1][1]+b[line]<startxpoint[n_check][2][1] || a[line]*startxpoint[n_check][1][2]+b[line]>startxpoint[n_check][2][2])) {'];
    l{end+1} =['                    los_tmp[line]=0; // line of sight blocked for line l!'];
    l{end+1} =['                }'];
    l{end+1} =['            }'];
    l{end+1} =['        }'];
    l{end+1} =['        if (los_tmp[1]+los_tmp[2]+los_tmp[3]+los_tmp[4]==0) {'];
    l{end+1} =['            // line of sight blocked for this specific n1 n2 combination'];
    l{end+1} =['            los_logic_single[n1][n2]=0;'];
    %l{end+1} =['            printf("los_logic_single[%%d][%%d]=0.\\n",n1,n2);'];
    l{end+1} =['        }'];
    l{end+1} =['        else {'];
    l{end+1} =['            // line of sight for this specific n1 n2 combination'];
    l{end+1} =['            los_logic_single[n1][n2]=1;'];
    %l{end+1} =['            printf("los_logic_single[%%d][%%d]=1.\\n",n1,n2);'];
    %l{end+1} =['            printf("kink rot = %lf.\\n\\n",rot' kinknumstr ');'];
    l{end+1} =['        }'];
    l{end+1} =['    }'];
    l{end+1} =['}'];
    l{end+1} =['los_check=0;'];
    l{end+1} =['for (n1=0;n1<' kinknumstr ';++n1) {'];
    l{end+1} =['    for (n2=' kinknumstr ';n2<' num2str(last+1) ';++n2) {'];
    l{end+1} =['        los_check=los_check+los_logic_single[n1][n2];'];
    l{end+1} =['    }'];
    l{end+1} =['}'];
    l{end+1} =['if (los_check==0) {'];
    l{end+1} =['    los_logic = 0;'];
    l{end+1} =['    printf("kink rot = %%lf.\\n\\n",rot' kinknumstr ');'];
    l{end+1} =['}'];
    l{end+1} =['}']; % to end the while loop
    
    
    % Write c function to write out the information nessecary to draw the
    % guide.
    
    
    
    end
    
    %l{end+1}=['printf("los_logic=%%d \\n",los_logic);'];
    
%                 % if los 
%             % globalt koordinatsystem uden approksimationer
%             los_logic=0; % antages lukket
%             for n1=1:n_kink
%                 for n2=last:n_kink
%                     % lav linjen: (4 muligheder -- ++ +- -+)
%                     punkt1=start_n1
%                     punkt2=end_n2
%                     % giv linjernes ligninger
%                     linje_1:4
% 
%                     for line_p (p=1:4)
%                     los_line(p)=1;
%                     if (linje_p < min_end_1 || linje_p>max_end_1); los_line1=0; end
%                     for i=1:last
%                         % tjek om linjen går igennem åbning
%                         if linje_p< min_start_i || linje_p>max_start_i && i ~= n1 && i ~=n2
%                            %los blocked for this line
%                            lons_linep=0;
%                         end
%                     end
%                     if los_linep==1; los_logic=1; end;
% 
%                 end
%             end
%    l{end+1} =['            // if (n_check != n1 &&  != n2 && (a[l]*endpointx[1][1]+b[l]<endpointx[2][1] || a[l]*endpointx[1][2]+b[l]>endpointx[2][2])) {'];
%    l{end+1} =['            //        los_tmp[l]=0; // line of sight blocked for line l!'];
%    l{end+1} =['            // }'];
    
    end %End index == 1 statements
    
    pcalcstring='';
    for i=1:length(l)
        pcalcstring=[pcalcstring l{i} '\n'];
    end
    clear l;
    McStasStr.pcalcstring=[McStasStr.pcalcstring '\n\n'  pcalcstring];    

%End kink treatment 
elseif (isfield(globalinfo,'curves') && globalinfo.kinkoverwrite(globalinfo.curves(1))==0 && globalinfo.raytracereq==1);
    % Will the above behave properly in case there is several curved guides
    % some of which is locked and one free?

    % declare needed c variables
    
    if sum(ismember(McStasStr.declare,['startxpoint[' num2str(last+1) '][3][3]']))<0.5 % not always needed
                McStasStr.declare{end+1}=['startxpoint[' num2str(last+1) '][3][3]'];
    end
    if sum(ismember(McStasStr.declare,['startypoint[' num2str(last+1) '][3][3]']))<0.5 % not always needed
                McStasStr.declare{end+1}=['startypoint[' num2str(last+1) '][3][3]'];
    end
    if sum(ismember(McStasStr.declare,['startXdirec[' num2str(last+1) '][3]']))<0.5 % not always needed
                McStasStr.declare{end+1}=['startXdirec[' num2str(last+1) '][3]'];
    end
    if sum(ismember(McStasStr.declare,['startYdirec[' num2str(last+1) '][3]']))<0.5 % not always needed
                McStasStr.declare{end+1}=['startYdirec[' num2str(last+1) '][3]'];
    end
    if sum(ismember(McStasStr.declare,['startXposition[' num2str(last+1) '][3]']))<0.5 % not always needed
                McStasStr.declare{end+1}=['startXposition[' num2str(last+1) '][3]'];
    end
    if sum(ismember(McStasStr.declare,['startYposition[' num2str(last+1) '][3]']))<0.5 % not always needed
                McStasStr.declare{end+1}=['startYposition[' num2str(last+1) '][3]'];
    end
    if sum(ismember(McStasStr.declare,'curve_small_radius'))<0.5 % not always needed
                McStasStr.declare{end+1}=['curve_small_radius'];
    end
    if sum(ismember(McStasStr.declare,'dist'))<0.5 % not always needed
                McStasStr.declare{end+1}=['dist'];
    end
    if globalinfo.rotlogic(globalinfo.curves(1)) > 0 % horizontal plane
       if sum(ismember(McStasStr.declare,['curveXcenter[3]']))<0.5 % not always needed
                McStasStr.declare{end+1}=['curveXcenter[3]'];
       end 
    else % vertical plane
       if sum(ismember(McStasStr.declare,['curveYcenter[3]']))<0.5 % not always needed
                McStasStr.declare{end+1}=['curveYcenter[3]'];
       end 
    end
    % initialize
    
    
    % what to do for each index
    if index == last
        % Initalize direction vector
        % Calculate corner positions
        % This uses the fact that the index == last can not be curved
        
        l{1}='';
        l{end+1} = ['startXdirec[' num '][1] = 1.0;']; % z
        l{end+1} = ['startXdirec[' num '][2] = 0.0;']; % x
        
        l{end+1} = ['startYdirec[' num '][1] = 1.0;']; % z
        l{end+1} = ['startYdirec[' num '][2] = 0.0;']; % y
        l{end+1}='';

        %l{end+1} = ['startXdirec[' numM1 '][1] = startXdirec[' num '][1];']; % z
        %l{end+1} = ['startXdirec[' numM1 '][2] = startXdirec[' num '][2];']; % x
        
        %l{end+1} = ['startYdirec[' numM1 '][1] = startYdirec[' num '][1];']; % z
        %l{end+1} = ['startYdirec[' numM1 '][2] = startYdirec[' num '][2];']; % y
        
        l{end+1} = ['startXposition[' num '][1] = endPoint' num ' - length' num ';']; % z
        l{end+1} = ['startXposition[' num '][2] = 0.0;']; % x
        
        l{end+1} = ['startYposition[' num '][1] = endPoint' num ' - length' num ';']; % z 
        l{end+1} = ['startYposition[' num '][2] = 0.0;']; % y
        l{end+1}='';
                
        
    elseif index + 1 == globalinfo.curves(1)
        curveindex = index + 1;
        % declare variables
        
        
        if sum(ismember(McStasStr.declare,['DeltaA' num ]))<0.5 % not always needed
                    McStasStr.declare{end+1}=['DeltaA' num ];
        end
        if sum(ismember(McStasStr.declare,['DeltaB' num ]))<0.5 % not always needed
                    McStasStr.declare{end+1}=['DeltaB' num ];
        end
        if sum(ismember(McStasStr.declare,['sinrot' numP1 ]))<0.5 % not always needed
                    McStasStr.declare{end+1}=['sinrot' numP1 ];
        end
        if sum(ismember(McStasStr.declare,['cosrot' numP1 ]))<0.5 % not always needed
                    McStasStr.declare{end+1}=['cosrot' numP1 ];
        end
        % calculate corner positions and new directions vector
        
        l{1}='';
        % Direction vector turned
        l{end+1} = ['cosrot' numP1 ' = cos(rot' numP1 '*DEG2RAD);'];
        l{end+1} = ['sinrot' numP1 ' = sin(rot' numP1 '*DEG2RAD);'];
        
        
        %l{end+1} = ['startXdirec[' num '][1] = startXdirec[' num '][1]*cosrot' numP1 ' - startXdirec[' num '][2]*sinrot' numP1 ';']; % z
        %l{end+1} = ['startXdirec[' num '][2] = startXdirec[' num '][1]*sinrot' numP1 ' + startXdirec[' num '][2]*cosrot' numP1 ';']; % x
        
        %l{end+1} = ['startYdirec[' num '][1] = startYdirec[' num '][1];']; % z
        %l{end+1} = ['startYdirec[' num '][2] = startYdirec[' num '][2];']; % y
        
        
        % Position vector
        l{end+1} = ['DeltaA' num ' = curve_radius' numP1 '*(1-cos(rot' numP1 '*DEG2RAD));'];
        l{end+1} = ['DeltaB' num ' = curve_radius' numP1 '* sin(rot' numP1 '*DEG2RAD);'];
        
        if globalinfo.rotlogic(curveindex) > 0 % Horizontal curve
        
        if globalinfo.rotsign(curveindex) > 0 % bends negative x direction
        % direc affected
        l{end+1} = ['startXdirec[' num '][1] = startXdirec[' numP1 '][1]*cosrot' numP1 ' + startXdirec[' numP1 '][2]*sinrot' numP1 ';']; % z
        l{end+1} = ['startXdirec[' num '][2] = - startXdirec[' numP1 '][1]*sinrot' numP1 ' + startXdirec[' numP1 '][2]*cosrot' numP1 ';']; % x
        % position affected
        l{end+1} = ['startXposition[' num '][1] = startXposition[' numP1 '][1] + startXdirec[' numP1 '][2]*DeltaA' num ' + startXdirec[' numP1 '][1]*DeltaB' num ';']; % z
        l{end+1} = ['startXposition[' num '][2] = startXposition[' numP1 '][2] - startXdirec[' numP1 '][1]*DeltaA' num ' + startXdirec[' numP1 '][2]*DeltaB' num ';']; % x
        else % bends positive x direction
        % direc affected
        l{end+1} = ['startXdirec[' num '][1] = startXdirec[' numP1 '][1]*cosrot' numP1 ' - startXdirec[' numP1 '][2]*sinrot' numP1 ';']; % z
        l{end+1} = ['startXdirec[' num '][2] = startXdirec[' numP1 '][1]*sinrot' numP1 ' + startXdirec[' numP1 '][2]*cosrot' numP1 ';']; % x
        % position affected
        l{end+1} = ['startXposition[' num '][1] = startXposition[' numP1 '][1] + startXdirec[' numP1 '][2]*DeltaA' num ' + startXdirec[' numP1 '][1]*DeltaB' num ';']; % z
        l{end+1} = ['startXposition[' num '][2] = startXposition[' numP1 '][2] + startXdirec[' numP1 '][1]*DeltaA' num ' - startXdirec[' numP1 '][2]*DeltaB' num ';'];
        end
        
        % direc Y not affected
        l{end+1} = ['startYdirec[' num '][1] = startYdirec[' numP1 '][1];']; % z
        l{end+1} = ['startYdirec[' num '][2] = startYdirec[' numP1 '][2];']; % y
        % position not affected
        l{end+1} = ['startYposition[' num '][1] = startYposition[' numP1 '][1] + startYdirec[' numP1 '][2]*DeltaA' num ' + startYdirec[' numP1 '][1]*DeltaB' num ';']; % z
        l{end+1} = ['startYposition[' num '][2] = startYposition[' numP1 '][2] + startYdirec[' numP1 '][2]*length' numP1 ';']; % x
        else % Vertical curve
        
        
        if globalinfo.rotsign(curveindex) > 0 % bends negative y direction
        % direc affected
        l{end+1} = ['startYdirec[' num '][1] = startYdirec[' numP1 '][1]*cosrot' numP1 ' + startYdirec[' numP1 '][2]*sinrot' numP1 ';']; % z
        l{end+1} = ['startYdirec[' num '][2] = - startYdirec[' numP1 '][1]*sinrot' numP1 ' + startYdirec[' numP1 '][2]*cosrot' numP1 ';']; % x
        % position affected
        l{end+1} = ['startYposition[' num '][1] = startYposition[' numP1 '][1] + startYdirec[' numP1 '][2]*DeltaA' num ' + startYdirec[' numP1 '][1]*DeltaB' num ';']; % z
        l{end+1} = ['startYposition[' num '][2] = startYposition[' numP1 '][2] - startYdirec[' numP1 '][1]*DeltaA' num ' + startYdirec[' numP1 '][2]*DeltaB' num ';']; % x
        else % bends positive y direction
        % direc affected
        l{end+1} = ['startYdirec[' num '][1] = startYdirec[' numP1 '][1]*cosrot' numP1 ' - startYdirec[' numP1 '][2]*sinrot' numP1 ';']; % z
        l{end+1} = ['startYdirec[' num '][2] = startYdirec[' numP1 '][1]*sinrot' numP1 ' + startYdirec[' numP1 '][2]*cosrot' numP1 ';']; % x
        % position affected
        l{end+1} = ['startYposition[' num '][1] = startYposition[' numP1 '][1] + startYdirec[' numP1 '][2]*DeltaA' num ' + startYdirec[' numP1 '][1]*DeltaB' num ';']; % z
        l{end+1} = ['startYposition[' num '][2] = startYposition[' numP1 '][2] + startYdirec[' numP1 '][1]*DeltaA' num ' - startYdirec[' numP1 '][2]*DeltaB' num ';']; % x
        end
        
        % direction not affected
        l{end+1} = ['startXdirec[' num '][1] = startXdirec[' numP1 '][1];']; % z
        l{end+1} = ['startXdirec[' num '][2] = startXdirec[' numP1 '][2];']; % y
        % position not affected
        l{end+1} = ['startXposition[' num '][1] = startXposition[' numP1 '][1] + startXdirec[' numP1 '][2]*DeltaA' num ' + startXdirec[' numP1 '][1]*DeltaB' num ';']; % z
        l{end+1} = ['startXposition[' num '][2] = startXposition[' numP1 '][2] + startXdirec[' numP1 '][2]*length' numP1 ';']; % x
        end
        
        l{end+1}='';
        % These makes it easier to calculate the next direction and
        % position vector.
        
        % Variables actually used in the los check, which can be calculated
        % from the position and direction vectors.
    else
        % calculate corner positions from direction vector
        
        l{1}='';
        % Direction vector
%         l{end+1} = ['startXdirec[' numM1 '][1] = startXdirec[' num '][1];']; % z
%         l{end+1} = ['startXdirec[' numM1 '][2] = startXdirec[' num '][2];']; % x
%         
%         l{end+1} = ['startYdirec[' numM1 '][1] = startYdirec[' num '][1];']; % z
%         l{end+1} = ['startYdirec[' numM1 '][2] = startYdirec[' num '][2];']; % y
        
        
        l{end+1} = ['startXdirec[' num '][1] = startXdirec[' numP1 '][1];']; % z
        l{end+1} = ['startXdirec[' num '][2] = startXdirec[' numP1 '][2];']; % x
        
        l{end+1} = ['startYdirec[' num '][1] = startYdirec[' numP1 '][1];']; % z
        l{end+1} = ['startYdirec[' num '][2] = startYdirec[' numP1 '][2];']; % y
        
        % Position vector
        l{end+1} = ['startXposition[' num '][1] = startXposition[' numP1 '][1] + startXdirec[' numP1 '][1]*length' numP1 ';']; % z
        l{end+1} = ['startXposition[' num '][2] = startXposition[' numP1 '][2] + startXdirec[' numP1 '][2]*length' numP1 ';']; % x
        
        l{end+1} = ['startYposition[' num '][1] = startYposition[' numP1 '][1] + startYdirec[' numP1 '][1]*length' numP1 ';']; % z
        l{end+1} = ['startYposition[' num '][2] = startYposition[' numP1 '][2] + startYdirec[' numP1 '][2]*length' numP1 ';']; % x
        l{end+1}='';
        % These makes it easier to calculate the next direction and
        % position vector.
        
        % Variables actually used in the los check, which can be calculated
        % from the position and direction vectors.
    end
    
        l{end+1} =['startxpoint[' num '][1][1] = startXposition[' num '][1]+startXdirec[' num '][2]*0.5*startx' num ';'];
        l{end+1} =['startxpoint[' num '][2][1] = startXposition[' num '][2]-startXdirec[' num '][1]*0.5*startx' num ';'];
        l{end+1} =['startxpoint[' num '][1][2] = startXposition[' num '][1]-startXdirec[' num '][2]*0.5*startx' num ';'];
        l{end+1} =['startxpoint[' num '][2][2] = startXposition[' num '][2]+startXdirec[' num '][1]*0.5*startx' num ';'];        

        l{end+1} =['startypoint[' num '][1][1] = startYposition[' num '][1]+startYdirec[' num '][2]*0.5*starty' num ';'];
        l{end+1} =['startypoint[' num '][2][1] = startYposition[' num '][2]-startYdirec[' num '][1]*0.5*starty' num ';'];
        l{end+1} =['startypoint[' num '][1][2] = startYposition[' num '][1]-startYdirec[' num '][2]*0.5*starty' num ';'];
        l{end+1} =['startypoint[' num '][2][2] = startYposition[' num '][2]+startYdirec[' num '][1]*0.5*starty' num ';'];
        l{end+1}='';
        
    if index == 1
        
        
        l{end+1} = ['startXdirec[' numM1 '][1] = startXdirec[' num '][1];']; % z
        l{end+1} = ['startXdirec[' numM1 '][2] = startXdirec[' num '][2];']; % x
        
        l{end+1} = ['startYdirec[' numM1 '][1] = startYdirec[' num '][1];']; % z
        l{end+1} = ['startYdirec[' numM1 '][2] = startYdirec[' num '][2];']; % y
        
        l{end+1} = ['startXposition[' numM1 '][1] = startXposition[' num '][1] + startXdirec[' num '][1]*length' num ';']; % z
        l{end+1} = ['startXposition[' numM1 '][2] = startXposition[' num '][2] + startXdirec[' num '][2]*length' num ';']; % x
        
        l{end+1} = ['startYposition[' numM1 '][1] = startYposition[' num '][1] + startYdirec[' num '][1]*length' num ';']; % z
        l{end+1} = ['startYposition[' numM1 '][2] = startYposition[' num '][2] + startYdirec[' num '][2]*length' num ';']; % x
        
        l{end+1} =['startxpoint[0][1][1] = startXposition[0][1]+startXdirec[0][2]*0.5*endx1;'];
        l{end+1} =['startxpoint[0][2][1] = startXposition[0][2]-startXdirec[0][1]*0.5*endx1;'];
        l{end+1} =['startxpoint[0][1][2] = startXposition[0][1]-startXdirec[0][2]*0.5*endx1;'];
        l{end+1} =['startxpoint[0][2][2] = startXposition[0][2]+startXdirec[0][1]*0.5*endx1;'];        

        l{end+1} =['startypoint[0][1][1] = startYposition[0][1]+startYdirec[0][2]*0.5*endy1;'];
        l{end+1} =['startypoint[0][2][1] = startYposition[0][2]-startYdirec[0][1]*0.5*endy1;'];
        l{end+1} =['startypoint[0][1][2] = startYposition[0][1]-startYdirec[0][2]*0.5*endy1;'];
        l{end+1} =['startypoint[0][2][2] = startYposition[0][2]+startYdirec[0][1]*0.5*endy1;'];
        l{end+1}='';
    end
    
    % Add info on the center and radius of the inner wall in the curved
    % guide
    if index == globalinfo.curves(1)
       % Allocate nessecary variables 
       
       if globalinfo.rotlogic(index) > 0
          % center in X Z plane 
          if globalinfo.rotsign(index) < 0 % bends in -x direction % That is not -x! REVERSE!
          l{end+1} =['curveXcenter[1] = startXposition[' num '][1] - startXdirec[' num '][2] * curve_radius' num ';'];
          l{end+1} =['curveXcenter[2] = startXposition[' num '][2] + startXdirec[' num '][1] * curve_radius' num ';'];
          else % bends in the x direction
          l{end+1} =['curveXcenter[1] = startXposition[' num '][1] + startXdirec[' num '][2] * curve_radius' num ';'];
          l{end+1} =['curveXcenter[2] = startXposition[' num '][2] - startXdirec[' num '][1] * curve_radius' num ';'];
          end
          l{end+1} =['curve_small_radius = curve_radius' num ' - 0.5*startx' num ';'];
       else % center in Y Z plane
          if globalinfo.rotsign(index) > 0 % bends in -y direction
          l{end+1} =['curveYcenter[1] = startYposition[' num '][1] - startYdirec[' num '][2] * curve_radius' num ';'];
          l{end+1} =['curveYcenter[2] = startYposition[' num '][2] + startYdirec[' num '][1] * curve_radius' num ';'];
          else % bend in the y direction
          l{end+1} =['curveYcenter[1] = startYposition[' num '][1] + startYdirec[' num '][2] * curve_radius' num ';'];
          l{end+1} =['curveYcenter[2] = startYposition[' num '][2] - startYdirec[' num '][1] * curve_radius' num ';'];     
          end
          l{end+1} =['curve_small_radius = curve_radius' num ' - 0.5*starty' num ';'];
       end
    end
    
    % if the first index, write the raytracer.
    if index==1 && isfield(globalinfo,'curves') && globalinfo.kinkoverwrite(globalinfo.curves(1))==0;   
    curvenumstr=num2str(globalinfo.curves(1));
    
    
    if sum(ismember(McStasStr.declareint,'n_check'))<0.5 % not always needed
             McStasStr.declareint{end+1}='n_check';
             McStasStr.declareint{end+1}='los_logic';
             McStasStr.declareint{end+1}='n1';
             McStasStr.declareint{end+1}='n2';
             McStasStr.declareint{end+1}='line';
             McStasStr.declareint{end+1}='los_tmp[5]';
             McStasStr.declareint{end+1}='n_check';
             McStasStr.declareint{end+1}='los_logic';
             McStasStr.declareint{end+1}='los_check';
             McStasStr.declareint{end+1}='ii';
             McStasStr.declareint{end+1}=['los_logic_single[' curvenumstr '][' num2str(last+1) ']'];      
    end
    
    if globalinfo.rotlogic(globalinfo.curves(1))>0
        if sum(ismember(McStasStr.declare,'X1[5]'))<0.5 % not always needed
                 McStasStr.declare{end+1}='X1[5]';
                 McStasStr.declare{end+1}='X2[5]';
                 McStasStr.declare{end+1}='Z1[5]';
                 McStasStr.declare{end+1}='Z2[5]';
                 McStasStr.declare{end+1}='a[5]';
                 McStasStr.declare{end+1}='b[5]';
                 McStasStr.declare{end+1}='tmp_double';
                 McStasStr.declare{end+1}='dx_circ';
                 McStasStr.declare{end+1}='dy_circ';
                 McStasStr.declare{end+1}='dr_circ';
                 McStasStr.declare{end+1}='D_circ';
                 McStasStr.declare{end+1}='x1_circ';
                 McStasStr.declare{end+1}='x2_circ';
                 McStasStr.declare{end+1}='y1_circ';
                 McStasStr.declare{end+1}='y2_circ';
                 McStasStr.declare{end+1}='sign_dy';
                 McStasStr.declare{end+1}='k_circ';
                 McStasStr.declare{end+1}='x_solution[2]';
                 McStasStr.declare{end+1}='y_solution[2]';

                 
        end
    else
        if sum(ismember(McStasStr.declare,'Y1[5]'))<0.5 % not always needed
                 McStasStr.declare{end+1}='Y1[5]';
                 McStasStr.declare{end+1}='Y2[5]';
                 McStasStr.declare{end+1}='Z1[5]';
                 McStasStr.declare{end+1}='Z2[5]';
                 McStasStr.declare{end+1}='a[5]';
                 McStasStr.declare{end+1}='b[5]';
                 McStasStr.declare{end+1}='tmp_double';
                 McStasStr.declare{end+1}='dx_circ';
                 McStasStr.declare{end+1}='dy_circ';
                 McStasStr.declare{end+1}='dr_circ';
                 McStasStr.declare{end+1}='D_circ';
                 McStasStr.declare{end+1}='x1_circ';
                 McStasStr.declare{end+1}='x2_circ';
                 McStasStr.declare{end+1}='y1_circ';
                 McStasStr.declare{end+1}='y2_circ';
                 McStasStr.declare{end+1}='sign_dy';
                 McStasStr.declare{end+1}='k_circ';
                 McStasStr.declare{end+1}='x_solution[2]';
                 McStasStr.declare{end+1}='y_solution[2]';                 
        end
    end
    
    
    if globalinfo.rotlogic(globalinfo.curves(1))>0
    l{end+1} =['los_logic = 1; // assume line of sight'];
    l{end+1} =['for (n1=0;n1<' curvenumstr ';++n1) {'];
    l{end+1} =['    for (n2=' curvenumstr ';n2<' num2str(last+1) ';++n2) {'];
    l{end+1} =['        Z1[1]=startxpoint[n1][1][1];'];
    l{end+1} =['        X1[1]=startxpoint[n1][2][1];'];
    l{end+1} =['        Z2[1]=startxpoint[n2][1][1];'];
    l{end+1} =['        X2[1]=startxpoint[n2][2][1];'];
    l{end+1} =['        Z1[2]=startxpoint[n1][1][2];'];
    l{end+1} =['        X1[2]=startxpoint[n1][2][2];'];
    l{end+1} =['        Z2[2]=startxpoint[n2][1][2];'];
    l{end+1} =['        X2[2]=startxpoint[n2][2][2];'];
    l{end+1} =['        Z1[3]=startxpoint[n1][1][1];'];
    l{end+1} =['        X1[3]=startxpoint[n1][2][1];'];
    l{end+1} =['        Z2[3]=startxpoint[n2][1][2];'];
    l{end+1} =['        X2[3]=startxpoint[n2][2][2];'];
    l{end+1} =['        Z1[4]=startxpoint[n1][1][2];'];
    l{end+1} =['        X1[4]=startxpoint[n1][2][2];'];
    l{end+1} =['        Z2[4]=startxpoint[n2][1][1];'];
    l{end+1} =['        X2[4]=startxpoint[n2][2][1];'];
    l{end+1} =['        for (line=1;line<5;++line) { // 4 lines need to be checked'];
    l{end+1} =['            tmp_double=(Z1[line]-Z2[line])/100;'];
    %l{end+1} =['            a[l]=(X1[l]-X2[l])/(Z1[l]-Z2[l]);'];
    l{end+1} =['            a[line]=(X1[line]-X2[line])/tmp_double;'];
    l{end+1} =['            tmp_double=a[line]*Z1[line];'];
    %l{end+1} =['            b[l]=X1[l]-a[l]*Z1[l];'];
    l{end+1} =['            b[line]=X1[line]-tmp_double/100;'];
    l{end+1} =['            los_tmp[line]=1;'];
    l{end+1} =['            for (n_check=0;n_check<' num2str(last+1) ';++n_check) {'];
    l{end+1} =['                // check if the line goes between the points'];
    l{end+1} =['                if (n_check != n1 && n_check != n2 && ((a[line]*startxpoint[n_check][1][1])/100+b[line]<startxpoint[n_check][2][1] || (a[line]*startxpoint[n_check][1][2])/100+b[line]>startxpoint[n_check][2][2])) {'];
    l{end+1} =['                    los_tmp[line]=0; // line of sight blocked for line l!'];
    l{end+1} =['                }'];
    l{end+1} =['            }'];
                                    % line intersects cylinder wall
                                    % 1 generer to punkter på linjen
                                    % In this z => x and x=> y.
    l{end+1} =['            // x1 and x2 chosen as start and end of curved section'];                                 
    l{end+1} =['            x1_circ = startXposition[' num2str(globalinfo.curves(1)) '][1] - curveXcenter[1];'];
    l{end+1} =['            x2_circ = startXposition[' num2str(globalinfo.curves(1)-1) '][1] - curveXcenter[1];'];
    l{end+1} =['            y1_circ = a[line]/100*startXposition[' num2str(globalinfo.curves(1)) '][1] + b[line] - curveXcenter[2];'];
    l{end+1} =['            y2_circ = a[line]/100*startXposition[' num2str(globalinfo.curves(1)-1) '][1] + b[line] - curveXcenter[2];'];
    l{end+1} =['            dx_circ = x2_circ - x1_circ;'];
    l{end+1} =['            dy_circ = y2_circ - y1_circ;'];
    l{end+1} =['            dr_circ = sqrt(dx_circ*dx_circ + dy_circ*dy_circ);'];
    l{end+1} =['            D_circ = x1_circ*y2_circ - y1_circ*x2_circ;'];
    l{end+1} =['            if (dy_circ >= 0)'];
    l{end+1} =['              sign_dy = 1;'];
    l{end+1} =['            else'];
    l{end+1} =['              sign_dy = -1;'];
    l{end+1} =['            k_circ = curve_small_radius*curve_small_radius * dr_circ * dr_circ - D_circ*D_circ;'];
    l{end+1} =[''];
    l{end+1} =['            if (k_circ > 0){'];
    l{end+1} =['             x_solution[0] = (D_circ*dy_circ + sign_dy*dx_circ*sqrt(k_circ))/(dr_circ*dr_circ);'];
    l{end+1} =['             x_solution[1] = (D_circ*dy_circ - sign_dy*dx_circ*sqrt(k_circ))/(dr_circ*dr_circ);'];
    l{end+1} =['             y_solution[0] = (-D_circ*dx_circ + fabs(dy_circ)*sqrt(k_circ))/(dr_circ*dr_circ);'];
    l{end+1} =['             y_solution[1] = (-D_circ*dx_circ - fabs(dy_circ)*sqrt(k_circ))/(dr_circ*dr_circ);'];
    l{end+1} =[''];         
    l{end+1} =['             for (ii=0;ii<2;++ii) {'];
    l{end+1} =['              if (x_solution[ii] > x1_circ && x_solution[ii] < x2_circ) {'];
    %l{end+1} =['                dist = sqrt((x_solution[ii]-curveXcenter[1])*(x_solution[ii]-curveXcenter[1]) + (y_solution[ii]-curveXcenter[2])*(y_solution[ii]-curveXcenter[2]));'];
    %l{end+1} =['                if ( dist < curve_small_radius ) {'];
    l{end+1} =['                  los_tmp[line]=0; // line of sight blocked for line l!'];
    %l{end+1} =['                }'];             
    l{end+1} =['              }'];
    l{end+1} =['             }'];
    l{end+1} =['            }'];
%     l{end+1} =['            dist = fabs(curveXcenter[1]*a[line]/100-curveXcenter[2]+b[line])/sqrt(a[line]*a[line]/10000+1);'];
%     l{end+1} =['            if ( dist < curve_small_radius ) {'];
%     l{end+1} =['                 los_tmp[line]=0; // line of sight blocked for line l!'];
%     l{end+1} =['            }'];
    l{end+1} =['        }'];
    l{end+1} =['        if (los_tmp[1]+los_tmp[2]+los_tmp[3]+los_tmp[4]==0) {'];
    l{end+1} =['            los_logic_single[n1][n2]=0;'];
    %l{end+1} =['            printf("los_logic_single[%%d][%%d]=0.\\n",n1,n2);'];
    l{end+1} =['        }'];
    l{end+1} =['        else {'];
    l{end+1} =['            los_logic_single[n1][n2]=1;'];
    %l{end+1} =['            printf("los_logic_single[%%d][%%d]=1.\\n",n1,n2);'];
    %l{end+1} =['            printf("kink rot = %%lf.\\n",rot' kinknumstr ');'];
    %l{end+1} =['            if (los_tmp[1]==1) printf("line=1 had los.\\n");'];
    %l{end+1} =['            if (los_tmp[2]==1) printf("line=2 had los.\\n");'];
    %l{end+1} =['            if (los_tmp[3]==1) printf("line=3 had los.\\n");'];
    %l{end+1} =['            if (los_tmp[4]==1) printf("line=4 had los.\\n");'];
    %l{end+1} =['            printf("\\n");'];
    l{end+1} =['        }'];
    l{end+1} =['    }'];
    l{end+1} =['}'];
    l{end+1} =['los_check=0;'];
    l{end+1} =['for (n1=0;n1<' curvenumstr ';++n1) {'];
    l{end+1} =['    for (n2=' curvenumstr ';n2<' num2str(last+1) ';++n2) {'];
    l{end+1} =['        los_check=los_check+los_logic_single[n1][n2];'];
    l{end+1} =['    }'];
    l{end+1} =['}'];
    l{end+1} =['if (los_check==0) {'];
    l{end+1} =['    los_logic = 0;'];
    l{end+1} =['    printf("kink rot = %%lf.\\n\\n",rot' curvenumstr ');'];
    l{end+1} =['}'];
    l{end+1} =['}']; % to end the while loop
    %l{end+1}=['printf("los_logic=%%d \\n",los_logic);'];
    else % Kink in the y direction
    l{end+1} =['los_logic = 1; // assume line of sight'];
    l{end+1} =['for (n1=0;n1<' curvenumstr ';++n1) {'];
    l{end+1} =['    for (n2=' curvenumstr ';n2<' num2str(last+1) ';++n2) {'];
    l{end+1} =['        Z1[1]=startypoint[n1][1][1];'];
    l{end+1} =['        Y1[1]=startypoint[n1][2][1];'];
    l{end+1} =['        Z2[1]=startypoint[n2][1][1];'];
    l{end+1} =['        Y2[1]=startypoint[n2][2][1];'];
    l{end+1} =['        Z1[2]=startypoint[n1][1][2];'];
    l{end+1} =['        Y1[2]=startypoint[n1][2][2];'];
    l{end+1} =['        Z2[2]=startypoint[n2][1][2];'];
    l{end+1} =['        Y2[2]=startypoint[n2][2][2];'];
    l{end+1} =['        Z1[3]=startypoint[n1][1][1];'];
    l{end+1} =['        Y1[3]=startypoint[n1][2][1];'];
    l{end+1} =['        Z2[3]=startypoint[n2][1][2];'];
    l{end+1} =['        Y2[3]=startypoint[n2][2][2];'];
    l{end+1} =['        Z1[4]=startypoint[n1][1][2];'];
    l{end+1} =['        Y1[4]=startypoint[n1][2][2];'];
    l{end+1} =['        Z2[4]=startypoint[n2][1][1];'];
    l{end+1} =['        Y2[4]=startypoint[n2][2][1];'];
    l{end+1} =['        for (line=1;line<5;++line) { // 4 lines need to be checked'];
    l{end+1} =['            a[line]=(Y1[line]-Y2[line])/(Z1[line]-Z2[line]);'];
    l{end+1} =['            b[line]=Y1[line]-a[line]*Z1[line];'];
    l{end+1} =['            los_tmp[line]=1;'];
    l{end+1} =['            for (n_check=0;n_check<' num2str(last+1) ';++n_check) {'];
    l{end+1} =['                // check if the line goes between the points'];
    l{end+1} =['                if (n_check != n1 && n_check != n2 && (a[line]*startxpoint[n_check][1][1]+b[line]<startxpoint[n_check][2][1] || a[line]*startxpoint[n_check][1][2]+b[line]>startxpoint[n_check][2][2])) {'];
    l{end+1} =['                    los_tmp[line]=0; // line of sight blocked for line l!'];
    l{end+1} =['                }'];
    l{end+1} =['            }'];
    l{end+1} =['            // x1 and x2 chosen as start and end of curved section in a coordinate system with the center in 0,0.'];                                 
    l{end+1} =['            x1_circ = startYposition[' num2str(globalinfo.curves(1)) '][1] - curveYcenter[1];']; % forgot -1
    l{end+1} =['            x2_circ = startYposition[' num2str(globalinfo.curves(1)-1) '][1] - curveYcenter[1];'];
    l{end+1} =['            y1_circ = a[line]/100*startYposition[' num2str(globalinfo.curves(1)) '][1] + b[line] - curveYcenter[2];'];
    l{end+1} =['            y2_circ = a[line]/100*startYposition[' num2str(globalinfo.curves(1)-1) '][1] + b[line] - curveYcenter[2];'];
    l{end+1} =['            dx_circ = x2_circ - x1_circ;'];
    l{end+1} =['            dy_circ = y2_circ - y1_circ;'];
    l{end+1} =['            dr_circ = sqrt(dx_circ*dx_circ + dy_circ*dy_circ);'];
    l{end+1} =['            D_circ = x1_circ*y2_circ - y1_circ*x2_circ;'];
    l{end+1} =['            if (dy_circ >= 0)'];
    l{end+1} =['              sign_dy = 1;'];
    l{end+1} =['            else'];
    l{end+1} =['              sign_dy = -1;'];
    l{end+1} =['            k_circ = curve_small_radius*curve_small_radius * dr_circ * dr_circ - D_circ*D_circ;'];
    l{end+1} =[''];
    l{end+1} =['            if (k_circ > 0){'];
    l{end+1} =['             x_solution[0] = (D_circ*dy_circ + sign_dy*dx_circ*sqrt(k_circ))/(dr_circ*dr_circ);'];
    l{end+1} =['             x_solution[1] = (D_circ*dy_circ - sign_dy*dx_circ*sqrt(k_circ))/(dr_circ*dr_circ);'];
    l{end+1} =['             y_solution[0] = (-D_circ*dx_circ + fabs(dy_circ)*sqrt(k_circ))/(dr_circ*dr_circ);'];
    l{end+1} =['             y_solution[1] = (-D_circ*dx_circ - fabs(dy_circ)*sqrt(k_circ))/(dr_circ*dr_circ);'];
    l{end+1} =[''];         
    l{end+1} =['             for (ii=0;ii<2;++ii) {'];
    l{end+1} =['              if (x_solution[ii] > x1_circ && x_solution[ii] < x2_circ) {'];
    %l{end+1} =['                dist = sqrt((x_solution[ii]-curveXcenter[1])*(x_solution[ii]-curveXcenter[1]) + (y_solution[ii]-curveXcenter[2])*(y_solution[ii]-curveXcenter[2]));'];
    %l{end+1} =['                if ( dist < curve_small_radius ) {'];
    l{end+1} =['                  los_tmp[line]=0; // line of sight blocked for line l!'];
    %l{end+1} =['                }'];             
    l{end+1} =['              }'];
    l{end+1} =['             }'];
    l{end+1} =['            }'];
    l{end+1} =['        }'];
    l{end+1} =['        if (los_tmp[1]+los_tmp[2]+los_tmp[3]+los_tmp[4]==0) {'];
    l{end+1} =['            // line of sight blocked for this specific n1 n2 combination'];
    l{end+1} =['            los_logic_single[n1][n2]=0;'];
    %l{end+1} =['            printf("los_logic_single[%%d][%%d]=0.\\n",n1,n2);'];
    l{end+1} =['        }'];
    l{end+1} =['        else {'];
    l{end+1} =['            // line of sight for this specific n1 n2 combination'];
    l{end+1} =['            los_logic_single[n1][n2]=1;'];
    %l{end+1} =['            printf("los_logic_single[%%d][%%d]=1.\\n",n1,n2);'];
    %l{end+1} =['            printf("kink rot = %lf.\\n\\n",rot' kinknumstr ');'];
    l{end+1} =['        }'];
    l{end+1} =['    }'];
    l{end+1} =['}'];
    l{end+1} =['los_check=0;'];
    l{end+1} =['for (n1=0;n1<' curvenumstr ';++n1) {'];
    l{end+1} =['    for (n2=' curvenumstr ';n2<' num2str(last+1) ';++n2) {'];
    l{end+1} =['        los_check=los_check+los_logic_single[n1][n2];'];
    l{end+1} =['    }'];
    l{end+1} =['}'];
    l{end+1} =['if (los_check==0) {'];
    l{end+1} =['    los_logic = 0;'];
    l{end+1} =['    printf("kink rot = %%lf.\\n\\n",rot' curvenumstr ');'];
    l{end+1} =['}'];
    l{end+1} =['}']; % to end the while loop 
    end
   
    
    l{end+1} =['fp = fopen("geometry.dat","w");'];
    l{end+1} =['fprintf(fp,"=x-z_position\\n");'];
    l{end+1} =['for (part=0;part<' num2str(last+1) ';++part) {'];
    l{end+1} =['fprintf(fp,"%%lf,%%lf\\n",startXposition[part][2],startXposition[part][1]);'];
    l{end+1} =['}'];
    l{end+1} =['fprintf(fp,"=x-z_direc\\n");'];
    l{end+1} =['for (part=0;part<' num2str(last+1) ';++part) {'];
    l{end+1} =['fprintf(fp,"%%lf,%%lf\\n",startXdirec[part][2],startXdirec[part][1]);'];
    l{end+1} =['}'];
    l{end+1} =['fprintf(fp,"=x-z_upper\\n");'];
    l{end+1} =['for (part=0;part<' num2str(last+1) ';++part) {'];
    l{end+1} =['fprintf(fp,"%%lf,%%lf\\n",startxpoint[part][2][1],startxpoint[part][1][1]);'];
    l{end+1} =['}'];
    l{end+1} =['fprintf(fp,"=x-z_lower\\n");'];
    l{end+1} =['for (part=0;part<' num2str(last+1) ';++part) {'];
    l{end+1} =['fprintf(fp,"%%lf,%%lf\\n",startxpoint[part][2][2],startxpoint[part][1][2]);'];
    l{end+1} =['}'];
    l{end+1} =['fprintf(fp,"=y-z_position\\n");'];
    l{end+1} =['for (part=0;part<' num2str(last+1) ';++part) {'];
    l{end+1} =['fprintf(fp,"%%lf,%%lf\\n",startYposition[part][2],startYposition[part][1]);'];
    l{end+1} =['}'];
    l{end+1} =['fprintf(fp,"=y-z_direc\\n");'];
    l{end+1} =['for (part=0;part<' num2str(last+1) ';++part) {'];
    l{end+1} =['fprintf(fp,"%%lf,%%lf\\n",startYdirec[part][2],startYdirec[part][1]);'];
    l{end+1} =['}'];
    l{end+1} =['fprintf(fp,"=y-z_upper\\n");'];
    l{end+1} =['for (part=0;part<' num2str(last+1) ';++part) {'];
    l{end+1} =['fprintf(fp,"%%lf,%%lf\\n",startypoint[part][2][1],startypoint[part][1][1]);'];
    l{end+1} =['}'];
    l{end+1} =['fprintf(fp,"=y-z_lower\\n");'];
    l{end+1} =['for (part=0;part<' num2str(last+1) ';++part) {'];
    l{end+1} =['fprintf(fp,"%%lf,%%lf\\n",startypoint[part][2][2],startypoint[part][1][2]);'];
    l{end+1} =['}'];
    if globalinfo.rotlogic(globalinfo.curves(1))>0
    l{end+1} =['fprintf(fp,"=curved\\n");'];
    l{end+1} =['fprintf(fp,"%%lf,%%lf,%%lf,%%lf,%%lf,%%lf,%%lf,%%lf,%%lf\\n",curve_small_radius,curveXcenter[2],curveXcenter[1],x1_circ,x2_circ,x_solution[0],x_solution[1],y_solution[0],y_solution[1]);'];
    l{end+1} =['fclose(fp);'];
    else
    l{end+1} =['fprintf(fp,"=curved\\n");'];
    l{end+1} =['fprintf(fp,"%%lf,%%lf,%%lf,%%lf,%%lf,%%lf,%%lf,%%lf,%%lf\\n",curve_small_radius,curveYcenter[2],curveYcenter[1],x1_circ,x2_circ,x_solution[0],x_solution[1],y_solution[0],y_solution[1]);'];
    l{end+1} =['fclose(fp);'];
    end

    end
    
    pcalcstring='';
    for i=1:length(l)
        pcalcstring=[pcalcstring l{i} '\n'];
    end
    clear l;
    McStasStr.pcalcstring=[McStasStr.pcalcstring '\n\n'  pcalcstring];    


end %End if statement which controlls which type of los break is used
    
else % New unified point_calc
   % This section needs to handle both kinks and curves, and not just from start to end of the guide, preferably in general.
   
   % Psudo code.
   % The resulting code should be writen directly to the initialize string,
   % which mean the order of module run and point_calc call should be
   % reversed.
   
   % if index == last
   %  Initialize position and direction
   % end
   
   
   % Write the points and directions for each module.
   
   % if index - 1 is a kink:
   %    write the special after kink section
   %    (the module will start a while loop)
   % elseif index - 1 is a curve:
   %    write the special after curve section
   %    (the module will start a while loop)
   % else
   %    write the standard concave guide section
   % end
   
   % if index == los_divider
   %    write the end of the while loop.
   % end
   
    
   % Testing the idears about multiple out of los.
   %
   % Have running list of los breakers.
   %
   % if los_breaker
   %    if los_list empty
   %       start while loop normally
   %    else
   %       creep in to while loop:
   %        add another free variable: ratio between direction changes
   %        add lines to control the second as a function of the first
   %    end
   % add module name to list (or just index number)
   % end
   
   % Additionally, needs to write the end of the loop so it can check
   % multiple non standard guides in an additive manner.
   % Lets play with that! Will make the basic task easier if succeseded
   
   
       % declare needed c variables
    
    if sum(ismember(McStasStr.declare,['startxpoint[' num2str(last+1) '][3][3]']))<0.5 % not always needed
                McStasStr.declare{end+1}=['startxpoint[' num2str(last+1) '][3][3]'];
    end
    if sum(ismember(McStasStr.declare,['startypoint[' num2str(last+1) '][3][3]']))<0.5 % not always needed
                McStasStr.declare{end+1}=['startypoint[' num2str(last+1) '][3][3]'];
    end
    if sum(ismember(McStasStr.declare,['startXdirec[' num2str(last+1) '][3]']))<0.5 % not always needed
                McStasStr.declare{end+1}=['startXdirec[' num2str(last+1) '][3]'];
    end
    if sum(ismember(McStasStr.declare,['startYdirec[' num2str(last+1) '][3]']))<0.5 % not always needed
                McStasStr.declare{end+1}=['startYdirec[' num2str(last+1) '][3]'];
    end
    if sum(ismember(McStasStr.declare,['startXposition[' num2str(last+1) '][3]']))<0.5 % not always needed
                McStasStr.declare{end+1}=['startXposition[' num2str(last+1) '][3]'];
    end
    if sum(ismember(McStasStr.declare,['startYposition[' num2str(last+1) '][3]']))<0.5 % not always needed
                McStasStr.declare{end+1}=['startYposition[' num2str(last+1) '][3]'];
    end
    if sum(ismember(McStasStr.declare,'dist'))<0.5 % not always needed
                McStasStr.declare{end+1}=['dist'];
    end
    % These three need num at the end. Use globalinfo.curve_index_list.
%     if sum(ismember(McStasStr.declare,'curve_small_radius'))<0.5 % not always needed
%                 McStasStr.declare{end+1}=['curve_small_radius'];
%     end
%     if globalinfo.rotlogic(globalinfo.curves(1)) > 0 % horizontal plane
%        if sum(ismember(McStasStr.declare,['curveXcenter[3]']))<0.5 % not always needed
%                 McStasStr.declare{end+1}=['curveXcenter[3]'];
%        end 
%     else % vertical plane
%        if sum(ismember(McStasStr.declare,['curveYcenter[3]']))<0.5 % not always needed
%                 McStasStr.declare{end+1}=['curveYcenter[3]'];
%        end 
%     end

    
    
    % First check for absolute starts, either fake or real.
    if index == last
        % REAL START
        % Initalize direction vector
        % Calculate corner positions
        % This uses the fact that the index == last can not be curved
        
        l{1}='';
        l{end+1} = ['startXdirec[' num '][1] = 1.0;']; % z
        l{end+1} = ['startXdirec[' num '][2] = 0.0;']; % x
        
        l{end+1} = ['startYdirec[' num '][1] = 1.0;']; % z
        l{end+1} = ['startYdirec[' num '][2] = 0.0;']; % y
        l{end+1}='';
        
        real_last_index = ghost_globalinfo.real_index(last);
        real_last_index_str = num2str(real_last_index);
        
        l{end+1} = ['startXposition[' num '][1] = endPoint' real_last_index_str ' - length' real_last_index_str ';']; % z
        l{end+1} = ['startXposition[' num '][2] = 0.0;']; % x
        
        l{end+1} = ['startYposition[' num '][1] = endPoint' real_last_index_str ' - length' real_last_index_str ';']; % z 
        l{end+1} = ['startYposition[' num '][2] = 0.0;']; % y
        l{end+1}='';
    elseif enable_raytracer
        if index == ghost_globalinfo.active_los{ghost_globalinfo.los_section}(1)
        % In this case one needs to write a fake start direction vector
        % It should not matter what this is.
        
        % All following point_calculations are releative, and the raytracer
        % will find if there is los or not, regardless of a
        % translation/rotation of the entire guide.
        
        % Needs to reffer to the start of the element.
        real_element_logic = ghost_globalinfo.real_index == ghost_globalinfo.real_index(index+1); % Unsure of which index to use in the last part. % Reused this in index = 1, but it malfunctioned.
        real_element = find(real_element_logic);
        ref_index = real_element(end);
        ref_index_str = num2str(ref_index);
        
        l{1}='';
        l{end+1} = '// Temporary fake start point, for this specific los section';
        l{end+1} = ['startXdirec[' ref_index_str '][1] = 1.0;']; % z
        l{end+1} = ['startXdirec[' ref_index_str '][2] = 0.0;']; % x
        
        l{end+1} = ['startYdirec[' ref_index_str '][1] = 1.0;']; % z
        l{end+1} = ['startYdirec[' ref_index_str '][2] = 0.0;']; % y
        l{end+1}='';

        %l{end+1} = ['startXposition[' num '][1] = endPoint' num ' - length' num ';']; % z
        l{end+1} = ['startXposition[' ref_index_str '][1] = 0.0;']; % x
        l{end+1} = ['startXposition[' ref_index_str '][2] = 0.0;']; % x
        
        %l{end+1} = ['startYposition[' num '][1] = endPoint' num ' - length' num ';']; % z 
        l{end+1} = ['startYposition[' ref_index_str '][1] = 0.0;']; % y
        l{end+1} = ['startYposition[' ref_index_str '][2] = 0.0;']; % y
        l{end+1}='';
        end
    end
       
   
    % Now add the relative change in position, different for each module
    if index ~= last
    if globalinfo.curve_logic(ghost_globalinfo.real_index(index+1)) % was real_index(index)+1
        curveindex = ghost_globalinfo.real_index(index+1);
        curveindex_str = num2str(curveindex);
        
        % Find refference_index, the real start of the curve
        real_curve_logic = ghost_globalinfo.real_index == curveindex;
        real_curve = find(real_curve_logic);
        ref_index = real_curve(end);
        ref_index_str = num2str(ref_index);
        
        % declare variables
        if sum(ismember(McStasStr.declare,['DeltaA' num ]))<0.5 % not always needed
                    McStasStr.declare{end+1}=['DeltaA' num ];
        end
        if sum(ismember(McStasStr.declare,['DeltaB' num ]))<0.5 % not always needed
                    McStasStr.declare{end+1}=['DeltaB' num ];
        end
        if sum(ismember(McStasStr.declare,['sinrot' numP1 ]))<0.5 % not always needed
                    McStasStr.declare{end+1}=['sinrot' numP1 ];
        end
        if sum(ismember(McStasStr.declare,['cosrot' numP1 ]))<0.5 % not always needed
                    McStasStr.declare{end+1}=['cosrot' numP1 ];
        end
        % calculate corner positions and new directions vector
        
        % Need to calculate the points at a certain point in the guide:
        % simply change rot accordingly? CHECK WHEN FINISHED: Worked fine.
        
        McStasStr.declare{end+1}=['rot_factor_g' num ];
        
        l{1}='';
        if ghost_globalinfo.los_logic(index)
            
            if ghost_globalinfo.los_mode(index) == 1
                % mode 1: percentage of guide
                l{end+1} = ['rot_factor_g' num ' = ' num2str(ghost_globalinfo.los_data(index)) ';'];
                % could make special cases for 0 and 1, but above will work
            elseif ghost_globalinfo.los_mode(index) == 2
                % mode 2: absolute length in meters.
                l{end+1} = ['rot_factor_g' num ' = ' num2str(ghost_globalinfo.los_data(index)) '/length' num2str(ghost_globalinfo.real_index(index)) ';'];
            elseif ghost_globalinfo.los_mode(index) == 3
                l{end+1} = ['rot_factor_g' num ' = (length' num2str(ghost_globalinfo.real_index(index)) ' - ' num2str(ghost_globalinfo.los_data(index)) ')/length' num2str(ghost_globalinfo.real_index(index)) ';'];
            else
                disp('ERROR, unrecognized mode!')
            end
            
            % Direction vector turned
            l{end+1} = ['cosrot' numP1 ' = cos(rot' curveindex_str '*rot_factor_g' num '*DEG2RAD);'];
            l{end+1} = ['sinrot' numP1 ' = sin(rot' curveindex_str '*rot_factor_g' num '*DEG2RAD);'];
        
            % Position vector
            l{end+1} = ['DeltaA' num ' = curve_radius' curveindex_str '*(1-cos(rot' curveindex_str '*rot_factor_g' num '*DEG2RAD));'];
            l{end+1} = ['DeltaB' num ' = curve_radius' curveindex_str '*sin(rot' curveindex_str '*rot_factor_g' num '*DEG2RAD);'];
            
        else
            l{end+1} = ['rot_factor_g' num ' = 1;'];
            
            % Direction vector turned
            l{end+1} = ['cosrot' numP1 ' = cos(rot' curveindex_str '*DEG2RAD);'];
            l{end+1} = ['sinrot' numP1 ' = sin(rot' curveindex_str '*DEG2RAD);'];
        
            % Position vector
            l{end+1} = ['DeltaA' num ' = curve_radius' curveindex_str '*(1-cos(rot' curveindex_str '*DEG2RAD));'];
            l{end+1} = ['DeltaB' num ' = curve_radius' curveindex_str '*sin(rot' curveindex_str '*DEG2RAD);'];
            
        end
        
        % Problem: If the all ghost_index build on the previous, then there
        % will be a problem when the last point on a curved adds an entire
        % curved guide after a segment of the curved guide. Will end up
        % with rot + rot*rot_factor_g angle diviation.
        % This is a problem for C,S,P and E.
        
        % Solution 1:
        % Start from scratch for each component. Needs to find the first
        % ghost_index which would have started the curved and select the
        % one just before as a starting point.
        % exchange numP1 with starting_ghost_index
        % determined starting_ghost_index from above, but is not
        % nessecaryly index+1.
        % Solution 2:
        % Could just add the rest of the guide using (1-rot_factor_g), the
        % problem with this solution is that C can be split into three
        % segments, making this difficult. Check if it is really relevant
        % to make C able to split into three segments.
        
        % Finding the first ghost_index which corresponds to the start of
        % the curved guide:
        
        % startXposition[index] is the start of module(index), which means
        % the task is to find the index which corresponds to the curve
        % section.
        
        % This is obviously curveindex in real_index, so need to find the
        % first ghost_index which have real_index = curve_index.
        
       
       
        
        
        if globalinfo.rotlogic(curveindex) > 0 % Horizontal curve
        
         if globalinfo.rotsign(curveindex) > 0 % bends negative x direction % DEBUG, normally >
          % direc affected
          l{end+1} = ['startXdirec[' num '][1] = startXdirec[' ref_index_str '][1]*cosrot' numP1 ' + startXdirec[' ref_index_str '][2]*sinrot' numP1 ';']; % z
          l{end+1} = ['startXdirec[' num '][2] = - startXdirec[' ref_index_str '][1]*sinrot' numP1 ' + startXdirec[' ref_index_str '][2]*cosrot' numP1 ';']; % x
          % position affected
          l{end+1} = ['startXposition[' num '][1] = startXposition[' ref_index_str '][1] + startXdirec[' ref_index_str '][2]*DeltaA' num ' + startXdirec[' ref_index_str '][1]*DeltaB' num ';']; % z
          l{end+1} = ['startXposition[' num '][2] = startXposition[' ref_index_str '][2] - startXdirec[' ref_index_str '][1]*DeltaA' num ' + startXdirec[' ref_index_str '][2]*DeltaB' num ';']; % x
         else % bends positive x direction
          % direc affected
          l{end+1} = ['startXdirec[' num '][1] = startXdirec[' ref_index_str '][1]*cosrot' numP1 ' - startXdirec[' ref_index_str '][2]*sinrot' numP1 ';']; % z
          l{end+1} = ['startXdirec[' num '][2] = startXdirec[' ref_index_str '][1]*sinrot' numP1 ' + startXdirec[' ref_index_str '][2]*cosrot' numP1 ';']; % x
          % position affected
          l{end+1} = ['startXposition[' num '][1] = startXposition[' ref_index_str '][1] + startXdirec[' ref_index_str '][2]*DeltaA' num ' + startXdirec[' ref_index_str '][1]*DeltaB' num ';']; % z
          %l{end+1} = ['startXposition[' num '][2] = startXposition['ref_index_str '][2] + startXdirec[' ref_index_str '][1]*DeltaA' num ' - startXdirec[' ref_index_str '][2]*DeltaB' num ';']; %BUG!
          l{end+1} = ['startXposition[' num '][2] = startXposition[' ref_index_str '][2] + startXdirec[' ref_index_str '][1]*DeltaA' num ' + startXdirec[' ref_index_str '][2]*DeltaB' num ';'];
         end
         
         % direc Y not affected
         l{end+1} = ['startYdirec[' num '][1] = startYdirec[' ref_index_str '][1];']; % z
         l{end+1} = ['startYdirec[' num '][2] = startYdirec[' ref_index_str '][2];']; % y
         % position not affected
         l{end+1} = ['startYposition[' num '][1] = startYposition[' ref_index_str '][1] + startYdirec[' ref_index_str '][2]*DeltaA' num ' + startYdirec[' ref_index_str '][1]*DeltaB' num ';']; % z
         l{end+1} = ['startYposition[' num '][2] = startYposition[' ref_index_str '][2] + startYdirec[' ref_index_str '][2]*length' curveindex_str '*rot_factor_g' num ';']; % x
        else % Vertical curve
        
        
         if globalinfo.rotsign(curveindex) > 0 % bends negative y direction
          % direc affected
          l{end+1} = ['startYdirec[' num '][1] = startYdirec[' ref_index_str '][1]*cosrot' numP1 ' + startYdirec[' ref_index_str '][2]*sinrot' numP1 ';']; % z
          l{end+1} = ['startYdirec[' num '][2] = - startYdirec[' ref_index_str '][1]*sinrot' numP1 ' + startYdirec[' ref_index_str '][2]*cosrot' numP1 ';']; % x
          % position affected
          l{end+1} = ['startYposition[' num '][1] = startYposition[' ref_index_str '][1] + startYdirec[' ref_index_str '][2]*DeltaA' num ' + startYdirec[' ref_index_str '][1]*DeltaB' num ';']; % z
          l{end+1} = ['startYposition[' num '][2] = startYposition[' ref_index_str '][2] - startYdirec[' ref_index_str '][1]*DeltaA' num ' + startYdirec[' ref_index_str '][2]*DeltaB' num ';']; % x
         else % bends positive y direction
          % direc affected
          l{end+1} = ['startYdirec[' num '][1] = startYdirec[' ref_index_str '][1]*cosrot' numP1 ' - startYdirec[' ref_index_str '][2]*sinrot' numP1 ';']; % z
          l{end+1} = ['startYdirec[' num '][2] = startYdirec[' ref_index_str '][1]*sinrot' numP1 ' + startYdirec[' ref_index_str '][2]*cosrot' numP1 ';']; % x
          % position affected
          l{end+1} = ['startYposition[' num '][1] = startYposition[' ref_index_str '][1] + startYdirec[' ref_index_str '][2]*DeltaA' num ' + startYdirec[' ref_index_str '][1]*DeltaB' num ';']; % z
          l{end+1} = ['startYposition[' num '][2] = startYposition[' ref_index_str '][2] + startYdirec[' ref_index_str '][1]*DeltaA' num ' + startYdirec[' ref_index_str '][2]*DeltaB' num ';']; % x % SIMILAR BUG?
         end
        
         % direction not affected
         l{end+1} = ['startXdirec[' num '][1] = startXdirec[' ref_index_str '][1];']; % z
         l{end+1} = ['startXdirec[' num '][2] = startXdirec[' ref_index_str '][2];']; % y
         % position not affected
         l{end+1} = ['startXposition[' num '][1] = startXposition[' ref_index_str '][1] + startXdirec[' ref_index_str '][2]*DeltaA' num ' + startXdirec[' ref_index_str '][1]*DeltaB' num ';']; % z
         l{end+1} = ['startXposition[' num '][2] = startXposition[' ref_index_str '][2] + startXdirec[' ref_index_str '][2]*length' curveindex_str '*rot_factor_g' num ';']; % x
        end
        l{end+1}='';
        
    elseif globalinfo.kink_logic(ghost_globalinfo.real_index(index+1))
        kinkindex = ghost_globalinfo.real_index(index+1);
        kinkindex_str = num2str(kinkindex);
        % Do the appropriate calculations for a kink
        
        % Only two allowed values for ghost_globalinfo.los_data, 0 or 1.
        % And ghost_globalinfo.los_mode = 1. (%).
        % Not a very likely senario, may not even allow it at all.
    
        if sum(ismember(McStasStr.declare,['sinrot' numP1 ]))<0.5 % not always needed
                    McStasStr.declare{end+1}=['sinrot' numP1 ];
        end
        if sum(ismember(McStasStr.declare,['cosrot' numP1 ]))<0.5 % not always needed
                    McStasStr.declare{end+1}=['cosrot' numP1 ];
        end
        
        l{1}='';
        % Direction vector turned
        l{end+1} = ['cosrot' numP1 ' = cos(rot' kinkindex_str '*DEG2RAD);'];
        l{end+1} = ['sinrot' numP1 ' = sin(rot' kinkindex_str '*DEG2RAD);'];
        
        % Position vector
        
        if globalinfo.rotlogic(kinkindex) > 0 % Horizontal kink
        
         if globalinfo.rotsign(kinkindex) > 0 % bends negative x direction
          % direc affected
          l{end+1} = ['startXdirec[' num '][1] = startXdirec[' numP1 '][1]*cosrot' numP1 ' + startXdirec[' numP1 '][2]*sinrot' numP1 ';']; % z
          l{end+1} = ['startXdirec[' num '][2] = - startXdirec[' numP1 '][1]*sinrot' numP1 ' + startXdirec[' numP1 '][2]*cosrot' numP1 ';']; % x
          % position affected
          % Needs to check wether the translation is in the correct direction
          l{end+1} = ['startXposition[' num '][1] = startXposition[' numP1 '][1] + startXdirec[' numP1 '][1]*length' kinkindex_str ' + startXdirec[' numP1 '][2]*translation' kinkindex_str ';']; % z
          l{end+1} = ['startXposition[' num '][2] = startXposition[' numP1 '][2] + startXdirec[' numP1 '][2]*length' kinkindex_str ' - startXdirec[' numP1 '][1]*translation' kinkindex_str ';']; % x
         else % bends positive x direction
          % direc affected
          l{end+1} = ['startXdirec[' num '][1] = startXdirec[' numP1 '][1]*cosrot' numP1 ' - startXdirec[' numP1 '][2]*sinrot' numP1 ';']; % z
          l{end+1} = ['startXdirec[' num '][2] = startXdirec[' numP1 '][1]*sinrot' numP1 ' + startXdirec[' numP1 '][2]*cosrot' numP1 ';']; % x
          % position affected
          % Need to be varified thorughly
          l{end+1} = ['startXposition[' num '][1] = startXposition[' numP1 '][1] + startXdirec[' numP1 '][1]*length' kinkindex_str ' - startXdirec[' numP1 '][2]*translation' kinkindex_str ';']; % z
          l{end+1} = ['startXposition[' num '][2] = startXposition[' numP1 '][2] + startXdirec[' numP1 '][2]*length' kinkindex_str ' + startXdirec[' numP1 '][1]*translation' kinkindex_str ';']; % x
         end
         
         % direc Y not affected
         l{end+1} = ['startYdirec[' num '][1] = startYdirec[' numP1 '][1];']; % z
         l{end+1} = ['startYdirec[' num '][2] = startYdirec[' numP1 '][2];']; % y
         % position not affected
         % The kink have the same length along the direction of startYdirec
         %  regardless of the angle being kinked and the translation
         l{end+1} = ['startYposition[' num '][1] = startYposition[' numP1 '][1] + startYdirec[' numP1 '][1]*length' kinkindex_str ';']; % z
         l{end+1} = ['startYposition[' num '][2] = startYposition[' numP1 '][2] + startYdirec[' numP1 '][2]*length' kinkindex_str ';']; % y
         
        else % Vertical kink
        
         if globalinfo.rotsign(kinkindex) > 0 % bends negative y direction
          % direc affected
          l{end+1} = ['startYdirec[' num '][1] = startYdirec[' numP1 '][1]*cosrot' numP1 ' + startYdirec[' numP1 '][2]*sinrot' numP1 ';']; % z
          l{end+1} = ['startYdirec[' num '][2] = - startYdirec[' numP1 '][1]*sinrot' numP1 ' + startYdirec[' numP1 '][2]*cosrot' numP1 ';']; % x
          % position affected
          % Needs to be verified
          l{end+1} = ['startYposition[' num '][1] = startYposition[' numP1 '][1] + startYdirec[' numP1 '][1]*length' kinkindex_str ' + startYdirec[' numP1 '][2]*translation' kinkindex_str ';']; % z
          l{end+1} = ['startYposition[' num '][2] = startYposition[' numP1 '][2] + startYdirec[' numP1 '][2]*length' kinkindex_str ' - startYdirec[' numP1 '][1]*translation' kinkindex_str ';']; % y
         else % bends positive y direction
          % direc affected
          l{end+1} = ['startYdirec[' num '][1] = startYdirec[' numP1 '][1]*cosrot' numP1 ' - startYdirec[' numP1 '][2]*sinrot' numP1 ';']; % z
          l{end+1} = ['startYdirec[' num '][2] = startYdirec[' numP1 '][1]*sinrot' numP1 ' + startYdirec[' numP1 '][2]*cosrot' numP1 ';']; % x
          % position affected
          % Needs to be verified
          l{end+1} = ['startYposition[' num '][1] = startYposition[' numP1 '][1] + startYdirec[' numP1 '][1]*length' kinkindex_str ' - startYdirec[' numP1 '][2]*translation' kinkindex_str  ';']; % z
          l{end+1} = ['startYposition[' num '][2] = startYposition[' numP1 '][2] + startYdirec[' numP1 '][2]*length' kinkindex_str ' - startYdirec[' numP1 '][1]*translation' kinkindex_str ';']; % y
         end
        
         % direction not affected
         l{end+1} = ['startXdirec[' num '][1] = startXdirec[' numP1 '][1];']; % z
         l{end+1} = ['startXdirec[' num '][2] = startXdirec[' numP1 '][2];']; % y
         % position not affected
         % The kink have the same length along the direction of startXdirec.
         l{end+1} = ['startXposition[' num '][1] = startXposition[' numP1 '][1] + startXdirec[' numP1 '][1]*length' kinkindex_str ';'];
         l{end+1} = ['startXposition[' num '][2] = startXposition[' numP1 '][2] + startXdirec[' numP1 '][2]*length' kinkindex_str ';'];
        end
        l{end+1}='';        
        
    elseif globalinfo.Selene_logic(ghost_globalinfo.real_index(index+1))
        selene_index = ghost_globalinfo.real_index(index+1);
        selene_index_str = num2str(selene_index);
        
        % Selene does not have any controls for direction, rots or rotd
        
        if sum(ismember(McStasStr.declare,['sinrot_x' numP1 ]))<0.5 % not always needed
                    McStasStr.declare{end+1}=['sinrot_x' numP1 ];
        end
        if sum(ismember(McStasStr.declare,['cosrot_x' numP1 ]))<0.5 % not always needed
                    McStasStr.declare{end+1}=['cosrot_x' numP1 ];
        end
        if sum(ismember(McStasStr.declare,['sinrot_y' numP1 ]))<0.5 % not always needed
                    McStasStr.declare{end+1}=['sinrot_y' numP1 ];
        end
        if sum(ismember(McStasStr.declare,['cosrot_y' numP1 ]))<0.5 % not always needed
                    McStasStr.declare{end+1}=['cosrot_y' numP1 ];
        end
        l{1}='';
        % Direction vector does not change, but the center is translated
        % both horizontally and vertically.
        l{end+1} = ['startXdirec[' num '][1] = startXdirec[' numP1 '][1];']; % z
        l{end+1} = ['startXdirec[' num '][2] = startXdirec[' numP1 '][2];']; % x
        l{end+1} = ['startYdirec[' num '][1] = startYdirec[' numP1 '][1];']; % z
        l{end+1} = ['startYdirec[' num '][2] = startYdirec[' numP1 '][2];']; % y
        
       
        l{end+1} = ['cosrot_x' numP1 ' = cos(selene_epsilon_y' selene_index_str '*DEG2RAD);'];  % named something else selene_epsilon_x
        l{end+1} = ['sinrot_x' numP1 ' = sin(selene_epsilon_y' selene_index_str '*DEG2RAD);']; % x and y swithed around. original naming correpond to the plane they change the angle in, here it is the axis around which they are rotated.
        l{end+1} = ['cosrot_y' numP1 ' = cos(selene_epsilon_x' selene_index_str '*DEG2RAD);'];
        l{end+1} = ['sinrot_y' numP1 ' = sin(selene_epsilon_x' selene_index_str '*DEG2RAD);'];
        
        % [0 0 1] * Rx Ry  = [ -Cox*Siy, Six, Cox*Coy]
        % [0 0 1] * Ry Rx  = [ -Siy, Coy*Six, Cox*Coy]
        % using the first one
        if sum(ismember(McStasStr.declare,['translation_x' numP1 ]))<0.5 % not always needed
             McStasStr.declare{end+1}=['translation_x' numP1 ];
             McStasStr.declare{end+1}=['translation_y' numP1 ];
             McStasStr.declare{end+1}=['translation_z' numP1 ];            
        end
        %l{end+1} = ['translation_z' numP1 ' =length' numP1 '*(cosrot_x' numP1 '*cosrot_y' numP1 ');'];
        %l{end+1} = ['translation_y' numP1 ' =length' numP1 '*sinrot_x' numP1 ';'];
        %l{end+1} = ['translation_x' numP1 ' =length' numP1 '*(-cosrot_x' numP1 '*sinrot_y' numP1 ');'];
        
        l{end+1} = ['translation_z' numP1 ' =length' numP1 '*(cosrot_x' numP1 '*cosrot_y' numP1 ');'];
        l{end+1} = ['translation_y' numP1 ' =length' numP1 '*(cosrot_y' numP1 '*sinrot_x' numP1 ');'];
        l{end+1} = ['translation_x' numP1 ' =-length' numP1 '*sinrot_y' numP1 ';'];
        
        l{end+1} = ['startXposition[' num '][1] = startXposition[' numP1 '][1] + startXdirec[' numP1 '][1]*translation_z' numP1 ' + startXdirec[' numP1 '][2]*translation_x' numP1 ';']; % z
        l{end+1} = ['startXposition[' num '][2] = startXposition[' numP1 '][2] + startXdirec[' numP1 '][2]*translation_z' numP1 ' - startXdirec[' numP1 '][1]*translation_x' numP1 ';']; % x
        l{end+1} = ['startYposition[' num '][1] = startYposition[' numP1 '][1] + startYdirec[' numP1 '][1]*translation_z' numP1 ' + startYdirec[' numP1 '][2]*translation_y' numP1 ';']; % z
        l{end+1} = ['startYposition[' num '][2] = startYposition[' numP1 '][2] + startYdirec[' numP1 '][2]*translation_z' numP1 ' - startYdirec[' numP1 '][1]*translation_y' numP1 ';']; % x
        
    else
        % Here in case of a straight concave guide, all other cases should
        % be taken care of.
        % Selene should not appear in active_los sections.
        real_index_P1 = ghost_globalinfo.real_index(index+1);
        real_index_str_P1 = num2str(real_index_P1);
        
        % Find refference_index, the real start of the element
        real_element_logic = ghost_globalinfo.real_index == ghost_globalinfo.real_index(index+1); % Unsure of which index to use in the last part.
        real_element = find(real_element_logic);
        ref_index = real_element(end);
        ref_index_str = num2str(ref_index);
        
        
        l{1}='';
        % Direction vector
        l{end+1} = ['startXdirec[' num '][1] = startXdirec[' ref_index_str '][1];']; % z
        l{end+1} = ['startXdirec[' num '][2] = startXdirec[' ref_index_str '][2];']; % x
        l{end+1} = ['startYdirec[' num '][1] = startYdirec[' ref_index_str '][1];']; % z
        l{end+1} = ['startYdirec[' num '][2] = startYdirec[' ref_index_str '][2];']; % y
        
        if ghost_globalinfo.los_logic(index)
            if ghost_globalinfo.los_mode(index) == 1
                % mode 1: percentage of guide
                factor_string = num2str(ghost_globalinfo.los_data(index));
                % could be done in a way so the user can affect this
                % afterwords, instead of hardcoding it into the instrument
                % file.
                
                % Position vector
                l{end+1} = ['startXposition[' num '][1] = startXposition[' ref_index_str '][1] + startXdirec[' ref_index_str '][1]*' factor_string '*length' real_index_str_P1 ';']; % z
                l{end+1} = ['startXposition[' num '][2] = startXposition[' ref_index_str '][2] + startXdirec[' ref_index_str '][2]*' factor_string '*length' real_index_str_P1 ';']; % x
                l{end+1} = ['startYposition[' num '][1] = startYposition[' ref_index_str '][1] + startYdirec[' ref_index_str '][1]*' factor_string '*length' real_index_str_P1 ';']; % z
                l{end+1} = ['startYposition[' num '][2] = startYposition[' ref_index_str '][2] + startYdirec[' ref_index_str '][2]*' factor_string '*length' real_index_str_P1 ';']; % y
                l{end+1}='';
            elseif ghost_globalinfo.los_mode(index) == 2
                % mode 2: absolute length in meters, from the start
                
                substitute_length = num2str(ghost_globalinfo.los_data(index));
                
                % Position vector
                l{end+1} = ['startXposition[' num '][1] = startXposition[' ref_index_str '][1] + startXdirec[' ref_index_str '][1]*' substitute_length ';']; % z
                l{end+1} = ['startXposition[' num '][2] = startXposition[' ref_index_str '][2] + startXdirec[' ref_index_str '][2]*' substitute_length ';']; % x
                l{end+1} = ['startYposition[' num '][1] = startYposition[' ref_index_str '][1] + startYdirec[' ref_index_str '][1]*' substitute_length ';']; % z
                l{end+1} = ['startYposition[' num '][2] = startYposition[' ref_index_str '][2] + startYdirec[' ref_index_str '][2]*' substitute_length ';']; % y
                l{end+1}='';
            elseif ghost_globalinfo.los_mode(index) == 3
                % mode 2: absolute length in meters, from the end
                
                substitute_length = ['(length' real_index_str_P1 ' - ' num2str(ghost_globalinfo.los_data(index)) ')'];
                
                % Position vector
                l{end+1} = ['startXposition[' num '][1] = startXposition[' ref_index_str '][1] + startXdirec[' ref_index_str '][1]*' substitute_length ';']; % z
                l{end+1} = ['startXposition[' num '][2] = startXposition[' ref_index_str '][2] + startXdirec[' ref_index_str '][2]*' substitute_length ';']; % x
                l{end+1} = ['startYposition[' num '][1] = startYposition[' ref_index_str '][1] + startYdirec[' ref_index_str '][1]*' substitute_length ';']; % z
                l{end+1} = ['startYposition[' num '][2] = startYposition[' ref_index_str '][2] + startYdirec[' ref_index_str '][2]*' substitute_length ';']; % y
                l{end+1}='';
            else
                % ERROR catch
            end
        else
            
            % Position vector
            l{end+1} = ['startXposition[' num '][1] = startXposition[' ref_index_str '][1] + startXdirec[' ref_index_str '][1]*length' real_index_str_P1 ';']; % z
            l{end+1} = ['startXposition[' num '][2] = startXposition[' ref_index_str '][2] + startXdirec[' ref_index_str '][2]*length' real_index_str_P1 ';']; % x
            l{end+1} = ['startYposition[' num '][1] = startYposition[' ref_index_str '][1] + startYdirec[' ref_index_str '][1]*length' real_index_str_P1 ';']; % z
            l{end+1} = ['startYposition[' num '][2] = startYposition[' ref_index_str '][2] + startYdirec[' ref_index_str '][2]*length' real_index_str_P1 ';']; % y
            l{end+1}='';
        end
        
        % Variables actually used in the los check, which can be calculated
        % from the position and direction vectors. 
    end
    end
    
   
    % Regardless of the module which index corresponds to, the following
    % code needs to be evaluated (outside of the main elseif).
    real_index = ghost_globalinfo.real_index(index);
    real_index_str = num2str(real_index);
    
    if ghost_globalinfo.los_logic(index) && index == 1 
        if ghost_globalinfo.los_end_data(index) == 1 && ghost_globalinfo.los_end_mode(index) == 1
            normal_case = 0;
        else
            normal_case = 1;
        end
    else
        normal_case = 1;
    end
    
    if ghost_globalinfo.los_logic(index) && normal_case == 1
            if strcmp(globalinfo.modules(globalinfo.modulelist(ghost_globalinfo.real_index(index))),'S')% corresponds to straight guide
                    % Calculates the corners of the module from the direction and position
                    % BUG! In case of los_end_logic = 1, los_end_data = 1,
                    % los_end_mode = 1 (e.g automatic guide end) this gives
                    % wrong results.
                    if ghost_globalinfo.los_mode(index) == 1
                      x_string =['0.5*(startx' real_index_str '+(endx' real_index_str ' - startx' real_index_str ')*' num2str(ghost_globalinfo.los_data(index)) ')'];
                      y_string =['0.5*(starty' real_index_str '+(endy' real_index_str ' - starty' real_index_str ')*' num2str(ghost_globalinfo.los_data(index)) ')'];
                    elseif ghost_globalinfo.los_mode(index) == 2
                      x_string =['0.5*(startx' real_index_str '+(endx' real_index_str ' - startx' real_index_str ')*' num2str(ghost_globalinfo.los_data(index)) '/length' real_index_str ')'];
                      y_string =['0.5*(starty' real_index_str '+(endy' real_index_str ' - starty' real_index_str ')*' num2str(ghost_globalinfo.los_data(index)) '/length' real_index_str ')'];
                    elseif ghost_globalinfo.los_mode(index) == 3
                      x_string =['0.5*(startx' real_index_str '+(endx' real_index_str ' - startx' real_index_str ')*(length' real_index_str '-' num2str(ghost_globalinfo.los_data(index)) ')/length' real_index_str ')'];
                      y_string =['0.5*(starty' real_index_str '+(endy' real_index_str ' - starty' real_index_str ')*(length' real_index_str '-' num2str(ghost_globalinfo.los_data(index)) ')/length' real_index_str ')'];
                    else
                        % ERROR catcher
                    end
                    
                    l{end+1} =['startxpoint[' num '][1][1] = startXposition[' num '][1]+startXdirec[' num '][2]*' x_string ';'];
                    l{end+1} =['startxpoint[' num '][2][1] = startXposition[' num '][2]-startXdirec[' num '][1]*' x_string ';'];
                    l{end+1} =['startxpoint[' num '][1][2] = startXposition[' num '][1]-startXdirec[' num '][2]*' x_string ';'];
                    l{end+1} =['startxpoint[' num '][2][2] = startXposition[' num '][2]+startXdirec[' num '][1]*' x_string ';'];

                    l{end+1} =['startypoint[' num '][1][1] = startYposition[' num '][1]+startYdirec[' num '][2]*' y_string ';'];
                    l{end+1} =['startypoint[' num '][2][1] = startYposition[' num '][2]-startYdirec[' num '][1]*' y_string ';'];
                    l{end+1} =['startypoint[' num '][1][2] = startYposition[' num '][1]-startYdirec[' num '][2]*' y_string ';'];
                    l{end+1} =['startypoint[' num '][2][2] = startYposition[' num '][2]+startYdirec[' num '][1]*' y_string ';'];
                    l{end+1}='';
            elseif strcmp(globalinfo.modules(globalinfo.modulelist(ghost_globalinfo.real_index(index))),'P') % corresponds to parabolic guide
                   % insperation from visualize (which was inspired from
                   % Kaspars guide_gravity ellipse)
                   %y_height(ii) = smallaxis_v*sqrt(1-(((xx(ii)-focus_s_v)-elength_v/2)/(elength_v/2))*(((xx(ii)-focus_s_v)-elength_v/2)/(elength_v/2)));
                   %x_width(ii)  = smallaxis_h*sqrt(1-(((xx(ii)-focus_s_h)-elength_h/2)/(elength_h/2))*(((xx(ii)-focus_s_h)-elength_h/2)/(elength_h/2)));                 
                   % declare variables
                   if sum(ismember(McStasStr.declare,'position'))<0.5 % not always needed
                    McStasStr.declare{end+1}='position';
                    McStasStr.declare{end+1}='focus_e';
                    McStasStr.declare{end+1}='elength';
                    McStasStr.declare{end+1}='width';
                    McStasStr.declare{end+1}='height';
                   end
                   
                   real_index = ghost_globalinfo.real_index(index);
                   real_index_str = num2str(real_index);
                   
                   if ghost_globalinfo.los_mode(index) == 1
                        l{end+1} =['position = ' num2str(ghost_globalinfo.los_data(index)) '*length' real_index_str ';'];
                   elseif ghost_globalinfo.los_mode(index) == 2
                        l{end+1} =['position = ' num2str(ghost_globalinfo.los_data(index)) ';'];
                   elseif ghost_globalinfo.los_mode(index) == 3
                        l{end+1} =['position = length' real_index_str ' - ' num2str(ghost_globalinfo.los_data(index)) ';'];
                   end
                   
                    l{end+1} =['focus_e = length' real_index_str ' + Loutx' real_index_str ';'];
                    l{end+1} =['elength = focus_e + Linx' real_index_str ';'];
                    l{end+1} =['width  = smallaxis_parabolic_x' real_index_str '*sqrt(1-(((position+Linx' real_index_str ')-elength/2)/(elength/2))*(((position+Linx' real_index_str ')-elength/2)/(elength/2)));'];
                    l{end+1} =['focus_e = length' real_index_str ' + Louty' real_index_str ';'];
                    l{end+1} =['elength = focus_e + Liny' real_index_str ';'];
                    l{end+1} =['height = smallaxis_parabolic_y' real_index_str '*sqrt(1-(((position+Liny' real_index_str ')-elength/2)/(elength/2))*(((position+Liny' real_index_str ')-elength/2)/(elength/2)));'];
                   
                    l{end+1} =['startxpoint[' num '][1][1] = startXposition[' num '][1]+startXdirec[' num '][2]*width;'];
                    l{end+1} =['startxpoint[' num '][2][1] = startXposition[' num '][2]-startXdirec[' num '][1]*width;'];
                    l{end+1} =['startxpoint[' num '][1][2] = startXposition[' num '][1]-startXdirec[' num '][2]*width;'];
                    l{end+1} =['startxpoint[' num '][2][2] = startXposition[' num '][2]+startXdirec[' num '][1]*width;'];

                    l{end+1} =['startypoint[' num '][1][1] = startYposition[' num '][1]+startYdirec[' num '][2]*height;'];
                    l{end+1} =['startypoint[' num '][2][1] = startYposition[' num '][2]-startYdirec[' num '][1]*height;'];
                    l{end+1} =['startypoint[' num '][1][2] = startYposition[' num '][1]-startYdirec[' num '][2]*height;'];
                    l{end+1} =['startypoint[' num '][2][2] = startYposition[' num '][2]+startYdirec[' num '][1]*height;'];
                   
            elseif strcmp(globalinfo.modules(globalinfo.modulelist(ghost_globalinfo.real_index(index))),'E') % corresponds to elliptic guide    
                   
                   % insperation from visualize (which was inspired from
                   % Kaspars guide_gravity ellipse)
                   %y_height(ii) = smallaxis_v*sqrt(1-(((xx(ii)-focus_s_v)-elength_v/2)/(elength_v/2))*(((xx(ii)-focus_s_v)-elength_v/2)/(elength_v/2)));
                   %x_width(ii)  = smallaxis_h*sqrt(1-(((xx(ii)-focus_s_h)-elength_h/2)/(elength_h/2))*(((xx(ii)-focus_s_h)-elength_h/2)/(elength_h/2)));                 
                   % declare variables
                   if sum(ismember(McStasStr.declare,'position'))<0.5 % not always needed
                    McStasStr.declare{end+1}='position';
                    McStasStr.declare{end+1}='focus_e';
                    McStasStr.declare{end+1}='elength';
                    McStasStr.declare{end+1}='width';
                    McStasStr.declare{end+1}='height';
                   end
                   
                   real_index = ghost_globalinfo.real_index(index);
                   real_index_str = num2str(real_index);
                   
                   if ghost_globalinfo.los_mode(index) == 1
                        l{end+1} =['position = ' num2str(ghost_globalinfo.los_data(index)) '*length' real_index_str ';'];
                   elseif ghost_globalinfo.los_mode(index) == 2
                        l{end+1} =['position = ' num2str(ghost_globalinfo.los_data(index)) ';'];
                   elseif ghost_globalinfo.los_mode(index) == 3
                        l{end+1} =['position = length' real_index_str ' - ' num2str(ghost_globalinfo.los_data(index)) ';'];
                   end
                   
                    l{end+1} =['focus_e = length' real_index_str ' + Loutx' real_index_str ';'];
                    l{end+1} =['elength = focus_e + Linx' real_index_str ';'];
                    l{end+1} =['width  = smallaxis_x' real_index_str '*sqrt(1-(((position+Linx' real_index_str ')-elength/2)/(elength/2))*(((position+Linx' real_index_str ')-elength/2)/(elength/2)));'];
                    l{end+1} =['focus_e = length' real_index_str ' + Louty' real_index_str ';'];
                    l{end+1} =['elength = focus_e + Liny' real_index_str ';'];
                    l{end+1} =['height = smallaxis_y' real_index_str '*sqrt(1-(((position+Liny' real_index_str ')-elength/2)/(elength/2))*(((position+Liny' real_index_str ')-elength/2)/(elength/2)));'];
                   
                    l{end+1} =['startxpoint[' num '][1][1] = startXposition[' num '][1]+startXdirec[' num '][2]*width;'];
                    l{end+1} =['startxpoint[' num '][2][1] = startXposition[' num '][2]-startXdirec[' num '][1]*width;'];
                    l{end+1} =['startxpoint[' num '][1][2] = startXposition[' num '][1]-startXdirec[' num '][2]*width;'];
                    l{end+1} =['startxpoint[' num '][2][2] = startXposition[' num '][2]+startXdirec[' num '][1]*width;'];

                    l{end+1} =['startypoint[' num '][1][1] = startYposition[' num '][1]+startYdirec[' num '][2]*height;'];
                    l{end+1} =['startypoint[' num '][2][1] = startYposition[' num '][2]-startYdirec[' num '][1]*height;'];
                    l{end+1} =['startypoint[' num '][1][2] = startYposition[' num '][1]-startYdirec[' num '][2]*height;'];
                    l{end+1} =['startypoint[' num '][2][2] = startYposition[' num '][2]+startYdirec[' num '][1]*height;'];
                   
            elseif strcmp(globalinfo.modules(globalinfo.modulelist(ghost_globalinfo.real_index(index))),'C') || strcmp(globalinfo.modules(globalinfo.modulelist(ghost_globalinfo.real_index(index))),'Cg')% corresponds to curved guide
                    % Calculates the corners of the module from the direction and position
                    l{end+1} =['startxpoint[' num '][1][1] = startXposition[' num '][1]+startXdirec[' num '][2]*0.5*startx' real_index_str ';'];
                    l{end+1} =['startxpoint[' num '][2][1] = startXposition[' num '][2]-startXdirec[' num '][1]*0.5*startx' real_index_str ';'];
                    l{end+1} =['startxpoint[' num '][1][2] = startXposition[' num '][1]-startXdirec[' num '][2]*0.5*startx' real_index_str ';'];
                    l{end+1} =['startxpoint[' num '][2][2] = startXposition[' num '][2]+startXdirec[' num '][1]*0.5*startx' real_index_str ';'];        

                    l{end+1} =['startypoint[' num '][1][1] = startYposition[' num '][1]+startYdirec[' num '][2]*0.5*starty' real_index_str ';'];
                    l{end+1} =['startypoint[' num '][2][1] = startYposition[' num '][2]-startYdirec[' num '][1]*0.5*starty' real_index_str ';'];
                    l{end+1} =['startypoint[' num '][1][2] = startYposition[' num '][1]-startYdirec[' num '][2]*0.5*starty' real_index_str ';'];
                    l{end+1} =['startypoint[' num '][2][2] = startYposition[' num '][2]+startYdirec[' num '][1]*0.5*starty' real_index_str ';'];
                    l{end+1}='';
            else
                disp(['ERROR, Module ' globalinfo.modules(globalinfo.modulelist(ghost_globalinfo.real_index(index))) ' can not begin or end the guide, or have los options'])
            end
            
    else            
        % Calculates the corners of the module from the direction and position
        l{end+1} =['startxpoint[' num '][1][1] = startXposition[' num '][1]+startXdirec[' num '][2]*0.5*startx' real_index_str ';'];
        l{end+1} =['startxpoint[' num '][2][1] = startXposition[' num '][2]-startXdirec[' num '][1]*0.5*startx' real_index_str ';'];
        l{end+1} =['startxpoint[' num '][1][2] = startXposition[' num '][1]-startXdirec[' num '][2]*0.5*startx' real_index_str ';'];
        l{end+1} =['startxpoint[' num '][2][2] = startXposition[' num '][2]+startXdirec[' num '][1]*0.5*startx' real_index_str ';'];        

        l{end+1} =['startypoint[' num '][1][1] = startYposition[' num '][1]+startYdirec[' num '][2]*0.5*starty' real_index_str ';'];
        l{end+1} =['startypoint[' num '][2][1] = startYposition[' num '][2]-startYdirec[' num '][1]*0.5*starty' real_index_str ';'];
        l{end+1} =['startypoint[' num '][1][2] = startYposition[' num '][1]-startYdirec[' num '][2]*0.5*starty' real_index_str ';'];
        l{end+1} =['startypoint[' num '][2][2] = startYposition[' num '][2]+startYdirec[' num '][1]*0.5*starty' real_index_str ';'];
        l{end+1}='';
    end
        
    
        
    if index == 1 && globalinfo.Selene_logic(ghost_globalinfo.real_index(index))
        selene_index = 1;
        selene_index_str = num2str(selene_index);
        % Selene does not have any controls for direction, rots or rotd
        
        if sum(ismember(McStasStr.declare,['sinrot_x1']))<0.5 % not always needed
                    McStasStr.declare{end+1}=['sinrot_x1'];
        end
        if sum(ismember(McStasStr.declare,['cosrot_x1']))<0.5 % not always needed
                    McStasStr.declare{end+1}=['cosrot_x1'];
        end
        if sum(ismember(McStasStr.declare,['sinrot_y1']))<0.5 % not always needed
                    McStasStr.declare{end+1}=['sinrot_y1'];
        end
        if sum(ismember(McStasStr.declare,['cosrot_y1']))<0.5 % not always needed
                    McStasStr.declare{end+1}=['cosrot_y1'];
        end
        l{1}='';
        % Direction vector does not change, but the center is translated
        % both horizontally and vertically.
        l{end+1} = ['startXdirec[0][1] = startXdirec[1][1];']; % z
        l{end+1} = ['startXdirec[0][2] = startXdirec[1][2];']; % x
        l{end+1} = ['startYdirec[0][1] = startYdirec[1][1];']; % z
        l{end+1} = ['startYdirec[0][2] = startYdirec[1][2];']; % y
        
       
        l{end+1} = ['cosrot_x1 = cos(selene_epsilon_y' selene_index_str '*DEG2RAD);'];  % named something else selene_epsilon_x
        l{end+1} = ['sinrot_x1 = sin(selene_epsilon_y' selene_index_str '*DEG2RAD);']; % x and y swithed around. original naming correpond to the plane they change the angle in, here it is the axis around which they are rotated.
        l{end+1} = ['cosrot_y1 = cos(selene_epsilon_x' selene_index_str '*DEG2RAD);'];
        l{end+1} = ['sinrot_y1 = sin(selene_epsilon_x' selene_index_str '*DEG2RAD);'];
        
        % [0 0 1] * Rx Ry  = [ -Cox*Siy, Six, Cox*Coy]
        % [0 0 1] * Ry Rx  = [ -Siy, Coy*Six, Cox*Coy]
        % using the first one
        if sum(ismember(McStasStr.declare,['translation_x1']))<0.5 % not always needed
             McStasStr.declare{end+1}=['translation_x1'];
             McStasStr.declare{end+1}=['translation_y1'];
             McStasStr.declare{end+1}=['translation_z1'];            
        end
        %l{end+1} = ['translation_z' numP1 ' =length' numP1 '*(cosrot_x' numP1 '*cosrot_y' numP1 ');'];
        %l{end+1} = ['translation_y' numP1 ' =length' numP1 '*sinrot_x' numP1 ';'];
        %l{end+1} = ['translation_x' numP1 ' =length' numP1 '*(-cosrot_x' numP1 '*sinrot_y' numP1 ');'];
        
        l{end+1} = ['translation_z1= length1*(cosrot_x1*cosrot_y1);'];
        l{end+1} = ['translation_y1= length1*(cosrot_y1*sinrot_x1);'];
        l{end+1} = ['translation_x1=-length1*sinrot_y1;'];
        
        l{end+1} = ['startXposition[0][1] = startXposition[1][1] + startXdirec[1][1]*translation_z1 + startXdirec[1][2]*translation_x1;']; % z
        l{end+1} = ['startXposition[0][2] = startXposition[1][2] + startXdirec[1][2]*translation_z1 - startXdirec[1][1]*translation_x1;']; % x
        l{end+1} = ['startYposition[0][1] = startYposition[1][1] + startYdirec[1][1]*translation_z1 + startYdirec[1][2]*translation_y1;']; % z
        l{end+1} = ['startYposition[0][2] = startYposition[1][2] + startYdirec[1][2]*translation_z1 - startYdirec[1][1]*translation_y1;']; % x
        
    elseif index == 1
        %ghost_globalinfo.real_index(index) == 1 
        % Needs to reffer to the start of the element.
        real_element_logic = ghost_globalinfo.real_index == ghost_globalinfo.real_index(index); % Unsure of which index to use in the last part. % BUG was wrong index+1 in this case
        real_element = find(real_element_logic);
        ref_index = real_element(end);
        ref_index_str = num2str(ref_index);
        
        l{end+1} = ['startXdirec[0][1] = startXdirec[' ref_index_str '][1];']; % z
        l{end+1} = ['startXdirec[0][2] = startXdirec[' ref_index_str '][2];']; % x
        l{end+1} = ['startYdirec[0][1] = startYdirec[' ref_index_str '][1];']; % z
        l{end+1} = ['startYdirec[0][2] = startYdirec[' ref_index_str '][2];']; % y
        
        l{end+1} = ['startXposition[0][1] = startXposition[' ref_index_str '][1] + startXdirec[' ref_index_str '][1]*length' real_index_str ';']; % z
        l{end+1} = ['startXposition[0][2] = startXposition[' ref_index_str '][2] + startXdirec[' ref_index_str '][2]*length' real_index_str ';']; % x
        l{end+1} = ['startYposition[0][1] = startYposition[' ref_index_str '][1] + startYdirec[' ref_index_str '][1]*length' real_index_str ';']; % z
        l{end+1} = ['startYposition[0][2] = startYposition[' ref_index_str '][2] + startYdirec[' ref_index_str '][2]*length' real_index_str ';']; % y
        
    end
    
    if index == 1
        % It is okay to hardcode the below numbers, because they are the
        % only possible case for this code.
        l{end+1} =['startxpoint[0][1][1] = startXposition[0][1]+startXdirec[0][2]*0.5*endx1;'];
        l{end+1} =['startxpoint[0][2][1] = startXposition[0][2]-startXdirec[0][1]*0.5*endx1;'];
        l{end+1} =['startxpoint[0][1][2] = startXposition[0][1]-startXdirec[0][2]*0.5*endx1;'];
        l{end+1} =['startxpoint[0][2][2] = startXposition[0][2]+startXdirec[0][1]*0.5*endx1;'];        

        l{end+1} =['startypoint[0][1][1] = startYposition[0][1]+startYdirec[0][2]*0.5*endy1;'];
        l{end+1} =['startypoint[0][2][1] = startYposition[0][2]-startYdirec[0][1]*0.5*endy1;'];
        l{end+1} =['startypoint[0][1][2] = startYposition[0][1]-startYdirec[0][2]*0.5*endy1;'];
        l{end+1} =['startypoint[0][2][2] = startYposition[0][2]+startYdirec[0][1]*0.5*endy1;'];
        l{end+1}='';
    end
        
    % Add info on the center and radius of the inner wall in the curved
    % guide
    if globalinfo.curve_logic(ghost_globalinfo.real_index(index))
        real_index = ghost_globalinfo.real_index(index);
        real_index_str = num2str(real_index);
        
       if sum(ismember(McStasStr.declare,['curve_small_radius' real_index_str ]))<0.5 % not always needed
                McStasStr.declare{end+1}=['curve_small_radius' real_index_str];
       end
       if globalinfo.rotlogic(ghost_globalinfo.real_index(index)) > 0 % horizontal plane
        if sum(ismember(McStasStr.declare,['curveXcenter' real_index_str '[3]']))<0.5 % not always needed
                McStasStr.declare{end+1}=['curveXcenter' real_index_str '[3]'];
        end 
       else % vertical plane
        if sum(ismember(McStasStr.declare,['curveYcenter' real_index_str '[3]']))<0.5 % not always needed
                McStasStr.declare{end+1}=['curveYcenter' real_index_str '[3]'];
        end 
       end
       % Allocate nessecary variables 
       if globalinfo.rotlogic(ghost_globalinfo.real_index(index)) > 0
          % center in X Z plane 
          if globalinfo.rotsign(ghost_globalinfo.real_index(index)) < 0 % bends in -x direction % That is not -x! REVERSE!
          l{end+1} =['curveXcenter' real_index_str '[1] = startXposition[' num '][1] - startXdirec[' num '][2] * curve_radius' real_index_str ';'];
          l{end+1} =['curveXcenter' real_index_str '[2] = startXposition[' num '][2] + startXdirec[' num '][1] * curve_radius' real_index_str ';'];
          else % bends in the x direction
          l{end+1} =['curveXcenter' real_index_str '[1] = startXposition[' num '][1] + startXdirec[' num '][2] * curve_radius' real_index_str ';'];
          l{end+1} =['curveXcenter' real_index_str '[2] = startXposition[' num '][2] - startXdirec[' num '][1] * curve_radius' real_index_str ';'];
          end
          l{end+1} =['curve_small_radius' real_index_str ' = curve_radius' real_index_str ' - 0.5*startx' real_index_str ';'];
       else % center in Y Z plane
          if globalinfo.rotsign(ghost_globalinfo.real_index(index)) < 0 % bends in -y direction
          l{end+1} =['curveYcenter' real_index_str '[1] = startYposition[' num '][1] - startYdirec[' num '][2] * curve_radius' real_index_str ';'];
          l{end+1} =['curveYcenter' real_index_str '[2] = startYposition[' num '][2] + startYdirec[' num '][1] * curve_radius' real_index_str ';'];
          else % bend in the y direction
          l{end+1} =['curveYcenter' real_index_str '[1] = startYposition[' num '][1] + startYdirec[' num '][2] * curve_radius' real_index_str ';'];
          l{end+1} =['curveYcenter' real_index_str '[2] = startYposition[' num '][2] - startYdirec[' num '][1] * curve_radius' real_index_str ';'];     
          end
          l{end+1} =['curve_small_radius' real_index_str ' = curve_radius' real_index_str ' - 0.5*starty' real_index_str ';'];
       end
    end
   
    
    % RAYTRACER SECTION (actually writing the raytracer using the points
    % calculated above)
    
    % condition for inserting a raytracer / while stop
    if enable_raytracer % needs to be in a seperate if, because active_los is not allocated in some cases
    if index == ghost_globalinfo.active_los{ghost_globalinfo.los_section}(end) 
        
    % Needed regardless of rotlogic.
    if sum(ismember(McStasStr.declareint,'n_check'))<0.5 % not always needed
         McStasStr.declareint{end+1}='n_check';
         McStasStr.declareint{end+1}='los_logic';
         McStasStr.declareint{end+1}='n1';
         McStasStr.declareint{end+1}='n2';
         McStasStr.declareint{end+1}='line';
         McStasStr.declareint{end+1}='los_tmp[5]';
         McStasStr.declareint{end+1}='n_check';
         McStasStr.declareint{end+1}='los_logic';
         McStasStr.declareint{end+1}='los_check';
         McStasStr.declareint{end+1}='ii';
         McStasStr.declareint{end+1}=['los_logic_single[' num2str(last) '][' num2str(last+1) ']']; % May need special attention
    end
   
    raytracer_start_index = ghost_globalinfo.active_los{ghost_globalinfo.los_section}(end); % -1 because it is the start of the earlier module.
    raytracer_end_index = ghost_globalinfo.active_los{ghost_globalinfo.los_section}(1); 
    
    if raytracer_start_index == 1
        % This fix does not work in case the last module have an los_end
        % different from the actual end of the guide! 
        if ghost_globalinfo.los_end_data(1) == 1 && ghost_globalinfo.los_end_mode(1) == 1
        raytracer_start_index = 0;
        end
    end
    
    % Time to check if what direction the overall raytracer should check.
    X_direc = 0 < sum(globalinfo.rotlogic(ghost_globalinfo.real_index(ghost_globalinfo.active_los{ghost_globalinfo.los_section})) == 1);
    Y_direc = 0 < sum(globalinfo.rotlogic(ghost_globalinfo.real_index(ghost_globalinfo.active_los{ghost_globalinfo.los_section})) == -1);
    
    if X_direc && ~Y_direc
        section_rotlogic = 1;
        % X direction
    elseif ~X_direc && Y_direc
        section_rotlogic = -1;
        % Y direction
    elseif X_direc && Y_direc
        % both
        disp('ERROR, only use rotd=''h'' or rotd=''v'' for each line-of-sight section! (Sorry for spamming this error)')
    else
        % none
        disp('ERROR, trying to write a raytracer without need for a raytracer? BUG or invalid input string (rotlogic)')
        globalinfo.rotlogic
        globalinfo.rotsign
    end
    
    
    % Time to check if what direction the overall raytracer should check.
    right_direc = 0 < sum(globalinfo.rotsign(ghost_globalinfo.real_index(ghost_globalinfo.active_los{ghost_globalinfo.los_section})) == 1);
    left_direc = 0 < sum(globalinfo.rotsign(ghost_globalinfo.real_index(ghost_globalinfo.active_los{ghost_globalinfo.los_section})) == -1);
    
    if right_direc && left_direc
        % both
        disp('ERROR, only use rots=1 or rots=-1 for each line-of-sight section, will not be out of los! (Sorry for spamming this error)')
    elseif right_direc && ~left_direc
        % Normal, ok!
    elseif ~right_direc && left_direc
        % Normal, ok!
    elseif ~right_direc && ~left_direc
        disp('ERROR, trying to write a raytracer without need for a raytracer? BUG or invalid input string (rots)')
        globalinfo.rotlogic
        globalinfo.rotsign
    end
    
        
    if section_rotlogic > 0
        % run rotlogic corresponding to bending in X direction
     if sum(ismember(McStasStr.declare,'X1[5]'))<0.5 % not always needed
         McStasStr.declare{end+1}='X1[5]';
         McStasStr.declare{end+1}='X2[5]';
         McStasStr.declare{end+1}='Z1[5]';
         McStasStr.declare{end+1}='Z2[5]';
         McStasStr.declare{end+1}='a[5]';
         McStasStr.declare{end+1}='b[5]';
         McStasStr.declare{end+1}='tmp_double';
     end
    
    
    l{end+1} =['los_logic = 1; // assume line of sight'];
    l{end+1} =['for (n1=' num2str(raytracer_start_index) ';n1<' num2str(raytracer_end_index) ';++n1) {'];
    l{end+1} =['    for (n2=n1+1;n2<' num2str(raytracer_end_index+1) ';++n2) {'];
    l{end+1} =['        Z1[1]=startxpoint[n1][1][1];'];
    l{end+1} =['        X1[1]=startxpoint[n1][2][1];'];
    l{end+1} =['        Z2[1]=startxpoint[n2][1][1];'];
    l{end+1} =['        X2[1]=startxpoint[n2][2][1];'];
    l{end+1} =['        Z1[2]=startxpoint[n1][1][2];'];
    l{end+1} =['        X1[2]=startxpoint[n1][2][2];'];
    l{end+1} =['        Z2[2]=startxpoint[n2][1][2];'];
    l{end+1} =['        X2[2]=startxpoint[n2][2][2];'];
    l{end+1} =['        Z1[3]=startxpoint[n1][1][1];'];
    l{end+1} =['        X1[3]=startxpoint[n1][2][1];'];
    l{end+1} =['        Z2[3]=startxpoint[n2][1][2];'];
    l{end+1} =['        X2[3]=startxpoint[n2][2][2];'];
    l{end+1} =['        Z1[4]=startxpoint[n1][1][2];'];
    l{end+1} =['        X1[4]=startxpoint[n1][2][2];'];
    l{end+1} =['        Z2[4]=startxpoint[n2][1][1];'];
    l{end+1} =['        X2[4]=startxpoint[n2][2][1];'];
    l{end+1} =['        for (line=1;line<5;++line) { // 4 lines need to be checked'];
    l{end+1} =['            tmp_double=(Z1[line]-Z2[line])/100;'];
    %l{end+1} =['            a[l]=(X1[l]-X2[l])/(Z1[l]-Z2[l]);'];
    l{end+1} =['            a[line]=(X1[line]-X2[line])/tmp_double;'];
    l{end+1} =['            tmp_double=a[line]*Z1[line];'];
    %l{end+1} =['            b[l]=X1[l]-a[l]*Z1[l];'];
    l{end+1} =['            b[line]=X1[line]-tmp_double/100;'];
    l{end+1} =['            los_tmp[line]=1;'];
    l{end+1} =['            for (n_check=' num2str(raytracer_start_index) ';n_check<' num2str(raytracer_end_index + 1) ';++n_check) {'];
    l{end+1} =['                // check if the line goes between the points'];
    l{end+1} =['                if (n_check != n1 && n_check != n2 && ((a[line]*startxpoint[n_check][1][1])/100+b[line]<startxpoint[n_check][2][1] || (a[line]*startxpoint[n_check][1][2])/100+b[line]>startxpoint[n_check][2][2])) {'];
    l{end+1} =['                    los_tmp[line]=0; // line of sight blocked for line l!'];
    l{end+1} =['                }'];
    l{end+1} =['            }'];
    for index_inner = ghost_globalinfo.active_los{ghost_globalinfo.los_section}
    
     if globalinfo.kink_logic(ghost_globalinfo.real_index(index_inner))
         % Nothing needs to be added for the kink, but is on the list to
         % demonstrate how to add modules to this list.
         
     elseif globalinfo.Selene_logic(ghost_globalinfo.real_index(index_inner))
         % Selene guides should have special care in this loop.
         % Not done yet.
         disp('ERRRO, should not use Selene in los breaker segments!')
        
     elseif globalinfo.curve_logic(ghost_globalinfo.real_index(index_inner))
         real_index_str = num2str(ghost_globalinfo.real_index(index_inner));
         
        if sum(ismember(McStasStr.declare,'dx_circ'))<0.5
                 McStasStr.declare{end+1}='dx_circ';
                 McStasStr.declare{end+1}='dy_circ';
                 McStasStr.declare{end+1}='dr_circ';
                 McStasStr.declare{end+1}='D_circ';
                 McStasStr.declare{end+1}='x1_circ';
                 McStasStr.declare{end+1}='x2_circ';
                 McStasStr.declare{end+1}='y1_circ';
                 McStasStr.declare{end+1}='y2_circ';
                 McStasStr.declare{end+1}='sign_dy';
                 McStasStr.declare{end+1}='k_circ';
                 McStasStr.declare{end+1}='x_solution[2]';
                 McStasStr.declare{end+1}='y_solution[2]';                 
        end
                                    % line intersects cylinder wall
                                    % 1 generer to punkter på linjen
                                    % In this z => x and x=> y.
    l{end+1} =['            // x1 and x2 chosen as start and end of curved section'];                                 
    l{end+1} =['            x1_circ = startXposition[' num2str(index_inner) '][1] - curveXcenter' real_index_str '[1];'];
    l{end+1} =['            x2_circ = startXposition[' num2str(index_inner-1) '][1] - curveXcenter' real_index_str '[1];'];
    l{end+1} =['            y1_circ = a[line]/100*startXposition[' num2str(index_inner) '][1] + b[line] - curveXcenter' real_index_str '[2];'];
    l{end+1} =['            y2_circ = a[line]/100*startXposition[' num2str(index_inner-1) '][1] + b[line] - curveXcenter' real_index_str '[2];'];
    l{end+1} =['            dx_circ = x2_circ - x1_circ;'];
    l{end+1} =['            dy_circ = y2_circ - y1_circ;'];
    l{end+1} =['            dr_circ = sqrt(dx_circ*dx_circ + dy_circ*dy_circ);'];
    l{end+1} =['            D_circ = x1_circ*y2_circ - y1_circ*x2_circ;'];
    l{end+1} =['            if (dy_circ >= 0)'];
    l{end+1} =['              sign_dy = 1;'];
    l{end+1} =['            else'];
    l{end+1} =['              sign_dy = -1;'];
    l{end+1} =['            k_circ = curve_small_radius' real_index_str '*curve_small_radius' real_index_str '* dr_circ * dr_circ - D_circ*D_circ;'];
    l{end+1} =[''];
    l{end+1} =['            if (k_circ > 0){'];
    l{end+1} =['             x_solution[0] = (D_circ*dy_circ + sign_dy*dx_circ*sqrt(k_circ))/(dr_circ*dr_circ);'];
    l{end+1} =['             x_solution[1] = (D_circ*dy_circ - sign_dy*dx_circ*sqrt(k_circ))/(dr_circ*dr_circ);'];
    l{end+1} =['             y_solution[0] = (-D_circ*dx_circ + fabs(dy_circ)*sqrt(k_circ))/(dr_circ*dr_circ);'];
    l{end+1} =['             y_solution[1] = (-D_circ*dx_circ - fabs(dy_circ)*sqrt(k_circ))/(dr_circ*dr_circ);'];
    l{end+1} =[''];         
    l{end+1} =['             for (ii=0;ii<2;++ii) {'];
    l{end+1} =['              if (x_solution[ii] > x1_circ && x_solution[ii] < x2_circ) {'];
    %l{end+1} =['                dist = sqrt((x_solution[ii]-curveXcenter[1])*(x_solution[ii]-curveXcenter[1]) + (y_solution[ii]-curveXcenter[2])*(y_solution[ii]-curveXcenter[2]));'];
    %l{end+1} =['                if ( dist < curve_small_radius ) {'];
    l{end+1} =['                  los_tmp[line]=0; // line of sight blocked for line l!'];
    %l{end+1} =['                }'];             
    l{end+1} =['              }'];
    l{end+1} =['             }'];
    l{end+1} =['            }'];
     end
    end
    
    
    l{end+1} =['        }'];
    l{end+1} =['        if (los_tmp[1]+los_tmp[2]+los_tmp[3]+los_tmp[4]==0) {'];
    l{end+1} =['            los_logic_single[n1][n2]=0;'];
    %l{end+1} =['            printf("los_logic_single[%%d][%%d]=0.\\n",n1,n2);'];
    l{end+1} =['        }'];
    l{end+1} =['        else {'];
    l{end+1} =['            los_logic_single[n1][n2]=1;'];
    %l{end+1} =['            printf("los_logic_single[%%d][%%d]=1.\\n",n1,n2);'];
    %l{end+1} =['            printf("kink rot = %%lf.\\n",rot' kinknumstr ');'];
    %l{end+1} =['            if (los_tmp[1]==1) printf("line=1 had los.\\n");'];
    %l{end+1} =['            if (los_tmp[2]==1) printf("line=2 had los.\\n");'];
    %l{end+1} =['            if (los_tmp[3]==1) printf("line=3 had los.\\n");'];
    %l{end+1} =['            if (los_tmp[4]==1) printf("line=4 had los.\\n");'];
    %l{end+1} =['            printf("\\n");'];
    l{end+1} =['        }'];
    l{end+1} =['    }'];
    l{end+1} =['}'];
    l{end+1} =['los_check=0;'];
    l{end+1} =['for (n1=' num2str(raytracer_start_index) ';n1<' num2str(raytracer_end_index) ';++n1) {'];
    l{end+1} =['    for (n2=n1+1;n2<' num2str(raytracer_end_index + 1) ';++n2) {'];
    l{end+1} =['        los_check=los_check+los_logic_single[n1][n2];'];
    l{end+1} =['    }'];
    l{end+1} =['}'];
    l{end+1} =['if (los_check==0) {'];
    l{end+1} =['    los_logic = 0;'];
    %l{end+1} =['    printf("kink rot = %%lf.\\n\\n",rot' curvenumstr ');'];
    l{end+1} =['}'];
    %l{end+1} =['if (rot' real_index_str ' > 2) los_logic = 0;']; % DELETE WHEN DEBUGGING IS DONE
    l{end+1} =['}']; % to end the while loop
        
    else
        % run rotlogic corresponding to bending in Y direction
    
    if sum(ismember(McStasStr.declare,'Y1[5]'))<0.5 % not always needed
            McStasStr.declare{end+1}='Y1[5]';
            McStasStr.declare{end+1}='Y2[5]';
            McStasStr.declare{end+1}='Z1[5]';
            McStasStr.declare{end+1}='Z2[5]';
            McStasStr.declare{end+1}='a[5]';
            McStasStr.declare{end+1}='b[5]';
            McStasStr.declare{end+1}='tmp_double';
    end
   
    l{end+1} =['los_logic = 1; // assume line of sight'];
    l{end+1} =['for (n1=' num2str(raytracer_start_index) ';n1<' num2str(raytracer_end_index) ';++n1) {'];
    l{end+1} =['    for (n2=n1+1;n2<' num2str(raytracer_end_index + 1) ';++n2) {'];
    l{end+1} =['        Z1[1]=startypoint[n1][1][1];'];
    l{end+1} =['        Y1[1]=startypoint[n1][2][1];'];
    l{end+1} =['        Z2[1]=startypoint[n2][1][1];'];
    l{end+1} =['        Y2[1]=startypoint[n2][2][1];'];
    l{end+1} =['        Z1[2]=startypoint[n1][1][2];'];
    l{end+1} =['        Y1[2]=startypoint[n1][2][2];'];
    l{end+1} =['        Z2[2]=startypoint[n2][1][2];'];
    l{end+1} =['        Y2[2]=startypoint[n2][2][2];'];
    l{end+1} =['        Z1[3]=startypoint[n1][1][1];'];
    l{end+1} =['        Y1[3]=startypoint[n1][2][1];'];
    l{end+1} =['        Z2[3]=startypoint[n2][1][2];'];
    l{end+1} =['        Y2[3]=startypoint[n2][2][2];'];
    l{end+1} =['        Z1[4]=startypoint[n1][1][2];'];
    l{end+1} =['        Y1[4]=startypoint[n1][2][2];'];
    l{end+1} =['        Z2[4]=startypoint[n2][1][1];'];
    l{end+1} =['        Y2[4]=startypoint[n2][2][1];'];
    l{end+1} =['        for (line=1;line<5;++line) { // 4 lines need to be checked'];
    %l{end+1} =['            a[line]=(Y1[line]-Y2[line])/(Z1[line]-Z2[line]);'];
    %l{end+1} =['            b[line]=Y1[line]-a[line]*Z1[line];'];
    l{end+1} =['            tmp_double=(Z1[line]-Z2[line])/100;'];
    %l{end+1} =['            a[l]=(X1[l]-X2[l])/(Z1[l]-Z2[l]);'];
    l{end+1} =['            a[line]=(Y1[line]-Y2[line])/tmp_double;'];
    l{end+1} =['            tmp_double=a[line]*Z1[line];'];
    %l{end+1} =['            b[l]=X1[l]-a[l]*Z1[l];'];
    l{end+1} =['            b[line]=Y1[line]-tmp_double/100;'];
    l{end+1} =['            los_tmp[line]=1;'];
    l{end+1} =['            for (n_check=' num2str(raytracer_start_index) ';n_check<' num2str(raytracer_end_index + 1) ';++n_check) {'];
    l{end+1} =['                // check if the line goes between the points'];
    l{end+1} =['                if (n_check != n1 && n_check != n2 && (a[line]*startypoint[n_check][1][1]/100+b[line]<startypoint[n_check][2][1] || a[line]*startypoint[n_check][1][2]/100+b[line]>startypoint[n_check][2][2])) {'];
    l{end+1} =['                    los_tmp[line]=0; // line of sight blocked for line l!'];
    l{end+1} =['                }'];
    l{end+1} =['            }'];
    % ------------ Special section, can be generated like so much other
    % code in this program when it is needed. Simply use the list of which
    % weird modules that was in this los section and loop over them in a
    % loop that adds the appropriate code.
    
    for index_inner = ghost_globalinfo.active_los{ghost_globalinfo.los_section}
    
     if globalinfo.kink_logic(ghost_globalinfo.real_index(index_inner))
         % Nothing needs to be added for the kink, but is on the list to
         % demonstrate how to add modules to this list.
         
     elseif globalinfo.Selene_logic(ghost_globalinfo.real_index(index_inner))
         % Selene guides should have special care in this loop.
         % Not done yet.
         disp('ERROR, should not use Selene in los breaker segments!')
        
     elseif globalinfo.curve_logic(ghost_globalinfo.real_index(index_inner))
         real_index_str = num2str(ghost_globalinfo.real_index(index_inner));
    % --------------------------------------------- Curve only!
    
        % Only declare these if never done before.
        % Make sure the name used to check is unique
        if sum(ismember(McStasStr.declare,'dx_circ'))<0.5
                 McStasStr.declare{end+1}='dx_circ';
                 McStasStr.declare{end+1}='dy_circ';
                 McStasStr.declare{end+1}='dr_circ';
                 McStasStr.declare{end+1}='D_circ';
                 McStasStr.declare{end+1}='x1_circ';
                 McStasStr.declare{end+1}='x2_circ';
                 McStasStr.declare{end+1}='y1_circ';
                 McStasStr.declare{end+1}='y2_circ';
                 McStasStr.declare{end+1}='sign_dy';
                 McStasStr.declare{end+1}='k_circ';
                 McStasStr.declare{end+1}='x_solution[2]';
                 McStasStr.declare{end+1}='y_solution[2]';                 
        end
    
    
     l{end+1} =['            // x1 and x2 chosen as start and end of curved section in a coordinate system with the center in 0,0.'];                                 
     l{end+1} =['            x1_circ = startYposition[' num2str(index_inner) '][1] - curveYcenter' real_index_str '[1];'];
     l{end+1} =['            x2_circ = startYposition[' num2str(index_inner-1) '][1] - curveYcenter' real_index_str '[1];'];
     l{end+1} =['            y1_circ = a[line]/100*startYposition[' num2str(index_inner) '][1] + b[line] - curveYcenter' real_index_str '[2];'];
     l{end+1} =['            y2_circ = a[line]/100*startYposition[' num2str(index_inner-1) '][1] + b[line] - curveYcenter' real_index_str '[2];'];
     l{end+1} =['            dx_circ = x2_circ - x1_circ;'];
     l{end+1} =['            dy_circ = y2_circ - y1_circ;'];
     l{end+1} =['            dr_circ = sqrt(dx_circ*dx_circ + dy_circ*dy_circ);'];
     l{end+1} =['            D_circ = x1_circ*y2_circ - y1_circ*x2_circ;'];
     l{end+1} =['            if (dy_circ >= 0)'];
     l{end+1} =['              sign_dy = 1;'];
     l{end+1} =['            else'];
     l{end+1} =['              sign_dy = -1;'];
     l{end+1} =['            k_circ = curve_small_radius' real_index_str '*curve_small_radius' real_index_str '* dr_circ * dr_circ - D_circ*D_circ;'];
     l{end+1} =[''];
     l{end+1} =['            if (k_circ > 0){'];
     l{end+1} =['             x_solution[0] = (D_circ*dy_circ + sign_dy*dx_circ*sqrt(k_circ))/(dr_circ*dr_circ);'];
     l{end+1} =['             x_solution[1] = (D_circ*dy_circ - sign_dy*dx_circ*sqrt(k_circ))/(dr_circ*dr_circ);'];
     l{end+1} =['             y_solution[0] = (-D_circ*dx_circ + fabs(dy_circ)*sqrt(k_circ))/(dr_circ*dr_circ);'];
     l{end+1} =['             y_solution[1] = (-D_circ*dx_circ - fabs(dy_circ)*sqrt(k_circ))/(dr_circ*dr_circ);'];
     l{end+1} =[''];         
     l{end+1} =['             for (ii=0;ii<2;++ii) {'];
     l{end+1} =['              if (x_solution[ii] > x1_circ && x_solution[ii] < x2_circ) {'];
     %l{end+1} =['                dist = sqrt((x_solution[ii]-curveXcenter[1])*(x_solution[ii]-curveXcenter[1]) + (y_solution[ii]-curveXcenter[2])*(y_solution[ii]-curveXcenter[2]));'];
     %l{end+1} =['                if ( dist < curve_small_radius ) {'];
     l{end+1} =['                  los_tmp[line]=0; // line of sight blocked for line l!'];
     %l{end+1} =['                }'];             
     l{end+1} =['              }'];
     l{end+1} =['             }'];
     l{end+1} =['            }'];
    % --------------------------------------------- Curve only!
     end
    end
    % ---------------
    l{end+1} =['        }'];
    l{end+1} =['        if (los_tmp[1]+los_tmp[2]+los_tmp[3]+los_tmp[4]==0) {'];
    l{end+1} =['            // line of sight blocked for this specific n1 n2 combination'];
    l{end+1} =['            los_logic_single[n1][n2]=0;'];
    %l{end+1} =['            printf("los_logic_single[%%d][%%d]=0.\\n",n1,n2);'];
    l{end+1} =['        }'];
    l{end+1} =['        else {'];
    l{end+1} =['            // line of sight for this specific n1 n2 combination'];
    l{end+1} =['            los_logic_single[n1][n2]=1;'];
    %l{end+1} =['            printf("los_logic_single[%%d][%%d]=1.\\n",n1,n2);'];
    %l{end+1} =['            printf("kink rot = %lf.\\n\\n",rot' kinknumstr ');'];
    l{end+1} =['        }'];
    l{end+1} =['    }'];
    l{end+1} =['}'];
    l{end+1} =['los_check=0;'];
    l{end+1} =['for (n1=' num2str(raytracer_start_index) ';n1<' num2str(raytracer_end_index) ';++n1) {'];
    l{end+1} =['    for (n2=n1+1;n2<' num2str(raytracer_end_index + 1) ';++n2) {'];
    l{end+1} =['        los_check=los_check+los_logic_single[n1][n2];'];
    l{end+1} =['    }'];
    l{end+1} =['}'];
    l{end+1} =['if (los_check==0) {'];
    l{end+1} =['    los_logic = 0;'];
    %l{end+1} =['    printf("kink rot = %%lf.\\n\\n",rot' curvenumstr ');'];
    l{end+1} =['}'];
    l{end+1} =['}']; % to end the while loop 
    
    end % end rotlogic 
    
    
    % In case curved guide have a second raytracer, the while loop needs to
    % be inserted here.
    if 1 < ghost_globalinfo.los_section  
        % Check if there is an overlap between this and the next active_los
        
        found=0;
        overlap=0;
        clear overlap_ghost_indexes;
        for ii = 1:length(ghost_globalinfo.real_index(ghost_globalinfo.active_los{ghost_globalinfo.los_section}))
            for jj = 1:length(ghost_globalinfo.real_index(ghost_globalinfo.active_los{ghost_globalinfo.los_section-1}))
                if ghost_globalinfo.real_index(ghost_globalinfo.active_los{ghost_globalinfo.los_section}(ii)) == ghost_globalinfo.real_index(ghost_globalinfo.active_los{ghost_globalinfo.los_section-1}(jj))
                   found = found + 1;
                   overlap(found) = ghost_globalinfo.real_index(ghost_globalinfo.active_los{ghost_globalinfo.los_section}(ii));
                   overlap_ghost_indexes{found} = [ ghost_globalinfo.active_los{ghost_globalinfo.los_section}(ii) ghost_globalinfo.active_los{ghost_globalinfo.los_section-1}(jj)];
                end
            end
        end
        if overlap(1) ~= 0
            %disp('overlap from point calc')
            %disp(overlap)
            %disp('With these indexes')
            
        
            if sum(globalinfo.curve_logic(overlap)) > 0.5
                
                abort = 0;
                if length(overlap) == 1
                   same_real_index_logic = ghost_globalinfo.real_index == overlap;
                   same_real_index = find(same_real_index_logic);
                   %disp('same_real_index')
                   %disp(same_real_index)
                   if same_real_index(1) == overlap_ghost_indexes{1}(1) || same_real_index(1) == overlap_ghost_indexes{1}(2)
                       % do nothing
                       abort = 1;
                       %disp('killed overlap correction --------------- !!!!!!!! ---------------')
                   end
                end
                
                
                if abort == 0;
                  % can only happen with curved guides, if more modules are added with the properties nessecary for this problem, simply make and elseif here.
                    
                    % This happens because curved guide can have several
                    % active_los sections at the same time, while being an
                    % active line of sight blocker.
                    
                    
                    % Should not count if one of the los sections is only
                    % using the end point of the curve which is not
                    % affected by the rotation.
                    %disp('overlap with curved! ---------------- * ---------------')
                    %disp(overlap)
                    
                    used_num = num2str(overlap(1));
                    
                    l{end+1} =['// Experimental code'];
                    l{end+1} =['//' num2str(overlap)];
                    l{end+1}=['var_divreq_x_protected = var_divreq_x;'];
                    l{end+1}=['var_divreq_y_protected = var_divreq_y;'];
                    l{end+1}=[''];
                    % If this is the controling rot for the next section,
                    % it should have rot = 0 first.
                    l{end+1}=['rot' used_num ' = rot' used_num ' - 0.0002;'];
                    l{end+1}=['los_logic = 1;'];
                    l{end+1}='while(los_logic==1) {';
                    l{end+1}=['rot' used_num ' = rot' used_num ' + 0.0002;'];
                    l{end+1}=['curve_radius' used_num ' = length' used_num '/(rot' used_num '*DEG2RAD);'];
                    l{end+1}=[''];
                    l{end+1}=['var_divreq_x = var_divreq_x_protected;'];
                    l{end+1}=['var_divreq_y = var_divreq_y_protected;'];
                    l{end+1}=[''];
            
                    % find C modules in current los and next los
                    %ghost_curve_modules_current = false(length(ghost_globalinfo.active_los{ghost_globalinfo.section_los}),1);
                    %for ii = 1:length(ghost_globalinfo.active_los{ghost_globalinfo.section_los})
                    %ghost_curve_modules_current(ii) = globalinfo.curve_logic(ghost_globalinfo.real_index(ghost_globalinfo.active_los{ghost_globalinfo.section_los}(ii)));
                    %end
                    %ghost_curve_modules_next
                end
            end
        end
    end    
    
    end % end raytracer section 
    end % end enable raytracer
    
         pcalcstring='';
        for i=1:length(l)
            pcalcstring=[pcalcstring l{i} '\n'];
        end
        clear l;
        %McStasStr.initialize=[pcalcstring '\n\n' McStasStr.initialize];
        McStasStr.pcalcstring=[McStasStr.pcalcstring '\n\n' pcalcstring ];
        
end

    % Will need to write a unified point_calc_writer function to prepare
    % for more times out of line of sight
    
    % Needs to use the vector approach as with the curved guide
    % Can almost be copy pasted to the kinked
    % Or the kink can be included in the logic behind the curved (easier)
    
    % Will need to check for not just one kink / curved section, but a more
    % general list of modules which will need special attention.
    % if the index of the special modules are given as a list of, one can
    % use sum(list == index) which will be 0 or 1.
    % if the index of the special modules are given as a logic, simply
    % check if the logic corresponding to the current index is true.
    
    % With more than one los breaker several while loops will be initiated.
    % They should not be nested, but one after the other.
    % In this process, the initialize code and point_calc / raytracer must
    % be mixed. Not seperated as they are now.
    

end

