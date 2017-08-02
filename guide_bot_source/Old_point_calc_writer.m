function [McStasStr] = point_calc_writer(McStasStr,index,last,options,globalinfo,defaults);
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
    if sum(ismember(McStasStr.declare,'curve_small_radius'))<0.5 % not always needed
                McStasStr.declare{end+1}=['curve_small_radius'];
    end
    if sum(ismember(McStasStr.declare,'dist'))<0.5 % not always needed
                McStasStr.declare{end+1}=['dist'];
    end
    % These two need num at the end. Use globalinfo.curve_index_list.
    if globalinfo.rotlogic(globalinfo.curves(1)) > 0 % horizontal plane
       if sum(ismember(McStasStr.declare,['curveXcenter[3]']))<0.5 % not always needed
                McStasStr.declare{end+1}=['curveXcenter[3]'];
       end 
    else % vertical plane
       if sum(ismember(McStasStr.declare,['curveYcenter[3]']))<0.5 % not always needed
                McStasStr.declare{end+1}=['curveYcenter[3]'];
       end 
    end

    
    
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
                
        
    elseif globalinfo.curve_logic(index + 1)
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
    elseif globalinfo.kink_logic(index + 1)
        % Do the appropriate calculations for a kink
        
    %elseif index==globalinfo.kinks(1) % need to be generalized to have more kinks
    % at kink (may need to be done in the kink module instead)
    
        l{1}='';
        % Direction vector turned
        l{end+1} = ['cosrot' numP1 ' = cos(rot' numP1 '*DEG2RAD);'];
        l{end+1} = ['sinrot' numP1 ' = sin(rot' numP1 '*DEG2RAD);'];
        
        % Position vector
        
        if globalinfo.rotlogic(index+1) > 0 % Horizontal kink
        
         if globalinfo.rotsign(index+1) > 0 % bends negative x direction
          % direc affected
          l{end+1} = ['startXdirec[' num '][1] = startXdirec[' numP1 '][1]*cosrot' numP1 ' + startXdirec[' numP1 '][2]*sinrot' numP1 ';']; % z
          l{end+1} = ['startXdirec[' num '][2] = - startXdirec[' numP1 '][1]*sinrot' numP1 ' + startXdirec[' numP1 '][2]*cosrot' numP1 ';']; % x
          % position affected
          % Needs to check wether the translation is in the correct direction
          l{end+1} = ['startXposition[' num '][1] = startXposition[' numP1 '][1] + startXdirec[' numP1 '][1]*length' numP1 ' + startXdirec[' numP1 '][2]*0.5*translation' numP1 ';']; % z
          l{end+1} = ['startXposition[' num '][2] = startXposition[' numP1 '][2] + startXdirec[' numP1 '][2]*length' numP1 ' - startXdirec[' numP1 '][1]*0.5*translation' numP1 ';']; % x
         else % bends positive x direction
          % direc affected
          l{end+1} = ['startXdirec[' num '][1] = startXdirec[' numP1 '][1]*cosrot' numP1 ' - startXdirec[' numP1 '][2]*sinrot' numP1 ';']; % z
          l{end+1} = ['startXdirec[' num '][2] = startXdirec[' numP1 '][1]*sinrot' numP1 ' + startXdirec[' numP1 '][2]*cosrot' numP1 ';']; % x
          % position affected
          % Need to be varified thorughly
          l{end+1} = ['startXposition[' num '][1] = startXposition[' numP1 '][1] + startXdirec[' numP1 '][1]*length' numP1 ' - startXdirec[' numP1 '][2]*0.5*translation' numP1 ';']; % z
          l{end+1} = ['startXposition[' num '][2] = startXposition[' numP1 '][2] + startXdirec[' numP1 '][2]*length' numP1 ' + startXdirec[' numP1 '][1]*0.5*translation' numP1 ';']; % x
         end
         
         % direc Y not affected
         l{end+1} = ['startYdirec[' num '][1] = startYdirec[' numP1 '][1];']; % z
         l{end+1} = ['startYdirec[' num '][2] = startYdirec[' numP1 '][2];']; % y
         % position not affected
         % The kink have the same length along the direction of startYdirec
         %  regardless of the angle being kinked and the translation
         l{end+1} = ['startYposition[' num '][1] = startYposition[' numP1 '][1] + startYdirec[' numP1 '][1]*length' numP1 ';']; % z
         l{end+1} = ['startYposition[' num '][2] = startYposition[' numP1 '][2] + startYdirec[' numP1 '][2]*length' numP1 ';']; % y
         
        else % Vertical curve
        
         if globalinfo.rotsign(curveindex) > 0 % bends negative y direction
          % direc affected
          l{end+1} = ['startYdirec[' num '][1] = startYdirec[' numP1 '][1]*cosrot' numP1 ' + startYdirec[' numP1 '][2]*sinrot' numP1 ';']; % z
          l{end+1} = ['startYdirec[' num '][2] = - startYdirec[' numP1 '][1]*sinrot' numP1 ' + startYdirec[' numP1 '][2]*cosrot' numP1 ';']; % x
          % position affected
          % Needs to be verified
          l{end+1} = ['startYposition[' num '][1] = startYposition[' numP1 '][1] + startYdirec[' numP1 '][1]*length' numP1 ' + startYdirec[' numP1 '][2]*0.5*translation' numP1 ';']; % z
          l{end+1} = ['startYposition[' num '][2] = startYposition[' numP1 '][2] + startYdirec[' numP1 '][2]*length' numP1 ' - startYdirec[' numP1 '][1]*0.5*trasnlation' numP1 ';']; % y
         else % bends positive y direction
          % direc affected
          l{end+1} = ['startYdirec[' num '][1] = startYdirec[' numP1 '][1]*cosrot' numP1 ' - startYdirec[' numP1 '][2]*sinrot' numP1 ';']; % z
          l{end+1} = ['startYdirec[' num '][2] = startYdirec[' numP1 '][1]*sinrot' numP1 ' + startYdirec[' numP1 '][2]*cosrot' numP1 ';']; % x
          % position affected
          % Needs to be verified
          l{end+1} = ['startYposition[' num '][1] = startYposition[' numP1 '][1] + startYdirec[' numP1 '][1]*length' numP1 ' - startYdirec[' numP1 '][2]*0.5*translation' numP1  ';']; % z
          l{end+1} = ['startYposition[' num '][2] = startYposition[' numP1 '][2] + startYdirec[' numP1 '][2]*length' numP1 ' - startYdirec[' numP1 '][1]*0.5*translation' numP1 ';']; % y
         end
        
         % direction not affected
         l{end+1} = ['startXdirec[' num '][1] = startXdirec[' numP1 '][1];']; % z
         l{end+1} = ['startXdirec[' num '][2] = startXdirec[' numP1 '][2];']; % y
         % position not affected
         % The kink have the same length along the direction of startXdirec.
         l{end+1} = ['startXposition[' num '][1] = startXposition[' numP1 '][1] + startXdirec[' numP1 '][1]*length' numP1 ';'];
         l{end+1} = ['startXposition[' num '][2] = startXposition[' numP1 '][2] + startXdirec[' numP1 '][2]*length' numP1 ';'];
        end
        l{end+1}='';        
        
    else
        l{1}='';
        % Direction vector
        l{end+1} = ['startXdirec[' num '][1] = startXdirec[' numP1 '][1];']; % z
        l{end+1} = ['startXdirec[' num '][2] = startXdirec[' numP1 '][2];']; % x
        l{end+1} = ['startYdirec[' num '][1] = startYdirec[' numP1 '][1];']; % z
        l{end+1} = ['startYdirec[' num '][2] = startYdirec[' numP1 '][2];']; % y
        
        % Position vector
        l{end+1} = ['startXposition[' num '][1] = startXposition[' numP1 '][1] + startXdirec[' numP1 '][1]*length' numP1 ';']; % z
        l{end+1} = ['startXposition[' num '][2] = startXposition[' numP1 '][2] + startXdirec[' numP1 '][2]*length' numP1 ';']; % x
        l{end+1} = ['startYposition[' num '][1] = startYposition[' numP1 '][1] + startYdirec[' numP1 '][1]*length' numP1 ';']; % z
        l{end+1} = ['startYposition[' num '][2] = startYposition[' numP1 '][2] + startYdirec[' numP1 '][2]*length' numP1 ';']; % y
        l{end+1}='';
        % Variables actually used in the los check, which can be calculated
        % from the position and direction vectors.
    end
    
    
    % Regardless of the module which index corresponds to, the following
    % code needs to be evaluated (outside of the main elseif).
    
        % Calculates the corners of the module from the direction and position
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
        l{end+1} = ['startYposition[' numM1 '][2] = startYposition[' num '][2] + startYdirec[' num '][2]*length' num ';']; % y
        
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
   
    if 1==1 % condition for inserting a raytracer / while stop
        
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
         McStasStr.declareint{end+1}=['los_logic_single[' curvenumstr '][' num2str(last+1) ']']; % May need special attention
    end
   
    if globalinfo.rotlogic(globalinfo.curves(1))>0 % Way to simple.
        % run rotlogic corresponding to bending in X direction
     if sum(ismember(McStasStr.declare,'X1[5]'))<0.5 % not always needed
         McStasStr.declare{end+1}='X1[5]';
         McStasStr.declare{end+1}='X2[5]';
         McStasStr.declare{end+1}='Z1[5]';
         McStasStr.declare{end+1}='Z2[5]';
         McStasStr.declare{end+1}='a[5]';
         McStasStr.declare{end+1}='b[5]';
     end
     
    l{end+1} =['los_logic = 1; // assume line of sight'];
    l{end+1} =['for (n1=' los_start_index ';n1<' los_end_index ';++n1) {'];
    l{end+1} =['    for (n2=n1;n2<' los_end_index + 1 ';++n2) {'];
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
    l{end+1} =['            for (n_check=' los_start_index ';n_check<' los_end_index + 1 ';++n_check) {'];
    l{end+1} =['                // check if the line goes between the points'];
    l{end+1} =['                if (n_check != n1 && n_check != n2 && ((a[line]*startxpoint[n_check][1][1])/100+b[line]<startxpoint[n_check][2][1] || (a[line]*startxpoint[n_check][1][2])/100+b[line]>startxpoint[n_check][2][2])) {'];
    l{end+1} =['                    los_tmp[line]=0; // line of sight blocked for line l!'];
    l{end+1} =['                }'];
    l{end+1} =['            }'];
    for index_inner = los_breakers.index_list
    
     if globalinfo.kinks(index_inner) == 1
         % Nothing needs to be added for the kink, but is on the list to
         % demonstrate how to add modules to this list.
         
     elseif globalinfo.Selene(index_inner) == 1
         % Selene guides should have special care in this loop.
         % Not done yet.
         disp('ERRRO, should not use Selene in los breaker segments!')
        
     elseif globalinfo.curves(index_inner) == 1
         
        if sum(ismember(McStasStr.declare,'tmp_double'))<0.5
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
                                    % line intersects cylinder wall
                                    % 1 generer to punkter på linjen
                                    % In this z => x and x=> y.
    l{end+1} =['            // x1 and x2 chosen as start and end of curved section'];                                 
    l{end+1} =['            x1_circ = startXposition[' num2str(inner_index) '][1] - curveXcenter[1];'];
    l{end+1} =['            x2_circ = startXposition[' num2str(inner_index-1) '][1] - curveXcenter[1];'];
    l{end+1} =['            y1_circ = a[line]/100*startXposition[' num2str(inner_index) '][1] + b[line] - curveXcenter[2];'];
    l{end+1} =['            y2_circ = a[line]/100*startXposition[' num2str(inner_index-1) '][1] + b[line] - curveXcenter[2];'];
    l{end+1} =['            dx_circ = x2_circ - x1_circ;'];
    l{end+1} =['            dy_circ = y2_circ - y1_circ;'];
    l{end+1} =['            dr_circ = sqrt(dx_circ*dx_circ + dy_circ*dy_circ);'];
    l{end+1} =['            D_circ = x1_circ*y2_circ - y1_circ*x2_circ;'];
    l{end+1} =['            if (dy_circ >= 0)'];
    l{end+1} =['              sign_dy = 1;'];
    l{end+1} =['            else'];
    l{end+1} =['              sign_dy = -1;'];
    l{end+1} =['            k_circ = curve_small_radius' num2str(inner_index) '*curve_small_radius' num2str(inner_index) ' * dr_circ * dr_circ - D_circ*D_circ;'];
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
    l{end+1} =['for (n1=' los_start_index ';n1<' los_end_index ';++n1) {'];
    l{end+1} =['    for (n2=n1;n2<' los_end_index + 1 ';++n2) {'];
    l{end+1} =['        los_check=los_check+los_logic_single[n1][n2];'];
    l{end+1} =['    }'];
    l{end+1} =['}'];
    l{end+1} =['if (los_check==0) {'];
    l{end+1} =['    los_logic = 0;'];
    l{end+1} =['    printf("kink rot = %%lf.\\n\\n",rot' curvenumstr ');'];
    l{end+1} =['}'];
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
    end
   
    l{end+1} =['los_logic = 1; // assume line of sight'];
    l{end+1} =['for (n1=' los_start_index ';n1<' los_end_index ';++n1) {'];
    l{end+1} =['    for (n2=n1;n2<' los_end_index + 1 ';++n2) {'];
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
    l{end+1} =['            for (n_check=' los_start_index ';n_check<' los_end_index + 1 ';++n_check) {'];
    l{end+1} =['                // check if the line goes between the points'];
    l{end+1} =['                if (n_check != n1 && n_check != n2 && (a[line]*startxpoint[n_check][1][1]+b[line]<startxpoint[n_check][2][1] || a[line]*startxpoint[n_check][1][2]+b[line]>startxpoint[n_check][2][2])) {'];
    l{end+1} =['                    los_tmp[line]=0; // line of sight blocked for line l!'];
    l{end+1} =['                }'];
    l{end+1} =['            }'];
    % ------------ Special section, can be generated like so much other
    % code in this program when it is needed. Simply use the list of which
    % weird modules that was in this los section and loop over them in a
    % loop that adds the appropriate code.
    
    for index_inner = los_breakers.index_list
    
     if globalinfo.kinks(index_inner) == 1
         % Nothing needs to be added for the kink, but is on the list to
         % demonstrate how to add modules to this list.
         
     elseif globalinfo.Selene(index_inner) == 1
         % Selene guides should have special care in this loop.
         % Not done yet.
         disp('ERRRO, should not use Selene in los breaker segments!')
        
     elseif globalinfo.curves(index_inner) == 1
    % --------------------------------------------- Curve only!
    
        % Only declare these if never done before.
        % Make sure the name used to check is unique
        if sum(ismember(McStasStr.declare,'tmp_double'))<0.5
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
    
    
     l{end+1} =['            // x1 and x2 chosen as start and end of curved section in a coordinate system with the center in 0,0.'];                                 
     l{end+1} =['            x1_circ = startYposition[' num2str(index_inner) '][1] - curveYcenter[1];'];
     l{end+1} =['            x2_circ = startYposition[' num2str(index_inner-1) '][1] - curveYcenter[1];'];
     l{end+1} =['            y1_circ = a[line]/100*startYposition[' num2str(index_inner) '][1] + b[line] - curveYcenter[2];'];
     l{end+1} =['            y2_circ = a[line]/100*startYposition[' num2str(index_inner-1) '][1] + b[line] - curveYcenter[2];'];
     l{end+1} =['            dx_circ = x2_circ - x1_circ;'];
     l{end+1} =['            dy_circ = y2_circ - y1_circ;'];
     l{end+1} =['            dr_circ = sqrt(dx_circ*dx_circ + dy_circ*dy_circ);'];
     l{end+1} =['            D_circ = x1_circ*y2_circ - y1_circ*x2_circ;'];
     l{end+1} =['            if (dy_circ >= 0)'];
     l{end+1} =['              sign_dy = 1;'];
     l{end+1} =['            else'];
     l{end+1} =['              sign_dy = -1;'];
     l{end+1} =['            k_circ = curve_small_radius' num2str(inner_index) '*curve_small_radius' num2str(inner_index) ' * dr_circ * dr_circ - D_circ*D_circ;'];
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
    los_breakers.index_list = 'empty';
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
    l{end+1} =['for (n1=' los_start_index ';n1<' los_end_index ';++n1) {'];
    l{end+1} =['    for (n2=n1;n2<' los_end_index + 1 ';++n2) {'];
    l{end+1} =['        los_check=los_check+los_logic_single[n1][n2];'];
    l{end+1} =['    }'];
    l{end+1} =['}'];
    l{end+1} =['if (los_check==0) {'];
    l{end+1} =['    los_logic = 0;'];
    %l{end+1} =['    printf("kink rot = %%lf.\\n\\n",rot' curvenumstr ');'];
    l{end+1} =['}'];
    l{end+1} =['}']; % to end the while loop 
    
    end % end rotlogic
    end % end raytracer section 
    
        pcalcstring='';
    for i=1:length(l)
        pcalcstring=[pcalcstring l{i} '\n'];
    end
    clear l;
    McStasStr.pcalcstring=[McStasStr.pcalcstring '\n\n'  pcalcstring];
   
    
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

