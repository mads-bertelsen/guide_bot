% 


% if los 
% globalt koordinatsystem uden approksimationer
los_logic=0; % antages lukket
for n1=1:n_kink
    for n2=last:n_kink
        % lav linjen: (4 muligheder -- ++ +- -+)
        punkt1=start_n1
        punkt2=end_n2
        % giv linjernes ligninger
        linje_1:4
        
        for line_p (p=1:4)
        los_line(p)=1;
        if (linje_p < min_end_1 || linje_p>max_end_1); los_line1=0; end
        for i=1:last
            % tjek om linjen går igennem åbning
            if linje_p< min_start_i || linje_p>max_start_i && i ~= n1 && i ~=n2
               %los blocked for this line
               lons_linep=0;
            end
        end
        if los_linep==1; los_logic=1; end;
        
    end
end