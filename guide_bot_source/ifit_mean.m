function [ mean_value ] = ifit_mean( iData_object )
%ifit_mean takes mean of iFit data set
%   Detailed explanation goes here


if isfield(iData_object.Data,'xlimits') 
    
    mean_value = mean(iData_object.signal);
    
elseif isfield(iData_object.Data,'xylimits') 
    
    mean_value = mean(mean(iData_object.signal));

end
end

