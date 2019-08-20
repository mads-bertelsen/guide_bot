function [ n_bins] = ifit_bins_1d( iData_object )
%ifit_max takes max of iFit data set
%   Detailed explanation goes here

if isfield(iData_object.Data,'array_1d') 
    n_bins = iData_object.Data.array_1d;
    
else
    n_bins = iData_object.Data.size(1);

end
end

