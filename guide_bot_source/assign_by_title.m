function [ output_iData ] = assign_by_title(title_string,iData_array)
%ASSIGN_BY_TITLE finds an object in a iData_array of objects with matching title
%   Searches for the title_string to match the first part of the title
%   Will also fix axis if they are missing
title_length = length(title_string);
found = 0;

for index=1:length(iData_array)
    
    if length(iData_array(index).title) > title_length
      if (strcmp(title_string,iData_array(index).title(1:title_length)))
       output_iData = iData_array(index);
       found = 1;
       break;
      end
    end
end
if found == 0
disp(['Did not find a matching object with the title ' title_string ])
end

if isfield(output_iData.Data,'xlimits') 
   % Check axis are assigned correctly for 1 dimensional monitor 
   if isempty(getaxis(output_iData))
       output_signal_size = size(output_iData.Signal);
       output_iData{1} = linspace(output_iData.Data.xlimits(1),output_iData.Data.xlimits(2),output_signal_size(1));
   end
   
   % Attempt to fix differences in versions of McStas, as later versions of
   % McStas does not return an x axis when monitor files are read
%    if isfield(output_iData,'x') == 0
%        %output_iData.x = output_iData{1};
%        setalias(output_iData,'x','this{1}');
%    end
   
elseif isfield(output_iData.Data,'xylimits') 
   % Check axis are assigned correctly for 2 dimensional monitor
   if isempty(getaxis(output_iData))
       output_signal_size = size(output_iData.Signal);
       output_iData{1} = linspace(output_iData.Data.xylimits(3),output_iData.Data.xylimits(4),output_signal_size(2));
       output_iData{2} = linspace(output_iData.Data.xylimits(1),output_iData.Data.xylimits(2),output_signal_size(1));
   end
   
%    if isfield(output_iData,'x') == 0
%        %output_iData.x = output_iData{1};
%        setalias(output_iData,'x','this{2}');
%    end
%    if isfield(output_iData,'y') == 0
%        %output_iData.x = output_iData{1};
%        setalias(output_iData,'y','this{1}');
%    end
   
end



end

