function [ output ] = str2num_safe(input)
%str2num_safe str2num, but checks if the input is a value first

if isnumeric(input) == 1
    output = input;
else
    output = str2num(input);
end 

end

