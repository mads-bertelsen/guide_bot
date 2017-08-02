% This function takes the cell form of the options associated with a component constructor
% and converts it into a structure.

% Leland Harriger NCNR January 2017.


function StructOut = Parse_options(options)
for ops=1:length(options)
    command = options{ops};
    try % Assume command is numeric
        eval(['StructArg.' command ';']);
    catch % If not numeric then store command as a string;
        for n = 1:length(command)
            if strcmp(command(n),'=');
                equals=n;
            end
        end
        try
            eval(['StructArg.' command(1:equals-1) '=''' command(equals+1:end) ''';']);
        end
    end
end
StructOut = StructArg;
end