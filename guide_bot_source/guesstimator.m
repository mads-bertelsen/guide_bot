function guesstimate = guesstimator(globalinfo,modules,modulelist)

% vector with the indexes of modules which break line of sight
list_los_breakers = [3];

% vector with indexes of modules which should be short ( gap, kink )
list_short = [2 3];

% vector with indexes of transport guide modules (can be long) (elliptic)
list_transporters = [4];

% vector with indexes of guide modules which are good for changing the
% cross section of the guide quickly. (Straight, parabolic)
list_optics = [1 5];


% The intention is to seperate the guide into parts. A part is a section of
% guide between fixed points or kinks and moderator / sample.

% PM(fixed)ESKSE - P ES SE

% For each part the guesses of length should be done from the above
% definitions of the characteristics of each module.

% P No length guesses needed
% ES Guess long ellipse short straight section
% SE Guess long ellipse short straight section
% M (if not set) guess very short, 10 cm or so
% K Guess the position of the kink to be in the center of the kinked part
% and the length of the kink to extremely short, 10-30cm

% Dimension guesses (start)
% P Start, no guesses needed other than guide start which is difficult
% E Guess narrow in horizontal because of kink, large in vertical
% S Guess large start (small end hopefully)
% K no guess neeeded, right? (check)
% S Guess narrow because of kink horizontal, vertical large
% E Guess large, middle of transporting guide

% Version for creating seperate parts based on fixed and splitters
% Algorithm sketch:

% Break into parts:
% Find splitters [M K]
% Find fixed (if they coincide with above, ignore, if not, add)

% for modules
%       if splitter 
%           Add to splitter list
%       elseif fixed
%           Add to fixed list
%       end
% end

% part_index = 1
% for indexer = modules
%    if indexer on splitter
%       part_index++;part_length=0
%    elseif indexer on fixed
%       part_index++;part_length=1
%       part{part_index} = modules(indexer)
%    else
%       if part_length=0
%           part{part_index} = modules(indexer);
%       else
%           part{part_index} = [part{indexer} modules(indexer)]
%       end
%    end
% end

% Version based on just splitting on splitters
% for modules
%       if splitter
%           add to splitter list
%       end
% end

% make parts of modules between the splitters
% add moderator (1) and length(modules)+1 sample to splitters
% for 1:splitters-1
%    part{i} = modules(splitters(i):splitters(i+1))
% end

% Each part must now be examined for "intent" in the form of length and
% start point options including maximum and minimum, as the guess must obey
% these requirements.

% for i = 1:parts
%   for j = 1:length(part{i})
%       if all maxlength
%       part maxlength = sum maxlength + length
%       end
%       part minlength = sum minlength + length
%
%       contained startpoints = startpoints(splitter(i):splitter(i+1))
%   end
% end

% convert minstartpoint and maxstartpoint to fixedstart point
% 3 cases: both max and min, max, min
% both max and min: take the average
% only min: naieve answer would be mean(min,end)

% if no fixed split the n parts into n equally long guides

% if guide == start of part || guide == end of part
%       if los_breaker
%          narrow in kink direction
%          large in other direction3
%       else
%          medium horizontal and vertical
%       end
% else
%       large horizontal and vertical
% end

% guestimate. structure

end