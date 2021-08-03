function sendmsg(msgstr, varargin)
%------------------------------------------------------------------------
% sendmsg(msgstr, sepchar)
%------------------------------------------------------------------------
% displays message between two separation strings of dashes (or other
% character given by sepchar
%------------------------------------------------------------------------
% Input Arguments:
% 	msgstr		string to display
%  
%  Optional Inputs:
%   sepchar    alternative "separation" character (default is '-')
%
% Output Arguments:
% 	none
%------------------------------------------------------------------------

if nargin > 1
	sepchar = varargin{1};
else
	sepchar = '-';
end

% sepstr = '----------------------------------------------------';

if iscell(msgstr)
   len = zeros(size(msgstr));
   for n = 1:length(msgstr)
      len(n) = length(msgstr{n});
   end
   sepstr = repmat(sepchar, 1, max(len) + 2); 

   fprintf('%s\n', sepstr);
   for n = 1:length(msgstr)
      fprintf('%s\n', msgstr{n});
   end
   fprintf('%s\n', sepstr);   
else
   sepstr = repmat(sepchar, 1, length(msgstr) + 2); 
   fprintf('%s\n%s\n%s\n', sepstr, msgstr, sepstr);
end


