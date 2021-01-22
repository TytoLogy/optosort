function sendmsg(msgstr, varargin)
%------------------------------------------------------------------------
% sendmsg(msgstr)
%------------------------------------------------------------------------
% displays message between two separation strings of dashes
%------------------------------------------------------------------------
% Input Arguments:
% 	msgstr		string to display
% Output Arguments:
% 	none
%------------------------------------------------------------------------
if nargin > 1
	sepchar = varargin{1};
else
	sepchar = '-';
end
% sepstr = '----------------------------------------------------';
sepstr = repmat(sepchar, 1, length(msgstr) + 2); 

fprintf('%s\n%s\n%s\n', sepstr, msgstr, sepstr);

