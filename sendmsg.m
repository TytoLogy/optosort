function sendmsg(varargin)
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
sepstr = '----------------------------------------------------';

fprintf('%s\n', sepstr);
if nargin
	fprintf('%s\n', varargin{1});
end
fprintf('%s\n', sepstr);

