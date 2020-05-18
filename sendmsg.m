function sendmsg(msgstr)
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
fprintf('%s\n%s\n%s\n', sepstr, msgstr, sepstr);

