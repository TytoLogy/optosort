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

if iscell(msgstr)
   nchar = zeros(length(msgstr), 1);
   for n = 1:length(msgstr)
      nchar(n) = length(msgstr{n});
   end
   sepstr = repmat(sepchar, 1, max(nchar) + 2);
   
   fprintf('%s\n', sepstr);
   for n = 1:length(msgstr)
      fprintf('%s\n', msgstr{n});
   end
   fprintf('%s\n', sepstr);
         
else
   sepstr = repmat(sepchar, 1, length(msgstr) + 2);
   fprintf('%s\n%s\n%s\n', sepstr, msgstr, sepstr);

end

