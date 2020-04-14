function checkstatus = check_sweeps(sweepbins)
%------------------------------------------------------------------------
% checkstatus = check_sweeps(sweepbins)
%------------------------------------------------------------------------
% TytoLogy:optosort
%------------------------------------------------------------------------
% Check consistency of sweep bins across channels
%------------------------------------------------------------------------
% Input Arguments:
% 	sweepbins		{nchannels, 1} cell array of vectors, where 
% 						each element is an array [1, nsweeps] of sweep bins
% Output Arguments:
% 	checkstatus		true if sweep times are inconsistent, false if all
%						are consistent (equal)
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%  Sharad J. Shanbhag
%	sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 14 January, 2020 (SJS)
%
% Revisions:
%	14 Apr, 2020 (SJS): moved out of export_for_plexon into separate file
%------------------------------------------------------------------------
% TO DO:
%------------------------------------------------------------------------


tmp = cell2mat(sweepbins);
[nr, ~] = size(tmp);
tmp2 = 0;
for r = 1:nr
	tmp2 = tmp2 + sum(tmp(r, :) - tmp(1, :));
end

if tmp2 ~= 0
	checkstatus = true;
else
	checkstatus = false;
end
