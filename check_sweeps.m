function checkstatus = check_sweeps(sweepbins)
%------------------------------------------------------------------------
% checkstatus = check_sweeps(sweepbins)
%------------------------------------------------------------------------
% <project>:<subproject>:<function_name>
%------------------------------------------------------------------------
% 
% Description
% 
%------------------------------------------------------------------------
% Input Arguments:
% 	sweepbins		{nchannels, 1} cell array of vectors, where 
% 						each element is an array [1, nsweeps] of sweep bins
% 
% Output Arguments:
% 	checkstatus		true if sweep times are inconsistent, false if all
%						are consistent (equal)
%
%------------------------------------------------------------------------
% See also: 
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% Sharad J. Shanbhag
% sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 13 January, 2020 (SJS)
%
% Revisions:
%------------------------------------------------------------------------
% TO DO:
%------------------------------------------------------------------------

tmp = cell2mat(sweepbins);

if sum((sum(tmp - tmp(1, :)))) ~= 0
	checkstatus = true;
else
	checkstatus = false;
end
