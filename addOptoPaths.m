%------------------------------------------------------------------------
% addOptoPaths.m
%------------------------------------------------------------------------
% TytoLogy:Experiments:optosort
%------------------------------------------------------------------------
% adds paths to opo related dirs
% 
%------------------------------------------------------------------------
% See also: 
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%  Sharad J. Shanbhag
%	sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 8 January, 2020 (SJS)
%
% Revisions:
%	10 Mar 2020 (SJS): streamlined code, added optosort and OptoObjects dirs
%------------------------------------------------------------------------
% TO DO: customize for each user's setup?
%------------------------------------------------------------------------


paths_to_add = { ...
	'~/Work/Code/Matlab/dev/TytoLogy/Experiments/Opto', ...
	'~/Work/Code/Matlab/dev/TytoLogy/Experiments/OptoAnalysis', ...
	'~/Work/Code/Matlab/dev/TytoLogy/Experiments/optosort', ...
	'~/Work/Code/Matlab/dev/TytoLogy/Experiments/optosort/OptoObjects' ...
};
	
fprintf('Adding path to opto functions:\n');

for p = 1:length(paths_to_add)
	fprintf('\t%s\n', paths_to_add{p});
	addpath(paths_to_add{p});
end
