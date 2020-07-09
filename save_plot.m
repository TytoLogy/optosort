function varargout = save_plot(figH, outputFormat, outputPath)
%------------------------------------------------------------------------
% figH = save_plot(figH, outputFormat, outputPath)
%------------------------------------------------------------------------
% optosorts
%------------------------------------------------------------------------
% 
% saves figure pointed to by figH. outputFormat determines format,
% outputPath is location of output file
% 
% name of file will be taken from name of figure 
% so, if 'Name' property of figH os 'aplot' and saveFormat is 'PNG', output
% file will by aplot.png
% 
%------------------------------------------------------------------------
% Input Arguments:
% 	figH	figure handle
%	outputFormat	output format for plot.
% 				'PNG'		.png, 300 dpi resolution
% 				'PDF'		.pdf
% 				'FIG'		.fig (saves MATLAB figure)
%  outputpath	output directory
%
% Output Arguments:
%	figH	figure handle
%
%------------------------------------------------------------------------
% See also: SpikeData.plotAllData
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% Sharad J. Shanbhag
% sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created:8 July, 2020 (SJS)
%
% Revisions:
%------------------------------------------------------------------------
% TO DO: should probably be in utilities...
%------------------------------------------------------------------------

% check inputs
if nargin ~= 3
	error('save_plot: need 3 input args, fig handle, format and path');
end


if ~any(strcmpi(outputFormat, {'PNG', 'PDF', 'FIG'}))
	% if format is invalid, throw error
	error('save_plot: unknown save format %s', outputFormat);
end

% get figure name from handle and append to output path
pname = fullfile(outputPath, get(figH, 'Name'));

% save plot
switch upper(outputFormat)
	case 'FIG'
		savefig(figH, pname, 'compact');
	case 'PDF'
		print(figH, pname, '-dpdf');
	case 'PNG'
		print(figH, pname, '-dpng', '-r300');
end

if nargout
	varargout{1} = figH;
end



