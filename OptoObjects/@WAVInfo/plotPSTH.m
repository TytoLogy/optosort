function H = plotPSTH(obj, ST, binSize, varargin)
%------------------------------------------------------------------------
% opto project
%------------------------------------------------------------------------
% WAVInfo.plotPSTH method
%------------------------------------------------------------------------
% interface between WAVInfo class and optoproc_plotPSTH_WAVbyLevel or
% optoproc_plotPSTH_WAV functions from OptoAnalysis
%------------------------------------------------------------------------
% Input Args
% ST	struct (usually returned by SpikeData.getSpikesByStim method)
%
% ST struct fields:
% FREQ_TUNING
%      spiketimes: {14×1 cell}
%       stimindex: {14×1 cell}
%         stimvar: [140×1 double]
%     unique_stim: [14×1 double]
%           nstim: 14
%      spiketable: {140×1 cell}
%       fileIndex: 2
%         channel: 5
%            unit: 1
%        fileName: '1407_20200309_03_01_1350_FREQ_TUNING.dat'
% FRA
%   struct with fields:
% 
%      spiketimes: {9×14 cell}
%       stimindex: {9×14 cell}
%         stimvar: {[14×1 double]  [9×1 double]}
%     unique_stim: {2×1 cell}
%           nstim: [2×1 double]
%      spiketable: {2520×1 cell}
%       fileIndex: 4
%         channel: 5
%            unit: 1
%        fileName: '1407_20200309_03_01_1350_FRA.dat'
%
% binSize	psth bin size, milliseconds
%
% optional input: 
%	'LEVEL'	separate plot pages by stimulus level (this is default)
%
% Output Args:
% H	cell array of plot handles
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%  Sharad J. Shanbhag
%	sshadnbhag@neomed.edu
%------------------------------------------------------------------------
% Created: June, 2020 (SJS)
% Revisions:
%------------------------------------------------------------------------
% TO DO:
%------------------------------------------------------------------------

%
% check if plots are to be separated by stimulus level
if ~isempty(varargin)
	switch(upper(varargin{1}))
		case {'LEVEL', 'BYLEVEL', 'PLOT_PSTH_BY_LEVEL'}
			plotPSTH_BY_LEVEL = true;
		otherwise
			warning('%s: unknown option %s', mfilename, varargin{1})
			plotPSTH_BY_LEVEL = false;
	end
else
	plotPSTH_BY_LEVEL = false;
end


% NOT NEEDED..... ???? 7 jul 2020
%{
% determine rows, cols for plots
if numel(ST.nstim) == 1
	% for data that are not "2D" (e.g., FRA), adjust # of columns based
	% on the number of variable levels or types (nvars)
	if ST.nstim <= 6
		prows = ST.nstim;
		pcols = 1;
	elseif iseven(ST.nstim)
		prows = ST.nstim/2;
		pcols = 2;
	else
		prows = ceil(ST.nstim/2);
		pcols = 2;
	end
else
	% otherwise, plot levels in rows, freqs in columns
	% rows = levels, cols = freqs
	prows = ST.nstim(2);
	pcols = ST.nstim(1);
end
%}

% get time limits - needed for input to optoproc_plotPSTH_WAVxxxxx
% functions 
timeLimits = [0 obj.Dinf.test.AcqDuration];


% get plot title string(s)
% titleStr = obj.getTitleString
% hack to have wavlist in Dinf
% might want to move this to constructor method of WAVInfo...?>???
if ~isfield(obj.Dinf.test, 'wavlist')
	obj.Dinf.test.wavlist = obj.getwavList;
elseif isempty(obj.Dinf.test.wavlist)
	obj.Dinf.test.wavlist = obj.getwavList;
end

% get stimulusTimes {onset, offset} by stimulus - need to do this AFTER
% getting obj.Dinf.test.wavlist values
stimulusTimes = obj.getstimulusTimes;


% plot using the optoproc_plotPSTH_WAVxxxxx functions
if plotPSTH_BY_LEVEL
	H = optoproc_plotPSTH_WAVbyLevel(ST.spiketimes, obj.Dinf, binSize, ...
		[], timeLimits, [], stimulusTimes);
% 	H = optoproc_plotPSTH_WAVbyLevel(ST.spiketimes, obj.Dinf, binSize, ...
% 		[prows pcols], timeLimits, []);
else
	H = optoproc_plotPSTH_byWAV(ST.spiketimes, obj.Dinf, binSize, ...
									timeLimits, []);
end

