function varargout = assign_spikes_to_sweeps(allSpikes, ...
													sweepStartBins, sweepEndBins, ...
													unit, Fs, ...
													varargin)
%------------------------------------------------------------------------
% 
%------------------------------------------------------------------------
% % TytoLogy:Experiments:optosort
%--------------------------------------------------------------------------
% 
%
%------------------------------------------------------------------------
% Input Arguments:
%	allSpikes	matrix of sorted spike information from Plexon Offline
%					Sorter
%		Column 1: channel (?)
%		Column 2: unit #
%		Column 3: timestamp (in seconds)
%		Column 4: PCA1 weight
%		Column 5: PCA2 weight
%		Column 6: PCA3 weight
%		Column 7-end : waveform	sweepStartBins	vector of sweep start samples (bins)
% 	sweepEndBins	vector of sweep end samples (bins)
% 		so, each sweep is defined by
% 			[sweepStartBins(1)  sweepEndBins(1)
% 			 sweepStartBins(2)  sweepEndBins(2)
% 			 .
% 			 .
% 			 .
% 			 sweepStartBins(n)  sweepEndBins(n) ]
%
%	unit
%		unit ID to return.  if empty, all units will be returned.
% 	Fs	sample rate for spike data (samples/second)
%	
%	Optional:
%	 <align_mode>	realign spike timestamps
% 			where align_mode is:
% 				'original'	leave as is (default)
% 				'file'		align to start time of file (first sweep)
% 				'sweep'		align to start of each sweep (useful for histograms,
% 								analysis, etc)
%
% Output Arguments:
%	spikesBySweep	cell vector of spike info per sweep
%
%------------------------------------------------------------------------
% See Also: import_from_plexon
%------------------------------------------------------------------------
%------------------------------------------------------------------------
%  Sharad Shanbhag
%   sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 11 February, 2020 (SJS)
%           - pulled code out from import_from_plexon.m
%				- removed ts as input
% Revisions:
%	14 Feb 2020 (SJS): added channel
%	deprecated - subsumed algorithm into SpikeData object methods
%------------------------------------------------------------------------
% TO DO:
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% some definitions
%--------------------------------------------------------------------------
[~, nc] = size(allSpikes);
CHAN_COL = 1; %#ok<NASGU>
UNIT_COL = 2; 
TS_COL = 3;
PCA_COL = 4:6; %#ok<NASGU>
WAV_COL = 7:nc; %#ok<NASGU>

%--------------------------------------------------------------------------
% process options and inputs
%--------------------------------------------------------------------------
if isempty(varargin)
	ALIGN = 'original';
else
	if ~any(strcmpi(varargin{1}, {'original', 'file', 'sweep'}))
		error('%s: unknown align option %s', varargin{1})
	else
		ALIGN = lower(varargin{1});
	end
end

% get # of sweeps
nsweeps = length(sweepStartBins);
% some checks
if nsweeps ~= length(sweepEndBins)
	error('%s: mismatch in sweepStartBins and sweepEndBins', mfilename);
elseif nsweeps == 0
	error('%s: no sweeps!', mfilename);	
end

% convert bins to start and end times for the sweeps
sweepStartTimes = (sweepStartBins - 1) * (1/Fs);
sweepEndTimes = (sweepEndBins - 1) * (1/Fs);

% alignment issues
switch ALIGN
	case 'original'
		% don't modify time stamp values
		alignval = zeros(nsweeps, 1);
	case 'file'
		% adjust timestamp value by first sweep start time
		alignval = sweepStartTimes(1) * ones(nsweeps, 1);
	case 'sweep'
		% adjust timestamp value by each sweep start time
		alignval = sweepStartTimes;
end

%--------------------------------------------------------------------------
% process data
%--------------------------------------------------------------------------

% allocate cell to store spike info for each sweep
spikesBySweep = cell(nsweeps, 1);
% isolate unit
if ~isempty(unit)
	allSpikes = allSpikes(allSpikes(:, UNIT_COL) == unit, :);
end
% local copy of ts
ts = allSpikes(:, TS_COL);

% loop through sweeps
for s = 1:nsweeps
	% find the valid time stampes (between sweepStartTimes and
	% sweepEndTimes)
	valid_ts = (ts >= sweepStartTimes(s)) & (ts < sweepEndTimes(s));
	spikesBySweep{s} = allSpikes(valid_ts, :);
	% apply offset correction
	spikesBySweep{s}(:, TS_COL) = spikesBySweep{s}(:, TS_COL) - alignval(s);
end

% assign outputs
varargout{1} = spikesBySweep;
if nargout > 1
	varargout{2} = nsweeps;
	varargout{3} = {sweepStartTimes, sweepEndTimes};
end