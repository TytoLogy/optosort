function varargout = assign_spikes_to_sweeps(allSpikes, ...
													sweepStartBins, sweepEndBins, Fs, ...
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
%
%------------------------------------------------------------------------
% TO DO:
%--------------------------------------------------------------------------

% process options
if isempty(varargin)
	ALIGN = 'original';
else
	if ~any(strcmpi(varargin{1}, {'original', 'file', 'sweep'}))
		error('%s: unknown align option %s', varargin{1})
	else
		ALIGN = lower(varargin{1});
	end;
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
		
% allocate cell to store spike info for each sweep
spikesBySweep = cell(nsweeps, 1);

% loop through sweeps
for s = 1:nsweeps
	% find the valid time stampes (between sweepStartTimes and
	% sweepEndTimes)
	valid_ts = (ts >= sweepStartTimes(s)) & (ts < sweepEndTimes(s));
	% apply offset correction
	valid_ts = valid_ts - alignval(s);
	spikesBySweep{s} = allSpikes(valid_ts, :);
end

% assign outputs
varargout{1} = spikesBySweep;
if nargout > 1
	varargout{2} = nsweeps;
	varargout{3} = {sweepStartTimes, sweepEndTimes};
end