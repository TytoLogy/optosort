function varargout = assign_spikes_to_sweeps(ts, allSpikes, ...
													sweepStartBins, sweepEndBins, Fs)
%------------------------------------------------------------------------
% 
%------------------------------------------------------------------------
% % TytoLogy:Experiments:optosort
%--------------------------------------------------------------------------
% make no assumptions about location of timestamps in allSpikes matrix, so
% get the timestamps (in seconds!) as a separate argument.
%
%------------------------------------------------------------------------
% Input Arguments:
%
%	allSpikes	matrix of sorted spike information from Plexon Offline
%					Sorter
%		Column 1: channel (?)
%		Column 2: unit #
%		Column 3: timestamp (in seconds)
%		Column 4: PCA1 weight
%		Column 5: PCA2 weight
%		Column 6: PCA3 weight
%		Column 7-end : waveform
%
% Output Arguments:
%
%
%------------------------------------------------------------------------
% See Also: import_from_plexon
%------------------------------------------------------------------------
%------------------------------------------------------------------------
%  Sharad Shanbhag
%   sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 10 June, 2016 (SJS)
%           - adapted from readHPData.m
%
% Revisions:
%	10 Jul 2017 (SJS): some minor tweaks
%	10 Oct 2017 (SJS): for some reason, this was in OptoAnalysis dir; 
%							 relocated to opto program dir
%	4 Feb 2019 (SJS): added FREQ+LEVEL and OPTO, not fully implemented for
%							finding tracesByStim...
%------------------------------------------------------------------------
% TO DO:
%   *Documentation!
%--------------------------------------------------------------------------


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

% allocate cell to store spike info for each sweep
spikesBySweep = cell(nsweeps, 1);

% loop through sweeps
for s = 1:nsweeps
	valid_ts = (ts >= sweepStartTimes(s)) & (ts < sweepEndTimes(s));
	spikesBySweep{s} = allSpikes(valid_ts, :);
end


varargout{1} = spikesBySweep;
