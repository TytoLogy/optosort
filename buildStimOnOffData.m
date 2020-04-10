% function [varargout] = buildStimOnOffData(Dinf)
%------------------------------------------------------------------------
% [cD, varargout] = buildChannelData(Dinf)
%------------------------------------------------------------------------
% TytoLogy
%------------------------------------------------------------------------
% 
% builds vector of stimulus onset and offset bins referenced to start
% of file
% 
%------------------------------------------------------------------------
% Input Arguments:
% 	Dinf			data information struct from opto program output file
%
% Output Arguments:
% 	startI	[# sweeps, 1] array with each element stimulus onset timestamp 
%					for each sweep
% 	endI		{# Channels, 1} cell array with each element being
% 					[# sweeps] vector holding end samples timestamp for each sweep
%	sweepLen	[# Channels, # sweeps] matrix of # of samples in each sweep
%------------------------------------------------------------------------
% See also: 
%------------------------------------------------------------------------


%------------------------------------------------------------------------
% Sharad J. Shanbhag
% sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 9 April, 2020 using some code from buildChannelData (SJS)
%
% Revisions:
%------------------------------------------------------------------------
% TO DO:
%	include option to have onset include/exclude stim delay???
%
%------------------------------------------------------------------------


%{

idea for algorithm:

loop through trials

get stimulus onset, offset info:
	onset = stimulus delay
	offset = stimulus delay + stimulus duration

use indev.Fs to convert to samples from milliseconds

To reference to absolute "time"/samples, add CurveInfo.startSweepBin(n) to
n'th onset and offset

some caveats:

will be different synthesized (tone, bbn) vs. wav stimuli, so will need to
handle those different cases - look at optoproc_plotPSTH_WAVbyLevel to get ideas!

10 Apr 2020 (SJS):
Roadmap:
 work here on more simple approach for BBN, freq tuning, etc. Then
 incorporate as method in CurveInfo object. Then override for use with WAV
 data in WAVInfo object
%}

%% read in test data
DataPath = '~/Work/Data/TestData/MT';
DataFile = '1382_20191212_02_02_3200_WAV.dat';

% use readOptoData to read in raw data. 
[D, tmpDinf] = readOptoData(fullfile(DataPath, DataFile));
% Fix test info
Dinf = correctTestType(tmpDinf);

%------------------------------------------------------------------------
%% settings
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% process data
%------------------------------------------------------------------------

% initialize startI and endI to store start and end sample bins
startI = cell(nChannelsToRead, 1);
endI = cell(nChannelsToRead, 1);
% samples in each sweep
sweepLen = zeros(nChannelsToRead, Dinf.nread);
% loop through channels
for c = 1:nChannelsToRead
	channel = channelIndx(c);
	% initialize startI and endI to store stop/start locations
	tmpStartI = zeros(1, Dinf.nread);
	tmpEndI = zeros(1, Dinf.nread);

	% loop through each sweep
	for s = 1:Dinf.nread
		% assign channel sweep data to cD after filtering
		% need to transpose to row vector
		cD{c, s} = D{s}.datatrace(:, channel)';
		if deOffset
			% remove initial offset
			cD{c, s} = cD{c, s} - cD{c, s}(1);
		end
		
		% filter data
		if filterData
			cD{c, s} = filtfilt(BPfilt.b, BPfilt.a, ...
										sin2array(cD{c, s}, ...
										Dinf.indev.Fs, BPfilt.ramp));
		end

		% store length of sweep
		sweepLen(c, s) = length(cD{c, s});
		
		% plot raw and filtered data
		if plotSweeps
			plot(D{s}.datatrace(:, channel)', 'k');
			hold on
				plot(cD{c, s}, 'b');
			hold off
			title(sprintf('Channel: %d  Sweep: %d(%d)', channel, s, Dinf.nread));
			drawnow
		end
		
		% build list of sweep start and end indices (in units of samples)
		% store index points
		if s ~= 1
			tmpStartI(s) = tmpEndI(s-1) + 1;
		else
			% start index is 1
			tmpStartI(s) = 1;
		end
		tmpEndI(s) = tmpStartI(s) + length(cD{c, s}) - 1;
	end
	% assign sweep indices to cell arrays
	startI{c} = tmpStartI;
	endI{c} = tmpEndI;
end

if DEBUG
	% loop through channels
	for c = 1:nChannelsToRead
		% loop through each sweep
		for s = 1:Dinf.nread
			% replace start and end value of each sweep with known value
			cD{c, s}(1) = exp(1);
			cD{c, s}(end) = -exp(1);
		end
	end
end

if nargout > 1
	varargout{1} = startI;
	varargout{2} = endI;
	varargout{3} = sweepLen;
end
