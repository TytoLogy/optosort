function [obj, varargout] = buildChannelData(obj, Channels, BPfilt, D, varargin)
%------------------------------------------------------------------------
% [obj, varargout] = CurveInfo.buildChannelData(Channels, BPfilt, D, varargin)
%------------------------------------------------------------------------
% TytoLogy:optosort:CurveInfo Object method
%------------------------------------------------------------------------
% 
% get data for each channel
% 
%------------------------------------------------------------------------
% Input Arguments:
% 	Channels		list of A/D (spike) channels to export
%	BPfilt		struct specifying bandpass filter settings
% 					if BPfilt is empty ([]), no filtering will be performed
% 			required fields:
% 			forder	filter order
% 			ramp		onset/offset signal ramp in milliseconds
% 			Fs			signal sampling rate in samples/sec
% 			b			bandpass filter coefficients (e.g., from butter())
% 			a			bandpass filter coefficients (e.g., from butter())
% 	D				raw data from opto program output file
% 	Dinf			data information struct from opto program output file
% 	
% 	Optional Input Args:
% 		'REMOVE_OFFSET'	<'Y'/'N'>	remove DC offset from sweeps?
% 		'PLOT_SWEEPS'		<'N'/'Y'>	plot sweeps?
%
% Output Arguments:
%	obj		updated copy of object
% 	cD			{# Channels, # sweeps} cell array of channel sweep data
% 				each element of cD will be a row vector of raw data for a sweep
% 	startI	{# Channels, 1} cell array with each element being
% 					[# sweeps] vector holding start sample timestamp for each sweep
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
% Created: 10 January, 2020 (SJS)
%
% Revisions:
%	13 Jan 2020 (SJS): added comments , optional args
%	13 Apr 2020 (SJS): moved into CurveInfo class
%	14 Apr 2020 (SJS): need to return obj...
%------------------------------------------------------------------------
% TO DO: 
% - probably should have filter as an option vs. default input arg
% - make Channels an option and have, by default, all channels returned??
%------------------------------------------------------------------------

% original call (from export_for_plexon):
%{
	% build into sweeps by channel format
	[cSweeps{f}, ...
		cInfo{f}.startSweepBin, cInfo{f}.endSweepBin, cInfo{f}.sweepLen] = ...
					buildChannelData(Channels, BPfilt, D, cInfo{f}.Dinf);
%}

%------------------------------------------------------------------------
% process options
%------------------------------------------------------------------------
% initialize options to default settings
deOffset = true;
plotSweeps = false;
filterData = true;
DEBUG = false;

if ~isempty(varargin)
	argIndx = 1;
	while argIndx <= length(varargin)
		fprintf('%s\n', upper(varargin{argIndx}))
		switch upper(varargin{argIndx})
			case {'REMOVE_OFFSET'}
				tmpArg = varargin{argIndx + 1};
				if strcmpi(tmpArg(1), 'Y')
					deOffset = true;
				else
					deOffset = false;
				end
				argIndx = argIndx + 2;
			case {'PLOT_SWEEPS'}
				tmpArg = varargin{argIndx + 1};
				if strcmpi(tmpArg(1), 'Y')
					plotSweeps = true;
				else
					plotSweeps = false;
				end
				argIndx = argIndx + 2;
			case {'DEBUG'}
				tmpArg = varargin{argIndx + 1};
				if strcmpi(tmpArg(1), 'Y')
					DEBUG = true;
				else
					DEBUG = false;
				end
				argIndx = argIndx + 2;
				
			otherwise
				error('%s: unknown input arg %s', mfilename, varargin{argIndx});
		end
	end
end

%------------------------------------------------------------------------
% settings
%------------------------------------------------------------------------
% if BPfilt is empty, don't filter data
if isempty(BPfilt)
	filterData = false;
end

% channel information
nChannelsToRead = length(Channels);
% check if list of channels to read is longer than # recorded
if nChannelsToRead > obj.Dinf.channels.nRecordChannels
	error('%s: mismatch in nChannelsToRead (%d) and nRecordChannels (%d)', ...
				mfilename, nChannelsToRead, obj.Dinf.channels.nRecordChannels);
end
% make sure requested channel was recorded
for c = 1:nChannelsToRead
	if ~any(Channels(c) == obj.Dinf.channels.RecordChannelList)
		error('%s: Channel %d not in RecordChannelList!', mfilename, Channels(c));
	end
end
% build channel index
if obj.Dinf.channels.nRecordChannels == 1
	% only one channel!
	channelIndx = 1;
	if Channels ~= obj.Dinf.channels.RecordChannelList
		warning('%s: requested channel %d not found in RecordChannelList', ...
						mfilename, Channels);
		fprintf('Using only available channel %d\n', ...
						obj.Dinf.channels.RecordChannelList);
		channelIndx = 1;
	end
elseif obj.Dinf.channels.nRecordChannels == 16
	channelIndx = Channels;
else
	channelIndx = zeros(nChannelsToRead);
	for c = 1:nChannelsToRead
		channelIndx(c) = find(Channels(c) == obj.Dinf.channels.RecordChannelList);
	end
end

%------------------------------------------------------------------------
% process data: get start and end of sweeps
%------------------------------------------------------------------------
% initialize cD to a store sweeps for each channel
cD = cell(nChannelsToRead, obj.Dinf.nread);
% initialize startI and endI to store start and end sample bins
obj.startSweepBin = cell(nChannelsToRead, 1);
obj.endSweepBin = cell(nChannelsToRead, 1);
% samples in each sweep
obj.sweepLen = zeros(nChannelsToRead, obj.Dinf.nread);
% loop through channels
for c = 1:nChannelsToRead
	channel = channelIndx(c);
	% initialize startI and endI to store stop/start locations
	tmpStartI = zeros(1, obj.Dinf.nread);
	tmpEndI = zeros(1, obj.Dinf.nread);

	% loop through each sweep
	for s = 1:obj.Dinf.nread
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
										obj.Dinf.indev.Fs, BPfilt.ramp));
		end

		% store length of sweep
		obj.sweepLen(c, s) = length(cD{c, s});
		
		% plot raw and filtered data
		if plotSweeps
			plot(D{s}.datatrace(:, channel)', 'k');
			hold on
				plot(cD{c, s}, 'b');
			hold off
			title(sprintf('Channel: %d  Sweep: %d(%d)', channel, s, obj.Dinf.nread));
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
	obj.startSweepBin{c} = tmpStartI;
	obj.endSweepBin{c} = tmpEndI;
end

% check the start and end sweep bin data for consistency
if check_sweeps(obj.startSweepBin)
	warning(['File %s: Inconsistent startSweepBin' ...
						'values across channels!!!!'], obj.F.file);
end
if check_sweeps(obj.endSweepBin)
	warning(['File %s: Inconsistent endSweepBin' ...
						'values across channels!!!!'], obj.F.file);
end

%------------------------------------------------------------------------
% if DEBUG specified, write a known value (+/- exp(1)) to the start and
% beginning to each sweep
%------------------------------------------------------------------------
if DEBUG
	% loop through channels
	for c = 1:nChannelsToRead
		% loop through each sweep
		for s = 1:obj.Dinf.nread
			% replace start and end value of each sweep with known value
			cD{c, s}(1) = exp(1);
			cD{c, s}(end) = -exp(1);
		end
	end
end

%------------------------------------------------------------------------
% assign output args
%------------------------------------------------------------------------
if nargout > 1
	varargout{1} = cD;
	varargout{1} = obj.startSweepBin;
	varargout{2} = obj.endSweepBin;
	varargout{3} = obj.sweepLen;
end
