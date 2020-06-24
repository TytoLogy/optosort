function varargout = getSpikesByStim(obj, fIndx, chanNum, unitNum)
%------------------------------------------------------------------------
% [data, datainfo, tracesByStim] = viewOptoData(varargin)
%------------------------------------------------------------------------
% TytoLogy:opto
%--------------------------------------------------------------------------
%	method: SpikeData.getSpikesByStim
%--------------------------------------------------------------------------
%
%	
%------------------------------------------------------------------------
% Output Arguments:
%
%
%------------------------------------------------------------------------
% See Also:
%------------------------------------------------------------------------
%------------------------------------------------------------------------
%  Sharad Shanbhag
%   sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 15 June, 2020 (SJS)
%	- adapted from OptoAnalysis:getFilteredOptoData
%
% Revisions:
%	18 Jun 2020 (SJS):
%		- modifying to return values for a single file/test, channel and unit
%		- using test fake data for now
%		- moved into method for SpikeData
%------------------------------------------------------------------------
% TO DO:
%   *Documentation!
%--------------------------------------------------------------------------

%{
--------------------------------------------------------------------------
--------------------------------------------------------------------------
ALGORITHM

Table in obj.Spikes has SpikeTimes in sequential order, aligned to start of
entire file

obj.spikesForAnalysis will return spikes for a specific file (can also
limit to specific channel and/or unit). for now, get spiketimes organized
into sweeps and across all channels
** NOT ACROSS ALL CHANNELS!

for each test, will need:
	list of stimuli and indices
		these are available via the CurveInfo.getStimulusIndices method


output will be in format that can be passed to computeRLF, computeFRA, etc:
% assumes spikeTimes is in format:
% 		spikeTimes{nLevels, 1}
% 			spikeTimes{n} = {nTrials, 1}
% 				spikeTimes{n}{t} = [spike1_ms spike2_ms spike3ms ...
%
--------------------------------------------------------------------------
--------------------------------------------------------------------------
%}

%------------------------------------------------------------------------
% definitions
%------------------------------------------------------------------------
% string used in error messages
mname = 'SpikeData.getSpikesByStim';

%------------------------------------------------------------------------
% parse inputs
%------------------------------------------------------------------------
% need to use nargin == 4 since obj counts as an input in 
% object methods!!!!!
if nargin < 4
	error('%s: need file index #, channel and unit', mname);
end

% file index
if (fIndx < 1) || (fIndx > obj.Info.nFiles)
	error('%s: fIndx out of bounds', mname);
end

% channels
if length(chanNum) > 1
	error('%s: only works for single channel at a time', mname);
elseif obj.check_channels(chanNum) == -1
	error('%s: invalid channel %d', mname, chanNum);
end

% unit
if length(unitNum)  > 1
	error('%s: only works for single unit at a time', mname);
elseif ~obj.check_units(chanNum, unitNum)
	error('%s: invalid unit %d', mname, unitNum);
end	

%------------------------------------------------------------------------
% get spikes for desired file, channel and unit
%------------------------------------------------------------------------
spiketable = obj.spikesForAnalysis(fIndx, 'channel', chanNum, ...
											'Unit', unitNum, 'Align', 'Sweep');

%------------------------------------------------------------------------
% get stim indices, varlist
%------------------------------------------------------------------------
%{
% THIS NEEDS TO BE REWORKED crashes on FRA data
%
%}
%{
% stimindex is a cell array with each element (corresponding to a different
% stimulus level/parameter) consisting of a list of indices into each data
% sweep.
% stimvar is a list of the variables in the sweeps
[stimindex, stimvar] = obj.Info.FileInfo{fIndx}.getStimulusIndices;


unique_stim = unique(stimvar);
nstim = length(unique_stim);

%------------------------------------------------------------------------
% convert to spiketimes format
%------------------------------------------------------------------------
% 		spikeTimes{nLevels, 1}
% 			spikeTimes{n} = {nTrials, 1}
% 				spikeTimes{n}{t} = [spike1_ms spike2_ms spike3ms ...
%

spiketimes = cell(nstim, 1);
% loop through stimuli
for s = 1:nstim
	fprintf('stimvar(%d) = %d\n', s, unique_stim(s));
	% allocate spiketimes storage
	spiketimes{s} = cell(length(stimindex{s}), 1);
	% loop through sweeps (aka trials, reps) for this stimulus
	for r = 1:length(stimindex{s})
		% get the proper index into spikeTable for this stimulus and sweep
		% combination
		rIndx = stimindex{s}(r);
		% get table for currrent sweep
		tmpT = spiketable{rIndx};
		% assign spike timestamps to spikeTimes, ...
		% converting to milliseconds
		spiketimes{s}{r} = force_row(1000 * tmpT.TS);
	end
end

%------------------------------------------------------------------------
% create output struct
%------------------------------------------------------------------------
str.spiketimes = spiketimes;
str.stimindex = stimindex;
str.stimvar = stimvar;
str.unique_stim = unique_stim;
str.nstim = nstim;
str.spiketable = spiketable;
%}

str = obj.Info.FileInfo{fIndx}.convertSpikeTableToSpikeTimes(spiketable)

% add fields to str
str.fileIndex = fIndx;
str.channel = chanNum;
str.unit = unitNum;
str.fileName = obj.Info.FileInfo{fIndx}.F.file;
% assign output
varargout{1} = str;

%{

%------------------------------------------------------------------------
%------------------------------------------------------------------------
%% Old attempt at all channels/units
%{
% specify file
fnum = 1;

% get list of stimuli and indices into sweep arrays
[stimindex, stimvar] = obj.Info.FileInfo{fnum}.getStimulusIndices;
unique_stim = unique(stimvar);
nstim = length(unique_stim);

% get the spikes for this file as a table
% spikeTable = obj.spikesForAnalysis(fnum, 'Channel', 7, 'align', 'sweep');
% spikeTable = obj.spikesForAnalysis(fnum, 'align', 'sweep');
channelList = obj.listChannels;
nchan = length(channelList);
unitList = obj.listUnits;
nUnits = obj.nUnits;

dataByChannel = repmat(struct(	'dataByUnit', [], ...
											'channel', []), ...
								nchan, 1 );
							
for c = 1:nchan
	dataByChannel(c).channel = channelList(c);
	dataByChannel(c).dataByUnit = ...
						repmat( struct(	'spikeTable', {}, ...
												'spikeTimes', {}, ...
												'Unit', [] ), ...
									nUnits(c), 1);

% convert to proper format for computeRLF
% 		spikeTimes{nLevels, 1}
% 			spikeTimes{n} = {nTrials, 1}
% 				spikeTimes{n}{t} = [spike1_ms spike2_ms spike3ms ...
%

	for u = 1:nUnits((c))
		fprintf('getting data for Channel %d, Unit %d\n', ...
					channelList(c), unitList{c}(u));
		% get spike table for current channel and unit
		dataByChannel(c).dataByUnit(u).Unit = unitList{c}(u);
		dataByChannel(c).dataByUnit(u).spikeTable = ...
					obj.spikesForAnalysis(	fnum, ...
												'Channel', channelList(c), ...
												'Unit', unitList{c}(u), ...
												'align', 'sweep');

	end
end

%%
for c = 1:nchan
	% loop through stimuli
	for s = 1:nstim
		fprintf('level = %d\n', unique_stim(s));

		% allocate spiketimes storage
		dataByChannel(c).dataByUnitspikeTimes{s} = cell(length(stimindex{s}), 1);

		% loop through sweeps (aka trials, reps) for this stimulus
		for r = 1:length(stimindex{s})
			% get the proper index into spikeTable for this stimulus and sweep
			% combination
			rIndx = stimindex{s}(r);
			% get table for currrent sweep
			tmpT = dataByChannel(c).spikeTable{rIndx};
			% assign spike timestamps to spikeTimes, converting to milliseconds
			dataByChannel(c).spikeTimes{s}{r} = force_row(1000 * tmpT.TS);
		end
	end
end

%% plot data
Dinf = obj.Info.FileInfo{fnum}.Dinf;
analysisWindow = [Dinf.audio.Delay (Dinf.audio.Delay + Dinf.audio.Duration)];
for c = 1:nchan
	RLF = computeRLF(dataByChannel(c).spikeTimes, ...
													unique_stim, analysisWindow);
	hRLF = plotCurveAndCI(RLF, 'mean');
	drawnow
end

%}





%------------------------------------------------------------------------
%------------------------------------------------------------------------
%------------------------------------------------------------------------
%% code for getFilteredOptoData
%------------------------------------------------------------------------
%------------------------------------------------------------------------
%------------------------------------------------------------------------
%{
%----------------------------------------------------------------------
%% settings for processing data
%----------------------------------------------------------------------


%----------------------------------------------------------------------
% process inputs
%----------------------------------------------------------------------
if nargin
	[datapath, datafile, dataext] = fileparts(varargin{1});
	datafile = [datafile dataext];

	if nargin > 1
		argIndx = 2;
		while argIndx <= nargin
			switch upper(varargin{argIndx})
				case 'FILTER'
					HPFreq = varargin{argIndx+1}(1);
					LPFreq = varargin{argIndx+1}(2);
					argIndx = argIndx + 2;
				
				case {'CHANNEL', 'CHAN'}
					channelNumber = varargin{argIndx+1};
					if numel(channelNumber) ~= 1
						error('%s: Channel must be a single value', mfilename);
					end
					argIndx = argIndx + 2;
				
				otherwise
					error('Unknown input argument to %s: %s', mfilename, ...
																	varargin{argIndx});
			end
		end		
	end
end


%----------------------------------------------------------------------
%% Read Data
%----------------------------------------------------------------------
% read in raw data
[D, Dinf] = readOptoData(fullfile(datapath, datafile));
% read in test data (if it exists)
[~, fbase, ~] = fileparts(datafile);
testfile = [fbase '_testdata.mat'];
if exist(fullfile(datapath, testfile), 'file')
	load(fullfile(datapath, testfile), 'testdata');
else
	testdata = [];
end
% read in wav info (if it exists)
wavinfofile = [fbase '_wavinfo.mat'];
if exist(fullfile(datapath, wavinfofile), 'file')
	fprintf('\tReading wav info from %s\n', wavinfofile);
	load(fullfile(datapath, wavinfofile));
end

% Get test info
Dinf = correctTestType(Dinf);
fprintf('Test type: %s\n', Dinf.test.Type);

%----------------------------------------------------------------------
%% Some test-specific things...
%----------------------------------------------------------------------
% for FREQ test, find indices of stimuli with same frequency
switch upper(Dinf.test.Type)
	case 'FREQ'
		fprintf('\t%s test, finding indices\n', Dinf.test.Type);
		% list of frequencies, and # of freqs tested
		freqlist = cell2mat(Dinf.test.stimcache.FREQ);
		nfreqs = length(Dinf.test.stimcache.vrange);
		% locate where trials for each frequency are located in the
		% stimulus cache list - this will be used to pull out trials of
		% same frequency
		stimindex = cell(nfreqs, 1);
		for f = 1:nfreqs
			stimindex{f} = find(Dinf.test.stimcache.vrange(f) == freqlist);
		end
       
% for LEVEL test, find indices of stimuli with same level (dB SPL)
	case 'LEVEL'
		fprintf('\t%s test, finding indices\n', Dinf.test.Type);
		% list of levels, and # of levels tested
		levellist = Dinf.test.stimcache.LEVEL;
		nlevels = length(Dinf.test.stimcache.vrange);
		% locate where trials for each frequency are located in the
		% stimulus cache list - this will be used to pull out trials of
		% same frequency
		stimindex = cell(nlevels, 1);
		for l = 1:nlevels
			stimindex{l} = find(Dinf.test.stimcache.vrange(l) == levellist);
		end

% for FRA (FREQ+LEVEL) test, find indices of stimuli with
% freq and same level (dB SPL)
	case 'FREQ+LEVEL'
		fprintf('\t%s test, finding freq and level indices\n', Dinf.test.Type);
		% convert stimtype and curvetype to strings
		if isnumeric(Dinf.test.stimcache.stimtype)
			Dinf.test.stimcache.stimtype = char(Dinf.test.stimcache.stimtype);
		end
		if isnumeric(Dinf.test.stimcache.curvetype)
			Dinf.test.stimcache.curvetype = char(Dinf.test.stimcache.curvetype);
		end
		% if necessary, convert cells to matrices
		testcell = {'splval', 'rmsval', 'atten', 'FREQ', 'LEVEL'};
		for c = 1:length(testcell)
			if iscell(Dinf.test.stimcache.(testcell{c}))
				Dinf.test.stimcache.(testcell{c}) = ...
										cell2mat(Dinf.test.stimcache.(testcell{c}));
			end
		end
		% list of stimulus freqs, # of freqs tested
		freqlist = unique(Dinf.test.stimcache.FREQ, 'sorted');
		nfreqs = length(freqlist);
		% list of stimulus levels, # of levels tested
		levellist = unique(Dinf.test.stimcache.LEVEL, 'sorted');
		nlevels = length(levellist);
		%{
			Raw data are in a vector of length nstims, in order of
			presentation.

			values used for the two variables (Freq. and Level) are stored in
			vrange matrix, which is of length (nfreq X nlevel) and holds
			values as row 1 = freq, row 2 = level

			e.g. Dinf.test.stimcache.vrange(:, 1:5) = 
				4000  4000  4000  4000  4000
				0     10    20    30    40

			trialRandomSequence holds randomized list of indices into vrange,
			has dimensions of [nreps, ntrials]

			To sort the data for FRA: (1)	for each freq and level combination,
			locate the indices for that combination in the respective FREQ and
			LEVEL list. (2)	These indices can then be used within the D{}
			array
		%}
		stimindex = cell(nlevels, nfreqs);
		for f = 1:nfreqs
			for l = 1:nlevels
				currentF = freqlist(f);
				currentL = levellist(l);
				stimindex{l, f} = find( (Dinf.test.stimcache.FREQ == currentF) & ...
												(Dinf.test.stimcache.LEVEL == currentL) );
			end
		end
		

% for OPTO test...
    case 'OPTO'
  		fprintf('\t%s test, finding indices\n', Dinf.test.Type);
 
% for WavFile, need to find indices with same filename.
	case 'WAVFILE'
		fprintf('\t%s test, finding indices\n', Dinf.test.Type);
		% get list of stimuli (wav file names)
		nwavs = length(Dinf.stimList);
		wavlist = cell(nwavs, 1);
		stimindex = cell(nwavs, 1);
		for w = 1:nwavs
			stype = Dinf.stimList(w).audio.signal.Type;
			if strcmpi(stype, 'null')
				wavlist{w} = 'null';
			elseif strcmpi(stype, 'noise')
				wavlist{w} = 'BBN';
			elseif strcmpi(stype, 'wav')
				[~, wavlist{w}] = fileparts(Dinf.stimList(w).audio.signal.WavFile);
			else
				error('%s: unknown type %s', mfilename, stype);
			end
			stimindex{w} = find(Dinf.test.stimIndices == w);
		end
		Dinf.test.wavlist = wavlist;
	otherwise
		error('%s: unsupported test type %s', mfilename, Dinf.test.Type);
end

% find channel data
% get list of channels
if isfield(Dinf.channels, 'nRecordChannels')
	nchan = Dinf.channels.nRecordChannels; %#ok<NASGU>
	channelList = Dinf.channels.RecordChannelList;
else
	nchan = Dinf.channels.nInputChannels; %#ok<NASGU>
	channelList = Dinf.channels.InputChannels;
end
% make sure specified channel is in the list
channelIndex = find(channelList == channelNumber);
if isempty(channelIndex)
	error('%s: Channel %d not recorded', mfilename, channelNumber);
else
	fprintf('\tChannels:')
	fprintf('%d ', channelList);
	fprintf('\n');
end


%----------------------------------------------------------------------
%% Pull out trials, apply filter, store in matrix
%----------------------------------------------------------------------
% define filter for data
% sampling rate
Fs = Dinf.indev.Fs;
% build bandpass filter, store coefficients in filtB, filtA
fband = [HPFreq LPFreq] ./ (0.5 * Fs);
[filtB, filtA] = butter(5, fband);

fprintf('\tFiltering and sorting data\n');

% build tracesByStim; process will vary depending on stimulus type
if strcmpi(Dinf.test.Type, 'FREQ')
	tracesByStim = cell(nfreqs, 1);
	for f = 1:nfreqs
		dlist = stimindex{f};
		ntrials = length(dlist);
		tracesByStim{f} = zeros(length(D{1}.datatrace(:, 1)), ntrials);
		for n = 1:ntrials
			tracesByStim{f}(:, n) = filtfilt(filtB, filtA, ...
												D{dlist(n)}.datatrace(:, channelIndex));
		end
	end
end
if strcmpi(Dinf.test.Type, 'LEVEL')
	tracesByStim = cell(nlevels, 1);
	for l = 1:nlevels
		dlist = stimindex{l};
		ntrials = length(dlist);
		tracesByStim{l} = zeros(length(D{1}.datatrace(:, 1)), ntrials);
		for n = 1:ntrials
			tracesByStim{l}(:, n) = filtfilt(filtB, filtA, ...
													D{dlist(n)}.datatrace(:, channelIndex));
		end
	end
end
if strcmpi(Dinf.test.Type, 'FREQ+LEVEL')
	% allocate cell array for traces, levels in rows, freqs in cols
	tracesByStim = cell(nlevels, nfreqs);
	% loop through freqs and levels
	for f = 1:nfreqs
		for l = 1:nlevels
			% get list of indices for this level + freq combination
			dlist = stimindex{l, f};
			% how many trials?
			ntrials = length(dlist);
			if ntrials ~= Dinf.test.Reps
				wstr = [ mfilename ': ' ...
							sprintf('only %d of desired ', ntrials) ...
							sprintf('%d reps', Dinf.test.Reps) ...
							sprintf('for F(%d) = %d, ', f, freqlist(f)) ...
							sprintf('L(%d) = %d', l, levellist(l)) ];
				warning(wstr);
			end
			tracesByStim{l, f} = zeros(length(D{1}.datatrace(:, 1)), ntrials);
			for n = 1:ntrials
				tracesByStim{l, f}(:, n) = ...
							filtfilt(filtB, filtA, ...
											D{dlist(n)}.datatrace(:, channelIndex));
			end
		end
	end
end
if strcmpi(Dinf.test.Type, 'OPTO')
	tracesByStim = {};
end
if strcmpi(Dinf.test.Type, 'WavFile')
	tracesByStim = cell(nwavs, 1);
	for w = 1:nwavs
		% create temporary array to hold data
		tracesByStim{w} = zeros(length(D{1}.datatrace(:, 1)), Dinf.test.Reps);
		for n = 1:Dinf.test.Reps
			dIndx = stimindex{w}(n);
			tracesByStim{w}(:, n) = filtfilt(filtB, filtA, ...
														D{dIndx}.datatrace(:, channelIndex));
		end
	end
end

%----------------------------------------------------------------------
%% assign outputs
%----------------------------------------------------------------------
varargout{1} = D;
varargout{2} = Dinf;
varargout{3} = tracesByStim;

% if no testdata, get outta here
if isempty(testdata)
	return
elseif ~isfield(testdata, 'SpikeTimes')
	return
else
	varargout{4} = testdata;
end

%}
	
	
	
%}
	