% function varargout = getSpikesByStim(varargin)
%------------------------------------------------------------------------
% [data, datainfo, tracesByStim] = viewOptoData(varargin)
%------------------------------------------------------------------------
% TytoLogy:Experiments:opto Application
%--------------------------------------------------------------------------
%
%	should probably be a method within SpikeData object...
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



%----------------------------------------------------------------------
%% for script, need to load data
%----------------------------------------------------------------------
%------------------------------------------------------------------------
% sorted data locations
%------------------------------------------------------------------------
sortedPath = '/Users/sshanbhag/Work/Data/TestData/MT/1407';
rawPath = sortedPath;
nexPath = sortedPath;
nexInfoFile = '1407_20200309_03_01_1350_BBN_nexinfo.mat';
nexFile = '1407_20200309_03_01_1350_BBN.nex';
plxFile = '1407_20200309_03_01_1350_BBN-sorted.ch4,5,7,15.plx';

%------------------------------------------------------------------------
%% load S object from file
%------------------------------------------------------------------------
% use OptoFileName object to make this easier
OFobj = OptoFileName(fullfile(nexPath, nexFile));
sfile = [OFobj.fileWithoutOther '_Sobj.mat'];
fprintf('\n%s\n', sepstr);
fprintf('loading Sobj from file\n\t%s\n', fullfile(nexPath, sfile));
fprintf('%s\n', sepstr);
load(fullfile(nexPath, sfile));

%% get stim indices, varlist
fnum = 1;

[stimindex, stimvar] = S.Info.FileInfo{fnum}.getStimulusIndices;

% get the spikes for this file as a table
% spikeTable = S.spikesForAnalysis(fnum, 'Channel', 7, 'align', 'sweep');
spikeTable = S.spikesForAnalysis(fnum, 'Channel', 7, 'align', 'sweep');


%%


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