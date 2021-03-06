% function varargout = import_from_plexon_nonObj(varargin)
%------------------------------------------------------------------------
% import_from_plexon.m
%------------------------------------------------------------------------
% TytoLogy:Experiments:optosort
%------------------------------------------------------------------------
% working script for importing sorted data from plexon
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
%------------------------------------------------------------------------
% TO DO:
%------------------------------------------------------------------------
% Initial things to define
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% Definitions
%------------------------------------------------------------------------
% string used to separate text 
sepstr = '----------------------------------------------------';

% data locations
sortedPath = '~/Work/Data/TestData/MT';
rawPath = '~/Work/Data/TestData/MT';
nexPath = '~/Work/Data/TestData/MT';

% sorted data file
% sortedFile = '1323_20190722_03_02_632_MERGE.mat';
sortedFile = '1382_20191212_02_02_3200.mat';

% nexinfo file
% nexInfoFile = '1323_20190722_03_02_632_MERGE_nexinfo.mat';
nexInfoFile = '1382_20191212_02_02_3200_MERGE_nexinfo.mat';

% nex file
nexFile = '1382_20191212_02_02_3200_MERGE.mat';

%------------------------------------------------------------------------
% Setup
%------------------------------------------------------------------------
fprintf('\n%s\n', sepstr);
fprintf('import_from_plexon running...\n');
fprintf('%s\n', sepstr);


%------------------------------------------------------------------------
%% load data
%------------------------------------------------------------------------
fprintf('\n%s\n', sepstr);
fprintf('Loading nexInfo from file\n\t%s\n', fullfile(nexPath, nexInfoFile));
fprintf('%s\n', sepstr);
% nexinfo
load(fullfile(nexPath, nexInfoFile), 'nexInfo');

% plexon sorted data
plxvars = who('-file', fullfile(sortedPath, sortedFile));
if isempty(plxvars)
	error('No variables in plexon output file %s', ...
					fullfile(sortedPath, sortedFile));
else
	fprintf('\n%s\n', sepstr);
	fprintf('Found %d channels in %s\n', length(plxvars), ...
							fullfile(sortedPath, sortedFile))
	for n = 1:length(plxvars)
		fprintf('\t%s\n', plxvars{n});
	end
	fprintf('%s\n', sepstr);
end

% for now, just load the first channel data - need to figure out a way to
% deal with multiple channels
tmp = load(fullfile(sortedPath, sortedFile), '-MAT', plxvars{1});
spikesAll = tmp.(plxvars{1});
clear tmp;

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% Plexon-exported sorted data:
%------------------------------------------------------------------------
%------------------------------------------------------------------------
% Column 1: unit number (where 0 is unsorted)
% Column 2: timestamp where spike crosses threshold (in seconds)
% Columns 3-34 (assuming waveform window of 1311us / 32 samples):
%   waveform snippet, with or without prewindow as set in Offline Sorter
%   (prewindow default: 494us / 12 samples)
%   (window default: 1311us / 32 samples)
%   This is in units of samples/sec of raw data file 
%   (24414.063 Hz based on settings in data acquisition program 
%        HPSearch or PresentStimCurve in RosenLab)

% may not be accurate.... might be
% Column 1: channel (?)
% Column 2: unit #
% Column 3: timestamp (in seconds)
% Column 4-... : waveform
%
% 22 Jan 2019
% if waveforms are clipped, look at pp 109 and 118 in offline sorter manual
% to fix gain
% 
% other issue: transients on first 1-4 samples of each waveform.... ????
% just noise????
%
% 24 Jan 19: cols 4-6 are weights of 3 principal components used in PCA 
%	spike sorting. so:
%		Column 1: channel (?)
%		Column 2: unit #
%		Column 3: timestamp (in seconds)
%		Column 4: PCA1 weight
%		Column 5: PCA2 weight
%		Column 6: PCA3 weight
%		Column 7-end : waveform
%------------------------------------------------------------------------
%{
nexInfo:
	NexFileName: '1382_20191212_02_02_3200_test.nex'
	fData: [1�3 struct]
	nFiles: 3
	Fs: 4.8828e+04
	sweepStartBin: {[1�140 double]  [1�160 double]  [1�180 double]}
	sweepEndBin: {[1�140 double]  [1�160 double]  [1�180 double]}
	fileStartBin: [1 1708981 3662101]
	fileEndBin: [1708980 3662100 5859360]
	startTimes: [1�480 double]
	endTimes: [1�480 double]
	fileStartTime: [0 34.9999 74.9998]
	fileEndTime: [34.9999 74.9998 119.9997]
	Channels: [4 5 11 14]


nexInfo.fData
	DataPath: '~/Work/Data/TestData/MT'
	DataFile: '1382_20191212_02_02_3200_FREQ_TUNING.dat'
	startSweepBin: {# channels � 1 cell}
	endSweepBin: {# channels � 1 cell}
	sweepLen: {# channels � 1 cell}
	fileStartBin: 1
	fileEndBin: 1708980
	Dinf: [1�1 struct] =
	F: [1�1 struct]

nexInfo.fData.F
	base: '1382_20191212_02_02_3200_FREQ_TUNING'
	animal: '1382'
	datecode: '20191212'
	unit: '02'
	penetration: '02'
	depth: '3200'
	other: 'FREQ_TUNING'
	path: '~/Work/Data/TestData/MT'
	file: '1382_20191212_02_02_3200_FREQ_TUNING.dat'
	testfile: '1382_20191212_02_02_3200_FREQ_TUNING_testdata.mat'
%}

% Define some handy values for indexing into spikesAll and similar arrays
%		Column 1: channel (?)
%		Column 2: unit #
%		Column 3: timestamp (in seconds)
%		Column 4: PCA1 weight
%		Column 5: PCA2 weight
%		Column 6: PCA3 weight
%		Column 7-end : waveform
% 
[~, nc] = size(spikesAll);
CHAN_COL = 1;
UNIT_COL = 2;
TS_COL = 3;
PCA_COL = 4:6;
WAV_COL = 7:nc;

%% break up spiketimes by data file
% use file Start/End time to do this

spikesByFile = cell(nexInfo.nFiles, 1);

for f = 1:nexInfo.nFiles
	% could use between() function, but using something explicit here for
	% clarity
	valid_ts = (spikesAll(:, TS_COL) >= nexInfo.fileStartTime(f)) & ...
					(spikesAll(:, TS_COL) <= nexInfo.fileEndTime(f));
	spikesByFile{f} = spikesAll(valid_ts, :);
end

%{
for f = 1:nexInfo.nFiles

	% shift spike times to correspond to start of each data file
	% make a local copy...
	tmpS = spikesByFile{f}(:, TS_COL);
	% ...and subtract fileStartTime in seconds from all time stamps
	tmpS = tmpS - nexInfo.fileStartTime(f);
	% put back into spikes ByFile
	spikesByFile{f}(:, TS_COL) = tmpS;
end
%}

%% N next step: assign spike times to appropriate sweeps/stimuli
% store spikes for each file in spikes ByStim
spikesByStim = cell(nexInfo.nFiles, 1);
% loop through files
for f = 1:nexInfo.nFiles
	spikesByStim{f} = assign_spikes_to_sweeps(spikesByFile{f}, ...
												nexInfo.sweepStartBin{f}, ...
												nexInfo.sweepEndBin{f}, ...
												nexInfo.Fs, ...
												'sweep');
end
spikesOrig{end}(:, TS_COL), spikesByStim{end}{end}(:, TS_COL)

spikesByStim{end}{3}(:, TS_COL)
%% need to adjust spike times to start of each sweep
%{
for f = 1:nexInfo.nFiles
	% shift spike times to correspond to start of each data file
	
	% get offset to apply to timestamps (in seconds)
	ts_offset = nexInfo.fileStartTime(f);
	
	for 
	% make a local copy...
	tmpS = spikesByFile{f}(:, TS_COL);
	% ...and subtract fileStartTime in seconds from all time stamps
	tmpS = tmpS - nexInfo.fileStartTime(f);
	% put back into spikes ByFile
	spikesByFile{f}(:, TS_COL) = tmpS;
end
%}

%% convert spikes to table

Channel = spikesAll(:, CHAN_COL);
Unit = spikesAll(:, UNIT_COL);
TS = spikesAll(:, TS_COL);
PCAmat = spikesAll(:, PCA_COL);
WAVmat = spikesAll(:, WAV_COL);
PCA = num2cell(PCAmat, 2);
WAV = num2cell(WAVmat, 2);

T = table(Channel, Unit, TS, PCA, WAV)
