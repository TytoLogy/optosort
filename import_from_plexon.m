% function varargout = import_from_plexon(varargin)
%------------------------------------------------------------------------
% import_from_plexon.m
%------------------------------------------------------------------------
% TytoLogy:Experiments:optosort
%------------------------------------------------------------------------
% working script for importing sorted data from plexon
% uses objects
%------------------------------------------------------------------------
% See also: 
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%  Sharad J. Shanbhag
%	sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 12 February, 2020 (SJS)
%	- adapted from import_from_plexon_nonObj
% Revisions:
%------------------------------------------------------------------------
% TO DO:
%------------------------------------------------------------------------
% Initial things to define
%------------------------------------------------------------------------

%% path to readPXFileC

if ~exist('readPLXFileC')
	fprintf('adding readPLXFile to path\n');
	addpath('/Users/sshanbhag/Work/Code/Matlab/stable/Toolbox/Plexon/readPLXFileC');
end

%------------------------------------------------------------------------
% Definitions
%------------------------------------------------------------------------
% string used to separate text 
sepstr = '----------------------------------------------------';

% data locations
% sortedPath = '~/Work/Data/TestData/MT';
% rawPath = sortedPath;
% nexPath = sortedPath;
sortedPath = '/Users/sshanbhag/Work/Data/TestData/MT/1407';
rawPath = sortedPath;
nexPath = sortedPath;

% nexinfo file
nexInfoFile = '1407_20200309_03_01_1350_BBN_nexinfo.mat';

% nex file
% nexFile = '1382_20191212_02_02_3200_MERGE.nex';
nexFile = '1407_20200309_03_01_1350_BBN.nex';

% plx file
plxFile = '1407_20200309_03_01_1350_BBN-sorted.ch4,5,7,15.plx';
%------------------------------------------------------------------------
% Setup
%------------------------------------------------------------------------
fprintf('\n%s\n', sepstr);
fprintf('import_from_plexon running...\n');
fprintf('%s\n', sepstr);

%------------------------------------------------------------------------
%{
nexInfo:
	NexFileName: '1382_20191212_02_02_3200_test.nex'
	fData: [1×3 struct]
	nFiles: 3
	Fs: 4.8828e+04
	sweepStartBin: {[1×140 double]  [1×160 double]  [1×180 double]}
	sweepEndBin: {[1×140 double]  [1×160 double]  [1×180 double]}
	fileStartBin: [1 1708981 3662101]
	fileEndBin: [1708980 3662100 5859360]
	startTimes: [1×480 double]
	endTimes: [1×480 double]
	fileStartTime: [0 34.9999 74.9998]
	fileEndTime: [34.9999 74.9998 119.9997]
	Channels: [4 5 11 14]


nexInfo.fData
	DataPath: '~/Work/Data/TestData/MT'
	DataFile: '1382_20191212_02_02_3200_FREQ_TUNING.dat'
	startSweepBin: {# channels × 1 cell}
	endSweepBin: {# channels × 1 cell}
	sweepLen: {# channels × 1 cell}
	fileStartBin: 1
	fileEndBin: 1708980
	Dinf: [1×1 struct] =
	F: [1×1 struct]

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

%------------------------------------------------------------------------
%% load data
%------------------------------------------------------------------------

%% create SpikeData object
S = SpikeData();
fprintf('\n%s\n', sepstr);
fprintf('Loading nexInfo from file\n\t%s\n', fullfile(nexPath, nexInfoFile));
fprintf('%s\n', sepstr);
% read nexinfo from file
S.Info = SpikeInfo('file', fullfile(nexPath, nexInfoFile));

% get plx data
P = PLXData(fullfile(sortedPath, plxFile));

% add spikes
S = S.addPlexonSpikesFromPLXObj(P);

%% save Sobject in file
% get base from one of the file objects in S.Info
sfile = [S.Info.FileData(1).F.fileWithoutOther '_Sobj.mat'];
fprintf('\n%s\n', sepstr);
fprintf('writing Sobj to file\n\t%s\n', fullfile(nexPath, sfile));
fprintf('%s\n', sepstr);
save(fullfile(nexPath, sfile), '-MAT', 'S');


%% break up spiketimes by data file
% use file Start/End time to do this
[~, nc] = size(spikesAll);
CHAN_COL = 1;
UNIT_COL = 2;
TS_COL = 3;
PCA_COL = 4:6;
WAV_COL = 7:nc;

sbf = cell(S.Info.nFiles, 1);

for f = 1:S.Info.nFiles
	% could use between() function, but using something explicit here for
	% clarity
	valid_ts = (spikesAll(:, TS_COL) >= S.Info.fileStartTime(f)) & ...
					(spikesAll(:, TS_COL) <= S.Info.fileEndTime(f));
	sbf{f} = spikesAll(valid_ts, :);
end

%% try with object
spikesByFile = cell(S.Info.nFiles, 1);
for f = 1:S.Info.nFiles
	spikesByFile{f} = S.spikesForFile(f);
end

%% N next step: assign spike times to appropriate sweeps/stimuli
% OBJ - by unit
unitID = S.listUnits;
nunits = S.nUnits;
% store spikes for each file in spikes ByStim, all units
spikesBySweepAndUnit = cell(S.Info.nFiles, nunits);
% loop through files
for f = 1:S.Info.nFiles
	% loop through units
	for u = 1:nunits
		% get spikes (as a table) for each file (rows) and unit (columns)
		spikesBySweepAndUnit{f, u} = S.spikesForAnalysisByUnit(f, unitID(u), 'sweep');
	end
end

%% OBJ - don't separate by unit - can do posthoc
% store spikes for each file in spikes ByStim, all units
spikesBySweep = cell(S.Info.nFiles, 1);
% loop through files
for f = 1:S.Info.nFiles
		spikesBySweep{f} = S.spikesForAnalysis(f, 'sweep');
end

%% extract timestamps for use in analysis, raster plots, etc., separated by unit
% convert to timestamps
spikesForAnalysis = cell(S.Info.nFiles, nunits);

% loop through files
for f = 1:S.Info.nFiles
	% loop through units
	for u = 1:nunits
		% extract timestamp data from each table, store in 
		spikesForAnalysis{f, u} = S.spikesForAnalysisByUnit(f, unitID(u), 'sweep');
	end
end

extractTimeStamps



%% plot waveforms

S.plotUnitWaveforms(S.listUnits);


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
% 
% Channel = spikesAll(:, CHAN_COL);
% Unit = spikesAll(:, UNIT_COL);
% TS = spikesAll(:, TS_COL);
% PCAmat = spikesAll(:, PCA_COL);
% WAVmat = spikesAll(:, WAV_COL);
% PCA = num2cell(PCAmat, 2);
% WAV = num2cell(WAVmat, 2);
% 
% T = table(Channel, Unit, TS, PCA, WAV)
