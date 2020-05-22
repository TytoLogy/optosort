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

if ~exist('readPLXFileC', 'file')
	fprintf('import_from_plexon: adding readPLXFile to path\n');
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
% nexFile = '1382_20191212_02_02_3200_MERGE.nex';

%------------------------------------------------------------------------
%{
sortedPath = '/Users/sshanbhag/Work/Data/TestData/MT/1407';
rawPath = sortedPath;
nexPath = sortedPath;
nexInfoFile = '1407_20200309_03_01_1350_BBN_nexinfo.mat';
nexFile = '1407_20200309_03_01_1350_BBN.nex';
plxFile = '1407_20200309_03_01_1350_BBN-sorted.ch4,5,7,15.plx';
%}
%------------------------------------------------------------------------

%------------------------------------------------------------------------
sortedPath = '/Users/sshanbhag/Work/Data/TestData/MT/1408';
rawPath = sortedPath;
nexPath = sortedPath;
nexInfoFile = '1408_20200319_02_01_950_WAV_nexinfo-ch1,2,5,12,15.mat';
nexFile = '1408_20200319_02_01_950_WAV-ch1,2,5,12,15.nex';
% file without continuous data
% plxFile = '1408_20200319_02_01_950_WAV-ch1,2,5,12,15-01-KmeansScan.plx';
% file with continuous data
plxFile = '1408_20200319_950_WAV_ch1,2,5,12,15-01-SORTEDContKmeansScan.plx';
%------------------------------------------------------------------------

fprintf('\n%s\n', sepstr);
fprintf('import_from_plexon running...\n');
fprintf('File: %s\n', plxFile);
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
sendmsg(sprintf('Loading nexInfo from file\n\t%s\n', fullfile(nexPath, nexInfoFile)));
% read nexinfo from file
S.Info = SpikeInfo('file', fullfile(nexPath, nexInfoFile));
% get plx data
Plx = PLXData(fullfile(sortedPath, plxFile), 'all', 'continuous');


%% add spikes
S = S.addPlexonSpikesFromPLXObj(Plx);
% add Continuous Data
if Plx.hasContinuousData
	sendmsg('Adding ContinuousChannels data to SpikeData');
	S = S.addContinuousDataFromPLXObj(Plx);
end

%% save Sobject in file
% get base from one of the file objects in S.Info
sfile = [S.Info.FileInfo{1}.F.fileWithoutOther '_Sobj.mat'];
fprintf('\n%s\n', sepstr);
fprintf('writing Sobj to file\n\t%s\n', fullfile(nexPath, sfile));
fprintf('%s\n', sepstr);
save(fullfile(nexPath, sfile), '-MAT', 'S');


%% break up spiketimes by data file
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
spikesByFile = cell(S.Info.nFiles, 1);
% loop through files
for f = 1:S.Info.nFiles
		spikesByFile{f} = S.spikesForAnalysis(f, 'align', 'sweep');
end

%% extract timestamps for use in analysis, raster plots, etc., separated by channel
% convert to timestamps
spikesByFileChannelSweep = cell(S.Info.nFiles, S.Info.nChannels);

% loop through files
for f = 1:S.Info.nFiles
	% loop through channels
	for c = 1:S.Info.nChannels
		fprintf('Getting spikes for file %d, channel %d\n', f, c);
	
		% extract timestamp data from each table, store in cell matrix
		spikesForAnalysis{f, c} = S.spikesForAnalysis(f, 'Channel', c, 'align', 'sweep');
	end
end

%%
spikesBySweep = S.spikesForAnalysis(f, 'Channel', S.Info.ADchannel, 'align', 'sweep');


%% plot waveforms

S.plotUnitWaveforms(S.listUnits);


%%

tmpS = S.spikesForFile(1);
channel_rows = false(size(tmpS.Channel));
for c = 3:5
	channel_rows = channel_rows | (tmpS.Channel == S.Info.ADchannel(c));
	cm(:, c) = channel_rows;
end
unique(tmpS.Channel(channel_rows, :))



