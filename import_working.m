%------------------------------------------------------------------------
% import_working.m
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
%	sshadnbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 16 June, 2020(SJS)
%	- split off from importTest.m
% Revisions:
%------------------------------------------------------------------------
% TO DO:
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%% add path to readPLXFileC if needed
%------------------------------------------------------------------------
% readPLXFileC is a function downloaded from MATLAB Central that allows
% direct reading of PLX file data in Matlab
if ~exist('readPLXFileC', 'file')
	fprintf('import_from_plexon: adding readPLXFile to path\n');
	addpath('readPLXFileC');
end

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% data locations
%------------------------------------------------------------------------
%------------------------------------------------------------------------
% sortedPath		path to sorted data
% rawPath			path to raw data (.dat), not necessary?
% nexPath			path to nex file (.nex), not necessary?
% nexFile			name of nex file used for sorting (not necessary?)
% plxFile			name of sorted .plx file from Plexon OfflineSorter
%------------------------------------------------------------------------

%{
sortedPath = '~/Work/Data/TestData/MT';
rawPath = sortedPath;
nexPath = sortedPath;
nexFile = '1382_20191212_02_02_3200_MERGE.nex';
%}

%{
%------------------------------------------------------------------------
% sorted data locations
%------------------------------------------------------------------------
sortedPath = '/Users/sshanbhag/Work/Data/TestData/MT/1407';
rawPath = sortedPath;
nexPath = sortedPath;
nexInfoFile = '1407_20200309_03_01_1350_BBN_nexinfo.mat';
nexFile = '1407_20200309_03_01_1350_BBN.nex';
plxFile = '1407_20200309_03_01_1350_BBN-sorted.ch4,5,7,15.plx';
%}

%{
%------------------------------------------------------------------------
% raw data locations
%------------------------------------------------------------------------
exportOpts.DataPath = {'/Users/sshanbhag/Work/Data/TestData/MT/1407'};
% write to D drive on PETROL for testing with Plexon OFS
% exportOpts.OutputPath = '/Volumes/D/1407';
% local working dir
exportOpts.OutputPath = '/Users/sshanbhag/Work/Data/TestData/working';
exportOpts.DataFile = {	'1407_20200309_03_01_1350_BBN.dat', ...
								'1407_20200309_03_01_1350_FREQ_TUNING.dat'};
exportOpts.TestFile = { '1407_20200309_03_01_1350_BBN_testdata.mat', ...
								'1407_20200309_03_01_1350_FREQ_TUNING_testdata.mat'};
exportOpts.OutputFile = '1407_20200309_03_01_1350_MERGE.nex';
exportOpts.Channels = [4, 5, 7, 15];
% filter parameters for raw neural data
exportOpts.BPfilt.Fc = [250 4000];
exportOpts.BPfilt.forder = 5;
exportOpts.BPfilt.ramp = 1;
exportOpts.BPfilt.type = 'butter';
% exportOpts.resampleData = [];
%------------------------------------------------------------------------
%}

%{
%------------------------------------------------------------------------
% "fake" data for testing sampling rate/timing in OFS .plx file
%------------------------------------------------------------------------
sortedPath = '/Users/sshanbhag/Work/Data/TestData/working/FakeData';
rawPath = sortedPath;
nexPath = sortedPath;
nexInfoFile = '1407_20200309_03_01_1350_TIMETESTDATA_nexinfo.mat';
nexFile = '1407_20200309_03_01_1350_TIMETESTDATA.nex';
plxFile = '1407_20200309_03_01_1350_TIMETESTDATA-02.plx';
%------------------------------------------------------------------------
%}

%{
%------------------------------------------------------------------------
% "fake" data for testing sorting OFS .plx file
%------------------------------------------------------------------------
sortedPath = '/Users/sshanbhag/Work/Data/TestData/working/FakeData/TestData';
rawPath = sortedPath;
nexPath = sortedPath;
nexInfoFile = '1407_20200309_03_01_1350_TESTDATA_nexinfo.mat';
nexFile = '1407_20200309_03_01_1350_TESTDATA.nex';
plxFile = '1407_20200309_03_01_1350_TESTDATA-Sort.plx';
%------------------------------------------------------------------------
%}

%------------------------------------------------------------------------
%{
sortedPath = '/Users/sshanbhag/Work/Data/TestData/MT/1408';
rawPath = sortedPath;
nexPath = sortedPath;
nexInfoFile = '1408_20200319_02_01_950_WAV_nexinfo-ch1,2,5,12,15.mat';
nexFile = '1408_20200319_02_01_950_WAV-ch1,2,5,12,15.nex';
% file without continuous data
% plxFile = '1408_20200319_02_01_950_WAV-ch1,2,5,12,15-01-KmeansScan.plx';
% file with continuous data
plxFile = '1408_20200319_950_WAV_ch1,2,5,12,15-01-SORTEDContKmeansScan.plx';
%}
%------------------------------------------------------------------------

sortedPath = '/Users/sshanbhag/Work/Data/TestData/working';
rawPath = sortedPath;
nexPath = sortedPath;
nexInfoFile = '1407_20200309_03_01_1350_MERGEEVENTS2_nexinfo.mat';
nexFile = '1407_20200309_03_01_1350_MERGEEVENTS2.nex';
plxFile = '1407_20200309_03_01_1350_MERGEEVENTS2.plx';

%------------------------------------------------------------------------
%------------------------------------------------------------------------

sendmsg(['import_from_plexon running for ' sprintf('File: %s\n', plxFile)]);

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
%------------------------------------------------------------------------
%% load sorted data
%------------------------------------------------------------------------
% How to use:
% import_from_plexon(<plx file name>, <nexinfo file name>)
%
% will return a SpikeData object containing data from plx file:
%	spike times/unit information if sorted
%	continuous data, if saved in plx file
%	stimulus info
%------------------------------------------------------------------------
S = import_from_plexon(fullfile(sortedPath, plxFile), ...
									fullfile(nexPath, nexInfoFile));
%{
%------------------------------------------------------------------------
%------------------------------------------------------------------------
%% load raw data in order to compare with plx data
%------------------------------------------------------------------------
%------------------------------------------------------------------------
F  = defineSampleData(exportOpts.DataPath, ...
											exportOpts.DataFile, exportOpts.TestFile);
% storage for raw data
D = cell(length(exportOpts.DataFile), 1);
N = cell(length(exportOpts.DataFile), 1);

for f = 1:length(F)
	[tmp, N{f}] = read_data_for_export(F(f), exportOpts.Channels, exportOpts.BPfilt, []);
	D{f} = tmp{1};
end

%------------------------------------------------------------------------
%------------------------------------------------------------------------
%% compare data from raw and sorted continuous data
%------------------------------------------------------------------------
%------------------------------------------------------------------------

% file (index into cell) to plot
fnum = 1;
% channel (index into channel list) to plot
chan = 1;
% sweep to plot
sweep = 3;

figure(1)
subplot(211)
rawD = D{fnum}{chan, sweep};
nD = length(rawD);
plxStart = (1 + (sweep-1)*nD)
plxD = S.Continuous(chan).Values(plxStart:(plxStart+nD));
plot(rawD)
title({F(1).file; sprintf('Channel %d Sweep %d', chan, sweep)}, 'Interpreter', 'none');
subplot(212)
plot(plxD)

figure(2)
plot(xcorr(rawD, plxD))
%}
									
%------------------------------------------------------------------------
%% save Sobject in file
%------------------------------------------------------------------------
% get base from one of the file objects in S.Info
sfile = [S.Info.FileInfo{1}.F.fileWithoutOther '_Sobj.mat'];
fprintf('\n%s\n', sepstr);
fprintf('writing Sobj to file\n\t%s\n', fullfile(nexPath, sfile));
fprintf('%s\n', sepstr);
save(fullfile(nexPath, sfile), '-MAT', 'S');

%------------------------------------------------------------------------
%% break up spiketimes by data file
%------------------------------------------------------------------------
spikesByFile = cell(S.Info.nFiles, 1);
for f = 1:S.Info.nFiles
	spikesByFile{f} = S.spikesForFile(f);
end

%------------------------------------------------------------------------
%% assign spike times to appropriate sweeps/stimuli
%------------------------------------------------------------------------
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

%------------------------------------------------------------------------
%% OBJ - don't separate by unit - can do posthoc
%------------------------------------------------------------------------
% store spikes for each file in spikes ByStim, all units
spikesByFile = cell(S.Info.nFiles, 1);
% loop through files
for f = 1:S.Info.nFiles
		spikesByFile{f} = S.spikesForAnalysis(f, 'align', 'sweep');
end

%------------------------------------------------------------------------
%% extract timestamps for use in analysis, raster plots, etc.
% in this case, organize data by file and channel in the spikesFC array
% data will be returned in a MATLAB table data structure.
% for more about accessing data in a table, please see matlab docs
%------------------------------------------------------------------------
% convert to timestamps
spikesFC = cell(S.Info.nFiles, S.Info.nChannels);

% loop through files
for f = 1:S.Info.nFiles
	% loop through channels
	for c = 1:S.Info.nChannels
		fprintf('Getting spikes for file %d, channel %d\n', f, S.Info.ADchannel(c));
		% extract timestamp data from each table, store in cell matrix
		spikesFC{f, c} = S.spikesForAnalysis(f, ...
									'Channel', S.Info.ADchannel(c), 'align', 'sweep');
	end
end


%% Get spike data for a specific file aligned by sweep for given channels
% get data for first file
fNum = 1;
spikesBySweep = S.spikesForAnalysis(fNum, 'Channel', S.Info.ADchannel, 'align', 'sweep');


%% plot waveforms

S.plotUnitWaveforms(4, S.listUnits(4));
 


%% check two data files

tmpS = S.spikesForFile(1);
channel_rows = false(size(tmpS.Channel));
for c = 3:5
	channel_rows = channel_rows | (tmpS.Channel == S.Info.ADchannel(c));
	cm(:, c) = channel_rows;
end
unique(tmpS.Channel(channel_rows, :))





