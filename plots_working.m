%------------------------------------------------------------------------
% plots_working.m
%------------------------------------------------------------------------
% TytoLogy:Experiments:optosort
%------------------------------------------------------------------------
% working script for importing sorted data from plexon
% amd plotting data
%------------------------------------------------------------------------
% See also: 
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%  Sharad J. Shanbhag
%	sshadnbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 18 June, 2020(SJS)
%	- split off from import_working.m
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

%------------------------------------------------------------------------
%------------------------------------------------------------------------
%------------------------------------------------------------------------
% sorted data locations
%------------------------------------------------------------------------
%------------------------------------------------------------------------
%------------------------------------------------------------------------

%{
%------------------------------------------------------------------------
% BBN data 
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
% merged data
%------------------------------------------------------------------------
%-------------------------------------------------
% PRE FRAInfo
%-------------------------------------------------
sortedPath = '/Users/sshanbhag/Work/Data/TestData/MT/1408/1408_20200319_02_01_950_MERGE-ch1,2,5,12,15-01-KmeansScan';
plxFile = '1408_20200319_02_01_950_MERGE.plx';
nexInfoFile = '1408_20200319_02_01_950_MERGE_nexinfo.mat';
%}


%-------------------------------------------------
% FRA wih FRAInfo
%-------------------------------------------------
sortedPath = '/Users/sshanbhag/Work/Data/TestData/working';
plxFile = '1407_20200309_03_01_1350_FRA.plx';
nexInfoFile = '1407_20200309_03_01_1350_FRA_nexinfo.mat';
SobjFile = fullfile(sortedPath, 'Sobj_1407_20200309_03_01_1350_FRA.mat');


%{
%-------------------------------------------------
% MERGED, FRA wih FRAInfo
%-------------------------------------------------
sortedPath = '/Users/sshanbhag/Work/Data/TestData/working';
plxFile = '1407_20200309_03_01_1350_MERGEALL.plx';
nexInfoFile = '1407_20200309_03_01_1350_MERGEALL.mat';
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

%{
%------------------------------------------------------------------------
% WAV data
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
%}

%------------------------------------------------------------------------
%------------------------------------------------------------------------
%------------------------------------------------------------------------

sendmsg(['plots_working running for ' sprintf('File: %s\n', plxFile)]);

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
if exist('SobjFile', 'var')
	if exist(SobjFile, 'file')
	sendmsg(['Loading SpikeData object from ' SobjFile]);
	load(SobjFile);
	else
		S = import_from_plexon(fullfile(sortedPath, plxFile), ...
									fullfile(sortedPath, nexInfoFile));
	end
else
	S = import_from_plexon(fullfile(sortedPath, plxFile), ...
								fullfile(sortedPath, nexInfoFile));
end
%{
%------------------------------------------------------------------------
%------------------------------------------------------------------------
% ALTERNATIVE
%------------------------------------------------------------------------
%% load S object from file
%------------------------------------------------------------------------
% use OptoFileName object to make this easier
% create object using nexfile
OFobj = OptoFileName(fullfile(nexPath, nexFile));
% create name of mat file containing SpikeData object
sfile = [OFobj.fileWithoutOther '_Sobj.mat'];
sendmsg(sprintf('loading Sobj from file\n\t%s', ...
						fullfile(nexPath, sfile)));
% load object from mat file
load(fullfile(nexPath, sfile));
%------------------------------------------------------------------------
%------------------------------------------------------------------------
%}
	
%------------------------------------------------------------------------
%% show file, channel, unit
%------------------------------------------------------------------------
% display file, channel, unit info
[fileList, channelList, unitList] = S.printInfo;

%------------------------------------------------------------------------
%% specify file, channel, unit
%------------------------------------------------------------------------

% What test data do you want to plot?
testToPlot = 'FRA';
% testToPlot = 'BBN';

% figure out file index for this test
% get list of test names
testNameList = S.listTestNames;
findx = find(strcmpi(testToPlot, testNameList));
if isempty(findx)
	error('%s: test %s not found in file %s', mfilename, testToPlot, ...
								plxFile);
end
% select Channel (technically , index into array of channel numbers)
channel = 4;
% select unit ID num
unit = 1;


%% get spikes times struct for 
fprintf('Getting data for file %d, channel, %d unit %d\n', ...
								findx, channel, unit);
st = S.getSpikesByStim(findx, channel, unit);
% make a local copy of Dinf for this file to make things a little simpler
Dinf = S.Info.FileInfo{findx}.Dinf;

%% plot data

if any(strcmpi(testToPlot, {'LEVEL', 'BBN'}))
	% set analysis window to [stimulus onset   stimulus offset]
	analysisWindow = [Dinf.audio.Delay ...
									(Dinf.audio.Delay + Dinf.audio.Duration)];
	% compute rate level function
	RLF = computeRLF(st.spiketimes, st.unique_stim, analysisWindow);
	% plot
	hRLF = plotCurveAndCI(RLF, 'mean');
	% create title string with 2 rows:
	%	filenametestToPlot = 'FRA';
	%	channel and unit
	tstr = {	st.fileName, ...
				sprintf('Channel %d Unit %d', channel, unit)};
	% add title to plot
	title(tstr, 'Interpreter', 'none');
end

if any(strcmpi(testToPlot, 'FREQ_TUNING'))
	% set analysis window to [stimulus onset   stimulus offset]
	analysisWindow = [Dinf.audio.Delay ...
									(Dinf.audio.Delay + Dinf.audio.Duration)];
	% compute rate level function
	FTC = computeFTC(st.spiketimes, st.unique_stim, analysisWindow);
	% plot
	hFTC = plotCurveAndCI(FTC, 'mean');
	% create title string with 2 rows:
	%	filename
	%	channel and unit
	tstr = {	st.fileName, ...
				sprintf('Channel %d Unit %d', channel, unit)};
	% add title to plot
	title(tstr, 'Interpreter', 'none');	
end

if any(strcmpi(testToPlot, 'FRA'))
	% window for spike count
	frawin = [Dinf.audio.Delay (Dinf.audio.Delay + Dinf.audio.Duration)];
	% calculate FRA stored in struct FRA
	FRA = computeFRA(st.spiketimes, st.unique_stim{1}, st.unique_stim{2}, frawin);
	% set fname to data file name
	FRA.fname = st.fileName;
	hFRA = plotFRA(FRA, 'dB');	
end


% plot waveforms for this channel and unit
S.plotUnitWaveforms(channel, unit);


