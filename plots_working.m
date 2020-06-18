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
S = import_from_plexon(fullfile(sortedPath, plxFile), ...
									fullfile(nexPath, nexInfoFile));

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
%% specify file, channel, unit
%------------------------------------------------------------------------
sendmsg(sprintf('Data in file %s', sfile));
% get and display list of files
fList = S.listFiles;
fprintf('Files:\n');
fprintf('\tIndex\t\tFilename\n');
for f = 1:S.Info.nFiles
	fprintf('\t%d:\t\t%s\n', f, fList{f});
end
fprintf('\n');
% get and display list of channels...
cList = S.listChannels;
% ...and units
% n.b.: could also get both using listUnits: [uList, cList] = S.listUnits
uList = S.listUnits;
fprintf('Channels and Unit ID #s:\n');
fprintf('\tIndex\tChannel\tUnits\n');
for c = 1:length(cList)
	fprintf('\t%d:\t%d\t', c, cList(c));
	fprintf('%d ', uList{c});
	fprintf('\n');
end
fprintf('\n');

% What file, channel unit to plot?
% file index
findx = 1;
% select Channel (technically , index into array of channel numbers)
cindx = 1;
% select unit ID num
uindx = 1;







