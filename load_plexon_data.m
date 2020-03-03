function varargout = load_plexon_data(varargin)
%------------------------------------------------------------------------
% load_plexon_data.m
%------------------------------------------------------------------------
% TytoLogy:Experiments:optosort
%------------------------------------------------------------------------
% loads sorted data, returns SpikeData objects
%------------------------------------------------------------------------
% See also: 
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%  Sharad J. Shanbhag
%	sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 3 <arch, 2020 (SJS)
%	- adapted from import_from_plexon
% Revisions:
%------------------------------------------------------------------------
% TO DO:
%------------------------------------------------------------------------
% Initial things to define
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% definitions
%------------------------------------------------------------------------

% string used to separate text 
sepstr = '----------------------------------------------------';


%------------------------------------------------------------------------
% process inputs
%------------------------------------------------------------------------

% if no input file provided, get file from user
if nargin == 0
	[sortedPath, sortedFile] = uigetfile('*.mat', ...
														'Select Plexon output .mat file');
	% if user cancelled, exit gracefully
	if isempty(sortedPath) || isempty(sortedFile)
		fprintf('%s: cancelled\n', mfilename);
		varargout{1} = [];
		return
	end
elseif nargin == 1
	[sortedPath, tmpFile, tmpExt] = fileparts(varargin{1});
	sortedFile = [tmpFile, tmpExt];
else
	error('%s: huh????', mfilename);
end

% build nexinfo file name from sortedFile
nexInfoFile = buildNexInfoFromPlexonMat(sortedFile);
	
	



%------------------------------------------------------------------------
% Setup
%------------------------------------------------------------------------
fprintf('\n%s\n', sepstr);
fprintf('%s running...\n', mfilename);
fprintf('%s\n', sepstr);

%------------------------------------------------------------------------
%------------------------------------------------------------------------
%% Plexon-exported sorted data:
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

% Define some handy values for indexing into spikesAll and similar arrays
%		Column 1: channel (?)
%		Column 2: unit #
%		Column 3: timestamp (in seconds)
%		Column 4: PCA1 weight
%		Column 5: PCA2 weight
%		Column 6: PCA3 weight
%		Column 7-end : waveform
% 

%------------------------------------------------------------------------
%% load data
%------------------------------------------------------------------------

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

% create SpikeData object
S = SpikeData();
fprintf('\n%s\n', sepstr);
fprintf('Loading nexInfo from file\n\t%s\n', fullfile(nexPath, nexInfoFile));
fprintf('%s\n', sepstr);
% nexinfo
S.Info = SpikeInfo('file', fullfile(nexPath, nexInfoFile));

% add spikes
% for now, just load the first channel data - need to figure out a way to
% deal with multiple channels (array of SpikeData objects?)
tmp = load(fullfile(sortedPath, sortedFile), '-MAT', plxvars{1});
spikesAll = tmp.(plxvars{1});
S = S.addPlexonSpikes(tmp.(plxvars{1}), plxvars{1});
clear tmp;

varargout{1} = S;

end

function nexfilename = buildNexInfoFromPlexonMat(plexfilename)
% from plexfilename, build nexfilename

	% first, init opto file obj
	F = OptoFileName(plexfilename);
	

end