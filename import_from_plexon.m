function varargout = import_from_plexon(varargin)
%------------------------------------------------------------------------
% [SpikeData] = import_from_plexon(plxfile, nexinfofile, <options>)
%------------------------------------------------------------------------
% TytoLogy:Experiments:optosort
%------------------------------------------------------------------------
% function for importing sorted data from Plexon OfflineSorter .plx files
% 
% Input Args:
%	<no args>		Will open panel to ask user to select .plx file
%	plxfile			.plx output file from Plexon OfflineSorter
%	nexinfofile		_nexInfo.mat file from export_for_plexon() function
%
%	Options:
%		'Continuous'	load continuous data from .plx file (default)
% 								can be slow/use lots of RAM
% 		'NoContinuous'	do not load continuous data from .plx file
% 
% Output Args:
%	SpikeData	SpikeData object
%	
%------------------------------------------------------------------------
% See also: SpikeData (class), export_for_plexon, OptoObjects
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%  Sharad J. Shanbhag
%	sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 12 February, 2020 (SJS)
%	- adapted from import_from_plexon_nonObj
%
% Revisions:
%	25 Jun 2020 (SJS):
%		- added comments
%		- added option to load/not load continuous data
%------------------------------------------------------------------------
% TO DO:
% - maybe have way to specify what to load...
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% Definitions
%------------------------------------------------------------------------
% string used to separate text 
sepstr = '----------------------------------------------------';

% load continuous data? default is true
load_continuous = true;

%------------------------------------------------------------------------
% path to readPXFileC
%------------------------------------------------------------------------
if ~exist('readPLXFileC', 'file')
	try
		fprintf('import_from_plexon: adding readPLXFile to path\n');	
		addpath('/Users/sshanbhag/Work/Code/Matlab/stable/Toolbox/Plexon/readPLXFileC');
	catch
		error('import_from_plexon: function readPLXFileC not in path');		
	end
end

%------------------------------------------------------------------------
% process args
%------------------------------------------------------------------------
if isempty(varargin)
	[plxFile, plxPath] = uigetfile({	'*.plx', 'Plexon file (*.plx)'; ...
												'*.*', 'All Files (*.*)' }, ...
												'Select .plx file');
	% return if cancelled
	if plxFile == 0
		fprintf('cancelled\n');
		varargout{nargout} = [];
		return
	end
	% assign name to OptoFileName object to assist name generation
	tmpf = OptoFileName(fullfile(plxPath, plxFile));
	% build nexinfo name skeleton to search for
	nexinfoFile = [tmpf.fileWithoutOther '*_nexinfo.mat'];
	% get neinfofile and path
	[nexinfoFile, nexinfoPath] = uigetfile(fullfile(plxPath, nexinfoFile), ...
														'Select nexinfo.mat file');
	if nexinfoFile == 0
		fprintf('cancelled\n');
		varargout{1} = [];
		return
	end
	
elseif nargin < 2
		error('%s: input arg error, need .plx and nexinfo.mat file', ...
						mfilename);

else
	% user provided plx file and nexinfo file
	[plxPath, plxFile, ext] = fileparts(varargin{1});
	plxFile = [plxFile ext];
	[nexinfoPath, nexinfoFile, ext] = fileparts(varargin{2});
	nexinfoFile = [nexinfoFile ext];
	
	% check options
	if nargin > 2
		for n = 3:nargin
			argV = upper(varargin{n});
			switch(argV)
				case {'CONTINUOUS', 'CONT'}
					load_continuous = true;
				case {'NOCONTINUOUS', 'NOCONT', 'NO_CONT'}
					load_continuous = false;
				otherwise
					error('%s: invalid option %s', mfilename, argV);
			end
		end
	end
end


%------------------------------------------------------------------------
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%------------------------------------------------------------------------
%------------------------------------------------------------------------

fprintf('\n%s\n', sepstr);
fprintf('import_from_plexon running...\n');
fprintf('PLX File: %s\n', plxFile);
fprintf('nexinfo File: %s\n', nexinfoFile);
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
%------------------------------------------------------------------------
% load data
%------------------------------------------------------------------------
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% create SpikeData object that will be used to hold and access data
%------------------------------------------------------------------------
S = SpikeData();
sendmsg(sprintf('Loading nexInfo from file\n\t%s\n', fullfile(nexinfoPath, nexinfoFile)));
% read nexinfo from file
S.Info = SpikeInfo('file', fullfile(nexinfoPath, nexinfoFile));
% get plx data with or without continuous data
if load_continuous
	Plx = PLXData(fullfile(plxPath, plxFile), 'all', 'continuous');
else
	Plx = PLXData(fullfile(plxPath, plxFile), 'all', 'nocontinuous');
end

%------------------------------------------------------------------------
% add spikes from PLX data to SpikeData object (S)
%------------------------------------------------------------------------
S = S.addPlexonSpikesFromPLXObj(Plx);
% add Continuous Data (if present)
if Plx.hasContinuousData
	sendmsg('Adding ContinuousChannels data to SpikeData');
	S = S.addContinuousDataFromPLXObj(Plx);
end

%------------------------------------------------------------------------
% assign outputs
%------------------------------------------------------------------------
varargout{1} = S;




