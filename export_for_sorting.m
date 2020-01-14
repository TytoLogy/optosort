function varargout = export_for_sorting(varargin)
%------------------------------------------------------------------------
% [] = export_for_sorting('DataList', <data file list>)
%------------------------------------------------------------------------
% TytoLogy:Experiments:optosort
%------------------------------------------------------------------------
% 
% exports raw sweep data for analysis using Plexon. Output is in .nex
% (Neural Explorere) format
% 
%------------------------------------------------------------------------
% Input Arguments:
%	With no inputs provided, a dialog will open to specify a file with a 
%	list of .dat files to process.
%
% 	'DataList'	<file with list of data files.m>
% 
% Output Arguments:
% 	Output	output info
%
%------------------------------------------------------------------------
% See also: 
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%  Sharad J. Shanbhag
%	sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 14 January, 2020 (SJS)
%
% Revisions:
%------------------------------------------------------------------------
% TO DO:
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% Initial things to define
%------------------------------------------------------------------------

sepstr = '----------------------------------------------------';
NEX_UTIL_PATH = ['~/Work/Code/Matlab/stable/Toolbox/NeuroExplorer' ...
					'/HowToReadAndWriteNexAndNex5FilesInMatlab'];
%------------------------------------------------------------------------
% Setup
%------------------------------------------------------------------------

fprintf('\n%s\n', sepstr);
fprintf('%s running...\n', mfilename);
fprintf('%s\n', sepstr);

% add path to .nex file utils
if ~exist('nexCreateFileData', 'file')
	if exist(NEX_UTIL_PATH, 'dir')
		fprintf('%s: Adding NEX file utility path\n', mfilename)
		addpath(NEX_UTIL_PATH);
	else
		fprintf('¡¡¡¡¡NEX utilities not found!!!!!!!\n');
		fprintf('Please add to path');
		error('%s: NEX utils not found', mfilename);
	end
end

% check for path to readOptoData
if ~exist('readOptoData', 'file')
	fprintf('¡¡¡¡¡readOptoData function not found!!!!!!!\n');
	fprintf('Please add to path');
	error('%s: readOptoData (in Opto project folder) not found', mfilename);
end

%------------------------------------------------------------------------
% Get Data File(s) information
%------------------------------------------------------------------------

% define path to data file and data file for testing
F = defineSampleData('1372_20191126_93_01_1500_exportlist.m');
fprintf('DataPath = %s\n', DataPath);
for f = 1:nFiles
	fprintf('DataFile{%d} = %s\n', f, DataFile{f});
end
fprintf('Animal: %s\n', F(1).animal);

% !!!!!!!
% in future, will need to have user specify all the files for this
% recording session and recording location!!!
% ¡¡¡¡¡¡¡

% channel(s) of data to obtain
% channel = 8;
Channels = [11 9 14];
nChannels = length(Channels);

% filter info
BPfilt.Fc = [300 4000];
BPfilt.forder = 5;
BPfilt.ramp = 1;

%% loop through files

% allocate some things
% bins for start and end of each file's data
fileStartBin = zeros(1, nFiles);
fileEndBin = zeros(1, nFiles);
% bins for all sweep starts and ends
sweepStartBin = cell(1, nFiles);
sweepEndBin = cell(1, nFiles);
% each file's sampling rate for neural data
tmpFs = zeros(nFiles, 1);
% struct to hold everything for each file
fData = repmat(	struct(		'DataPath', '', ...
										'DataFile', '', ...
										'cSweeps', {}, ...
										'startSweepBin', {}, ...
										'endSweepBin', {}, ...
										'fileStartBin', [], ...
										'fileEndBin', [], ...
										'Dinf', [] ...
								), ...
						1, nFiles);
for f = 1:nFiles
	% save parse file info
	fData(f).Finfo = F(f);
	% get data for each file and channel and convert to row vector format
	% algorithm:
	%		(1) put each sweep for this channel in a {1, # sweeps} cell array
	%				cSweeps
	%		(2) make a note of the length of each sweep to use for
	%				markers/timestamps
	%		(3) after cSweeps is built, convert to a row vector using cell2mat
	% 
	% read in data
	% use readOptoData to read in raw data. 
	[D, Dinf] = readOptoData(fullfile(DataPath, DataFile{f}));
	% Fix test info
	Dinf = correctTestType(Dinf);
	% store file info
	fData(f).DataPath = DataPath;
	fData(f).DataFile = DataFile{f};
	% build filter for neural data
	BPfilt.Fs = Dinf.indev.Fs;
	BPfilt.Fnyq = Dinf.indev.Fs / 2;
	BPfilt.cutoff = BPfilt.Fc / BPfilt.Fnyq;
	[BPfilt.b, BPfilt.a] = butter(BPfilt.forder, BPfilt.cutoff, 'bandpass');

	% check to make sure consistent # of sweeps (aka trials)
	if Dinf.nread ~= length(D)
		error('%s: mismatch in Dinf.nread (%d) and length(D) (%d)', ...
					mfilename, Dinf.nread, length(D));
	end
	
	% build into sweeps by channel format
	fprintf('Test type: %s\n', Dinf.test.Type);
	[fData(f).cSweeps, fData(f).startSweepBin, fData(f).endSweepBin] = ...
					buildChannelData(Channels, BPfilt, D, Dinf);
	% check the start and end sweep bin data for consistency
	if check_sweeps(fData(f).startSweepBin)
		warning(['File %s: Inconsistent startSweepBin' ...
							'values across channels!!!!'], fData(f).DataFile);
	end
	if check_sweeps(fData(f).endSweepBin)
		warning('Inconsistent endSweepBin values across channels!!!!');
	end
	% store sample for start of this file (should be 1); use channel 1 value
	fData(f).fileStartBin = fData(f).startSweepBin{1}(1);
	% store sample for end of this file
	fData(f).fileEndBin = fData(f).endSweepBin{1}(end);
	% store file information struct
	fData(f).Dinf = Dinf;
	% to avoid any issues, should make sure sample rates are consistent
	% need in implement a check somehow.  this code doesn't work:
	% if any(Dinf.indev.Fs ~= [fData(:).Dinf.indev.Fs])
	%	error('Sample Rate Mismatch found!');
	% end
	tmpFs(f) = fData(f).Dinf.indev.Fs;
	
	% calculate overall start and end bins for each file's data
	if f == 1
		fileStartBin(f) = fData(f).fileStartBin;
		fileEndBin(f) = fData(f).fileEndBin;
	else
		% add 1 to prior end bin for start
		fileStartBin(f) = fileEndBin(f-1) + 1;
		fileEndBin(f) = fileStartBin(f) + fData(f).fileEndBin - 1;
	end
end

% assign values for bins
for f = 1:nFiles
	% calculate start and end sweep bins for each file's data
	if f == 1
		sweepStartBin{f} = fData(f).startSweepBin{1};
		sweepEndBin{f} = fData(f).endSweepBin{1};
	else
		% add previous file's final endSweepBin value as offset
		sweepStartBin{f} = fData(f).startSweepBin{1} + sweepEndBin{f-1}(end);
		sweepEndBin{f} = fData(f).endSweepBin{1} + sweepEndBin{f-1}(end);
	end
end

%% Create file and sweep start/end bin and timestamp vectors

%check sampling rates
if ~all(tmpFs(1) == tmpFs)
	error('Sample Rate Mismatch!!!');
else
	% store overall sample rate
	Fs = tmpFs(1);
end

% convert file start/end bins to times
fileStartTime = (fileStartBin - 1) ./ Fs;
fileEndTime = (fileEndBin - 1) ./ Fs;

% convert sweep bin cells to vectors...
startBins = [sweepStartBin{:}];
endBins = [sweepEndBin{:}];
% ... and then to times
startTimes = (startBins - 1) ./ Fs;
endTimes = (endBins - 1) ./ Fs;


%% convert (concatenate) cSweeps to vector for each channel, add to nex struct
% to save on memory requirements, clear channel data after adding to nex
% data struct.

% create output .nex file name 
%	assume data from first file is consistent with others!!!!!!!!!
NexFileName = [	fData(1).Finfo.animal '_' ...
						fData(1).Finfo.datecode '_' ...
						fData(1).Finfo.unit '_' ...
						fData(1).Finfo.penetration '_' ...
						fData(1).Finfo.depth ...
						'.nex'];
					
% start new nex file data struct
nD = nexCreateFileData(Fs);

% loop through channels
for c = 1:nChannels
	% this will be a matrix of format
	% 	[# channels, (# sweeps) * (# samples per sweep)
	cVector = cell(1, nFiles);
	for f = 1:nFiles
		cVector{1, f} = fData(f).cSweeps(c, :);
	end
	% concatenate cell array, convert to vector, add to nex struct
	% steps:
	%	concatenate: tmp = [cVector{:}];
	%	tmpVector = cell2mat(tmp);
	%  add to nex struct:
	%	[nexFile] = nexAddContinuous( nexFile, startTime, adFreq, values, name)
	nD = nexAddContinuous(nD, fileStartTime(1), Fs, ...
									cell2mat([cVector{:}]), ...
									sprintf('spikechan_%d', Channels(c)));
	% clear cVector to save memory
	clear cVector
end

% add start sweep time stamps as event - assume consistent across channels!
%  [nexFile] = nexAddEvent( nexFile, timestamps, name )
% events must be in column format...?
nD = nexAddEvent(nD, force_col(startTimes), 'startsweep');
% add end sweep time stamps as event - assume consistent across channels!
nD = nexAddEvent(nD, force_col(endTimes), 'endsweep');
% add file times
nD = nexAddEvent(nD, force_col(fileStartTime), 'filestart');
nD = nexAddEvent(nD, force_col(fileEndTime), 'fileend');
% write to nexfile
writeNexFile(nD, NexFileName);

end

function varargout = defineSampleData(varargin)
%------------------------------------------------------------------------
%------------------------------------------------------------------------
% defineSampleData.m
%------------------------------------------------------------------------
% defines path to sample data for testing
%------------------------------------------------------------------------
%------------------------------------------------------------------------
	DataPath = {};
	DataFile = {};
	TestFile = {};
	if nargin
		if exist(varargin{1}, 'file')
			eval(varargin{1})
		end
	end
	if any([isempty(DataPath) isempty(DataFile)])
		% do something to load data file names
	end

	for f = 1:length(DataFile)
		F(f) = parse_opto_filename(DataFile{f}); %#ok<AGROW>
		F(f).path = DataPath; %#ok<AGROW>
		F(f).file = DataFile{f}; %#ok<AGROW>
		F(f).testfile = TestFile{f}; %#ok<AGROW>
	end
	varargout{1} = F;
end
