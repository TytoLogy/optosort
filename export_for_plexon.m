function varargout = export_for_plexon(varargin)
%------------------------------------------------------------------------
% [nD, nexInfo] = export_for_plexon(DataFileStruct)
%------------------------------------------------------------------------
% TytoLogy:Experiments:optosort
%------------------------------------------------------------------------
% 
% exports raw sweep data for analysis using Plexon. Output is in .nex
% (Neural Explorere) format
%
% Will require toolbox from Neural Explorer:
% https://www.neuroexplorer.com/downloads/HowToReadAndWriteNexAndNex5FilesInMatlab.zip
% 
%------------------------------------------------------------------------
% Input Arguments:
%	With no inputs provided, a dialog will open to specify a file with a 
%	list of .dat files to process (not yet implemented!!!!)
%
% DataFileStruct	<file with list of data files.m>
%	see exportTest.m for example
%
% Output Arguments:
% 	nD		NeuroExplorer nex data struct written to output file
%	nexInfo	struct with information about data written to .nex file:
% 		NexFileName			name of _nexinfo.mat file
% 		fData					struct array with info about data files
% 		sweepStartBin		sample index for sweep start
% 		sweepEndBin			sample index for sweep end
% 		fileStartBin		sample index for data start for each data file
% 		fileEndBin			sample index for data end for each data file
% 		startTimes			start times (seconds) for each sweep
% 		endTimes				end times (seconds) for each sweep
% 		fileStartTime		start times (seconds) for data from each file
% 		fileEndTime			end times (seconds) for data from each file
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
%	15 Jan 2020 (SJS): added documentation, writing info to .mat file
%------------------------------------------------------------------------
% TO DO:
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% Initial things to define
%------------------------------------------------------------------------
sepstr = '----------------------------------------------------';
NEX_UTIL_PATH = ['~/Work/Code/Matlab/stable/Toolbox/NeuroExplorer' ...
					'/HowToReadAndWriteNexAndNex5FilesInMatlab'];
% filter info
defaultFilter.Fc = [300 4000];
defaultFilter.forder = 5;
defaultFilter.ramp = 1;
%------------------------------------------------------------------------
% Setup
%------------------------------------------------------------------------
fprintf('\n%s\n%s\n', sepstr, sepstr);
fprintf('%s running...\n', mfilename);
fprintf('\n%s\n%s\n', sepstr, sepstr);
sendmsg('Checking paths');
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
% Check inputs
%------------------------------------------------------------------------
if nargin == 1
	tmp = varargin{1};
	if ~isstruct(tmp)
		error('%s: input must be a valid sample data struct!');
	end
	% assign values, fixing some things as necessary
	% if file elements are not cells (i.e., just strings), convert to cell.
	if ~iscell(tmp.DataPath)
		tmp.DataPath = {tmp.DataPath};
	else
		tmp.DataPath = tmp.DataPath;
	end
	if ~iscell(tmp.DataFile)
		tmp.DataFile = {tmp.DataFile};
	else
		tmp.DataFile = tmp.DataFile;
	end
	if ~iscell(tmp.TestFile)
		tmp.TestFile = {tmp.TestFile};
	else
		tmp.TestFile = tmp.TestFile;
	end
	% check Channels
	if ~isnumeric(tmp.Channels)
		error('%s: Channels must be a numeric array!', mfilename);
	else
		Channels = tmp.Channels;
	end
	% filter options
	if isfield(tmp, 'BPfilt')
		BPfilt = tmp.BPfilt;
	else
		% use default
		BPfilt = defaultFilter;
	end
	% define path to data file and data file for testing
	F = defineSampleData(tmp.DataPath, tmp.DataFile, tmp.TestFile);
	clear tmp
else
	% define path to data file and data file for testing
	[F, Channels] = defineSampleData();
	% for now use default filter - probably want to have UI for user to
	% specify
	BPfilt = defaultFilter;
end

%------------------------------------------------------------------------
% Get Data File(s) information
%------------------------------------------------------------------------
sendmsg('Data Files:');
% determine # of files
nFiles = length(F);
for f = 1:nFiles
	fprintf('DataFile{%d} = %s\n', f, F(f).file);
end
fprintf('Animal: %s\n', F(1).animal);

% # channel(s) of data to obtain
nChannels = length(Channels);

%------------------------------------------------------------------------
% Read data
%------------------------------------------------------------------------
sendmsg('Reading data');

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
										'sweepLen', [], ...
										'fileStartBin', [], ...
										'fileEndBin', [], ...
										'Dinf', [] ...
								), ...
						1, nFiles);

% loop through files
for f = 1:nFiles
	% save parse file info
	fData(f).F = F(f);
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
	[D, Dinf] = readOptoData(fullfile(F(f).path, F(f).file));
	% Fix test info
	Dinf = correctTestType(Dinf);
% 	% store file info
% 	fData(f).DataPath = DataPath;
% 	fData(f).DataFile = DataFile{f};
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
	[fData(f).cSweeps, ...
		fData(f).startSweepBin, fData(f).endSweepBin, fData(f).sweepLen] = ...
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

sendmsg('Building start and end sweep indices:');
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

%------------------------------------------------------------------------
% Create file and sweep start/end bin and timestamp vectors
%------------------------------------------------------------------------
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

%------------------------------------------------------------------------
% convert (concatenate) cSweeps to vector for each channel, 
% add to nex struct, create event times, add to nex struct, write to 
% .nex file
%------------------------------------------------------------------------
% to save on memory requirements, clear channel data after adding to nex
% data struct.
%------------------------------------------------------------------------
sendmsg('Adding continuous and event data to nex struct:');

% create output .nex file name 
%	assume data from first file is consistent with others!!!!!!!!!
NexFileName = [	fData(1).F.animal '_' ...
						fData(1).F.datecode '_' ...
						fData(1).F.unit '_' ...
						fData(1).F.penetration '_' ...
						fData(1).F.depth ...
						'.nex'];
fprintf('Exporting data to %s\n', NexFileName);
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

sendmsg(sprintf('Writing nex file %s:', NexFileName));
% write to nexfile
writeNexFile(nD, NexFileName);

%------------------------------------------------------------------------
% write useful information to _nexinfo.mat file
%------------------------------------------------------------------------

% create output _nexinfo.mat file name - base is same as .nex file
NexinfoFileName = [	fData(1).F.animal '_' ...
							fData(1).F.datecode '_' ...
							fData(1).F.unit '_' ...
							fData(1).F.penetration '_' ...
							fData(1).F.depth ...
							'_nexinfo.mat'];

% create nexInfo struct to hold sweep/file bin and time data
nexInfo.NexFileName = NexFileName;
% need to remove cSweeps from nexInfo copy of fData to save memory
nexInfo.fData = rmfield(fData, 'cSweeps');
nexInfo.nFiles = nFiles;
nexInfo.Fs = Fs;
nexInfo.sweepStartBin = sweepStartBin;
nexInfo.sweepEndBin = sweepEndBin;
nexInfo.fileStartBin = fileStartBin;
nexInfo.fileEndBin = fileEndBin;
nexInfo.startTimes = startTimes;
nexInfo.endTimes = endTimes;
nexInfo.fileStartTime = fileStartTime;
nexInfo.fileEndTime = fileEndTime;

% save to matfile
sendmsg(sprintf('Writing _nexinfo.mat file %s:', NexinfoFileName));
save(NexinfoFileName, 'nexInfo', '-MAT');

%------------------------------------------------------------------------
% output
%------------------------------------------------------------------------
if nargout
	varargout{1} = nD;
	if nargout > 1
		varargout{2} = nexInfo;
	end
end
%------------------------------------------------------------------------
% END OF MAIN FUNCTION DEFINITION
%------------------------------------------------------------------------
end

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% NESTED FUNCTIONS
%------------------------------------------------------------------------
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%------------------------------------------------------------------------
function varargout = defineSampleData(varargin)
%------------------------------------------------------------------------
%------------------------------------------------------------------------
% defineSampleData.m
%------------------------------------------------------------------------
% defines path to sample data for testing
%------------------------------------------------------------------------
%------------------------------------------------------------------------
	Channels = [];
	
	% check inputs 
	if nargin == 0
		% do something to get data files from user
		DataPath = {};
		DataFile = {};
		TestFile = {};
		Channels = [];
	elseif nargin == 3
		% assign values to DataPath, DataFile and TestFile
		DataPath = varargin{1};
		DataFile = varargin{2};
		TestFile = varargin{3};
	else
		error('%s->defineSampleData: invalid inputs', mfilename)
	end

	% loop through # of data files
	for f = 1:length(DataFile)
		tmpF = parse_opto_filename(DataFile{f});
		if length(DataPath) == 1
			% only 1 element in DataPath so assume all data files are on this
			% path
			tmpF.path = DataPath{1};
		else
			tmpF.path = DataPath{f};
		end
		tmpF.file = DataFile{f}; 
		tmpF.testfile = TestFile{f}; 
		F(f) = tmpF; %#ok<AGROW>
	end
	varargout{1} = F;
	if nargout == 2
		varargout{2} = Channels;
	end
	
end

%------------------------------------------------------------------------
%------------------------------------------------------------------------
function checkstatus = check_sweeps(sweepbins)
%------------------------------------------------------------------------
% checkstatus = check_sweeps(sweepbins)
%------------------------------------------------------------------------
% Check consistency of sweep bins across channels
%------------------------------------------------------------------------
% Input Arguments:
% 	sweepbins		{nchannels, 1} cell array of vectors, where 
% 						each element is an array [1, nsweeps] of sweep bins
% Output Arguments:
% 	checkstatus		true if sweep times are inconsistent, false if all
%						are consistent (equal)
%------------------------------------------------------------------------
	tmp = cell2mat(sweepbins);

	if sum((sum(tmp - tmp(1, :)))) ~= 0
		checkstatus = true;
	else
		checkstatus = false;
	end
end

%------------------------------------------------------------------------
%------------------------------------------------------------------------
function sendmsg(msgstr)
%------------------------------------------------------------------------
% sendmsg(msgstr)
%------------------------------------------------------------------------
% displays message between two separation strings of dashes
%------------------------------------------------------------------------
% Input Arguments:
% 	msgstr		string to display
% Output Arguments:
% 	none
%------------------------------------------------------------------------
	sepstr = '----------------------------------------------------';
	fprintf('%s\n%s\n%s\n', sepstr, msgstr, sepstr);
end
