function varargout = export_for_plexon(varargin)
%------------------------------------------------------------------------
% [nD, nexInfo] = export_for_plexon(exportInfo)
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
% exportInfo	data file and options struct
%	see exportTest.m for examples
%
%----------------------
% 	Fields:
%----------------------
% 	exportInfo.DataPath: 
% 	 Path(s) to data files
% 		can specify individual file paths in a cell array...
% 			exportInfo.DataPath = {'~/Work/Data/TestData/MT_IC';
% 							'~/Work/Data/TestData/MT_IC'; ...
% 							'~/Work/Data/TestData/MT_IC'; };
% 		or single path that applies to all files specified:
% 			exportInfo.DataPath = '~/Work/Data/TestData/MT_IC';
% 			exportInfo.DataPath = 
% 								'/Volumes/Wenstrup Laboratory/By User/SJS/Data/SpikeSort';
% 
% 	exportInfo.DataFile:
% 	 list (cell array) of data files to export (and merge)
% 		exportInfo.DataFile = {'1372_20191126_03_01_1500_FREQ_TUNING.dat'; ...
% 						'1372_20191126_03_01_1500_BBN.dat'; ...
% 						'1372_20191126_03_01_1500_FRA.dat'; ...
% 						'1372_20191126_03_01_1500_WAV.dat'; };
% 
% 	exportInfo.TestFile:
% 	 list (cell array) of test data files corresponding to DataFiles;
%		exportInfo.TestFile = {'1372_20191126_03_01_1500_FREQ_TUNING_testdata.mat'; ...
% 										'1372_20191126_03_01_1500_BBN_testdata.mat'; ...
% 										'1372_20191126_03_01_1500_FRA_testdata.mat'; ...
% 										''; };
% 
%	exportInfo.OutputPath, exportInfo.OutputFile
% 	 you can specify an output path and nex file name, or just leave blank
% 	 and export_plexon_data will create one in current directory
% 		exportInfo.OutputPath = exportOpts.DataPath;
% 		exportInfo.OutputFile = '1372_20191126_03_01_1500_test.nex';
% 
%	exportInfo.Channels:
% 	 neural A/D channels to export from data file(s). 
% 	 note that channels must be consistent across all the data files to be
% 	 merged in the exported file!
% 	 If blank, all channels present in file will be included
% 		exportInfo.Channels = [11 9 14];
% 
% 	exportInfo.BPfilt:
% 	 specifies filter for processing output data.
% 	 if blank/unspecified, no filter will be applied to data
% 		[highpass lowpass] cutoff frequencies in Hz:
% 		 exportOpts.BPfilt.Fc = [300 4000];
% 		order of filter. note that the filtfilt() function in MATLAB is used,
% 		so the effective order is doubled. typically use 5:
% 		 exportOpts.BPfilt.forder = 5;
% 		ramp time (ms) to apply to each sweep in order to cutdown on onset/offset
% 		transients from filtering:
% 		 exportOpts.BPfilt.ramp = 1;
% 
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
%	12 Feb 2020 (SJS): revising for object oriented storage
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
defaultFilter = [];

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
	fprintf('Please add to MATLAB path\n');
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
	
	% checks for output path and file
	if isfield(tmp, 'OutputPath')
		if ~isempty(tmp.OutputPath)
			% if not empty, make sure path exists
			if ~exist(tmp.OutputPath, 'dir')
				error('%s: output directory %s not found', mfilename, tmp.OutputPath);
			end
		end
		NexFilePath = tmp.OutputPath;
	else
		% create empty field
		NexFilePath = '';
	end
	
	if isfield(tmp, 'OutputFile')
		NexFileName = tmp.OutputFile;
	else
		% if no OutputFile field, create empty field
		NexFileName = '';
	end
	
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
% If not provided, create output .nex file name - adjust depending on # of files
%	assume data from first file is consistent with others!!!!!!!!!
%------------------------------------------------------------------------
if isempty(NexFileName)
	if nFiles > 1
		% append MERGE to filename
		NexFileName = [	F(1).fileWithoutOther '_' ...
								'MERGE.nex'];
	else
		NexFileName = [	F(1).base '.nex'];
	end
end

%------------------------------------------------------------------------
% create nexInfo object (SpikeInfo) to hold sweep/file bin and time data
%------------------------------------------------------------------------
nexInfo = SpikeInfo();
nexInfo.FileName = NexFileName;
nexInfo.ADchannel = Channels;

%------------------------------------------------------------------------
% pre-allocate some things
%------------------------------------------------------------------------
% bins for start and end of each file's data
nexInfo.fileStartBin = zeros(1, nFiles);
nexInfo.fileEndBin = zeros(1, nFiles);
% bins for all sweep starts and ends
nexInfo.sweepStartBin = cell(1, nFiles);
nexInfo.sweepEndBin = cell(1, nFiles);
% each file's sampling rate for neural data
tmpFs = zeros(nFiles, 1);
% cell array to hold sweep data
cSweeps = cell(nFiles, 1);
% struct to hold everything for each file
fData = repmat(	struct(		'startSweepBin', {}, ...
										'endSweepBin', {}, ...
										'sweepLen', [], ...
										'fileStartBin', [], ...
										'fileEndBin', [], ...
										'Dinf', [] ...
								), ...
						1, nFiles);


%------------------------------------------------------------------------
% Read data
%------------------------------------------------------------------------
sendmsg('Reading data');

% loop through files
for f = 1:nFiles
	% save file info object for current data file
	fData(f).F = F(f);
	
	% get data for each file and channel and convert to row vector format
	% algorithm:
	%		(1) put each sweep for this channel in a {1, # sweeps} cell array
	%				cSweeps
	%		(2) make a note of the length of each sweep to use for
	%				markers/timestamps
	%		(3) after cSweeps is built, convert to a row vector using cell2mat

	% use readOptoData to read in raw data. 
	[D, Dinf] = readOptoData(fullfile(F(f).path, F(f).file));
	% Fix test info
	Dinf = correctTestType(Dinf);

	% build filter for neural data
	if ~isempty(BPfilt)
		BPfilt.Fs = Dinf.indev.Fs;
		BPfilt.Fnyq = Dinf.indev.Fs / 2;
		BPfilt.cutoff = BPfilt.Fc / BPfilt.Fnyq;
		[BPfilt.b, BPfilt.a] = butter(BPfilt.forder, BPfilt.cutoff, 'bandpass');
	end
	
	% check to make sure consistent # of sweeps (aka trials)
	if Dinf.nread ~= length(D)
		error('%s: mismatch in Dinf.nread (%d) and length(D) (%d)', ...
					mfilename, Dinf.nread, length(D));
	end
	
	% build into sweeps by channel format
	fprintf('Test type: %s\n', Dinf.test.Type);
	[cSweeps{f}, ...
		fData(f).startSweepBin, fData(f).endSweepBin, fData(f).sweepLen] = ...
					buildChannelData(Channels, BPfilt, D, Dinf);
	% check the start and end sweep bin data for consistency
	if check_sweeps(fData(f).startSweepBin)
		warning(['File %s: Inconsistent startSweepBin' ...
							'values across channels!!!!'], fData(f).F.file);
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
	% to do this, store list of sample rates and check once out of this loop
	tmpFs(f) = fData(f).Dinf.indev.Fs;
	
	% calculate overall start and end bins for each file's data
	if f == 1
		nexInfo.fileStartBin(f) = fData(f).fileStartBin;
		nexInfo.fileEndBin(f) = fData(f).fileEndBin;
	else
		% add 1 to prior end bin for start
		nexInfo.fileStartBin(f) = nexInfo.fileEndBin(f-1) + 1;
		nexInfo.fileEndBin(f) = nexInfo.fileStartBin(f) + fData(f).fileEndBin - 1;
	end
end

sendmsg('Building start and end sweep indices:');
% assign values for bins
for f = 1:nFiles
	% calculate start and end sweep bins for each file's data
	if f == 1
		nexInfo.sweepStartBin{f} = fData(f).startSweepBin{1};
		nexInfo.sweepEndBin{f} = fData(f).endSweepBin{1};
	else
		% add previous file's final endSweepBin value as offset
		nexInfo.sweepStartBin{f} = fData(f).startSweepBin{1} + ...
												nexInfo.sweepEndBin{f-1}(end);
		nexInfo.sweepEndBin{f} = fData(f).endSweepBin{1} + ...
												nexInfo.sweepEndBin{f-1}(end);
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
	nexInfo.Fs = Fs;
end

% % convert file start/end bins to times
% fileStartTime = (fileStartBin - 1) ./ Fs;
% fileEndTime = (fileEndBin - 1) ./ Fs;

% % convert sweep bin cells to vectors... 
% startBinVector = [sweepStartBin{:}];
% endBinVector = [sweepEndBin{:}];
% % ... and then to times... these will be written to the nex
% % file as event timestamps
% startTimes = (startBinVector - 1) ./ Fs;
% endTimes = (endBinVector - 1) ./ Fs;

%------------------------------------------------------------------------
% convert (concatenate) cSweeps to vector for each channel, 
% add to nex struct, create event times, add to nex struct, write to 
% .nex file
%------------------------------------------------------------------------
% to save on memory requirements, clear channel data after adding to nex
% data struct.
%------------------------------------------------------------------------
sendmsg('Adding continuous and event data to nex struct:');


fprintf('Exporting data to %s\n', fullfile(NexFilePath, NexFileName));

% start new nex file data struct
nD = nexCreateFileData(Fs);

% loop through channels
for c = 1:nChannels
	% this will be a matrix of format
	% 	[# channels, (# sweeps) * (# samples per sweep)
	cVector = cell(1, nFiles);
	for f = 1:nFiles
		cVector{1, f} = cSweeps{f}(c, :);
	end
	% concatenate cell array, convert to vector, add to nex struct
	% steps:
	%	concatenate: tmp = [cVector{:}];
	%	tmpVector = cell2mat(tmp);
	%  add to nex struct:
	%	[nexFile] = nexAddContinuous( nexFile, startTime, adFreq, values, name)
	nD = nexAddContinuous(nD, nexInfo.fileStartTime(1), Fs, ...
									cell2mat([cVector{:}]), ...
									sprintf('spikechan_%d', Channels(c)));
	% clear cVector to save memory
	clear cVector
end

% add start sweep time stamps as event - assume consistent across channels!
%  [nexFile] = nexAddEvent( nexFile, timestamps, name )
% events must be in column format...?
nD = nexAddEvent(nD, force_col(nexInfo.startTimeVector), 'startsweep');
% add end sweep time stamps as event - assume consistent across channels!
nD = nexAddEvent(nD, force_col(nexInfo.endTimeVector), 'endsweep');
% add file times
nD = nexAddEvent(nD, force_col(nexInfo.fileStartTime), 'filestart');
nD = nexAddEvent(nD, force_col(nexInfo.fileEndTime), 'fileend');

sendmsg(sprintf('Writing nex file %s:', fullfile(NexFilePath, NexFileName)));
% write to nexfile
writeNexFile(nD, fullfile(NexFilePath, NexFileName));

%------------------------------------------------------------------------
% write useful information to _nexinfo.mat file 
%------------------------------------------------------------------------

% create output _nexinfo.mat file name - base is same as .nex file
[~, nibase] = fileparts(NexFileName);
NexinfoFileName = [nibase '_nexinfo.mat'];

%{
%%%% pre OOP
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
nexInfo.Channels = Channels;
%}

% create nexInfo object (SpikeInfo) to hold sweep/file bin and time data
% nexInfo = SpikeInfo();
% nexInfo.FileName = NexFileName;
% need to remove cSweeps from nexInfo copy of fData to save memory
nexInfo.FileData = fData;
% nexInfo.Fs = Fs;
% nexInfo.sweepStartBin = sweepStartBin;
% nexInfo.sweepEndBin = sweepEndBin;
% nexInfo.fileStartBin = fileStartBin;
% nexInfo.fileEndBin = fileEndBin;
% nexInfo.startBinVector = startBinVector;
% nexInfo.endBinVector = endBinVector;

nexInfo.dataFilter = BPfilt;
% save to matfile
sendmsg(sprintf('Writing _nexinfo.mat file %s:', ...
											fullfile(NexFilePath, NexinfoFileName)));
save(fullfile(NexFilePath, NexinfoFileName), 'nexInfo', '-MAT');

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
		% do something to get data files from user !!!!
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

	% loop through # of data files, create file objects
	for f = 1:length(DataFile)
		if length(DataPath) == 1
			% only 1 element in DataPath so assume all data files are on this
			% path
			F(f) = OptoFileName(fullfile(DataPath{1}, DataFile{f})); %#ok<AGROW>
		else
			F(f) = OptoFileName(fullfile(DataPath{f}, DataFile{f})); %#ok<AGROW>
		end
		F(f).testfile = TestFile{f};  %#ok<AGROW>
	end	
	
	% assign outputs
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
	[nr, ~] = size(tmp);
	tmp2 = 0;
	for r = 1:nr
		tmp2 = tmp2 + sum(tmp(r, :) - tmp(1, :));
	end
	
	if tmp2 ~= 0
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
