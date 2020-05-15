function varargout = export_raw(varargin)
%------------------------------------------------------------------------
% [nD, nexInfo] = export_raw(exportInfo)
%------------------------------------------------------------------------
% TytoLogy:Experiments:optosort
%------------------------------------------------------------------------
% 
% exports raw sweep data for analysis. output is in raw binary format
%
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
%		filter type is either 'bessel' or 'butter' (butterworth)
%		 exportOpts.BPfilet.type = 'bessel';
% 
%
% Output Arguments:
% 	nD		NeuroExplorer nex data struct written to output file
%	nexInfo	SpikeInfo object

% OLD:
%     nexInfo struct with information about data written to .nex file:
% 		NexFileName			name of _nexinfo.mat file
% 		fData					array of CurveInfo objects with info about data files
% 		sweepStartBin		sample index for sweep start
% 		sweepEndBin			sample index for sweep end
% 		fileStartBin		sample index for data start for each data file
% 		fileEndBin			sample index for data end for each data file
% 		startTimes			start times (seconds) for each sweep
% 		endTimes				end times (seconds) for each sweep
% 		fileStartTime		start times (seconds) for data from each file
% 		fileEndTime			end times (seconds) for data from each file
%------------------------------------------------------------------------
% See also: SpikeInfo, CurveInfo, WAVInfo classes
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%  Sharad J. Shanbhag
%	sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 14 May, 2020  from export_for_plexon (SJS)
%
% Revisions:
%------------------------------------------------------------------------
% TO DO:
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% Initial things to define
%------------------------------------------------------------------------
sepstr = '----------------------------------------------------';
% filter info
defaultFilter = [];

%------------------------------------------------------------------------
% Setup
%------------------------------------------------------------------------
fprintf('\n%s\n%s\n', sepstr, sepstr);
fprintf('%s running...\n', mfilename);
fprintf('\n%s\n%s\n', sepstr, sepstr);
sendmsg('Checking paths');
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
		RawFileName = tmp.OutputFile;
	else
		% if no OutputFile field, create empty field
		RawFileName = '';
	end
	
	clear tmp
else
	% define path to data file and data file for testing
	[F, Channels] = defineSampleData();
	% for now use default filter - probably want to have UI for user to
	% specify
	BPfilt = defaultFilter;
end

% # channel(s) of data to obtain
nChannels = length(Channels);
% Get Data File(s) information
sendmsg('Data Files:');
% determine # of files
nFiles = length(F);
% let user know about files
for f = 1:nFiles
	fprintf('DataFile{%d} = %s\n', f, F(f).file);
end
fprintf('Animal: %s\n', F(1).animal);

%------------------------------------------------------------------------
% If not provided, create output .nex file name - adjust depending on # of files
%	assume data from first file is consistent with others!!!!!!!!!
%------------------------------------------------------------------------
if isempty(RawFileName)
	if nFiles > 1
		% append MERGE to filename
		RawFileName = [	F(1).fileWithoutOther '_' ...
								'MERGE.nex'];
	else
		RawFileName = [	F(1).base '.bin'];
	end
end

%------------------------------------------------------------------------
% create nexInfo object (SpikeInfo) to hold sweep/file bin and time data
%------------------------------------------------------------------------
nexInfo = SpikeInfo();
nexInfo.FileName = fullfile(NexFilePath, RawFileName);
% create output _nexinfo.mat file name - base is same as .nex file
[~, nibase] = fileparts(nexInfo.FileName);
nexInfo.InfoFileName = fullfile(NexFilePath, [nibase '_nexinfo.mat']);
clear nibase
% store channel information
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
% bins for stim onset, offset
nexInfo.stimStartBin = cell(1, nFiles);
nexInfo.stimEndBin = cell(1, nFiles);

% each file's sampling rate for neural data
tmpFs = zeros(nFiles, 1);
% cell array to hold sweep data - this will be converted to a single
% "vector" of values per channel that will be added to the .nex file
% (algorithm will be outlined in the next section)
cSweeps = cell(nFiles, 1);
% allocate cell array to hold information for each file
% we'll assign/allocate/initialize objects (CurveInfo, WAVINfo) as we loop
% through the files in the next section
cInfo = cell(nFiles, 1);

%------------------------------------------------------------------------
% Read and process data from raw files
%------------------------------------------------------------------------

% loop through files
for f = 1:nFiles
	sendmsg(sprintf('Reading data for file %s', F(f).file));
	
	% get data for each file and channel and convert to row vector format
	% algorithm:
	%		(1) put each sweep for this channel in a {1, # sweeps} cell array
	%				cSweeps
	%		(2) make a note of the length of each sweep to use for
	%				markers/timestamps
	%		(3) after cSweeps is built, convert to a row vector using cell2mat

	% use readOptoData to read in raw data. 
	[D, tmpDinf] = readOptoData(fullfile(F(f).path, F(f).file));
	% Fix test info
	tmpDinf = correctTestType(tmpDinf);
	
	% CurveInfo (or WavInfo) to hold everything for each file
	switch(upper(tmpDinf.test.Type))
		case {'FREQ', 'LEVEL', 'FREQ+LEVEL', 'OPTO'}
			cInfo{f} = CurveInfo(tmpDinf);
			
		case {'WAV', 'WAVFILE'}
			% build wavinfo file and load it
			wavinfo_filename = [F(f).base '_wavinfo.mat'];
			if ~exist(fullfile(F(f).path,wavinfo_filename), 'file')
				warning('%s: wavinfo file %s not found', mfilename, ...
																wavinfo_filename);
				cInfo{f} = WAVInfo(tmpDinf);
			else
				tmpW = load(fullfile(F(f).path,wavinfo_filename));
				cInfo{f} = WAVInfo(tmpDinf, tmpW);
			end
	
		otherwise
			error('%s: unknown test.Type %s', mfilename, tmpDinf.test.Type);
	end
	
	% build filter for neural data
	if ~isempty(BPfilt)
		BPfilt.Fs = cInfo{f}.Dinf.indev.Fs;
		BPfilt.Fnyq = cInfo{f}.Dinf.indev.Fs / 2;
		BPfilt.cutoff = BPfilt.Fc / BPfilt.Fnyq;
		if strcmpi(BPfilt.type, 'bessel')
			[BPfilt.b, BPfilt.a] = besself(BPfilt.forder, BPfilt.cutoff, 'bandpass');
		elseif strcmpi(BPfilt.type, 'butter')
			[BPfilt.b, BPfilt.a] = butter(BPfilt.forder, BPfilt.cutoff, 'bandpass');
		else
			error('%s: unknown filter type %s', mfilename, BPfilt.type)
		end
	end
	
	% check to make sure consistent # of sweeps (aka trials)
	if cInfo{f}.Dinf.nread ~= length(D)
		error('%s: mismatch in Dinf.nread (%d) and length(D) (%d)', ...
					mfilename, cInfo{f}.Dinf.nread, length(D));
	end
	
	% build into sweeps by channel format
	fprintf('Test type: %s\n', cInfo{f}.testtype);
	[cInfo{f}, cSweeps{f}] = cInfo{f}.buildChannelData(Channels, BPfilt, D);
	
	% store sample for start of this file (should be 1); use channel 1 value
	 cInfo{f}.fileStartBin = cInfo{f}.startSweepBin{1}(1);
	% store sample for end of this file
	 cInfo{f}.fileEndBin = cInfo{f}.endSweepBin{1}(end);

	% to avoid any issues, should make sure sample rates are consistent
	% to do this, store list of sample rates and check once out of this loop
	tmpFs(f) = cInfo{f}.Dinf.indev.Fs;
	
	% calculate overall start and end bins for each file's data
	if f == 1
		nexInfo.fileStartBin(f) = cInfo{f}.fileStartBin;
		nexInfo.fileEndBin(f) = cInfo{f}.fileEndBin;
	else
		% add 1 to prior end bin for start
		nexInfo.fileStartBin(f) = nexInfo.fileEndBin(f-1) + 1;
		nexInfo.fileEndBin(f) = nexInfo.fileStartBin(f) + cInfo{f}.fileEndBin - 1;
	end
	
	% compute stimulus onset/offset bins for this file
	cInfo{f} = cInfo{f}.buildStimOnOffData;
		
end

% assign cInfo to nexInfo.FileData
nexInfo.FileInfo = cInfo;
% store filter info
nexInfo.dataFilter = BPfilt;

sendmsg('Checking sample rates across files');
% check sampling rates
if ~all(tmpFs(1) == tmpFs)
	error('%s: Sample Rate Mismatch!!!', mfilename);
else
	% store overall sample rate
	Fs = tmpFs(1);
	nexInfo.Fs = Fs;
	sendmsg(sprintf('Neural A/D Fs = %.4f', Fs));
end

%------------------------------------------------------------------------
% Create file, sweep, stimulus start/end bin and timestamp vectors
%------------------------------------------------------------------------
sendmsg('Building start and end sweep indices:');
% assign values for bins
for f = 1:nFiles
	% calculate start and end sweep bins for each file's data
	if f == 1
		nexInfo.sweepStartBin{f} = cInfo{f}.startSweepBin{1};
		nexInfo.sweepEndBin{f} = cInfo{f}.endSweepBin{1};
	else
		% add previous file's final endSweepBin value as offset
		nexInfo.sweepStartBin{f} = cInfo{f}.startSweepBin{1} + ...
												nexInfo.sweepEndBin{f-1}(end);
		nexInfo.sweepEndBin{f} = cInfo{f}.endSweepBin{1} + ...
												nexInfo.sweepEndBin{f-1}(end);
	end
	
	% for stim onset/offset, align to file start bin *** would this work for
	% the sweep start and end computation above???? need to test!!!!!
	nexInfo.stimStartBin{f} = nexInfo.fileStartBin(f) + cInfo{f}.stimStartBin;
	nexInfo.stimEndBin{f} = nexInfo.fileStartBin(f) + cInfo{f}.stimEndBin;
end

%------------------------------------------------------------------------
% convert (concatenate) cSweeps to vector for each channel, 
% add to nex struct, create event times, add to nex struct, write to 
% .nex file
%------------------------------------------------------------------------
% to save on memory requirements, clear channel data after adding to nex
% data struct.
%------------------------------------------------------------------------
sendmsg(sprintf('Exporting data to %s', nexInfo.FileName));

% open raw file
fp = fopen(nexInfo.InfoFileName, 'wb');

% loop through files
for f = 1:nFiles
	% write data for this file
	fwrite(fp, cell2mat([cSweeps{f}]), 'float64')
	% could also:
	%	concatenate: tmp = [cVector{:}];
	%	fwrite(fp, cell2mat(tmp), 'float64');
end
fclose(fp);

%------------------------------------------------------------------------
% write useful information to _nexinfo.mat file 
%------------------------------------------------------------------------
% save to matfile
sendmsg(sprintf('Writing _nexinfo.mat file %s:', ...
											nexInfo.InfoFileName));
save(nexInfo.InfoFileName, 'nexInfo', '-MAT');

%------------------------------------------------------------------------
% output
%------------------------------------------------------------------------
if nargout
	varargout{1} = cSweeps;
	if nargout > 1
		varargout{2} = nexInfo;
		varargout{3} = cInfo;
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
