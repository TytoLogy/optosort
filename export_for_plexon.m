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
%		filter type is either 'bessel' or 'butter' (butterworth)
%		 exportOpts.BPfilet.type = 'bessel';
% 
%	exportInfo.resampleData
% 		resample data to nearest-lower integer value?
%		Default: true
% 		For exporting to many software packages for sorting, integer values are
% 		needed.
% 		if true, A/D data will be resampled so that sampling rate is an integer
% 		value. i.e., 48828.125 will be converted to 48828 samples/second
% 		if false, no change will be made to sampling rate.
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
% Created: 14 January, 2020 (SJS)
%
% Revisions:
%	15 Jan 2020 (SJS): added documentation, writing info to .mat file
%	12 Feb 2020 (SJS): revising for object oriented storage
%	3 Mar 2020 (SJS): converted fData to CurveInfo array
%	26 Mar 2020 (SJS): reworking for real world use. 
%		- added filter 'type' to exportOpts.BPfilt struct
%	10 Apr, 2020 (SJS): cInfo converted to cell array instead of array of
%	CurveInfo objects - this was done in order to include WAVInfo objects
%	(which is subclass of CurveInfo) to be included in the array. will need
%	to update downstream code.
% 9 June 2020 (SJS): added resampleData field to exportOptions
%------------------------------------------------------------------------
% TO DO:
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% Initial things to define
%------------------------------------------------------------------------
sepstr = '----------------------------------------------------';
NEX_UTIL_PATH = ['~/Work/Code/Matlab/stable/Toolbox/NeuroExplorer' ...
					'/NexTools'];
% filter info
defaultFilter = [];
% resample data? default is true
resampleData = true;

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
	% resample?
	if isfield(tmp, 'resampleData')
		resampleData = tmp.resampleData;
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
% get the data and information from the raw files
%------------------------------------------------------------------------
[cSweeps, nexInfo] = read_data_for_export(F, Channels, BPfilt, resampleData);

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
nexInfo.FileName = fullfile(NexFilePath, NexFileName);
% create output _nexinfo.mat file name - base is same as .nex file
[~, baseName] = fileparts(nexInfo.FileName);
nexInfo.InfoFileName = fullfile(NexFilePath, [baseName '_nexinfo.mat']);


%------------------------------------------------------------------------
% convert (concatenate) cSweeps to vector for each channel, 
% add to nex struct, create event times, add to nex struct, write to 
% .nex file
%------------------------------------------------------------------------
% to save on memory requirements, clear channel data after adding to nex
% data struct.
%------------------------------------------------------------------------
sendmsg('Adding continuous and event data to nex struct:');
fprintf('Exporting data to %s\n', nexInfo.FileName);

% start new nex file data struct
nD = nexCreateFileData(nexInfo.Fs);

% loop through channels
for c = 1:nChannels
	% this will be a matrix of format
	% 	[# channels, (# sweeps) * (# samples per sweep)
	cVector = cell(1, nFiles);
	for f = 1:nFiles
		cVector{1, f} = cSweeps{f}(c, :);
	end
	% for each channels's data, concatenate cell array, convert to vector, 
	% add to nex struct
	% steps:
	%	concatenate: tmp = [cVector{:}];
	%	tmpVector = cell2mat(tmp);
	%  add to nex struct:
	%	[nexFile] = nexAddContinuous( nexFile, startTime, adFreq, values, name)
	nD = nexAddContinuous(nD, nexInfo.fileStartTime(1), nexInfo.Fs, ...
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
% add stimulus onset and offset times
nD = nexAddEvent(nD, force_col(nexInfo.stimStartTimeVector), 'stimstart');
nD = nexAddEvent(nD, force_col(nexInfo.stimEndTimeVector), 'stimend');
% write to nexfile
sendmsg(sprintf('Writing nex file %s:', nexInfo.FileName));
writeNexFile(nD, nexInfo.FileName);

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
	varargout{1} = nD;
	if nargout > 1
		varargout{2} = nexInfo;
% 		varargout{3} = cInfo;
	end
end
%------------------------------------------------------------------------
% END OF MAIN FUNCTION DEFINITION
%------------------------------------------------------------------------

