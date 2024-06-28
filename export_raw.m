function varargout = export_raw(varargin)
%------------------------------------------------------------------------
% [nD, expInfo] = export_raw(exportInfo)
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
% 		Order of filter. note that the filtfilt() function in MATLAB is used,
% 		so the effective order is doubled. typically use 5:
% 		 exportOpts.BPfilt.forder = 5;
% 		ramp time (ms) to apply to each sweep in order to cutdown on onset/offset
% 		transients from filtering:
% 		 exportOpts.BPfilt.ramp = 1;
%		Filter type is either 'bessel' or 'butter' (butterworth)
%		 exportOpts.BPfilt.type = 'bessel';
% 
%	exportInfo.resampleData
% 		change data sampling rate to resampleData rate (samples/second)
%		Default: no resampling
% 		For exporting to many software packages for sorting, integer values 
% 		are needed.
% 		If value is provided, A/D data will be resampled i.e., original
%		rate of 48828.125 samples/second will be converted to 
%     resampleData value 
% 		If empty or not specified, no change will be made to sampling rate.
% 
%  exportInfo.referenceData
%     Default: 'raw'
%     Options:
%        'raw'    no common reference applied (but data will be filtered 
%                 as set by BPfilt option)
%        'avg'    common average reference - at each bin, computes 
%                 average across channels, subtracts this value from 
%                 each channel
%        'med'    common median reference - at each bin, computes 
%                 median across channels, subtracts this value from 
%                 each channel
% 
%  exportInfo.OutputShape
%     specify output data "shape" (i.e., orientation of channels 
%     and samples in output data matrix) 
%          'ChannelsSamples':
%              [channels, samples]  default, used by SpykingCircus
%          'SamplesChannels':
%              [samples, channels]  used by SpikeInterface, maybe
%                                   kilosort(?)
%     https://spikeinterface.readthedocs.io/en/latest/how_to/load_matlab_data.html
% 
% Output Arguments:
% 	nD		NeuroExplorer nex data struct written to output file
%	expInfo	SpikeInfo object=
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
% 9 June 2020 (SJS): added resampleData field to exportOptions
% 1 May 2024 (SJS): adding option OutputShape to specify output data in 
%                   format:
%                    [channels, samples]   default, used by SpykingCircus
%                    [samples, channels]   used by SpikeInterface, maybe
%                                          kilosort(?)
% https://spikeinterface.readthedocs.io/en/latest/how_to/load_matlab_data.html
% 28 Jun 2024 (SJS): adding common avg/median referencing option (code used
%                    from export_for_plexon
%------------------------------------------------------------------------
% TO DO:
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% Initial things to define
%------------------------------------------------------------------------
%------------------------------------------------------------------------
sepstr = '----------------------------------------------------';
% filter info
defaultFilter = [];
% resample data? default is no
defaultResampleRate = [];
% valid output formats (binary)
ValidFormats =	{	'int16', 'uint16', 'int8', ...
						'float32', 'float64', 'single', 'double'};
OutputShape = 'ChannelsSamples';
ValidOutputShapes = {'ChannelsSamples', 'SamplesChannels'};
% default for common referencing is none (use raw data)
defaultreferenceData = 'raw';

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% Setup
%------------------------------------------------------------------------
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
%------------------------------------------------------------------------
% Check inputs
%------------------------------------------------------------------------
%------------------------------------------------------------------------
if nargin == 1
	tmp = varargin{1};
	if ~isstruct(tmp)
		error('%s: input must be a valid sample data struct!');
	end
   %---------------------------------------------------------------------
   %---------------------------------------------------------------------
	% assign values, fixing some things as necessary
   %---------------------------------------------------------------------
   %---------------------------------------------------------------------

   %---------------------------------------------------------------------
   % if file elements are not cells (i.e., just strings), convert to cell.
   %---------------------------------------------------------------------
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
   %---------------------------------------------------------------------
	% check Channels
   %---------------------------------------------------------------------
	if ~isnumeric(tmp.Channels)
		error('%s: Channels must be a numeric array!', mfilename);
	else
		Channels = tmp.Channels;
	end
   %---------------------------------------------------------------------
	% check output format
   %---------------------------------------------------------------------
	if isfield(tmp, 'OutputFormat')
		if ~any(strcmpi(tmp.OutputFormat, ValidFormats))
			error('%s: invalid output format %s\n', mfilename, tmp.OutputFormat);
		else
			OutputFormat = tmp.OutputFormat;
		end
	else
		OutputFormat = 'float32';
	end
   %---------------------------------------------------------------------
	% filter options
   %---------------------------------------------------------------------
	if isfield(tmp, 'BPfilt')
		BPfilt = tmp.BPfilt;
	else
		% use default
		BPfilt = defaultFilter;
	end
   %---------------------------------------------------------------------
	% resample?
   %---------------------------------------------------------------------
	if isfield(tmp, 'resampleData')
		resampleData = tmp.resampleData;
	else
		resampleData = defaultResampleRate;
	end
   %---------------------------------------------------------------------	
   % apply reference?
   %---------------------------------------------------------------------	
   if isfield(tmp, 'referenceData')
      referenceData = tmp.referenceData;
   else
      referenceData = defaultreferenceData;
   end
   %---------------------------------------------------------------------	
	% define path to data file and data file for testing
   %---------------------------------------------------------------------
	F = defineSampleData(tmp.DataPath, tmp.DataFile, tmp.TestFile);
   %---------------------------------------------------------------------	
	% check for output path and file
   %---------------------------------------------------------------------
	if isfield(tmp, 'OutputPath')
		if ~isempty(tmp.OutputPath)
			% if not empty, make sure path exists
			if ~exist(tmp.OutputPath, 'dir')
            % create directory
				fprintf('%s: output directory\n %s\nnot found\n', ...
                           mfilename, ...
                           tmp.OutputPath);
            fprintf('Creating directory\n');
            mkdir(tmp.OutputPath);
			end
		end
		OutputPath = tmp.OutputPath;
	else
		% create empty field (use current directory)
		OutputPath = '';
	end
   %---------------------------------------------------------------------
	% check for output file
   %---------------------------------------------------------------------
   if isfield(tmp, 'OutputFile')
      RawFileName = tmp.OutputFile;
   else
      % if no OutputFile field, create empty field
      RawFileName = '';
   end
   %---------------------------------------------------------------------
   % output data shape ([channels, samples] or [samples, channels])
   %---------------------------------------------------------------------
   if isfield(tmp, 'OutputShape')
      if any(strcmpi(tmp.OutputShape, ValidOutputShapes))
         OutputShape = tmp.OutputShape;
      else
         error('%s: invalid OutputShape: %s', mfilename, tmp.OutputShape);
      end
   end
   %---------------------------------------------------------------------
   % clear tmp var
   %---------------------------------------------------------------------
	clear tmp
else
   %---------------------------------------------------------------------
	% define path to data file and data file for testing
   %---------------------------------------------------------------------
	[F, Channels] = defineSampleData();
   %---------------------------------------------------------------------
	% for now use default filter - probably want to have UI for user to
	% specify
   %---------------------------------------------------------------------
	BPfilt = defaultFilter;
   %---------------------------------------------------------------------
	% default output format is 32 bit floating point
   %---------------------------------------------------------------------
	OutputFormat = 'float32';
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
[cSweeps, expInfo] = read_data_for_export(F, Channels, ...
                                                BPfilt, resampleData);

%------------------------------------------------------------------------
% Apply common reference here, save reference mode in expInfo struct
%------------------------------------------------------------------------
switch referenceData
   case {'raw', ''}
      % do nothing
      expInfo.referenceMode = 'raw';
   case 'avg'
      expInfo.referenceMode = 'average';
      cSweeps = applyCommonReference(cSweeps, @common_avg_ref);
   case 'med'
      expInfo.referenceMode = 'median';
      cSweeps = applyCommonReference(cSweeps, @common_med_ref);
   otherwise
      error('%s: invalid referenceData option: %s', ...
                  mfilename, referenceData);
end

%------------------------------------------------------------------------
% If not provided, create output .bin file name - adjust depending 
% on # of files
% Assume data from first file is consistent with others!!!!!!!!!
%------------------------------------------------------------------------
if isempty(RawFileName)
	if nFiles > 1
		% append MERGE to filename
		RawFileName = [	F(1).fileWithoutOther '_' ...
								'MERGE.bin'];
	else
		RawFileName = [	F(1).base '.bin'];
	end
end
expInfo.FileName = fullfile(OutputPath, RawFileName);
% create output _nexinfo.mat file name - base is same as .nex file
[~, baseName] = fileparts(expInfo.FileName);
expInfo.InfoFileName = fullfile(OutputPath, [baseName '_info.mat']);

%------------------------------------------------------------------------
% create spyking params file
%------------------------------------------------------------------------
% get default params struct
% filename is <output path>/<filename>.params
params = default_spyking_params( ...
                              fullfile(OutputPath, [baseName '.params']));
params.Fs = expInfo.Fs;
params.OutputFormat = OutputFormat;
params.nChannels = nChannels;
params.DataOffset = 0;
params.DataOffsetFormat = 'auto';
params.Gain = 1;

%------------------------------------------------------------------------
% need to write parameters to config file
%------------------------------------------------------------------------
sendmsg(sprintf('Writing configuration data to %s', params.Filename));
write_spyking_params(params);

%------------------------------------------------------------------------
% convert (concatenate) cSweeps to [nchannels, nsamples] matrix and write
% to binary output file
%------------------------------------------------------------------------
sendmsg(sprintf('Exporting raw binary data to %s', expInfo.FileName));

% open raw file
fp = fopen(expInfo.FileName, 'wb');

% loop through files
for f = 1:nFiles
	fprintf('Writing data for %s\n', F(f).file);
	% write data for this file, specify shape based on outputShape setting
   % 1 May 2024: this may not have any effect since binary files store
   % everything in 1D format!!!!
   switch OutputShape
      case 'ChannelsSamples'
      	nw = fwrite(fp, cell2mat([cSweeps{f}]), OutputFormat);
      case 'SamplesChannels'
         nw = fwrite(fp, cell2mat([cSweeps{f}])', OutputFormat);
   end

	% could also:
	%	concatenate: tmp = [cVector{:}];
	%	fwrite(fp, cell2mat(tmp), 'float64');
	fprintf('\t%d values written\n', nw);
end
fclose(fp);

%------------------------------------------------------------------------
% write useful information to _nexinfo.mat file 
%------------------------------------------------------------------------
% save to matfile
sendmsg(sprintf('Writing _nexinfo.mat file %s:', ...
											expInfo.InfoFileName));
save(expInfo.InfoFileName, 'expInfo', '-MAT');

%------------------------------------------------------------------------
% output
%------------------------------------------------------------------------
if nargout
	varargout{1} = cSweeps;
	if nargout > 1
		varargout{2} = expInfo;
      varargout{3} = nw;
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
function varargout = write_spyking_params(varargin)
%------------------------------------------------------------------------
%------------------------------------------------------------------------
% write_spyking_params
%------------------------------------------------------------------------
% writes params file for Spyking Circus sorting
% RAW_BINARY (read/parallel write)
% 
% | The parameters for RAW_BINARY file format are:
% |
% | -- sampling_rate -- <type 'float'> [** mandatory **]
% | -- data_dtype -- <type 'str'> [** mandatory **]
% | -- nb_channels -- <type 'int'> [** mandatory **]
% |
% | -- data_offset -- <type 'int'> [default is 0]
% | -- dtype_offset -- <type 'str'> [default is auto]
% | -- gain -- <type 'int'> [default is 1]
%------------------------------------------------------------------------
%------------------------------------------------------------------------
	if isempty(varargin)
		params = default_spyking_params;
	elseif isstruct(varargin{1})
		params = varargin{1};
	else
		error('%s: input to write_spyking_params must be a struct or empty', ...
					mfilename);
	end
	if isempty(params.Filename)
		% write to stdout
		fp = 1;
	else
		fp = fopen(params.Filename, 'w');
	end
	
	fprintf(fp, '[data]\n');
	fprintf(fp, 'file_format   = raw_binary\n');
	% samples/second
	fprintf(fp, 'sampling_rate = %.4f\n', params.Fs);
	% should be int16,uint16,float32,...
	fprintf(fp, 'data_dtype    = %s\n', params.OutputFormat);
	% # of channels of output data (for demux the data)
	fprintf(fp, 'nb_channels   = %d\n', params.nChannels);
	% Optional, if a header with a fixed size is present
	fprintf(fp, 'data_offset   = %d\n', params.DataOffset);
	% Optional, if a header with a fixed size is present
	fprintf(fp, 'dtype_offset  = %s\n', params.DataOffsetFormat);
	% Optional, if you want a non unitary gain for the channels
	fprintf(fp, 'gain          = %d\n', params.Gain);
	if fp ~= 1
		fclose(fp);
	end
	varargout{1} = params;
end

%------------------------------------------------------------------------
%------------------------------------------------------------------------
function params = default_spyking_params(varargin)
%------------------------------------------------------------------------
% params = default_spyking_params(varargin)
%------------------------------------------------------------------------
% create default params file for Spyking Circus sorting
%------------------------------------------------------------------------
% Input: (optional) params filename
%------------------------------------------------------------------------
	if isempty(varargin)
		params.Filename = 'config_default.params';
	else
		params.Filename = varargin{1};
	end
	params.Fs = 48828.1250;
	params.OutputFormat = 'float32';
	params.nChannels = 16;
	params.DataOffset = 0;
	params.DataOffsetFormat = 'auto';
	params.Gain = 1;
end

