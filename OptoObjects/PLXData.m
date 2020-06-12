classdef PLXData
%------------------------------------------------------------------------
% TytoLogy:Experiments:optoproc:PLXDdata class
%------------------------------------------------------------------------
% Class to encapsulate data and code that uses the readPLXFileC()
% function by Benjamin Kraus (bkraus@bu.edu, ben@benkraus.com)
%------------------------------------------------------------------------


%------------------------------------------------------------------------
%  Sharad J. Shanbhag
%	sshanbhag@neomed.edu
%
%	readPLXFileC:
% 	Author: Benjamin Kraus (bkraus@bu.edu, ben@benkraus.com)
% 	Copyright (c) 2012-2013
% 	Last Modified: 2016-10-07 22:29:18
% 	Revision: 4886
%------------------------------------------------------------------------
% Created: 22 April 2020 (SJS)
%
% Revisions:
%------------------------------------------------------------------------
% TO DO: 
%------------------------------------------------------------------------

	properties
		pathname = '';
		filename = '';
		P
	end
	properties (Dependent)
		hasContinuousData
		plxfile
	end
	
	methods

		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		% Constructor
		%------------------------------------------------------------------------
		function obj = PLXData(varargin)
		% PLXDataObject = PLXData()
		% 	with no input args, creates empty PLXData object
		% PLXDataObject = PLXData(plx_file_name, <options>)
		% 	given a plx filename, data from the plx file will be read using 
		% 	readPLXFileC and stored	in the P as a struct. options to be passed
		% 	to readPLXFileC can be given after the filename
		% PLXDataObject = PLXData(struct_from_readPLXFileC)
		% 	assigns the input data to P
			if isempty(varargin)
				return;
			end
			
			% try different options depending on input
			if ischar(varargin{1})
				% if character string, assume a filename was provided
			
				% get parts of the provided file
				[tmpPath, tmpFile, tmpExt] = fileparts(varargin{1});
				% check if it is a .plx file
				if ~strcmpi(tmpExt, '.plx')
					% ...if not, warn user
					warning('PLXData: file %s does not have .plx extension!', ...
																					varargin{1});
				end
				% assign path and filename.
				obj.pathname = tmpPath;
				obj.filename = [tmpFile tmpExt];
				
				% use default if no arguments provided
				if nargin == 1
					% load data from filename with defaults ('all',
					% 'nocontinuous')
					obj.P = obj.read_plx_file(obj.plxfile);
				else
					% use provided options
					obj.P = obj.read_plx_file(obj.plxfile, varargin{2:end});
				end
				
			elseif isstruct(varargin{1})
				% if struct, assume it is output from readPLXFileC
				obj.P = varargin{1};
				% check if PLXFile field is defined
				if isfield(obj.P, 'PLXFile')
					if isempty(obj.P.PLXFile)
						obj.pathname = '';
						obj.filename = '';
					else
						[tmpPath, tmpFile, tmpExt] = fileparts(obj.P.PLXFile);
						% assign path and filename.
						obj.pathname = tmpPath;
						obj.filename = [tmpFile tmpExt];
					end
				else
					obj.pathname = '';
					obj.filename = '';
				end
			else
				error('PLXData: unknown input to constructor');
			end
			
		end
		%------------------------------------------------------------------------
		
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		% load data from file, without continuous data, store in P property
		%------------------------------------------------------------------------
		function obj = loadPLX_nocontinuous(obj)
			plxfile = obj.plxfile;
			if isempty(plxfile)
				error('PLXData: filename not defined');
			end
			% read in data using readPLXFileC
			obj.P = readPLXFileC(plxfile, 'all', 'nocontinuous');
			obj.P.PLXFile = plxfile;
		end
		%------------------------------------------------------------------------

		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		% load data from file, with continuous data, store in P property
		%------------------------------------------------------------------------
		function obj = loadPLX_continuous(obj)
			plxfile = obj.plxfile;
			if isempty(plxfile)
				error('PLXData: filename not defined');
			end
			% read in data using readPLXFileC
			obj.P = readPLXFileC(plxfile, 'all', 'continuous');
			obj.P.PLXFile = plxfile;
		end
		%------------------------------------------------------------------------

		
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		% get A/D channel number from the Channel name
		%------------------------------------------------------------------------
		function val = getADChannel(obj, varargin)
			% get all channel ids if channel is not provided
			if isempty(varargin)
				val = zeros(obj.P.NumSpikeChannels, 1);
				for c = 1:obj.P.NumSpikeChannels
					val(c) = obj.getADChannel(c);
				end
				return
			else
				cNumber = varargin{1};
			end
			% if requested number is out of bounds, throw error
			if (cNumber < 1) || (cNumber > obj.P.NumSpikeChannels)
				error('PLXData.getADChannel: requested channel # out of bounds');
			end
			% need to extract number from the P.SpikeChannels.name string
			val = sscanf(obj.P.SpikeChannels(cNumber).Name, 'spikechan_%d');
			% warning if val is empty
			if isempty(val)
				warning( ['PLXData.getADChannel: channel not found' ...
								'or misformed channel name in plx file']);
			end
		end
		%------------------------------------------------------------------------

		
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		% export data in OFS mat format
		%------------------------------------------------------------------------
		function [val, varargout] = export_as_mat(obj, varargin)
			% create cell array
			Carr = cell(obj.P.NumSpikeChannels, 1);
			% loop through channels
			for c = 1:obj.P.NumSpikeChannels
				Carr{c} = export_channel_as_mat(obj, c);
			end
			
			if ~isempty(varargin)
				if strcmpi(varargin{1}, 'sort_by_timestamp')
					val = sortrows(cell2mat(Carr), 3);
				end
			else
				val = cell2mat(Carr);
			end
			if nargout > 1
				varargout{1} = Carr;
			end
		end
		
		function val = export_channel_as_mat(obj, channel)
		% Match format expected by SpikeData object:
		%		Column 1: channel is AD Channel from
		%									P.SpikeChannels(cNumber).Name)
		%		Column 2: unit #
		%		Column 3: timestamp (in seconds)
		%		Column 4: PCA1 weight
		%		Column 5: PCA2 weight
		%		Column 6: PCA3 weight
		%		Column 7-end : waveform
		
			% convert Timestamps to seconds (need to do conversion 
			% to double from uint32
			ts_sec = double(obj.P.SpikeChannels(channel).Timestamps) / ...
								double(obj.P.ADFrequency);
			% unit #
			unit = double(obj.P.SpikeChannels(channel).Units);
			% AD channel
			ADchannel = double(obj.getADChannel(channel))*ones(size(ts_sec));
			% dummy data for PCA 
			PCA = zeros(length(ts_sec), 3);
			% need to convert waves from int16 to double and scale
			snips = double(obj.P.SpikeChannels(channel).Waves') / ...
												double(obj.P.SpikeMaxMagnitudeMV);
			% glue everything together
			val = [ADchannel unit ts_sec PCA snips];
		end
		%------------------------------------------------------------------------
		
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		% returns PLXstruct.ContinuousChannels struct
		%------------------------------------------------------------------------
		function val = getContinuousData(obj)
		%------------------------------------------------------------------------
		% val = getContinuousData
		%------------------------------------------------------------------------
		% returns ContinuousChannls struct array (if loaded from plx file),
		% otherwise returns empty array
		%	ContinuousChannels struct fields:
		% 		Name: 'spikechan_4'
		% 		Channel: 0
		% 		SpikeChannel: 1
		% 		SourceID: 101
		% 		ChannelID: 0
		% 		Comment: ''
		% 		Enabled: 1
		% 		ADFrequency: 48828
		% 		ADGain: 580
		% 		PreAmpGain: 1
		%------------------------------------------------------------------------
			if obj.hasContinuousData
				val = obj.P.ContinuousChannels;
			else
				val = [];
			end
		end

		%------------------------------------------------------------------------
		function val = getScaledContinousChannel(obj, varargin)
		%------------------------------------------------------------------------
		% val = getScaledContinousChannel(obj, <channels to get>)
		%------------------------------------------------------------------------
		% returns cell array of scaled continuous data
		% if channels are provided, only data from those channels (indexed as
		% by PLX data - not by A/D channel!) from specified channels will be
		% retrieved. 
		%
		% Data are scaled by PLXData.P.ContMaxMAgnitudeMV
		%------------------------------------------------------------------------
			if obj.hasContinuousData
				if isempty(varargin)
					channelList = 1:obj.P.NumContChannels;
				else
					channelList = varargin{1};
				end
				nchannels = length(channelList);
				val = cell(1, nchannels);
				for c = 1:nchannels
					val{c} = double(obj.P.ContinuousChannels(channelList(c)).Values) / ...
											obj.P.ContMaxMagnitudeMV;
				end
			else
				val = [];
			end
		end
		
		
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		% get/set access for dependent properties
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------

		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		% returns test.Type
		%------------------------------------------------------------------------
		function val = get.plxfile(obj)
			% check to make sure file is specified
			if isempty(obj.filename)
				warning('PLXData: empty filename');
				val = '';
				return
			end
			% create full path and file
			val = fullfile(obj.pathname, obj.filename);
			if ~exist(val, 'file')
				warning('PLXData: file %s not found', val);
				return
			end
		end			
		%------------------------------------------------------------------------

		function val = get.hasContinuousData(obj)
			if any(obj.P.ContSampleCounts > 0)
				val = true;
			else
				val = false;
			end
		end
		
	end % END OF METHODS (general)
	
	methods (Static)
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		% internal method to act as interface to readPLXFileC 
		% function - this will add filename to the struct returned
		% by readPLXFileC
		%
		% NEED TO INCLUDE OTHER OPTIONS (i.e. continuous)
		%------------------------------------------------------------------------
		function pStruct = read_plx_file(varargin)
			if isempty(varargin)
				% get file from user if file wasn't provided
				[fname, fpath] = uigetfile('*.plx', 'Select .plx file');
				if fname == 0
					fprintf('User Cancelled\n');
					pStruct = [];
					return
				else
					plxfile = fullfile(fpath, fname);
				end
				
			elseif ischar(varargin{1})
				% If user provided a filename, check if it exists
				if ~exist(varargin{1}, 'file')
					error('read_plx_file: File not found %s', varargin{1})
				else
					% if so, use it
					plxfile = varargin{1};
				end
			end
			
			% see if options were provided
			if length(varargin) > 1
				plxreadArgs = varargin(2:end);							
			else
				% otherwise, use defaults
				plxreadArgs = {'all', 'nocontinuous'};
			end

			% read in data
			fprintf('read_plx_file: Reading from %s\n', plxfile);
			fprintf('\treadPLXFileC options: %s\n', plxreadArgs{:});
			pStruct = readPLXFileC(plxfile, plxreadArgs{:});
			fprintf('read_plx_file: Done\n');
			pStruct.PLXFile = plxfile;
		end
		%------------------------------------------------------------------------

	end
end