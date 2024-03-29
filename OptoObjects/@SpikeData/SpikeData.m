classdef SpikeData
%------------------------------------------------------------------------
% SpikeData Class
%------------------------------------------------------------------------
% TytoLogy:Experiments:opto...
%------------------------------------------------------------------------
% Class Properties
%  Info    SpikeInfo object
%  Spikes  sorted spikes in table object
% 	        Table Variable Names:
% 					Channel		AD channel
% 					Unit			Unit ID (for given channel! note that units might
% 									 not have unique IDs across channels)
% 					TS				timestamp (seconds)
% 					PCA			PCA values (not valid for data imported 
% 									 directly from plx file
% 					Wave			Wave snippet data
%	Continuous	Continuous Data
%	plxvar		variable names from Plexon 
% 					(used when imported exported MAT from OfflineSorter)
%	hasContinuousData
% 					1 if continuous data are loaded, 
%
%------------------------------------------------------------------------
% Methods for class SpikeData:
% 
%    SpikeData
%        constructor method
%    addContinuousDataFromPLXObj
%        add continuous data from PLXdata object
%    addPlexonInfo
%    addPlexonSpikesFromMat
%    addPlexonSpikesFromPLXObj
%        add sorted spike data from PLXdata object
%    check_channels
%    check_units
%    getSpikesByStim
%    get_default_channel
%    get_default_channel_and_unit
%    get_default_unit
%    indexForChannel
%    indexForFile
%    indexForTestName
%    listChannels
%    listFiles
%    listInfo
%    listTestNames
%    listTestTypes
%    listUnits
%        display list of units and channels
%    nUnits
%    plotAllData
%    plotUnitWaveforms
%    printInfo
%    selectSpikes
%    spikesAsMatrix
%    spikesForAnalysis
%    spikesForAnalysisByUnit
%    spikesForChannel
%    spikesForFile
%    validTestNames
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%  Sharad J. Shanbhag
%	sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: ? 2020 (SJS)
%
% Revisions:
%	30 Apr 2020 (SJS): changing addPlexonSpikes to addPlexonSpikesFromMat
%	various dates (SJS): altered spikesforanalysis method
%	22 May, 2020 (SJS): added Continuous to hold continuous data
%	11 Jun 2020 (SJS): scaling continuous data
%	7 Jul 2020 (SJS): updated documentation
%  21 Jul 2021 (SJS):
%     - added code to scrub unsorted (putative) spike channels from
%       data exported via getSpikesBySt
%     - updated documentation
%------------------------------------------------------------------------
% TO DO: 
%------------------------------------------------------------------------

	%-------------------------------------------------------
	%-------------------------------------------------------
	%-------------------------------------------------------
	%-------------------------------------------------------
	% Properties
	%-------------------------------------------------------
	%-------------------------------------------------------
	properties
	% Info			SpikeInfo object
	% Spikes			sorted spikes in table object
	% Continuous	Continuous Data
	% plxvar			variable names from Plexon
		Info
		Spikes
		Continuous
		plxvar
	end
	
	%-------------------------------------------------------
	% dependent props
	%-------------------------------------------------------
	properties (Dependent)
		hasContinuousData
	end
	
	%-------------------------------------------------------
	% protected contants
	%-------------------------------------------------------
	properties (Access = protected, Constant = true)
		VALID_TESTNAMES =  {'BBN', 'FREQ_TUNING', 'WAV', 'FRA', ...
                           'OPTO-AMP', 'CLICK'};
	end
	
	%-------------------------------------------------------
	%-------------------------------------------------------
	%-------------------------------------------------------
	%-------------------------------------------------------
	% METHODS
	%-------------------------------------------------------
	%-------------------------------------------------------	
	methods
		
		%-------------------------------------------------------
		%-------------------------------------------------------
		% Constructor
		%-------------------------------------------------------
		% [SpikeDataObj] = SpikeData(SpikeInfo, SpikesTable, ContinuousData,
		%																			plxvars)
		%-------------------------------------------------------
		function obj = SpikeData(varargin)
			if length(varargin) ~= 2
				return;
			end
			obj.Info = varargin{1};
			obj.Spikes = varargin{2};
			if length(varargin) > 2
				obj.Continuous = varargin{3};
			end
			if length(varargin) == 4
				obj.plxvar = varargin{4};
			end
		end
		%-------------------------------------------------------
		
		%-------------------------------------------------------
		%-------------------------------------------------------
		function obj = addPlexonSpikesFromMat(obj, plxSpikes, varargin)
		%-------------------------------------------------------
		% addPlexonSpikes(plxSpikes, (optional) plxvar)
		%	plxSpikes is nspikes X ??? matrix exported from Plexon Offline
		%	Sorter (OFS)
		%
		%	optional plxvar is a string(s) from the matfile indicating the adc
		%	channel(s)
		%
		%-------------------------------------------------------
		
		% Define some handy values for indexing into plxSpikes 
		% and similar arrays
		%		Column 1: channel (? need to figure out how 
		%									this maps to original data)
		%		Column 2: unit #
		%		Column 3: timestamp (in seconds)
		%		Column 4: PCA1 weight
		%		Column 5: PCA2 weight
		%		Column 6: PCA3 weight
		%		Column 7-end : waveform
			[~, nc] = size(plxSpikes);
			CHAN_COL = 1;
			UNIT_COL = 2;
			TS_COL = 3;
			PCA_COL = 4:6;
			WAV_COL = 7:nc;
			% convert matrix to table
			% define variable names
			vNames = {'Channel', 'Unit', 'TS', 'PCA', 'Wave'};
			vUnits = {'', '', 'seconds', '', 'mV'};
			obj.Spikes = table(	plxSpikes(:, CHAN_COL), ...
										plxSpikes(:, UNIT_COL), ...
										plxSpikes(:, TS_COL), ...
										plxSpikes(:, PCA_COL), ...
										plxSpikes(:, WAV_COL), ...
										'VariableNames', vNames );
			obj.Spikes.Properties.VariableUnits = vUnits;
			if ~isempty(varargin)
				obj.plxvar = varargin{1};
			end
		end
		%-------------------------------------------------------
		
		%-------------------------------------------------------
		%-------------------------------------------------------
		function obj = addPlexonSpikesFromPLXObj(obj, Pobj)
		%-------------------------------------------------------
		% SpikeData.addPlexonSpikesFromPLXObj
		%	given a plexon data obect (PLXData), extract the sorted spike data
		%	as a matrix (similar to that generated by the Plexon offline
		%	sorter "Export to MAT") and then convert to table as required by
		%	the SpikeData class
		%-------------------------------------------------------
			% get spikes from object as a matrix
			plxSpikes = Pobj.export_as_mat('sort_by_timestamp');
			% number of columns in plxSpikes
			[~, nc] = size(plxSpikes);
			% some shorthand for the different columns
			CHAN_COL = 1;
			UNIT_COL = 2;
			TS_COL = 3;
			PCA_COL = 4:6;
			WAV_COL = 7:nc;
			% convert matrix to table, store as Spikes property
			% define variable names
			vNames = {'Channel', 'Unit', 'TS', 'PCA', 'Wave'};
			vUnits = {'', '', 'seconds', '', 'mV'};
			obj.Spikes = table(	plxSpikes(:, CHAN_COL), ...
										plxSpikes(:, UNIT_COL), ...
										plxSpikes(:, TS_COL), ...
										plxSpikes(:, PCA_COL), ...
										plxSpikes(:, WAV_COL), ...
										'VariableNames', vNames );			
			obj.Spikes.Properties.VariableUnits = vUnits;			
		end
		%-------------------------------------------------------
		
		%-------------------------------------------------------
		%-------------------------------------------------------
		function obj = addContinuousDataFromPLXObj(obj, Pobj)
		%-------------------------------------------------------
		% adds continuousChannels data from Pobj to SpikeData.Continuous
		%  - converts data to double and scales by ContMaxMagnitudeMV
		%-------------------------------------------------------
			obj.Continuous = Pobj.getContinuousData;
			% scale Continuous data
			for c = 1:length(obj.Continuous)
				obj.Continuous(c).Values = double(obj.Continuous(c).Values) / ...
											Pobj.P.ContMaxMagnitudeMV;
			end
		end
		%-------------------------------------------------------
		
				
		%-------------------------------------------------------
		%-------------------------------------------------------
		function obj = addPlexonInfo(obj, SpikeInfo)
		%-------------------------------------------------------
		% assign plexon information obj (SpikeInfo obj) 
		%-------------------------------------------------------
			if ~isstruct(SpikeInfo)
				error('need struct as input');
			end
			% initialize info object
			obj.Info = SpikeInfo;
		end
		%-------------------------------------------------------
		
		%-------------------------------------------------------
		function val = listFiles(obj)
		%-------------------------------------------------------
		% return a list of files in the source data
		% list will be a cell array of strings
		%-------------------------------------------------------
			val = cell(obj.Info.nFiles, 1);
			for f = 1:obj.Info.nFiles
				val{f} = obj.Info.FileInfo{f}.F.file;
			end
		end
		%-------------------------------------------------------
		
		%-------------------------------------------------------
		function val = listTestNames(obj)
		%-------------------------------------------------------
		% return a list of test names in the source data
		% list will be a cell array of strings
		%-------------------------------------------------------
			val = cell(obj.Info.nFiles, 1);
			for f = 1:obj.Info.nFiles
				val{f} = obj.Info.FileInfo{f}.testname;
			end
		end
		%-------------------------------------------------------
					
		%-------------------------------------------------------
		function val = listTestTypes(obj)
		%-------------------------------------------------------
		% return a list of test types in the source data
		% test type will refer to stimulus type (usually)
		% list will be a cell array of strings
		%-------------------------------------------------------
			val = cell(obj.Info.nFiles, 1);
			for f = 1:obj.Info.nFiles
				val{f} = obj.Info.FileInfo{f}.testtype;
			end
		end
		%-------------------------------------------------------
		
		%-------------------------------------------------------
		%-------------------------------------------------------
		function val = listChannels(obj)
		%-------------------------------------------------------
		% return a list of unique channels in the Spikes Table
		% �NOTE: channels will match the AD channel from  TDT!
		%
		% two ways to do this:
		%	val = unique(obj.Spikes{:, 'Channel'});
		% or
		%	val = unique(obj.Spikes.Channel);
		%-------------------------------------------------------
			val = unique(obj.Spikes.Channel);
		end
		%-------------------------------------------------------
		
		%-------------------------------------------------------
		%-------------------------------------------------------
		function [val, varargout] = listUnits(obj, varargin)
		%-------------------------------------------------------
		% return vector of unique unit numbers for each channel
		% stored in a cell array; data are tacken from the Spikes Table
		% �NOTE: channels will match the AD channel from  TDT!
		%
		%	val = obj.listUnits
		%		if no channel provided, all channels used, will return
		%		a {nchannels, 1} cell array, where each element is a 
		%		[# units on channel, 1] vector of unit numbers
		%	val = obj.listUnits(<channel #>)
		%		will return a cell containing array with list of 
		%		units for that channel
		%	val = obj.listUnits([channels])
		%		will return a cell vector with list of units on each of the
		%		requested channels
		%-------------------------------------------------------
		
			% check inputs (channels)
			[clist, nchan] = obj.check_channels(varargin);
			
			% allocate val cell array
			val = cell(nchan, 1);
			% find units on each channel. do this by 
			%	(1) getting indices of current channel in loop:
			%		obj.Spikes.Channel == clist(c)
			%	(2) getting the units for this channel:
			%		obj.Spikes.Unit(obj.Spikes.Channel == clist(c))
			%	(3) finding unique values
			%		unique(obj.Spikes.Unit(obj.Spikes.Channel == clist(c)));
			for c = 1:nchan
				val{c} = unique(obj.Spikes.Unit(obj.Spikes.Channel == clist(c)));
			end
			if nargout > 1
				varargout{1} = clist;
			end
		end
		%-------------------------------------------------------
		
		%-------------------------------------------------------
		%-------------------------------------------------------
		function [val, varargout] = nUnits(obj, varargin)
		%-------------------------------------------------------
		% [# units] = SpikeData.nUnits([channel numbers])
		%-------------------------------------------------------
		% return # of unique units for each channel
		%	if no channel(s) provided, list for all channels
		% �NOTE: channels will match the AD channel from  TDT!
		%-------------------------------------------------------
		
			% check inputs (channels)
			[clist, nchan] = obj.check_channels(varargin);

			val = zeros(nchan, 1);
			for c = 1:nchan
				tmplist = obj.listUnits(clist(c));
				val(c) = length(tmplist{1});
			end
			if nargout > 1
				varargout{1} = clist;
			end
		end
		%-------------------------------------------------------

		%-------------------------------------------------------
		%-------------------------------------------------------
		function varargout = listInfo(obj)
		%-------------------------------------------------------
		% returns cell array lists of Files, Channels, Units
		%-------------------------------------------------------
			varargout{1} = obj.listFiles;
			varargout{2} = obj.listChannels;
			varargout{3} = obj.listUnits;
		end
		%-------------------------------------------------------

		%-------------------------------------------------------
		%-------------------------------------------------------
		function varargout = printInfo(obj)
		%-------------------------------------------------------
		% [fileList, channelList, unitList] = SpikeData.printInfo
		%	Displays information about data
		%-------------------------------------------------------
			sendmsg(sprintf('Data in file %s', obj.Info.FileName));
			% get info
			[fList, cList, uList] = obj.listInfo;
			% continuous data?
			if obj.hasContinuousData
				sendmsg('Has continuous data: Yes');
				for s = 1:length(obj.Continuous)
					fprintf('Continuous Channel %d:\n', s)
					fprintf('\tName: %s\n', obj.Continuous(s).Name);
					fprintf('\t# samples: %d\n', ...
									length(obj.Continuous(s).Values));
				end
					
			else
				sendmsg('Has continuous data: No');
			end
			% display list of files
			sendmsg('Input Data Files:\n');
			fprintf('\tIndex\t\tFilename\n');
			for f = 1:obj.Info.nFiles
				fprintf('\t%d:\t\t%s\n', f, fList{f});
			end
			% display list of channels...
			% ...and units
			% n.b.: could also get both using listUnits: [uList, cList] = obj.listUnits
			sendmsg('Channels and Unit ID #s:\n');
			fprintf('\tIndex\tChannel\tUnits\n');
			for c = 1:length(cList)
				fprintf('\t%d:\t%d\t', c, cList(c));
				fprintf('%d ', uList{c});
				fprintf('\n');
			end
			% assign outputs
			if nargout
				varargout{1} = fList;
				varargout{2} = cList;
				varargout{3} = uList;				
			end
		end
		%-------------------------------------------------------
		
		
		%-------------------------------------------------------
		%-------------------------------------------------------
		function tbl = spikesForChannel(obj, chanNum, varargin)
		% get table of spikes for a specific unit and channel
		% �NOTE: channels will match the AD channel from  TDT!

		% check inputs
			if length(chanNum) ~= 1
				error(['SpikeData.spikesForChannel: requires single,' ...
							'valid channel number']);
			elseif obj.check_channels(chanNum) == -1
				error('SpikeData.spikesForChannel: invalid channel %d', ...
								chanNum);
			end
			% spikes for the desired channel
			tbl = obj.Spikes(obj.Spikes.Channel == chanNum, :);
			
			% if no unit number provided, we're done
			if isempty(varargin)
				return
			else
				unitNum = varargin{1};
			end
			
			% otherwise, get spikes for given unit number
			tbl = tbl( (tbl.Unit == unitNum), :);
		end
		%-------------------------------------------------------
		
		%-------------------------------------------------------
		%-------------------------------------------------------
		function [tbl, varargout] = spikesForFile(obj, fileNum)
		% get table of spikes for a specific file
			if ~between(fileNum, 1, obj.Info.nFiles)
				error('requested file %d out of range [1 %d]', ...
										fileNum, obj.Info.nFiles);
			else
				% look for all timestamps between start of file and end of file
				valid_rows = (obj.Spikes.TS >= obj.Info.fileStartTime(fileNum)) & ...
									(obj.Spikes.TS <= obj.Info.fileEndTime(fileNum));
				tbl = obj.Spikes(valid_rows, :);
				varargout{1} = valid_rows;
			end
		end
		%-------------------------------------------------------


		%-------------------------------------------------------
		%-------------------------------------------------------
		function arr = spikesAsMatrix(obj)
		% get spike information as original matrix (as exported from Plexon)
		
			% return Spikes table in original form
			arr = table2array(obj.Spikes);
		end
		%-------------------------------------------------------		

		%{
		%-------------------------------------------------------
		function events = stimEventTimesForFile(obj, fileNum, varargin)
		%-------------------------------------------------------
		% get array of event structs for given file
		%	events		struct array with fields:
		%		name
		%		samples
		%		timestamps
		%-------------------------------------------------------
			alignMode = 'FILE';
			
			if ~isempty(varargin)
				if any(strcmpi(varargin{1}, {'FILE', 'ORIG'}))
					alignMode = varargin{1};
				else
					error('SpikeData.stimEventTimesForFile: invalid mode %s', ...
									varargin{1});
				end
			end
			
			% get events for the given file/test
			events = obj.Info.FileInfo{fileNum}.geteventList;
			
			switch upper(alignMode)
				case 'FILE'
					% get onset of file
					offsetbin = obj.Info.fileStartBin(fileNum);
					
				case 'ORIG'
					offsetbin = 1;
			end
			% convert samples to timestamps, aligned
			% 1) loop through nevents
			for n = 1:length(events)
				% 2) event times  are ebins + (filestartbin -1) / sampling rate
				events(n).timestamps = ((offsetbin - 1) + events(n).samples) ...
							./ obj.Info.Fs;
			end
		end
		%-------------------------------------------------------
		%}
		
		% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% % % THESE SHOULD PROBABLY BE PROTECTED METHODS
		%-------------------------------------------------------
		%-------------------------------------------------------		
		function [channelNum, unitNum] = get_default_channel_and_unit(obj)
		% [channelNum, unitNum] = SpikeData.get_default_channel_and_unit
		%	looks for first non-zero unit within lowest channel number
		%	returns empty value ([]) if channel and/or unit is not found.
	
		% get list of channels and units
			channelList = obj.listChannels;
			unitList = obj.listUnits;
			% check lists
			if isempty(channelList)
				error('No channels in SpikeData object');
			end
			if isempty(unitList)
				error('No units in SpikeDataObject');
			end
			% get default channel
			channelNum = obj.get_default_channel;			
			% get default unit
			unitNum = obj.get_default_unit(channelNum);
			% raise alerts if needed
			if isempty(channelNum)
				warning('No channel found')
			elseif isempty(unitNum)
				warning('No unit found');
			end
		end
		
		%-------------------------------------------------------
		%-------------------------------------------------------		
		function channelNum = get_default_channel(obj)
		% [channelNum] = SpikeData.get_default_channel
		%	looks for first non-zero unit within lowest channel number
		%	returns empty value ([]) if channel and/or unit is not found.

			% setup: set channel and unit to empty
			channelNum = [];
			% get list of channels and units
			channelList = obj.listChannels;
			unitList = obj.listUnits;
			% check lists
			if isempty(channelList)
				error('No channels in SpikeData object');
			end
			if isempty(unitList)
				error('No units in SpikeDataObject');
			end
			% loop through channels to find first non-zero unit in lowest
			% channel number
			cI = 1;
			while isempty(channelNum) && cI <= length(channelList)
				% check if unitList is empty for current channel
				if ~isempty(unitList{cI})
					% nonzero units?
					nz = unitList{cI}(unitList{cI} ~= 0);
					if ~isempty(nz)
						% if there is a non-zero unit in this channel, we're done
						channelNum = channelList(cI);
					else
						% no non-zero unit, so go to next channel
						cI = cI + 1;
					end
				else
					% no units thos channel, go to next one
					cI = cI + 1;
				end
			end
		end
		%-------------------------------------------------------
		
		%-------------------------------------------------------
		%-------------------------------------------------------		
		function unitNum = get_default_unit(obj, channelNum)
		% [channelNum] = SpikeData.get_default_unit(channelNum)
		%	looks for first non-zero unit within given channelNum
		%	returns empty value ([]) if channel and/or unit is not found.
		
			% setup: set unit to empty
			unitNum = [];
			% need index of channel in channelList
			cI = find(channelNum == obj.listChannels);
			if isempty(cI)
				error('channel %d not found', channelNum);
			end
			% nonzero units?
			unitList = obj.listUnits;
			nz = unitList{cI}(unitList{cI} ~= 0);
			if ~isempty(nz)
				unitNum = nz(1);
			end
		end
		%-------------------------------------------------------
		%-------------------------------------------------------
		
		
		%-------------------------------------------------
		%-------------------------------------------------
		function indx = indexForFile(obj, filename)
		% index_number = SpikeData.indexForFile(filename)
		%-------------------------------------------------
		% given a filename (with or without path), returns index to that
		% file (into SpikeData.Info.FileInfo{}).
		% If not found, empty array returned.
		%-------------------------------------------------
		
			% strip path, save file and extension
			% convert to unix filesep if need be
			[~, fbase, fext] = fileparts(path_unix(filename));
			indx = find(strcmpi([fbase fext], obj.listFiles));
		end
		%-------------------------------------------------
		%-------------------------------------------------
		
		
		%-------------------------------------------------
		%-------------------------------------------------
		function indx = indexForChannel(obj, channel_number)
		% index_number = indexForChannel(obj, channel_number)
		%-------------------------------------------------
		% given a channel_number, returns index_number for accessing
		% that channel in various SpikeData arrays
		%-------------------------------------------------
			if obj.check_channels(channel_number) == -1
				indx = [];
			else
				cList = obj.listChannels;
				indx = find(channel_number == cList);
			end
		end
		%-------------------------------------------------
		%-------------------------------------------------

		%-------------------------------------------------
		%-------------------------------------------------
		function indx = indexForTestName(obj, test_name)
		%-------------------------------------------------
		% index_number = indexForTestName(obj, test_name)
		%-------------------------------------------------
		% given a test name, returns index_number for accessing
		% that file in various SpikeData arrays (e.g.,
		% SpikeData.Info.FileInfo{})
		% if file is not found, empty value will be returned
		%-------------------------------------------------		
			indx = find(strcmpi(test_name, obj.listTestNames));			
		end
		%-------------------------------------------------
		%-------------------------------------------------		
		
		
		%-------------------------------------------------
		%-------------------------------------------------		
		function val = validTestNames(obj)
			val = obj.VALID_TESTNAMES;
		end		
		%-------------------------------------------------

		% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		
		% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		% get/set access for dependent properties
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		% continuous data?
		function val = get.hasContinuousData(obj)
			val = ~isempty(obj.Continuous);
		end

		
		% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		% TEMPORARY (???)
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		function [clist, varargout] = check_channels(obj, channel_arg)
			if isempty(channel_arg)
				% if no channels provided, use all existing channels
				clist = obj.listChannels;
			else
				% otherwise, use provided channels
				% note that channel_arg might be a cell
				if iscell(channel_arg)
					clist = channel_arg{1};
				else
					clist = channel_arg;
				end
				% check channels existence
				cchk = ismember(clist, obj.listChannels);
				% if all channels are not in list of channels, find mismatch
				% and report error. set clist element for that channel to -1 to
				% indicate error
				if ~all(cchk)
					for c = 1:length(clist)
						if cchk(c) == false
							warning('unknown channel: %d\n', clist(c))
							clist(c) = -1;
						end
					end
				end
			end
			% number of channels
			varargout{1} = length(clist);
			% index or indices for channel(s)
			tmp = zeros(length(clist), 1);
			for c = 1:length(clist)
				tmpindx = find(clist(c) == obj.listChannels);
				if isempty(tmpindx)
					tmp(c) = -1;
				else
					tmp(c) = tmpindx;
				end
			end
			varargout{2} = tmp;
		end
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		
 		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		function [val, varargout] = check_units(obj, channels, units)
			% need to use nargin == 3 since obj counts as an input in 
			% object methods!!!!!
			if nargin ~= 3
				error('SpikeData.check_units: requires channel(s) and unit(s)');
			end
			
			% check channel provided
			cchk = obj.check_channels(channels);
			if any(cchk == -1)
				fprintf('SpikeData.check_units: invalid channelNum %d\n', ...
									channels(cchk == -1));
				error('SpikeData.check_units: invalid channel');
			end
			% check if units is empty
			if isempty(units)
				error(['SpikeData.check_units: ' ...
							'list of units to check is empty']);
			end

			% some checks depending on length of channels input
			if length(channels) > 1
				% if more than one channel is specified, units needs to be a
				% cell
				if ~iscell(units)
					% throw error
					error(['SpikeData.check_units: if multiple ' ...
								'channels are specified, cell vector of units' ...
								'to match must be provided to spikesForAnalysis']);
				else
					% no error - assign units cell to uCell;
					uCell = units;
				end
			elseif ~iscell(units)
				% if length channels is 1, put units into a cell to simplify
				% search/check code
				uCell = {units};
			end

			% check for length mismatch in channels and uCell
			if length(channels) ~= length(uCell)
				error(['SpikeData.check_units: mismatch in length of ' ...
								'channels and units']);
			end

			% allocate output
			uList = cell(size(uCell));
			val = true(size(uCell));

			% get list of units and the corresponding list of channels for the
			% units
			[allUnits, allChannels] = obj.listUnits;
			% loop through channels
			for c = 1:length(channels)
				% get list of units for this channel
				units_for_chan = allUnits{channels(c) == allChannels};
				% check if input units for this channel match
				uchk = ismember(uCell{c}, units_for_chan);
				% add list of zeros to output uList for this channel
				uList{c} = zeros(size(uchk));
				if ~all(uchk)
					% output for val for this channel is false
					val(c) = false;
					for u = 1:length(uCell{c})
						% assign -1 to bad units
						if uchk(u) == false
							uList{c}(u) = -1;
							warning(['SpikeData.check_units: unit %d not ' ...
										'found for channel %d'], uCell{c}(u), ...
										channels(c))
						end
					end % END u loop
				end % END if ~all(uchk)
			end % END c loop
			
			% assign output
			if nargout > 1
				varargout{1} = uList;
			end
		end % END check_units()
		
		
		
		%-------------------------------------------------
		%-------------------------------------------------
		% methods defined in separate files
		%-------------------------------------------------
		%-------------------------------------------------
		% get data  and sweep info for each channel
		varargout = getSpikesByStim(obj, fileNum, channel, unit)
		% plot waveforms
		H = plotUnitWaveforms(obj, channel, varargin)
		% get spikes (in matlab table) for file, channel and unit
		[tbl, varargout] = selectSpikes(obj, fileNum, channelNum, unitNum)
		% get spikes for analysis
		spikesBySweep = spikesForAnalysisByUnit(obj, fileNum, ...
																	unitNum, varargin)
		% plots all data for given channel and unit
		varargout = plotAllData(obj, channel, unit, varargin)
	end	% END OF METHODS (general)
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------

	
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	%{
	methods (Access = protected)
	end	% END METHODS (PROTECTED)
	%}
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
		
end	% END OF CLASS DEF
