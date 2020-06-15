classdef SpikeData
%------------------------------------------------------------------------
% TytoLogy:Experiments:opto...
%------------------------------------------------------------------------
% Info		SpikeInfo object
% Spikes		sorted spikes in table object
% 				Table Variable Names:
% 					Channel		AD channel
% 					Unit			Unit ID (for given channel! note that units might
% 									 not have unique IDs across channels)
% 					TS				timestamp (seconds)
% 					PCA			PCA values (not valid for data imported 
% 									 directly from plx file
% 					Wave			Wave snippet data
% Continuous	Continuous Data
% plxvar		variable names from Plexon 
% 					(used when imported exported MAT from OfflineSorter)
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
%------------------------------------------------------------------------
% TO DO: how to handle multiple channels?????
%------------------------------------------------------------------------

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
	properties (Dependent)
		hasContinuousData
	end
	
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
		%-------------------------------------------------------
% 			if Pobj.hasContinuousData
% 				obj.Continuous = Pobj.P.ContinuousChannels;
% 			else
% 				warning(['SpikeData.addContinuousDataFromPLXObj: ' ...
% 								'Pobj has no continupis channel data'])
% 				obj.Continuous = [];
% 			end
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
		function tbl = spikesForChannel(obj, chanNum, varargin)
		%-------------------------------------------------------
		% tbl = spikesForChannel(channelNum, <unitNum>)
		% get table of spikes for a specific unit and channel
		% �NOTE: channels will match the AD channel from  TDT!
		%-------------------------------------------------------
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
		% get table of spikes for a specific file
		%-------------------------------------------------------
		function tbl = spikesForFile(obj, fileNum)
			if ~between(fileNum, 1, obj.Info.nFiles)
				error('requested file %d out of range [1 %d]', ...
										fileNum, obj.Info.nFiles);
			else
				% look for all timestamps between start of file and end of file
				valid_rows = (obj.Spikes.TS >= obj.Info.fileStartTime(fileNum)) & ...
									(obj.Spikes.TS <= obj.Info.fileEndTime(fileNum));
				tbl = obj.Spikes(valid_rows, :);
			end
		end
		%-------------------------------------------------------

		%-------------------------------------------------------
		%-------------------------------------------------------
		function spikesBySweep = spikesForAnalysis(obj, fileNum, varargin)
		%-------------------------------------------------------
		% get table of spikes for a specific file aligned to file start,
		% sweep, or as-is (original)
		%
		%	spikesBySweep = obj.spikesForAnalysis(fileNum)
		%
		%	Input args:
		%		fileNum	file index for data 
		%-------------------------------------------------------
			%--------------------------------------
			% set defaults
			%--------------------------------------
			% default is no shift
			ALIGN = 'original';
			% set channel and unit to empty
			channelNum = [];
			unitNum = [];
					
			%--------------------------------------
			% process options and inputs
			%--------------------------------------			
			% process options
			argI = 1;
			while argI <= length(varargin)
				switch upper(varargin{argI})
					case 'ALIGN'
						% check alignment mode
						if ~any(strcmpi(varargin{argI+1}, ...
										{'original', 'file', 'sweep'}))
							% unknown mode
							error(['SpikeData.spikesForAnalysis:' ...
										'unknown align option %s'], varargin{argI+1});
						else
							ALIGN = lower(varargin{argI+1});
						end
						argI = argI + 2;
					case 'CHANNEL'
						% user specified channel option, so get desired list
						channelNum = varargin{argI + 1};
						argI = argI + 2;
					case 'UNIT'
						% user specified unit(s) so get them from input
						unitNum = varargin{argI + 1};
						argI = argI + 2;
					otherwise
						% unknown option provided by user
						error(['SpikeData.spikesForAnalysis:' ...
										'unknown option %s'], varargin{argI});
				end
			end

			%--------------------------------------
			% select valid timestamps/table entries
			%--------------------------------------
			% check that file is in range
			if ~between(fileNum, 1, obj.Info.nFiles)
				error('requested file %d out of range [1 %d]', ...
										fileNum, obj.Info.nFiles);
			end
			
			% create temp table of data for desired file
			tmpS = obj.spikesForFile(fileNum);
			
			% if channelNum is empty, use all channels and units
			if isempty(channelNum)
				fprintf(['SpikeData.spikesForAnalysis:' ...
									'using all channels and units\n']);
				channel_rows = true(size(tmpS.Channel));
				unit_rows = true(size(tmpS.Channel)); 
% 				channelNum = -1;
				unitNum = -1;
			else
				% check channel provided
				cchk = obj.check_channels(channelNum);
				if any(cchk == -1)
					fprintf('SpikeData.spikesForAnalysis: invalid channelNum %d\n', ...
									channelNum(cchk == -1));
					error('SpikeData.spikesForAnalysis: invalid channel');
				end
				% get indices for channel(s)
				channel_rows = false(size(tmpS.Channel));
				for c = 1:length(channelNum)
					channel_rows = channel_rows | (tmpS.Channel == channelNum(c));
				end
			end
			% if unitNum is empty, find all units
			if isempty(unitNum)
				fprintf(['SpikeData.spikesForAnalysis:' ...
									'using all units for channel %d\n'], channelNum);
				% set unit_rows to ones, size of tmpS.Channel
				unit_rows = true(size(tmpS.Channel));
			elseif unitNum ~= -1
				% get indices for unit
				unit_rows = tmpS.Unit == unitNum;				
			end
			% reduce table to valid channel and unit
			vS = tmpS( (channel_rows & unit_rows), :);
			% clear tmpS
			clear tmpS
			
			%--------------------------------------
			% process sweeps
			%--------------------------------------
			% number of sweeps for this file
			nsweeps = length(obj.Info.sweepStartTime{fileNum});
			% some checks
			if nsweeps ~= length(obj.Info.sweepEndTime{fileNum})
				error('mismatch in lengths of sweepStartTime and sweepEndTime');
			elseif nsweeps == 0
				error('no sweeps!');	
			end
			% alignment issues (reference timestamps to start of datafile,
			% start of sweep, or leave as is  - referenced to start of merged
			% file used for spike sorting)
			switch ALIGN
				case 'original'
					% don't modify time stamp values
					fprintf('SpikeData:spikesForAnalysis: no timestamp realignment\n');
					alignval = zeros(nsweeps, 1);
				case 'file'
					% adjust timestamp value by first sweep start time for each
					% file in the overall merged file used for sorting
					fprintf('SpikeData:spikesForAnalysis: align timestamp to file start\n');
					alignval = obj.Info.sweepStartTime{fileNum}(1) ...
										* ones(nsweeps, 1);
				case 'sweep'
					% adjust timestamp value by each sweep start time - this is
					% most useful when doing things like analysis and plots of
					% data 
					fprintf('SpikeData:spikesForAnalysis: align timestamp to sweep start\n');
					alignval = obj.Info.sweepStartTime{fileNum};
			end
			
			% allocate cell to store spike info for each sweep
			spikesBySweep = cell(nsweeps, 1);

			% loop through sweeps
			for s = 1:nsweeps
				% find the valid time stamps (between sweepStartTimes and
				% sweepEndTimes)
				valid_rows = (vS.TS >= obj.Info.sweepStartTime{fileNum}(s)) ...
								& (vS.TS < obj.Info.sweepEndTime{fileNum}(s));
				spikesBySweep{s} = vS(valid_rows, :);
				% apply offset correction
				spikesBySweep{s}{:, 'TS'} = spikesBySweep{s}{:, 'TS'} - alignval(s);
			end
		end
		%-------------------------------------------------------
		

		%-------------------------------------------------------
		%-------------------------------------------------------
		function spikesBySweep = spikesForAnalysisByUnit(obj, fileNum, ...
																			unitNum, varargin)
		%-------------------------------------------------------
		% get table of spikes for a specific file, unit and by sweep
		% this should be deprecated since it doesn't account for channel
		%-------------------------------------------------------
			%--------------------------------------
			% process options and inputs
			%--------------------------------------
			% check that file is in range
			if ~between(fileNum, 1, obj.Info.nFiles)
				error('requested file %d out of range [1 %d]', ...
										fileNum, obj.Info.nFiles);
			elseif ~any(unitNum == obj.listUnits)
				% and that unitnum is valid
				error('unit %d not in Spikes table', unitNum);
			end
			% check alignment mode
			if isempty(varargin)
				% default is no shift
				ALIGN = 'original';
			else
				if ~any(strcmpi(varargin{1}, {'original', 'file', 'sweep'}))
					% unknown mode
					error('%s: unknown align option %s', varargin{1})
				else
				ALIGN = lower(varargin{1});
				end
			end
			
			%--------------------------------------
			% process data
			%--------------------------------------
			% get valid rows for unit
			unit_rows = (obj.Spikes.Unit == unitNum);
			% get valid rows for file - this is done by finding the spike
			% times that occurred between the specified file's start and end
			% times
			file_rows = (obj.Spikes.TS >= obj.Info.fileStartTime(fileNum)) & ...
									(obj.Spikes.TS <= obj.Info.fileEndTime(fileNum));
			% get reduced Spikes table for specified unit and data file
			vS = obj.Spikes( (unit_rows & file_rows), :);
			
			% number of sweeps for this file
			nsweeps = length(obj.Info.sweepStartTime{fileNum});
			% some checks
			if nsweeps ~= length(obj.Info.sweepEndTime{fileNum})
				error('mismatch in lengths of sweepStartTime and sweepEndTime');
			elseif nsweeps == 0
				error('no sweeps!');	
			end
			% alignment issues
			switch ALIGN
				case 'original'
					% don't modify time stamp values
					alignval = zeros(nsweeps, 1);
				case 'file'
					% adjust timestamp value by first sweep start time
					alignval = obj.Info.sweepStartTime{fileNum}(1) ...
										* ones(nsweeps, 1);
				case 'sweep'
					% adjust timestamp value by each sweep start time
					alignval = obj.Info.sweepStartTime{fileNum};
			end
			
			% allocate cell to store spike info for each sweep
			spikesBySweep = cell(nsweeps, 1);
			
			% loop through sweeps
			for s = 1:nsweeps
				% find the valid time stamps (between sweepStartTimes and
				% sweepEndTimes), store the row indices...
				valid_rows = (vS.TS >= obj.Info.sweepStartTime{fileNum}(s)) ...
								& (vS.TS < obj.Info.sweepEndTime{fileNum}(s));
				% get Spikes that match
				spikesBySweep{s} = vS(valid_rows, :);
				% apply offset correction to timestamps (TS variable)
				spikesBySweep{s}{:, 'TS'} = spikesBySweep{s}{:, 'TS'} - alignval(s);
			end
		end
		%-------------------------------------------------------
		
		%-------------------------------------------------------
		function arr = spikesAsMatrix(obj)
		%-------------------------------------------------------
		% get spike information as original matrix (as exported from Plexon)
		%-------------------------------------------------------
			% return Spikes table in original form
			arr = table2array(obj.Spikes);
		end
		%-------------------------------------------------------
		
		%-------------------------------------------------------
		%-------------------------------------------------------
		function H = plotUnitWaveforms(obj, channel, varargin)
		%-------------------------------------------------------
		% [plot handles] = obj.plotUnitWaveforms([channels], 
		% Plot sorted waveforms for each identified unit for a given channel
		% This will work for individual channels and either all units for the
		% channel (if unit list is not provided) or a specified unit(s)
		%-------------------------------------------------------
			% check channel provided
			cchk = obj.check_channels(channel);
			if cchk == -1
				error('SpikeData.plotUnitWaveforms: invalid channel %d', channel);
			end
			if isempty(varargin)
				tmp = obj.listUnits(channel);
				unitList = tmp{1};
				nU = length(unitList);
			else
				unitList = varargin{1};
				% check units
				nU = length(unitList);
				if nU == 0
					error('SpikeData.plotUnitWaveforms: no units for channel %d', ...
												channel);
				else
					for u = 1:nU
						if ~any(unitList(u) == obj.listUnits)
							error('unit %d not found', unitList(u));
						end
					end
				end
			end			
			% allocate gobjects array to hold figure handles
			H = gobjects(nU, 1);
			% loop through units
			for u = 1:nU
				fprintf('Plotting unit %d waveforms\n', unitList(u));
				% create figure and store handle in H array
				H(u) = figure;
				% get spike waveforms for this unit
				W = obj.Spikes{obj.Spikes.Unit == unitList(u), 'Wave'};
				if ~isempty(W)
					[~, nBins] = size(W);
					ms = (1000/obj.Info.Fs) * (0:(nBins - 1));
					plot(ms, W', 'k');
				end
				[~, fstr] = fileparts(obj.Info.FileName);
				tstr = sprintf('File: %s   Unit: %d', [fstr '.mat'], unitList(u));
				title(tstr, 'Interpreter', 'none');
				xlabel('ms');
				grid on
				H(u).Name = sprintf('%s_unit%d', fstr, unitList(u));
			end
		end
		%-------------------------------------------------------

		
		% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%		
		% % % THESE SHOULD PROBABLY BE PROTECTED METHODS
		%-------------------------------------------------------
		%-------------------------------------------------------		
		function [channelNum, unitNum] = get_default_channel_and_unit(obj)
		%-------------------------------------------------------
		% [channelNum, unitNum] = SpikeData.get_default_channel_and_unit
		%	looks for first non-zero unit within lowest channel number
		%	returns empty value ([]) if channel and/or unit is not found.
		%-------------------------------------------------------
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
		%-------------------------------------------------------
		% [channelNum] = SpikeData.get_default_channel
		%	looks for first non-zero unit within lowest channel number
		%	returns empty value ([]) if channel and/or unit is not found.
		%-------------------------------------------------------
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
		%-------------------------------------------------------
		% [channelNum] = SpikeData.get_default_unit(channelNum)
		%	looks for first non-zero unit within given channelNum
		%	returns empty value ([]) if channel and/or unit is not found.
		%-------------------------------------------------------
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
		
	end	% END OF METHODS (general)
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------

	
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	methods (Access = protected)
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
		end
		%------------------------------------------------------------------------
	end	% END METHODS (PROTECTED)
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
		
end	% END OF CLASS DEF
