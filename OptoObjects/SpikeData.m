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
% plxvar		variable names from Plexon 
% 					(used when imported exported MAT)
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%  Sharad J. Shanbhag
%	sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: ? 2020 (SJS)
%
% Revisions:
%	30 Apr 2020 (SJS): changing addPlexonSpikes to addPlexonSpikesFromMat
%------------------------------------------------------------------------
% TO DO: how to handle multiple channels?????
%------------------------------------------------------------------------

	properties
	% Info		SpikeInfo object
	% Spikes		sorted spikes in table object
	% plxvar		variable names from Plexon
		Info
		Spikes
		plxvar
	end
	
	methods
		function obj = SpikeData(varargin)
			if length(varargin) ~= 2
				return;
			end
			obj.Info = varargin{1};
			obj.Spikes = varargin{2};
			if length(varargin) == 3
				plxvars = varargin{3}; %#ok<NASGU>
			end
		end
		
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
		% 
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
			plxSpikes = Pobj.export_as_mat('sort_by_timestamp');
			[~, nc] = size(plxSpikes);
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
		function obj = addPlexonInfo(obj, plxInfo)
		%-------------------------------------------------------
		% assign plexon information obj (SpikeInfo obj) 
		%-------------------------------------------------------
			if ~isstruct(plxInfo)
				error('need struct as input');
			end
			% initialize info object
			obj.Info = plxInfo;
		end
		%-------------------------------------------------------
		
		%-------------------------------------------------------
		%-------------------------------------------------------
		function val = listChannels(obj)
		%-------------------------------------------------------
		% return a list of unique channels in the Spikes Table
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
		function val = listUnits(obj, varargin)
		%-------------------------------------------------------
		% return vector of unique unit numbers for each channel
		% stored in a cell array
		%	if no channel provided, all channels used
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
		end
		%-------------------------------------------------------
		

		%-------------------------------------------------------
		%-------------------------------------------------------
		function val = nUnits(obj, varargin)
		%-------------------------------------------------------
		% [# units] = SpikeData.nUnits([channel numbers])
		% return # of unique units for each channel
		%	if no channel(s) provided, list for all channels
		%-------------------------------------------------------
			% check inputs (channels)
			[clist, nchan] = obj.check_channels(varargin);

			val = zeros(nchan, 1);
			for c = 1:nchan
				tmplist = obj.listUnits(clist(c));
				val(c) = length(tmplist{1});
			end
		end
		%-------------------------------------------------------
		

		%-------------------------------------------------------
		%-------------------------------------------------------
		function tbl = spikesForChannel(obj, chanNum, varargin)
		%-------------------------------------------------------
		% tbl = spikesForChannel(channelNum, <unitNum>)
		% get table of spikes for a specific unit and channel
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
		% get table of spikes for a specific file aligned to file start,
		% sweep, or as-is (original)
		%-------------------------------------------------------
		function spikesBySweep = spikesForAnalysis(obj, fileNum, varargin)
			%--------------------------------------
			% process options and inputs
			%--------------------------------------
			% check that file is in range
			if ~between(fileNum, 1, obj.Info.nFiles)
				error('requested file %d out of range [1 %d]', ...
										fileNum, obj.Info.nFiles);
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
			% get valid rows for file
			file_rows = (obj.Spikes.TS >= obj.Info.fileStartTime(fileNum)) & ...
									(obj.Spikes.TS <= obj.Info.fileEndTime(fileNum));
			% get reduced table
			vS = obj.Spikes( file_rows, :);
			
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
					alignval = zeros(nsweeps, 1);
				case 'file'
					% adjust timestamp value by first sweep start time for each
					% file in the overall merged file used for sorting
					alignval = obj.Info.sweepStartTime{fileNum}(1) ...
										* ones(nsweeps, 1);
				case 'sweep'
					% adjust timestamp value by each sweep start time - this is
					% most useful when doing things like analysis and plots of
					% data 
					alignval = obj.Info.sweepStartTime{fileNum};
			end
			
			% allocate cell to store spike info for each sweep
			spikesBySweep = cell(nsweeps, 1);
			
			% tmp2{:, 'TS'} = tmp{:, 'TS'} * 2

			% loop through sweeps
			for s = 1:nsweeps
				% find the valid time stampes (between sweepStartTimes and
				% sweepEndTimes)
% 				valid_rows = (vS(:, 'TS') >= obj.Info.sweepStartTime{fileNum}(s)) ...
% 								& (vS(:, 'TS') < obj.Info.sweepEndTime{fileNum}(s));
				valid_rows = (vS.TS >= obj.Info.sweepStartTime{fileNum}(s)) ...
								& (vS.TS < obj.Info.sweepEndTime{fileNum}(s));
				spikesBySweep{s} = vS(valid_rows, :);
				% apply offset correction
				spikesBySweep{s}{:, 'TS'} = spikesBySweep{s}{:, 'TS'} - alignval(s);
			end
		end
		

		%-------------------------------------------------------
		% get table of spikes for a specific file, unit and by sweep
		%-------------------------------------------------------
		function spikesBySweep = spikesForAnalysisByUnit(obj, fileNum, ...
																			unitNum, varargin)
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
			
			% tmp2{:, 'TS'} = tmp{:, 'TS'} * 2

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
		% get spike information as original matrix (as exported from Plexon)
		%-------------------------------------------------------
		function arr = spikesAsMatrix(obj)
			% return Spikes table in original form
			arr = table2array(obj.Spikes);
		end
		
		%-------------------------------------------------------
		% Plot sorted waveforms for each identified unit
		%-------------------------------------------------------
		function H = plotUnitWaveforms(obj, unitList)
			% check units
			nU = length(unitList);
			if nU == 0
				unitList = obj.listUnits;
				nU = obj.nUnits;
			else
				for u = 1:nU
					if ~any(unitList(u) == obj.listUnits)
						error('unit %d not found', unitList(u));
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
		
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		% get/set access for dependent properties
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------

		% returns spikes table
		function tbl = get.Spikes(obj)
			tbl = obj.Spikes;
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
