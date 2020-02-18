classdef SpikeData
	properties
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
		
		function obj = addPlexonSpikes(obj, plxSpikes, varargin)
			% addPlexonSpikes(plxSpikes, (optional) plxvar)
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
		% assign plexon information obj (SpikeInfo obj) 
		%-------------------------------------------------------
		function obj = addPlexonInfo(obj, plxInfo)
			if ~isstruct(plxInfo)
				error('need struct as input');
			end
			% initialize info object
			obj.Info = plxInfo;
		end
		
		%-------------------------------------------------------
		% return a list of unique channels in the Spikes Table
		%-------------------------------------------------------
		function val = listChannels(obj)
			% two ways to do this:
			%	val = unique(obj.Spikes{:, 'Channel'});
			% or
			%	val = unique(obj.Spikes.Channel);
			val = unique(obj.Spikes.Channel);
		end
		
		%-------------------------------------------------------
		% return vector of unique unit numbers
		%-------------------------------------------------------
		function val = listUnits(obj)
			val = unique(obj.Spikes.Unit);
		end
		function val = nUnits(obj)
			val = length(obj.listUnits);
		end
		
		function tbl = get.Spikes(obj)
			tbl = obj.Spikes;
		end
		
		%-------------------------------------------------------
		% get table of spikes for a specific unit
		%-------------------------------------------------------
		function tbl = spikesForUnit(obj, unitNum)
			% check unit_num
			if ~any(unitNum == obj.listUnits)
				warning('unit %d not in Spikes table', unitNum);
				tbl = [];
			else
				unit_rows = (obj.Spikes.Unit == unitNum);
				tbl = obj.Spikes(unit_rows, :);
			end
		end
		
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
		% get table of spikes for a specific file, unit and by sweep
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
			% get valid rows for file
			file_rows = (obj.Spikes.TS >= obj.Info.fileStartTime(fileNum)) & ...
									(obj.Spikes.TS <= obj.Info.fileEndTime(fileNum));
			% get reduced table
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
				H(u).Name = fstr;
			end
			
		end
		
		
		
	end
end
