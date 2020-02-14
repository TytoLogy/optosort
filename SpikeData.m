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
		
		
		function tbl = get.Spikes(obj)
			tbl = obj.Spikes;
		end
		
		%-------------------------------------------------------
		% get table of spikes for a specific unit
		%-------------------------------------------------------
		function tbl = spikesForUnit(obj, unit_num)
			% check unit_num
			if ~any(unit_num == obj.listUnits)
				warning('unit %d not in Spikes table', unit_num);
				tbl = [];
			else
				unit_rows = (obj.Spikes.Unit == unit_num);
				tbl = obj.Spikes(unit_rows, :);
			end
		end
		
		%-------------------------------------------------------
		% get spikes for a specific file
		%-------------------------------------------------------
		function tbl = spikesForFile(obj, fIndx)
			if ~between(fIndx, 1, obj.Info.nFiles)
				error('requested file %d out of range [1 %d]', ...
										fIndx, obj.Info.nFiles);
			else
				valid_rows = (obj.Spikes.TS >= obj.Info.fileStartTime(fIndx)) & ...
									(obj.Spikes.TS <= obj.Info.fileEndTime(fIndx));
				tbl = obj.Spikes(valid_rows, :);
			end
		end
		
		function arr = spikesAsMatrix(obj)
			% return Spikes table in original form
			arr = table2array(obj.Spikes);
		end
		
		
	end
end
