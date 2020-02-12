classdef SpikeData
	properties
		Info
		Spikes
	end
	
	methods
		function obj = SpikeData(varargin)
			if length(varargin) ~= 2
				return;
			end
			obj.Info = varargin{1};
			obj.Spikes = varargin{2};
		end
		
		function obj = addPlexonSpikes(obj, plxSpikes)
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
			obj.Spikes = table(	plxSpikes(:, CHAN_COL), ...
										plxSpikes(:, UNIT_COL), ...
										plxSpikes(:, TS_COL), ...
										num2cell(plxSpikes(:, PCA_COL), 2), ...
										num2cell(plxSpikes(:, WAV_COL), 2), ...
										'VariableNames', vNames );
	
		end
		
		function obj = addPlexonInfo(obj, plxInfo)

		end
		
		function obj = listChannels(obj)
			
		end
	end
end
