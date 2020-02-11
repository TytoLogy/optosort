classdef SpikeData
	properties
		Info
		Spikes
	end
	
	methods
		function obj = SpikeData(varargin)
			if length(varargin) ~= 2
				error('Need Info and Spikes');
			end
			
			obj.Info = varargin{1};
			obj.Spikes = varargin{2};
		end
		
	end
end
			