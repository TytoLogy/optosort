classdef LEVELcache < StimulusCache
	properties

		
	end	% END properties (main)
	properties (Dependent)
	end	% END properties(Dependent)

	methods
		%-------------------------------------------------
		%-------------------------------------------------
		% Constructor
		%-------------------------------------------------
		%-------------------------------------------------
		function obj = LEVELcache(varargin)
			if isempty(varargin)
				return
			end
			if isstruct(varargin{1})
				if length(varargin) == 1
					obj = initObject(obj, varargin{1});
				elseif length(varargin) == 3
					obj = initObject(obj, varargin{1}, varargin{2}, varargin{3});
				else
					error('incorrect number of input arguments');
				end
			else
				error('Unknown input type %s', varargin{1});
			end
		end

		%-------------------------------------------------
		%-------------------------------------------------
		% returns stimulus Indices and list of stim variables
		%-------------------------------------------------
		%-------------------------------------------------
		function [stimindex, levellist] = getStimulusIndices(obj)
			% for LEVEL test, find indices of stimuli with same level (dB SPL)
			fprintf('\t%s test, finding indices\n', obj.Type);
			% list of levels, and # of levels tested
			levellist = obj.LEVEL;
			nlevels = length(obj.vrange);
			% locate where trials for each frequency are located in the
			% stimulus cache list - this will be used to pull out trials of
			% same frequency
			stimindex = cell(nlevels, 1);
			for l = 1:nlevels
				stimindex{l} = find(obj.vrange(l) == levellist);
			end
		end
	end	% END methods
	
end	% END classdef