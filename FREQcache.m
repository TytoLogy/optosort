classdef FREQcache < StimulusCache
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
		function obj = FREQcache(varargin)
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
		function [stimindex, freqlist] = getStimulusIndices(obj)
			% for FREQ test, find indices of stimuli with same frequency
			fprintf('\t%s test, finding indices\n', obj.testtype);
			% list of frequencies, and # of freqs tested
			freqlist = cell2mat(obj.FREQ);
			nfreqs = length(obj.vrange);
			% locate where trials for each frequency are located in the
			% stimulus cache list - this will be used to pull out trials of
			% same frequency
			stimindex = cell(nfreqs, 1);
			for f = 1:nfreqs
				stimindex{f} = find(obj.vrange(f) == freqlist);
			end
		end
	end	% END methods
	
end	% END classdef