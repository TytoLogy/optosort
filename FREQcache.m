classdef FREQcache < StimulusCache
	properties

		
	end	% END properties (main)
	properties (Dependent)
	end	% END properties(Dependent)

	methods
		function [stimindex, freqlist] = getStimulusIndices(obj)
			fprintf('\t%s test, finding indices\n', obj.curvetype);
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