classdef StimulusCache
	properties
		%{
			stimtype: [116 111 110 101]
			curvetype: [70 82 69 81]
			freezeStim: 0
			nreps: 10
			saveStim: 0
			ntrials: 14
			nstims: 140
			repnum: [140×1 double]
			trialnum: [140×1 double]
			splval: {140×1 cell}
			rmsval: {140×1 cell}
			atten: {140×1 cell}
			FREQ: {140×1 cell}
			LEVEL: [140×1 double]
			opto: {140×1 cell}
			radvary: 1
			trialRandomSequence: [10×14 double]
			vname: [70 82 69 81]
			vrange: [1×14 double]
			stimvar: {1×140 cell}
		%}
			stimtype
			curvetype
			freezeStim
			nreps
			saveStim
			ntrials
			nstims
			repnum
			trialnum
			splval = {}
			rmsval =  {}
			atten =  {}
			FREQ =  {}
			LEVEL
			opto =  {}
			radvary
			trialRandomSequence
			vname
			vrange
			stimvar =  {}
	end	% END properties (main)
	properties (Dependent)
	end	% END properties(Dependent)
	
	methods
		%-------------------------------------------------
		%-------------------------------------------------
		% Constructor
		%-------------------------------------------------
		%-------------------------------------------------
		function obj = StimulusCache(varargin)
			if isempty(varargin)
				return
			end
			if isstruct(varargin{1})				
				obj = initObjectFromStruct(obj, varargin{1});
			else
				error('Unknown input type %s', varargin{1});
			end
		end

		
		function obj = initObjectFromStruct(obj, cacheStruct)
			% get list of properties
			proplist = properties(obj);
			% check
			if ~all(isfield(cacheStruct, proplist))
				error('field/property mismatch');
			end
			% check to see if character fields are stored as chars or need to
			% be converted from ASCII integers
			charfields = {'stimtype', 'curvetype', 'vname'};
			for c = 1:length(charfields)
				if isfield(cacheStruct, charfields{c})
					if ~ischar(cacheStruct.(charfields{c}))
						cacheStruct.(charfields{c}) = ...
															char(cacheStruct.(charfields{c}));
					end
				end
			end
			% assign properties from struct
			for p = 1:length(proplist)
				% check if property name is one of the charfields
				ccheck = strcmpi(proplist{p}, charfields);
				if any(ccheck)
					% if so, convert to char
					obj.(proplist{p}) = char(cacheStruct.(proplist{p}));
				else
					% otherwise, simple assignment
					obj.(proplist{p}) = cacheStruct.(proplist{p});
				end
			end
		end
		%-------------------------------------------------
		%-------------------------------------------------
		% utility methods
		%-------------------------------------------------
		%-------------------------------------------------
		% superclass implementation is trivial
		function val = stimulusIndex(obj) %#ok<MANU>
			val = {};
		end
		%-------------------------------------------------
		%-------------------------------------------------
		% access for dependent properties
		%-------------------------------------------------
		%-------------------------------------------------
		
	
		
	end	% END methods
	
end	% END classdef
	
		