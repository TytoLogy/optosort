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
			testtype
			testname
	end	% END properties (main)
	properties (Dependent)
	end	% END properties(Dependent)
	properties (Constant, Access = 'protected')
		coreProperties = {'stimtype', 'curvetype', 'freezeStim', 'nreps', ...
							'saveStim', 'ntrials', 'nstims', 'repnum', ...
							'trialnum', 'splval', 'rmsval', 'atten', ...
							'FREQ', 'LEVEL', 'opto', 'radvary', ...
							'trialRandomSequence', 'vname', 'vrange', 'stimvar'};
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

		
		function obj = initObject(obj, cacheStruct, varargin)
			% check
			if ~all(isfield(cacheStruct, obj.coreProperties))
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
			for p = 1:length(obj.coreProperties)
				obj.(obj.coreProperties{p}) = cacheStruct.(obj.coreProperties{p});
			end
			if ~isempty(varargin)
				if length(varargin) == 2
					obj.testtype = varargin{1};
					obj.testname = varargin{2};
				else
					error('incorrect input args');
				end
			else
				obj.testtype = '';
				obj.testname = '';
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
	
		