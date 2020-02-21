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
			if length(varargin) == 2
				if strcmpi(varargin{1}, 'struct')
					% this might be deprecated/removed in future
					obj = obj.initFromNexInfoStruct(varargin{2});
				elseif strcmpi(varargin{1}, 'file')
					obj = obj.initFromNexInfoFile(varargin{2});
				else
					error('Unknown input type %s', varargin{1});
				end
			else
				error('need both input mode and input value');
			end
		end

		%-------------------------------------------------
		%-------------------------------------------------
		% utility methods
		%-------------------------------------------------
		%-------------------------------------------------

	
		%-------------------------------------------------
		%-------------------------------------------------
		% access for dependent properties
		%-------------------------------------------------
		%-------------------------------------------------
		
		%-------------------------------------------------
		% # of files merged into nex file
		%-------------------------------------------------
		function val = get.nFiles(obj)
			val = length(obj.FileData);
		end
	
		
	end	% END methods
	
end	% END classdef
	
		