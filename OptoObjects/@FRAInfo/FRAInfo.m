classdef FRAInfo < CurveInfo
%------------------------------------------------------------------------
% Class: FRAInfo
%------------------------------------------------------------------------
% TytoLogy:Experiments:optosort
%------------------------------------------------------------------------
% implements and encapsulates some utilities for dealing with
% opto data -> Dinf struct
%
% code pulled from opto:getFilteredOptoData
% 
%------------------------------------------------------------------------
% class properties
%------------------------------------------------------------------------
%------------------------------------------------------------------------
% See also: CurveInfo, WAVInfo
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%  Sharad J. Shanbhag
%	sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 18 June 2020, 2020 (SJS)
%	- subclassed from CurveInfo
% Revisions:
%------------------------------------------------------------------------
% TO DO:
%------------------------------------------------------------------------

	%-------------------------------------------------
	% class properties
	%-------------------------------------------------
	properties

	end	% END properties (main)
	properties (Dependent)

	end	% END properties(Dependent)
	
	%-------------------------------------------------
	%-------------------------------------------------
	% methods defined here
	%-------------------------------------------------
	%-------------------------------------------------
	methods
		%-------------------------------------------------
		%-------------------------------------------------
		% Constructor
		%-------------------------------------------------
		%-------------------------------------------------
		function obj = FRAInfo(varargin)
			% invoke superclass constructor - this initializes obj.Dinf if it
			% was provided as input
			obj@CurveInfo(varargin{1})
			% return if nothing else to do (no Dinf)
			if isempty(varargin)
				return
			end
		end
		
		%-------------------------------------------------
		%-------------------------------------------------
		% returns stimulus Indices and list of stim variables
		%-------------------------------------------------
		%-------------------------------------------------
		function varargout = getStimulusIndices(obj)
		%-------------------------------------------------
		% returns stimindex{} and stimvar() lists
		%	stimindex	{# unique stimuli, 1} cell array where each
		%					element is a [nreps, 1] vector of indices into
		%					the total list of stimulus sweeps (trials). so, 
		%					if there are 8 sound levels and 20 reps of each 
		%					stimulus in a rate-level curve, then indices will be
		%					values from 1-160
		%					For data like FRA, this will be {nlevels, nfreqs}
		%					cell array (2 variables)
		%					For WAV data see WAVInfo subclass of CurveInfo
		%	stimvar		[total # of sweeps, 1] vector of varied values
		%					e.g., db SPL levels for rate-level curve
		%-------------------------------------------------		
			% make sure Dinf is initialized
			if isempty(obj.Dinf)
				error(['FRAInfo.getStimulusIndices: ' ...
							'Dinf not defined/is empty']);
			end

			% for FRA (FREQ+LEVEL) test, find indices of stimuli with
			% freq and same level (dB SPL)
			fprintf('\t%s test, finding freq and level indices\n', ...
																	obj.testtype);

			% if necessary, convert cells to matrices
			testcell = {'splval', 'rmsval', 'atten', 'FREQ', 'LEVEL'};
			for c = 1:length(testcell)
				if iscell(obj.Dinf.test.stimcache.(testcell{c}))
					obj.Dinf.test.stimcache.(testcell{c}) = ...
							cell2mat(obj.Dinf.test.stimcache.(testcell{c}));
				end
			end
			% list of stimulus freqs, # of freqs tested
			freqlist = unique(obj.freqs_bysweep, 'sorted');
			nfreqs = length(freqlist);
			% list of stimulus levels, # of levels tested
			levellist = unique(obj.levels_bysweep, 'sorted');
			nlevels = length(levellist);
			%{
			Raw data are in a vector of length nstims, in order of
			presentation.

			values used for the two variables (Freq. and Level) are
			stored in vrange matrix, which is of length (nfreq X nlevel)
			and holds values as row 1 = freq, row 2 = level

			e.g. obj.Dinf.test.stimcache.vrange(:, 1:5) = 4000  4000
			4000  4000  4000 0     10    20    30    40

			trialRandomSequence holds randomized list of indices into
			vrange, has dimensions of [nreps, ntrials]

			To sort the data for FRA: (1)	for each freq and level
			combination, locate the indices for that combination in the
			respective FREQ and LEVEL list. (2)	These indices can then
			be used within the D{} array
			%}
			stimindex = cell(nlevels, nfreqs);
			for f = 1:nfreqs
				for l = 1:nlevels
					currentF = freqlist(f);
					currentL = levellist(l);
					stimindex{l, f} = ...
						find( (obj.freqs_bysweep == currentF) & ...
								(obj.levels_bysweep == currentL) );
				end
			end
			% assign outputs
			varargout{1} = stimindex;
			varargout{2} = {freqlist levellist};

		end	% END getStimulusIndices method
		

		%-------------------------------------------------
		function titleString = getCurveTitleString(obj)
		%-------------------------------------------------
		% returns title string for curve type
		%-------------------------------------------------
			[~, fname, fext] = fileparts(obj.Dinf.filename);
			fname = [fname '.' fext];
			% list of freq, levels
			varlist = cell(2, 1);
			% # of freqs in nvars(1), # of levels in nvars(2)
			nvars = zeros(2, 1);
			for v = 1:2
				varlist{v} = unique(obj.varied_values(v, :), 'sorted');
				nvars(v) = length(varlist{v});
			end
			titleString = fname;
		end
		%-------------------------------------------------

		%-------------------------------------------------
		function [varlist, nvars] = varlist(obj)
		%-------------------------------------------------
		% returns list of variable value and # of vars..
		% use overloaded methods depending on curve/test type
		%-------------------------------------------------
			switch upper(obj.testtype)
				case {'FREQ', 'LEVEL'}
					error('use CurveInfo')
				case 'WAVFILE'
					error('use WAVInfo')

				case 'FREQ+LEVEL'
					% list of freq, levels
					varlist = cell(2, 1);
					% # of freqs in nvars(1), # of levels in nvars(2)
					nvars = zeros(2, 1);
					tmprange = obj.varied_values;
					for v = 1:2
						varlist{v} = unique(tmprange(v, :), 'sorted');
						nvars(v) = length(varlist{v});
					end

				case 'OPTO'
					warning('FRAInfo.varlist: OPTO not yet implemented');
					varlist = {obj.varied_values};
					nvars = length(varlist);

				otherwise
					error('%s: unsupported test type %s', mfilename, cInfo.testtype);
			end
		end
		%-------------------------------------------------
		%-------------------------------------------------

		
		%-------------------------------------------------
		%-------------------------------------------------
		% shortcut methods to stimcache values
		%-------------------------------------------------
		%-------------------------------------------------
	
		
		%-------------------------------------------------
		%-------------------------------------------------
		% get/set access for dependent properties
		%-------------------------------------------------
		%-------------------------------------------------
	
		
		%-------------------------------------------------
		%-------------------------------------------------
		% methods defined in separate files
		%-------------------------------------------------
		%-------------------------------------------------
		
	end	% END methods
	
end	% END classdef
	
		