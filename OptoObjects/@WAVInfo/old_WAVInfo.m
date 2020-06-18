classdef WAVInfo < CurveInfo
%------------------------------------------------------------------------
% Class: WAVInfo
%------------------------------------------------------------------------
% TytoLogy:Experiments:optosort
%------------------------------------------------------------------------
% implements and encapsulates some utilities for dealing with
% opto data -> Dinf struct
%------------------------------------------------------------------------
% class properties
%------------------------------------------------------------------------
% 
% wavOnsetBin		sample for onset of wav stimulus playback
%						this is effectively the same as stimStartBin for
%						CurveInfo. For WAV stimuli, stimStart and stimEnd refer
%						to the actual onset of stimulus waveform. This is done
%						because stimStart/End are typically used for analysis of
%						the neural response to a stimulus and we want these
%						onset/offset bin values to be consistent across stimuli
%
% Dependent properties:
% 
%------------------------------------------------------------------------
% See also: 
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%  Sharad J. Shanbhag
%	sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 24 February, 2020 (SJS)
%	- adapted from import_from_plexon_nonObj
% Revisions:
%	3 Mar 2020 (SJS): adding elements from fData struct in the 
%		export_for_plexon.m function to avoid future duplications and
%		streamline curve/test information handling
%	15 Apr 2020 (SJS): adding code to determine stim onset, offset and wav
%	onset
%------------------------------------------------------------------------
% TO DO:
%------------------------------------------------------------------------

	%-------------------------------------------------
	% class properties
	%-------------------------------------------------
	properties
		wavInfo
		wavOnsetBin = [];
	end	% END properties (main)
	properties (Dependent)

	end	% END properties(Dependent)
	
	methods
		%-------------------------------------------------
		%-------------------------------------------------
		% Constructor
		%-------------------------------------------------
		%-------------------------------------------------
		function obj = WAVInfo(varargin)
			% invoke superclass constructor - this initializes obj.Dinf if it
			% was provided as input
			obj@CurveInfo(varargin{1})
			% return if nothing else to do (no Dinf)
			if isempty(varargin)
				return
			end
			if isstruct(varargin{1})
				testfields = {'ScriptType', 'optovar_name', 'audiovar_name', ...
										'audiovar', 'curvetype'};
				for n = 1:length(testfields)
					obj.Dinf.test.(testfields{n}) = ...
								convert_to_text(obj.Dinf.test.(testfields{n}));
				end
			else
				error('Unknown input type %s', varargin{1});
			end
			% if wavinfo is provided, use it
			if nargin == 2
				if isstruct(varargin{2})
					obj.wavInfo = obj.init_wavInfo(varargin{2});
				else
					error('Unknown type for wavInfo input');
				end
			end
		end

		%-------------------------------------------------
		function wInf = init_wavInfo(obj, W) %#ok<INUSL>
		%-------------------------------------------------
		% takes variables read in from _wavinfo.mat file, stores in new
		% struct
		% don't need to store stimList, stimIndices!
		%-------------------------------------------------		
			% wavInfo needs to be remapped to differently named field
			% nwavs struct array with information about wav stimulus
			wInf.wavs = W.wavInfo;
			% store waveform inside .wavs field
			if length(W.wavInfo) ~= length(W.wavS0)
				error('mismatch in wavInfo and wavS0 lengths');
			else
				for n = length(W.wavInfo)
					wInf.wavs(n).wavS0 = W.wavS0{n};
				end
			end
			% audio has information about general audio properties, might be
			% redundant or not totally accurate (i.e. level info)
			wInf.audio = W.audio;
			% nullstim has information about the null (no) stimulus
			wInf.nullstim = W.nullstim;
			% noise has information about noise stimulus (BBN)
			wInf.noise = W.noise;
		end
		
		%-------------------------------------------------
		%-------------------------------------------------
		% overload the CurveInfo method
		% returns stimulus Indices and list of stim variables
		%-------------------------------------------------
		%-------------------------------------------------
		function varargout = getStimulusIndices(obj)
		%-------------------------------------------------
		% returns stimindex{} and stimvar{} lists
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
		% how to deal with different levels????
		%-------------------------------------------------
			% make sure Dinf is initialized
			if isempty(obj.Dinf)
				error('%s: Dinf not defined/is empty', mfilename)
			end
			
			% for WavFile, need to find indices with same filename.
			fprintf('\t%s test, finding indices\n', obj.testtype);
			% get list of stimuli (wav file names)
			nwavs = length(obj.Dinf.stimList);
			wavlist = cell(nwavs, 1);
			stimindex = cell(nwavs, 1);
			% loop through 
			for w = 1:nwavs
				stype = obj.Dinf.stimList(w).audio.signal.Type;
				if strcmpi(stype, 'null')
					wavlist{w} = 'null';
				elseif strcmpi(stype, 'noise')
					wavlist{w} = 'BBN';
				elseif strcmpi(stype, 'wav')
					[~, wavlist{w}] = ...
						fileparts(obj.Dinf.stimList(w).audio.signal.WavFile);
				else
					error('%s: unknown type %s', mfilename, stype);
				end
				stimindex{w} = find(obj.Dinf.test.stimIndices == w);
			end
			% assign outputs
			varargout{1} = stimindex;
			varargout{2} = wavlist;
		end	% END getStimulusIndices method
		
		
		%{
		function titleString = getCurveTitleString(obj)
		% returns title string for curve type
		
			[~, fname, fext] = fileparts(obj.Dinf.filename);
			fname = [fname '.' fext];
			switch obj.testtype
				case 'FREQ'
					% list of frequencies, and # of freqs tested
					varlist = obj.varied_values;
					nvars = length(varlist);
					titleString = cell(nvars, 1);
					for v = 1:nvars
						if v == 1
							titleString{v} = {fname, ...
													sprintf('Frequency = %.0f kHz', ...
																			0.001*varlist(v))};
						else
							titleString{v} = sprintf('Frequency = %.0f kHz', ...
													0.001*varlist(v));
						end
					end
				case 'LEVEL'
					% list of levels, and # of levels tested
					varlist = obj.varied_values;
					nvars = length(varlist);
					titleString = cell(nvars, 1);
					for v = 1:nvars
						if v == 1
							titleString{v} = {fname, sprintf('Level = %d dB SPL', ...
																					varlist(v))};
						else
							titleString{v} = sprintf('Level = %d dB SPL', varlist(v));
						end
					end
				case 'FREQ+LEVEL'
					% list of freq, levels
					varlist = cell(2, 1);
					% # of freqs in nvars(1), # of levels in nvars(2)
					nvars = zeros(2, 1);
					for v = 1:2
						varlist{v} = unique(obj.varied_values(v, :), 'sorted');
						nvars(v) = length(varlist{v});
					end
					titleString = fname;

				case 'OPTO'
					% not yet implemented
					
				case 'WAVFILE'
					% get list of stimuli (wav file names)
					varlist = obj.Dinf.test.wavlist;
					nvars = length(varlist);
					titleString = cell(nvars, 1);
					for v = 1:nvars
						if v == 1 
							titleString{v} = {fname, sprintf('wav name: %s', varlist{v})};
						else
							titleString{v} = sprintf('wav name: %s', varlist{v});
						end
					end
				otherwise
					error('%s: unsupported test type %s', mfilename, obj.testtype);
			end			
		end
		%}
		

		function [varlist, nvars] = varlist(obj)
		%---------------------------------------------------------------------
		% returns list of variable value and # of vars..
		% overloaded for WAV
		%---------------------------------------------------------------------
			switch upper(obj.testtype)
				case {'FREQ', 'LEVEL', 'FREQ+LEVEL'}
					% wrong class
					error('class mismatch in varlist')

				case 'OPTO'
					warning('WAVInfo.varlist: OPTO not yet implemented');
					varlist = {obj.varied_values};
					nvars = length(varlist);

				case 'WAVFILE'
					% get list of stimuli (wav file names)
					varlist = obj.Dinf.test.wavlist;
					nvars = length(varlist);

				otherwise
					error('%s: unsupported test type %s', mfilename, cInfo.testtype);
			end
		end


%{		
		%-------------------------------------------------
		%-------------------------------------------------
		% shortcut methods to values
		%-------------------------------------------------
		%-------------------------------------------------
		% returns test.stimcache.LEVELS, which is a list of db SPL 
		% stimulus levels used for each stimulus sweep
		%	this is a numerical array [nsweeps, 1]
		function val = levels_bysweep(obj)
			if obj.has_stimList
				levels = zeros(obj.ntrials, 1);
				for l = 1:obj.ntrials
					levels(l) = obj.Dinf.stimList(l).audio.Level;
				end
				val = levels(obj.Dinf.test.stimIndices);
			else
				val = [];
			end
		end
		% returns test.stimcache.vname, char string identifying
		% variable(s) for curve (similar to test.Name
		function val = varied_parameter(obj)
			if obj.has_stimList
				val = {'Level'; 'WAV'};
			else
				val = [];
			end
		end
		% returns test.stimcache.vrange, values of varied parameter
		function val = varied_values(obj)
			if obj.has_stimList
				val = {obj.Dinf.test.Level'; ...
							unique(obj.Dinf.test.wavlist, 'stable')};
			else
				val = [];
			end
		end
		% returns test.Reps: # of reps for each stimulus
		function val = nreps(obj)
			if obj.has_stimList
				val = obj.Dinf.test.Reps;
			else
				val = [];
			end
		end
		% returns test.stimcache.ntrials: # of stimulus types
		function val = ntrials(obj)
			if obj.has_stimList
				val = obj.Dinf.test.nCombinations;
			else
				val = [];
			end
		end
		% returns test.stimcache.nstims: total # of stimulus presentations
		% (usually equal to nreps * ntrials
		function val = nstims(obj)
			if obj.has_stimList
				val = length(obj.Dinf.test.stimIndices);
			else
				val = [];
			end
		end
%}		
		
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
		% builds vector of stimulus onset and offset bins referenced to start
		% of file
		[obj, varargout] = buildStimOnOffData(obj)		
	end	% END methods
	
end	% END classdef
	
function out = convert_to_text(in)
% convert double ascii arrays or cell array of ascii arrays to char
% this should really be moved to Utilities toolbox...
	if iscell(in)
		n_in = numel(in);
		out = cell(size(in));
		for n = 1:n_in
			out{n} = char(in{n});
		end
	else
		out = char(in);
	end
end