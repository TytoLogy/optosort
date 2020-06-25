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
% See also: CurveInfo, FRAInfo
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
%	18 Jun 2020 (SJS): broadened scope of subclass - overloaded methods
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
		%-------------------------------------------------		
	
		
		%-------------------------------------------------
		% overload method
		%-------------------------------------------------
		function varargout = convertSpikeTableToSpikeTimes(obj, spiketable)
		%-------------------------------------------------
		% = convertSpikeTableToSpikeTimes(obj, spiketable)
		%-------------------------------------------------
			
			%-----------------------------------------------------------
			% get stim indices, varlist
			%-----------------------------------------------------------
			% stimindex is a cell array with each element (corresponding to a
			% different stimulus level/parameter) consisting of a list of
			% indices into each data sweep. stimvar is a list of the variables
			% in the sweeps
			[stimindex, stimvar] = obj.getStimulusIndices;
% this worked for other curves, but won't for WAV
% 			unique_stim = unique(stimvar, 'stable');
% 			nstim = length(unique_stim);
			% use stimvar for unique_stim
			unique_stim = stimvar;
			nstim = length(unique_stim);
			%-----------------------------------------------------------
			% convert to spiketimes format (for 1-D data)
			%-----------------------------------------------------------
			% 		spikeTimes{nLevels, 1}
			% 			spikeTimes{n} = {nTrials, 1}
			% 				spikeTimes{n}{t} = [spike1_ms spike2_ms spike3ms ...
			%
			spiketimes = cell(nstim, 1);
			% loop through stimuli
			for s = 1:nstim
% 				fprintf('stimvar(%d) = %d\n', s, unique_stim(s));
				% allocate spiketimes storage
				spiketimes{s} = cell(length(stimindex{s}), 1);
				% loop through sweeps (aka trials, reps) for this stimulus
				for r = 1:length(stimindex{s})
					% get the proper index into spikeTable for this stimulus and
					% sweep combination
					rIndx = stimindex{s}(r);
					% get table for currrent sweep
					tmpT = spiketable{rIndx};
					% assign spike timestamps to spikeTimes, ...
					% converting to milliseconds
					spiketimes{s}{r} = force_row(1000 * tmpT.TS);
				end
			end			
			%-----------------------------------------------------------
			% create output struct
			%-----------------------------------------------------------
			str.spiketimes = spiketimes;
			str.stimindex = stimindex;
			str.stimvar = stimvar;
			str.unique_stim = unique_stim;
			str.nstim = nstim;
			str.spiketable = spiketable;
			varargout{1} = str;
		end
		%-------------------------------------------------
		
		%-------------------------------------------------
		function titleString = getCurveTitleString(obj)
		%-------------------------------------------------
		% returns title string for curve type
		%-------------------------------------------------
			[~, fname, fext] = fileparts(obj.Dinf.filename);
			fname = [fname '.' fext];
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
		end

		%-------------------------------------------------
		function [varlist, nvars] = varlist(obj)
		%-------------------------------------------------
		% returns list of variable value and # of vars..
		%-------------------------------------------------
			switch upper(obj.testtype)
				case {'FREQ', 'LEVEL'}
					error('use CurveInfo')
					
				case 'FREQ+LEVEL'
					error('Use FRAInfo');

				case 'OPTO'
					warning('WAVInfo.varlist: OPTO not yet implemented');
					varlist = {obj.varied_values};
					nvars = length(varlist);

				case 'WAVFILE'
					% get list of stimuli (wav file names)
					varlist = obj.getwavList;
					nvars = length(varlist);

				otherwise
					error('%s: unsupported test type %s', mfilename, cInfo.testtype);
			end
		end


		%-------------------------------------------------
		%-------------------------------------------------
		% shortcut methods to values
		%-------------------------------------------------
		%-------------------------------------------------
		% returns wavList
		function wavlist = getwavList(obj)
			% get list of stimuli (wav file names)
			nwavs = length(obj.Dinf.stimList);
			wavlist = cell(nwavs, 1);
			stimindex = cell(nwavs, 1);
			for w = 1:nwavs
				stype = obj.Dinf.stimList(w).audio.signal.Type;
				if strcmpi(stype, 'null')
					wavlist{w} = 'null';
				elseif strcmpi(stype, 'noise')
					wavlist{w} = 'BBN';
				elseif strcmpi(stype, 'wav')
					[~, wavlist{w}] = fileparts(obj.Dinf.stimList(w).audio.signal.WavFile);
				else
					error('%s: unknown type %s', mfilename, stype);
				end
				stimindex{w} = find(obj.Dinf.test.stimIndices == w);
			end
		end

		
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