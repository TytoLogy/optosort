classdef CurveInfo
%------------------------------------------------------------------------
% Class: CurveInfo
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
% Dinf				Data information struct from opto .dat files
% F					opto file object
% startSweepBin	sample for start of each sweep (for each channel)
%							cell array
% endSweepBin		sample for end of each sweep (for each channel)
%							cell array
% sweepLen			length (# of samples) for each sweep (for each channel)
%							vector
% stimStartBin		sample for stimulus onset (vector)
% stimEndBin		sample for stimulus offset (vector)
% fileStartBin		sample for start of file in merged file
% fileEndBin		sample for end of file in merged data file
% 
% validSweep		logical vector [nsweeps, 1], true for sweeps with valid
%						(no artifact, etc) data
%
% Dependent properties:
% 	testtype
% 	testname
% 	freqs_bysweep
% 	levels_bysweep
% 	varied_parameter
% 	varied_values
% 	analysis_window
% 	nreps
% 	ntrials
% 	nstims
% 	ADFs
% 	DAFs
%------------------------------------------------------------------------
% See also: Subclasses
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
%	9 Apr 2020 (SJS): adding stimStartBin, stimEndBin.
%	9 Jun 2020 (SJS): modifying to create test data for checking timing
%	15 Jun 2020 (SJS): fixed issue with trying to directly index into 
%							varied_values dependent property
%	2 Jul 2020 (SJS): pulled out FREQTUNING specific code into subclass
%  16 Jun 2021 (SJS): subclass OptoInfo to handle OPTO and OPTO-AMP,
%                     altered code here to deal with this
%------------------------------------------------------------------------
% TO DO:
%------------------------------------------------------------------------

	%-------------------------------------------------
	% class properties
	%-------------------------------------------------
	properties
		Dinf
		F
		startSweepBin = {}
		endSweepBin = {}
		stimStartBin = [];
		stimEndBin = [];
		sweepLen
		fileStartBin
		fileEndBin
		% validSweep = true(0); NOT SURE WHERE THIS SHOULD LIVE - or at what
		% level....??????? try here for now
		validSweep = true(0);
	end	% END properties (main)
	properties (Dependent)
		testtype
		testname
% 		freqs_bysweep
% 		levels_bysweep
% 		varied_parameter
% 		varied_values
		analysis_window
% 		nreps
% 		ntrials
% 		nstims
		ADFs
		DAFs
	end	% END properties(Dependent)
	properties (Access = protected)
		has_stimcache = 0;
		has_stimList = 0;
	end
	
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
		function obj = CurveInfo(varargin)
			if isempty(varargin)
				return
			end
			if isstruct(varargin{1})
				% save Dinf
				obj.Dinf = varargin{1};
				% extract filename, convert to optofilename object
				obj.F = OptoFileName(obj.Dinf.filename);
				% initialize validSweep
				obj.validSweep = true(obj.Dinf.nread, 1);
				% if necessary, convert stimtype and curvetype to strings
				% not all tests (WAV) have stimcache...
				if isfield(obj.Dinf.test, 'stimcache')
					if ~isempty(obj.Dinf.test.stimcache)
						obj.has_stimcache = 1;
						obj.Dinf.test.stimcache.stimtype = ...
													char(obj.Dinf.test.stimcache.stimtype);
						obj.Dinf.test.stimcache.curvetype = ...
													char(obj.Dinf.test.stimcache.curvetype);
					else
						obj.has_stimcache = 0;
					end
				else
					obj.has_stimcache = 0;
				end
				% not all tests (WAV) have stimList...
				if isfield(obj.Dinf, 'stimList')
					if ~isempty(obj.Dinf.stimList)
						obj.has_stimList = 1;
					else
						obj.has_stimList = 0;
					end
				else
					obj.has_stimList = 0;
				end
				% signal name should be a char
				if isfield(obj.Dinf, 'audio')
					if isfield(obj.Dinf.audio, 'signal')
						if isfield(obj.Dinf.audio.signal, 'Type')
							obj.Dinf.audio.signal.Type = ...
											char(obj.Dinf.audio.signal.Type);
						end
					end
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
				error(['CurveInfo.getStimulusIndices: ' ...
							'Dinf not defined/is empty']);
			end
			
			% for FREQ test, find indices of stimuli with same frequency
			switch upper(obj.testtype)
				case 'FREQ'
					error('%s: Please use FreqTuningInfo class for these data', ...
								mfilename);

			% for LEVEL test, find indices of stimuli with same level (dB SPL)
				case 'LEVEL'
					fprintf('\t%s test, finding indices\n', obj.testtype);
					% list of levels by sweep, levels and # of levels tested
					levellist = obj.levels_bysweep;
					levels = obj.varied_values;
					nlevels = length(levels);
					% locate where trials for each level are located in the
					% stimulus cache list - this will be used to pull out trials
					% of same level
					stimindex = cell(nlevels, 1);
					for l = 1:nlevels
						stimindex{l} = find(levels(l) == levellist);
					end
					% assign outputs
					varargout{1} = stimindex;
					varargout{2} = levellist;

			% for FRA (FREQ+LEVEL) test, need to use subclass FRAInfo
				case 'FREQ+LEVEL'
					error('%s: Please use FRAInfo class for these data', ...
								mfilename);

			% for WavFile, tell user to use WAVInfo class.
				case 'WAVFILE'
					error('%s: Please use WAVInfo class for these data', ...
								mfilename);
                     
			% for OPTO-AMP test...
				case {'OPTO', 'OPTO-AMP'}
					error('%s: Please use OptoInfo class for these data', ...
								mfilename);
				
			% unknown type
				otherwise
					error('%s: unsupported test type %s', mfilename, ...
								obj.testtype);
			end
		end	% END getStimulusIndices method
		%-------------------------------------------------

		
		%-------------------------------------------------
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
			unique_stim = unique(stimvar, 'sorted');
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

		%-------------------------------------------------
		%-------------------------------------------------
		function titleString = getCurveTitleString(obj)
		%-------------------------------------------------
		% returns title string for curve type
		%-------------------------------------------------
			[~, fname, fext] = fileparts(path_unix(obj.Dinf.filename));
			fname = [fname '.' fext];
			switch obj.testtype
				case 'FREQ'
						error(['CurveInfo.getCurveTitleString: '...
									'Use FreqTuningInfo subclass']);

				case 'LEVEL'
					% list of levels, and # of levels tested
					varlist = obj.varied_values;
					nvars = length(varlist);
					titleString = cell(nvars, 1);
					for v = 1:nvars
						if v == 1
							titleString{v} = {	fname, ...
														sprintf('Level = %d dB SPL', ...
																			varlist(v))};
						else
							titleString{v} = sprintf('Level = %d dB SPL', ...
																varlist(v));
						end
					end
					
				case 'FREQ+LEVEL'
					error(['CurveInfo.getCurveTitleString: '...
									'Use FRAInfo subclass']);

				case 'OPTO'
					% not yet implemented
					
				case 'OPTO-AMP'
					error(['CurveInfo.getCurveTitleString: '...
									'Use OptoAmpInfo subclass']);
               
				case 'WAVFILE'
						error(['CurveInfo.getCurveTitleString: '...
									'Use WAVInfo subclass']);
					
				otherwise
					error('%s: unsupported test type %s', ...
								'CurveInfo.getCurveTitleString', ...
								obj.testtype);
			end
		end
		%-------------------------------------------------
		%-------------------------------------------------

		%-------------------------------------------------
		%-------------------------------------------------
		function [varlist, nvars] = varlist(obj)
		%---------------------------------------------------------------------
		% returns list of variable value and # of vars..
		%---------------------------------------------------------------------
			switch upper(obj.testtype)
				case {'FREQ'}
					error('CurveInfo.varlist: Use FreqTuningInfo subclass');
					
				case {'LEVEL'}
					% list of levels, and # of levels tested
					varlist = {obj.varied_values};
					nvars = length(varlist);

				case 'FREQ+LEVEL'
					error('CurveInfo.varlist: Use FRAINfo subclass');

				case 'OPTO'
					warning('CurveInfo.varlist: OPTO not yet implemented');
					varlist = {obj.varied_values};
					nvars = length(varlist);

				case 'WAVFILE'
					error('CurveInfo.varlist: Use WAVINfo subclass');

				otherwise
					error('%s: unsupported test type %s', ...
								'CurveInfo.varlist', ...
								obj.testtype);
			end
		end
		%-------------------------------------------------
		%-------------------------------------------------

		%-------------------------------------------------
		% returns eventlist (stimulus type, levels)
		%-------------------------------------------------
		function events = geteventList(obj)
		%-------------------------------------------------
		% get list of stimuli (wav file names)
		%-------------------------------------------------
			% get stimulus indices
			[stimindex, ~] = obj.getStimulusIndices;
			% get varied values (wav files)
			varied_values = obj.varied_values;
			% get stimulus onset offset bins
			% to align to appended/merged file, will need to add
			% SpikeInfo.fileStartBin(findx) - 1
			onsetbins = obj.stimStartBin;
% 			offsetbins = obj.stimEndBin;

			%------------------------------------------------------------------------
			% create eventList as struct array
			%------------------------------------------------------------------------
			% get # of events
			nevents = length(varied_values);
			% init events struct
			events = repmat(	struct(	'name', '', ...
												'samples', [], ...
												'timestamps', [] ), ...
									nevents, 1);
			% format string depends on test type
			switch upper(obj.testtype)
				case 'LEVEL'
					formatstr = '%s_%s_%ddB';
					for n = 1:nevents
						events(n).name = sprintf(formatstr, obj.testtype, ...
																obj.testname, varied_values(n));
						events(n).samples = onsetbins(stimindex{n});
					end

				otherwise
					error('CurveInfo: use subclass or unsupported %s', obj.testtype);
			end
		end
		%-------------------------------------------------
		
		%-------------------------------------------------
		%-------------------------------------------------
		% shortcut methods to stimcache values
		%-------------------------------------------------
		%-------------------------------------------------
	
		% returns test.stimcache.FREQS, which is a list of frequencies (or
		% freq ranges for BBN) used for each stimulus sweep
		%	this is a cell array {nsweeps, 1}
		function val = freqs_bysweep(obj)
			if obj.has_stimcache
				val = obj.Dinf.test.stimcache.FREQ;
			else
				val = [];
			end
		end
		% returns test.stimcache.LEVELS, which is a list of db SPL 
		% stimulus levels used for each stimulus sweep
		%	this is a numerical array [nsweeps, 1]
		function val = levels_bysweep(obj)
			if obj.has_stimcache
				val = obj.Dinf.test.stimcache.LEVEL;
			else
				val = [];
			end
		end
		% returns test.stimcache.vname, char string identifying
		% variable(s) for curve (similar to test.Name)
		function val = varied_parameter(obj)
			if obj.has_stimcache
				val = char(obj.Dinf.test.stimcache.vname);
			else
				val = [];
			end
		end
		% returns test.stimcache.vrange, values of varied parameter
		function val = varied_values(obj)
			if obj.has_stimcache
				val = obj.Dinf.test.stimcache.vrange;
			else
				val = [];
			end
		end
		% returns test.stimcache.nreps: # of reps for each stimulus
		function val = nreps(obj)
			if obj.has_stimcache
				val = obj.Dinf.test.stimcache.nreps;
			else
				val = [];
			end
		end
		% returns test.stimcache.ntrials: # of stimulus types
		function val = ntrials(obj)
			if obj.has_stimcache
				val = obj.Dinf.test.stimcache.ntrials;
			else
				val = [];
			end
		end
		% returns test.stimcache.nstims: total # of stimulus presentations
		% (usually equal to nreps * ntrials
		function val = nstims(obj)
			if obj.has_stimcache
				val = obj.Dinf.test.stimcache.nstims;
			else
				val = [];
			end
		end
		
		%-------------------------------------------------
		%-------------------------------------------------
		% get/set access for dependent properties
		%-------------------------------------------------
		%-------------------------------------------------
		% returns test.Type
		function val = get.testtype(obj)
			val = obj.Dinf.test.Type;
		end
		% returns test.Name
		function val = get.testname(obj)
			val = obj.Dinf.test.Name;
		end
		% returns analysis_window : audio Delay to Delay+Duration interval
		function val = get.analysis_window(obj)
			val = [obj.Dinf.audio.Delay ...
							(obj.Dinf.audio.Delay + obj.Dinf.audio.Duration)];
		end
		% returns Dinf.indev.Fs
		function val = get.ADFs(obj)
			val = obj.Dinf.indev.Fs;
		end
		% returns Dinf.outdev.Fs
		function val = get.DAFs(obj)
			val = obj.Dinf.outdev.Fs;
		end

		%{
		% returns test.stimcache.FREQS, which is a list of frequencies (or
		% freq ranges for BBN) used for each stimulus sweep
		%	this is a cell array {nsweeps, 1}
		function val = freqs_bysweep(obj)
			val = obj.Dinf.test.stimcache.FREQ;
		end
		% returns test.stimcache.LEVELS, which is a list of db SPL 
		% stimulus levels used for each stimulus sweep
		%	this is a numerical array [nsweeps, 1]
		function val = levels_bysweep(obj)
			val = obj.Dinf.test.stimcache.LEVEL;
		end
		% returns test.stimcache.vname, char string identifying
		% variable(s) for curve (similar to test.Name
		function val = varied_parameter(obj)
			val = char(obj.Dinf.test.stimcache.vname);
		end
		% returns test.stimcache.vrange, values of varied parameter
		function val = varied_values(obj)
			val = obj.Dinf.test.stimcache.vrange;
		end
		% returns analysis_window : audio Delay to Delay+Duration interval
		function val = analysis_window(obj)
			val = [obj.Dinf.audio.Delay (obj.Dinf.audio.Delay + obj.Dinf.audio.Duration)];
		end
		% returns test.stimcache.nreps: # of reps for each stimulus
		function val = nreps(obj)
			if obj.has_stimcache
				val = obj.Dinf.test.stimcache.nreps;
			elseif obj.has_stimList
				val = obj.Dinf.test.Reps;
			end
		end
		% returns test.stimcache.ntrials: # of stimulus types
		function val = ntrials(obj)
			val = obj.Dinf.test.stimcache.ntrials;
		end
		% returns test.stimcache.nstims: total # of stimulus presentations
		% (usually equal to nreps * ntrials
		function val = nstims(obj)
			val = obj.Dinf.test.stimcache.nstims;
		end		
		% returns Dinf.indev.Fs
		function val = ADFs(obj)
			val = obj.Dinf.indev.Fs;
		end
		% returns Dinf.outdev.Fs
		function val = DAFs(obj)
			val = obj.Dinf.outdev.Fs;
		end
		%}
		
		%-------------------------------------------------
		%-------------------------------------------------
		% methods defined in separate files
		%-------------------------------------------------
		%-------------------------------------------------
		% get data  and sweep info for each channel
		[obj, varargout] = buildChannelData(obj, Channels, BPfilt, ...
															D, varargin)
		% create dummy/test data  and sweep info for each channel
		[obj, varargout] = buildTimingTestData(obj, Channels, BPfilt, ...
															D, varargin)
		% builds vector of stimulus onset and offset bins referenced to start
		% of file
		[obj, varargout] = buildStimOnOffData(obj)
		% plot PSTH
		H = plotPSTH(obj, ST, binSize)
	end	% END methods
	
end	% END classdef
	
		