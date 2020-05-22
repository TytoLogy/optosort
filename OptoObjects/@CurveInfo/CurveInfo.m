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
%	9 Apr 2020 (SJS): adding stimStartBin, stimEndBin.
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
		% returns stimindex{} and stimvar{} lists
		%-------------------------------------------------		
			% make sure Dinf is initialized
			if isempty(obj.Dinf)
				error('Dinf not defined/is empty')
			end
			
			% for FREQ test, find indices of stimuli with same frequency
			switch upper(obj.testtype)
				case 'FREQ'
					fprintf('\t%s test, finding indices\n', obj.testtype);
					% list of frequencies, and # of freqs tested
					freqlist = cell2mat(obj.freqs_bysweep);
					nfreqs = length(obj.varied_values);
					% locate where trials for each frequency are located in the
					% stimulus cache list - this will be used to pull out trials of
					% same frequency
					stimindex = cell(nfreqs, 1);
					for f = 1:nfreqs
						stimindex{f} = find(obj.varied_values(f) == freqlist);
					end
					% assign outputs
					varargout{1} = stimindex;
					varargout{2} = freqlist;

			% for LEVEL test, find indices of stimuli with same level (dB SPL)
				case 'LEVEL'
					fprintf('\t%s test, finding indices\n', obj.testtype);
					% list of levels, and # of levels tested
					levellist = obj.levels_bysweep;
					nlevels = length(obj.varied_values);
					% locate where trials for each frequency are located in the
					% stimulus cache list - this will be used to pull out trials of
					% same frequency
					stimindex = cell(nlevels, 1);
					for l = 1:nlevels
						stimindex{l} = find(obj.varied_values(l) == levellist);
					end
					% assign outputs
					varargout{1} = stimindex;
					varargout{2} = levellist;

			% for FRA (FREQ+LEVEL) test, find indices of stimuli with
			% freq and same level (dB SPL)
				case 'FREQ+LEVEL'
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

						values used for the two variables (Freq. and Level) are stored in
						vrange matrix, which is of length (nfreq X nlevel) and holds
						values as row 1 = freq, row 2 = level

						e.g. obj.Dinf.test.stimcache.vrange(:, 1:5) = 
							4000  4000  4000  4000  4000
							0     10    20    30    40

						trialRandomSequence holds randomized list of indices into vrange,
						has dimensions of [nreps, ntrials]

						To sort the data for FRA: (1)	for each freq and level combination,
						locate the indices for that combination in the respective FREQ and
						LEVEL list. (2)	These indices can then be used within the D{}
						array
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

			% for OPTO test...
				case 'OPTO'
					fprintf('\t%s test, finding indices\n', obj.testtype);

			% for WavFile, tell user to use WAVInfo class.
				case 'WAVFILE'
					error('%s: Please use WAVInfo class for these data', mfilename);
				
				% unknown type
				otherwise
					error('%s: unsupported test type %s', mfilename, obj.testtype);
			end
		end	% END getStimulusIndices method
		

		%-------------------------------------------------
		%-------------------------------------------------
		function titleString = getCurveTitleString(obj)
		%-------------------------------------------------
		% returns title string for curve type
		%-------------------------------------------------
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
		%-------------------------------------------------
		%-------------------------------------------------

		%-------------------------------------------------
		%-------------------------------------------------
		function [varlist, nvars] = varlist(obj)
		%---------------------------------------------------------------------
		% returns list of variable value and # of vars..
		%---------------------------------------------------------------------
			switch upper(obj.testtype)
				case {'FREQ', 'LEVEL'}
					% list of frequencies, and # of freqs tested
					% list of levels, and # of levels tested
					varlist = {obj.varied_values};
					nvars = length(varlist);

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
					warning('CurveInfo.varlist: OPTO not yet implemented');
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
		%-------------------------------------------------
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
		% variable(s) for curve (similar to test.Name
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
		% returns test.Type
		function val = testtype(obj)
			val = obj.Dinf.test.Type;
		end
		% returns test.Name
		function val = testname(obj)
			val = obj.Dinf.test.Name;
		end
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
		% builds vector of stimulus onset and offset bins referenced to start
		% of file
		[obj, varargout] = buildStimOnOffData(obj)
	end	% END methods
	
end	% END classdef
	
		