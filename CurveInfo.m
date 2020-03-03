classdef CurveInfo
%------------------------------------------------------------------------
% Class: CurveInfo
%------------------------------------------------------------------------
% TytoLogy:Experiments:optosort
%------------------------------------------------------------------------
% implements and encapsulates some utilities for dealing with
% optodata
% code pulled from opto:getFilteredOptoData
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
%------------------------------------------------------------------------
% TO DO:
%------------------------------------------------------------------------

	%-------------------------------------------------
	% class properties
	%-------------------------------------------------
	properties
		%{
			Dinf		Data information struct from opto .dat files
			startSweepBin		sample for start of each sweep 
			endSweepBin		sample for end of each sweep
			sweepLen			length (# of samples) for each sweep
			fileStartBin	sample for start of file in merged file
			fileEndBin		sample for end of file in merged data file
		%}
		Dinf
		startSweepBin = {}
		endSweepBin = {}
		sweepLen
		fileStartBin
		fileEndBin
	end	% END properties (main)
	properties (Dependent)
		testtype
		testname
		freqs_bysweep
		levels_bysweep
		varied_parameter
		varied_values
		analysis_window
	end	% END properties(Dependent)

	
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
				obj.Dinf = varargin{1};
				% if necessary, convert stimtype and curvetype to strings
				if isnumeric(obj.Dinf.test.stimcache.stimtype)
					obj.Dinf.test.stimcache.stimtype = ...
												char(obj.Dinf.test.stimcache.stimtype);
				end
				if isnumeric(obj.Dinf.test.stimcache.curvetype)
					obj.Dinf.test.stimcache.curvetype = ...
												char(obj.Dinf.test.stimcache.curvetype);
				end
% 			elseif isnumeric(varargin{1})
% 				% handles case where varargin is numeric i.e. an array is being
% 				% built
% 				return
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

			% for WavFile, need to find indices with same filename.
				case 'WAVFILE'
					fprintf('\t%s test, finding indices\n', obj.testtype);
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
				
				% unknown type
				otherwise
					error('%s: unsupported test type %s', mfilename, obj.testtype);
			end
		
	
		end	% END getStimulusIndices method
		
		
		
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
		
		
		%-------------------------------------------------
		%-------------------------------------------------
		% access for dependent properties
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
		% returns test.stimcache.FREQS, which is a list of frequencies (or
		% freq ranges for BBN) used for each stimulus sweep
		%	this is a cell array {nsweeps, 1}
		function val = get.freqs_bysweep(obj)
			val = obj.Dinf.test.stimcache.FREQ;
		end
		% returns test.stimcache.LEVELS, which is a list of db SPL 
		% stimulus levels used for each stimulus sweep
		%	this is a numerical array [nsweeps, 1]
		function val = get.levels_bysweep(obj)
			val = obj.Dinf.test.stimcache.LEVEL;
		end
		% returns test.stimcache.vname, char string identifying
		% variable(s) for curve (similar to test.Name
		function val = get.varied_parameter(obj)
			val = char(obj.Dinf.test.stimcache.vname);
		end
		% returns test.stimcache.vrange, values of varied parameter
		function val = get.varied_values(obj)
			val = obj.Dinf.test.stimcache.vrange;
		end
		% returns analysis_window : audio Delay to Delay+Duration interval
		function val = get.analysis_window(obj)
			val = [obj.Dinf.audio.Delay (obj.Dinf.audio.Delay + obj.Dinf.audio.Duration)];
		end
	end	% END methods
	
end	% END classdef
	
		