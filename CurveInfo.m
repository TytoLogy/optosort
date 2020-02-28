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
%------------------------------------------------------------------------
% TO DO:
%------------------------------------------------------------------------

	%-------------------------------------------------
	% class properties
	%-------------------------------------------------
	properties
		Dinf
	end	% END properties (main)
	properties (Dependent)
		testtype
		testname
		freqs_bysweep
		levels_bysweep
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
					nfreqs = length(obj.Dinf.test.stimcache.vrange);
					% locate where trials for each frequency are located in the
					% stimulus cache list - this will be used to pull out trials of
					% same frequency
					stimindex = cell(nfreqs, 1);
					for f = 1:nfreqs
						stimindex{f} = find(obj.Dinf.test.stimcache.vrange(f) ...
																					== freqlist);
					end
					% assign outputs
					varargout{1} = stimindex;
					varargout{2} = freqlist;

			% for LEVEL test, find indices of stimuli with same level (dB SPL)
				case 'LEVEL'
					fprintf('\t%s test, finding indices\n', obj.testtype);
					% list of levels, and # of levels tested
					levellist = obj.levels_bysweep;
					nlevels = length(obj.Dinf.test.stimcache.vrange);
					% locate where trials for each frequency are located in the
					% stimulus cache list - this will be used to pull out trials of
					% same frequency
					stimindex = cell(nlevels, 1);
					for l = 1:nlevels
						stimindex{l} = find(obj.Dinf.test.stimcache.vrange(l) ...
																					== levellist);
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
		
	end	% END methods
	
end	% END classdef
	
		