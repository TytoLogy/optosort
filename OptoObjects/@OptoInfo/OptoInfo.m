classdef OptoInfo < CurveInfo
%------------------------------------------------------------------------
% Class: OptoInfo < CurveInfo
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
% Created: 16 June, 2021 (SJS)
%	
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
		function obj = OptoInfo(varargin)
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
            error(['OptoInfo.getStimulusIndices: ' ...
                   'Dinf not defined/is empty']);
         end
			
         % look for levels, indices
			switch upper(obj.testtype)
			% for OPTO-AMP test, find indices of stimuli with same opto
			% amplitude
				case 'OPTO-AMP'
					fprintf('\t%s test, finding indices\n', obj.testtype);
					% list of amplitudes by sweep, amplitudes and 
               % # of amplitudes tested
					levellist = obj.optoAmp_bysweep;
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
				
			% unknown type
				otherwise
					error('%s: unsupported test type %s', mfilename, ...
								obj.testtype);
			end
		end	% END getStimulusIndices method
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

				case {'OPTO', 'OPTO-AMP'}
               % lumping opto and opto-amp together. this might
               % be stupid or brilliant...
					% list of levels, and # of levels tested
					varlist = obj.varied_values;
					nvars = length(varlist);
					titleString = cell(nvars, 1);
               for v = 1:nvars
                  if v == 1
                     titleString{v} = {   fname, ...
                                          sprintf('Light = %d mV', ...
                                          varlist(v))};
                  else
                     titleString{v} = sprintf('Light = %d mV', ...
                                                varlist(v));
                  end
               end
               					
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
  				case {'OPTO', 'OPTO-AMP'}
					warning('CurveInfo.varlist: OPTO not yet implemented');
					varlist = {obj.varied_values};
					nvars = length(varlist);

				otherwise
					error('%s: unsupported test type %s', ...
								'CurveInfo.varlist', ...
								cInfo.testtype);
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

%{
some notes:

opto-amp will run ntrials * nlevels 



%}
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
				case {'OPTO', 'OPTO-AMP'}
					formatstr = '%s_%dmV';
					for n = 1:nevents
						events(n).name = sprintf(formatstr, obj.testtype, ...
																varied_values(n));
						events(n).samples = onsetbins(stimindex{n});
					end

				otherwise
					error('OptoInfo: use subclass or unsupported %s', ...
                       obj.testtype);
			end
		end
		%-------------------------------------------------
		
		%-------------------------------------------------
		%-------------------------------------------------
		% shortcut methods to stimcache values
		%-------------------------------------------------
		%-------------------------------------------------
		
		%-------------------------------------------------
      % returns all values of test.stimcache.opto.amp
		%-------------------------------------------------
      function val = optoAmp_bysweep(obj)
         if obj.has_stimcache
            % opto is a cell array of structs... 
            %    why??? I don't remember...
            % so, will have to do a little machination
            % convert cell array of structs to an array 
            % of structs
            tmp = cell2mat(obj.Dinf.test.stimcache.opto);
            % then get all amplitude values
            val = [tmp.Amp]';
         else
            val = [];
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
      %{
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
      %}
	end	% END methods
	
end	% END classdef
	
		