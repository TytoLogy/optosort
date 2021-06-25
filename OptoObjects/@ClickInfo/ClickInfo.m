classdef ClickInfo < CurveInfo
%------------------------------------------------------------------------
% Class: ClickInfo < CurveInfo
%------------------------------------------------------------------------
% TytoLogy:Experiments:optosort
%------------------------------------------------------------------------
% implements and encapsulates some utilities for dealing with
% click train data -> Dinf struct
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
		function obj = ClickInfo(varargin)
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
			% for click test, find indices of stimuli with same opto
			% amplitude
				case 'LEVEL'
    				fprintf('\t%s test, finding indices\n', obj.testtype);

               % from opto_build_clickstimList.m :
               % stim properties are stored in stimList for newer 
               % data (post 23 Jun 2021!)
               
               % for data with empty stimList, need to take a different
               % approach
               
					% list of amplitudes by sweep, amplitudes and 
               % # of amplitudes tested
					levellist = obj.getLevels;
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

      function levellist = getLevels(obj)
         % get varied values
         levels = obj.varied_values;
         % check if NULL stimulus was played
         if obj.Dinf.test.NullStim
            % since null stim was played, ASSUME that 
            % stimindex == 1 is null,
            % stimindex == 2 is click
            % create list of levels by sweep
            levellist = zeros(size(obj.Dinf.test.stimIndices));
            % set levels for stimIndices == 1(null) to 0;
            levellist(obj.Dinf.test.stimIndices == 1) = levels(1);
            % set levels for stimIndices == 2 (click) to level
            levellist(obj.Dinf.test.stimIndices == 2) = levels(2);
         else
            % since no null stim was played, ASSUME that 
            % stimindex == 1 is click
            % create list of levels by sweep
            levellist = zeros(size(obj.Dinf.test.stimIndices));
            % set levels for stimIndices == 1 (click) to levels(1);
            levellist(obj.Dinf.test.stimIndices == 1) = levels(1);
         end
      end
      
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
  				case {'CLICK'}
					varlist = {obj.varied_values};
					nvars = length(varlist);

				otherwise
					error('%s: unsupported test type %s', ...
								'ClickInfo.varlist', ...
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
			% get stimulus onset bins
			% to align to appended/merged file, will need to add
			% SpikeInfo.fileStartBin(findx) - 1
			onsetbins = obj.stimStartBin;

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
				case {'LEVEL'}
               % make sure we have click test
               if ~strcmpi(obj.testname, 'CLICK')
                  error('ClickInfo: unsupported testname %s', ...
                           obj.testname);
               end
               
               formatstr = '%s_%dmV';
               for n = 1:nevents
                  events(n).name = sprintf(formatstr, obj.testtype, ...
                                                varied_values(n));
                  events(n).samples = onsetbins(stimindex{n});
               end
                  
				otherwise
					error('ClickInfo: unsupported testtype %s (LEVEL)', ...
                       obj.testtype);
			end
		end
		%-------------------------------------------------
		
		%-------------------------------------------------
		%-------------------------------------------------
		% shortcut methods to stimcache values
		%-------------------------------------------------
		%-------------------------------------------------
%{
		ci.Dinf.test

  struct with fields:

             Type: 'LEVEL'
             Name: 'CLICK'
       ScriptType: [83 84 65 78 68 65 76 79 78 69]
             Reps: 20
        Randomize: 1
            Block: 0
         saveStim: 0
            Level: 70
         NullStim: 1
        NoiseStim: 0
      AcqDuration: 250
      SweepPeriod: 255
      stimIndices: [40×1 double]
    nCombinations: 2
     optovar_name: [65 109 112]
          optovar: 0
    audiovar_name: [67 108 105 99 107 76 101 118 101 108]
         audiovar: [76 101 118 101 108]
        curvetype: [76 69 86 69 76 43 79 112 116 111 79 70 70]
      %}

		% returns test.stimcache.vname, char string identifying
		% variable(s) for curve (similar to test.Name)
		function val = varied_parameter(obj)
			if obj.has_stimcache
				val = char(obj.Dinf.test.stimcache.vname);
         else
            % might not need char() conversion, but do it anyway just in
            % case...
				val = char(obj.Dinf.test.Name);
			end
		end
		% returns test.stimcache.vrange, values of varied parameter
		function val = varied_values(obj)
			if obj.has_stimcache
				val = obj.Dinf.test.stimcache.vrange;
         else
            % return level(s) of stimulus
            val = obj.Dinf.test.Level;
            % if null stimulus was played, add it as level == 0 
            if obj.Dinf.test.NullStim
               val = [0 val];
            end
			end
		end      
		% returns test.stimcache.nreps: # of reps for each stimulus
		function val = nreps(obj)
			if obj.has_stimcache
				val = obj.Dinf.test.stimcache.nreps;
			else
				val = obj.Dinf.test.Reps;
			end
      end
		% returns test.stimcache.ntrials: # of stimulus types
		function val = ntrials(obj)
			if obj.has_stimcache
				val = obj.Dinf.test.stimcache.ntrials;
         else
            val = 1;
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
	
		