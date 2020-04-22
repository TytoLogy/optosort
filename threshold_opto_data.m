function varargout = threshold_opto_data(cInfo, tracesByStim, varargin)
%------------------------------------------------------------------------
% [spikedata, CurveInfo, ThresholdInfo, tracesByStim] = 
% 															  threshold_opto_data(see help)
%------------------------------------------------------------------------
% TytoLogy:Experiments:OptoAnalysis
%--------------------------------------------------------------------------
% Processes data collected by the opto program
%
%------------------------------------------------------------------------
% Inputs:
%	cInfo				CurveData object
%	tracesByStim	cell array of recording sweeps organized by stimulus
%
%  Input Options:
%	 FILTER_DATA			if given, data will be bandpass filtered. default
%									is no filter applied
%   HPFREQ 					high pass filter cutoff frequency for neural data (Hz)
%   LPFREQ					low pass filter cutoff frequency for neural data (Hz)
%	 FORDER					filter order (typically 5 or lower)
%	 METHOD					default is RMS, optional: absolute (ABS)
%   THRESHOLD				RMS spike threshold (# RMS, default: 3)
%	 TSETTINGS				apply previously computed threshold settings
%								stored in ThresholdInfo struct (see outputs)
%	 SPIKE_WINDOW			[pre_ts post_ts] window for grabbing spike
%								waveform snippets, in milliseconds
%	 SHOW_DEFAULTS			show default values for options
% 
% 	Outputs:
% 	 spikes	struct with fields:
%					ts  cell array of spiketimes (in milliseconds)
% 					snips		cell array of spike waveform snippets
% 					tset	thresholding information struct:
% 							netrmsvals
% 							maxvals
% 							mean_rms
% 							global_max
% 							ThresholdMethod
% 							Threshold
% 							SpikeWindow
% 							BPfilt
% 	 tracesByStim	cell array of raw traces - usually only useful if 
% 						FILTER_DATA was requested!
%------------------------------------------------------------------------
% See Also: opto, plotFRA, plotRLF, plotFTC
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%  Sharad Shanbhag
%   sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 3 March 2020 (SJS) from optoproc
%
% Revisions:
%--------------------------------------------------------------------------


%---------------------------------------------------------------------
% settings for processing data
%---------------------------------------------------------------------
% filter data?
FilterData = false;
% filter
HPFreq = 300;
LPFreq = 4000;
Forder = 5;

ThresholdMethod = 'RMS';
% RMS spike threshold
Threshold = 3;
% ThresholdSettings
USER_THRESHOLD = false;
TSet = [];
% Spike Window [preDetectTime postDetectTime] (ms)
SpikeWindow = [1 2];

%---------------------------------------------------------------------
% Parse inputs
%---------------------------------------------------------------------
nvararg = length(varargin);
if nvararg > 0
	argIndx = 1;
	while argIndx <= nvararg
		fprintf('%s\n', upper(varargin{argIndx}))
		switch upper(varargin{argIndx})
			case 'FILTER_DATA'
				FilterData = true;
				argIndx = argIndx + 1;
			case 'HPFREQ'
				HPFreq = varargin{argIndx + 1};
				argIndx = argIndx + 2;
			case 'LPFREQ'
				LPFreq = varargin{argIndx + 1};
				argIndx = argIndx + 2;
			case 'FORDER'
				Forder = varargin{argIndx + 1};
				argIndx = argIndx + 2;
			case 'METHOD'
				tmp = varargin{argIndx + 1};
				if strcmpi(tmp, 'RMS')
					ThresholdMethod = 'RMS';
				elseif strcmpi(tmp, 'ABS')
					ThresholdMethod = 'ABSOLUTE';
				else
					error('%s: invalid threshold method: %s', mfilename, tmp);
				end
				argIndx = argIndx + 2;
			case 'THRESHOLD'
				tmp = varargin{argIndx + 1};
				if ischar(tmp)
					if strcmpi(tmp, 'DEFAULT')
						fprintf('%s: using default threshold: %d\n', ...
															mfilename, Threshold);
					else
						error('%s: unknown threshold command: %s', mfilename, tmp);
					end
				elseif isnumeric(tmp)
					Threshold = tmp;
				else
					error('%s: invalid threshold value: %s', mfilename, tmp)
				end
				argIndx = argIndx + 2;
			case 'TSETTINGS'
				tmp = varargin{argIndx + 1};
				if ~isstruct(tmp)
					error('%s: TSETTINGS input must be a struct', mfilename);
				else
					TSet = varargin{argIndx + 1};
					USER_THRESHOLD = true;
				end
				argIndx = argIndx + 2;
					
			case 'SPIKE_WINDOW'
				tmp = varargin{argIndx + 1};
				if ~isnumeric(tmp)
					error('%s: spike window must be 2 element number vector', ...
									mfilename);
				end
				SpikeWindow = tmp;
				argIndx = argIndx + 2;
				
			case 'SHOW_DEFAULTS'
				% display default values for settings & options
				fprintf('%s: Default values:\n', mfilename)
				fprintf('\FILTER_DATA: %d\n', FilterData);
				fprintf('\tHPFREQ: %d\n', HPFreq);
				fprintf('\tLPFREQ: %d\n', LPFreq);
				fprintf('\tFORDER: %d\n', Forder);
				fprintf('\tTHRESHOLD: %d\n', Threshold);
				fprintf('\tSPIKE_WINDOW: %d\n', SpikeWindow);
				return
			otherwise
				error('%s: unknown input arg %s', mfilename, varargin{argIndx});
		end
	end
end

% check threshold method
if strcmpi(ThresholdMethod, 'RMS')
	if Threshold < 0
		error('%s: RMS threshold value must be > 0! : %.4f', mfilename, Threshold)
	end
end

% %---------------------------------------------------------------------
% % create CurveInfo from Dinf
% %---------------------------------------------------------------------
% cInfo = CurveInfo(Dinf);

%---------------------------------------------------------------------
% determine global RMS and max - used for thresholding
%---------------------------------------------------------------------
% first, get  # of stimuli (called ntrials by opto) as well as # of reps

% use size of tracesByStim to get number of stimuli
% check reps
nstim = numel(tracesByStim);
[trRows, trCols] = size(tracesByStim);

reps_by_stim = zeros(trRows, trCols);
for r = 1:trRows
	for c = 1:trCols
		[reps_by_stim(r, c), ~] = size(tracesByStim{r, c}');
% 		fprintf('tracesByStim{%d, %d} has %d reps\n', r, c, reps_by_stim(r, c));
	end
end
if length(unique(reshape(reps_by_stim, numel(reps_by_stim), 1))) > 1
	warning('%s: unequal number of reps in tracesByStim', mfilename);
	nreps = max((reshape(reps_by_stim, numel(reps_by_stim), 1)));
	fprintf('Using max nreps for allocation: %d\n', nreps);
else
	nreps = reps_by_stim(1);
end

%---------------------------------------------------------------------
% get RMS values, calculate threshold if needed
%---------------------------------------------------------------------

if ~USER_THRESHOLD
	% get threshold settings
	netrmsvals = zeros(nstim, nreps);
	maxvals = zeros(nstim, nreps);
	% find rms, max vals for each stim
	for s = 1:nstim
		netrmsvals(s, :) = rms(tracesByStim{s});
		maxvals(s, :) = max(abs(tracesByStim{s}));
	end
	% compute overall mean rms for threshold
	fprintf('Calculating mean RMS for data...\n');
	mean_rms = mean(reshape(netrmsvals, numel(netrmsvals), 1));
	fprintf('\tMean rms: %.4f\n', mean_rms);
	% find global max value (will be used for plotting)
	global_max = max(max(maxvals));
	fprintf('\tGlobal max abs value: %.4f\n', global_max);
	
else
	% user-supplied threshold settings
	fprintf('Using pre-calculated threshold settings\n');
	fprintf('Updating netrmsvals, maxvals\n');
	netrmsvals = zeros(nstim, nreps);
	maxvals = zeros(nstim, nreps);
	% find rms, max vals for each stim
	for s = 1:nstim
		netrmsvals(s, :) = rms(tracesByStim{s});
		maxvals(s, :) = max(abs(tracesByStim{s}));
	end
	% use mean_rms, max
	fprintf('Using user-specified mean RMS and max val for data...\n');
	mean_rms = TSet.mean_rms;
	fprintf('\tMean rms: %.4f\n', mean_rms);
	% find global max value (will be used for plotting)
	global_max = TSet.global_max;
	fprintf('\tGlobal max abs value: %.4f\n', global_max);
	% using user specified Threshold
	ThresholdMethod = TSet.ThresholdMethod;
	Threshold = TSet.Threshold;
	% spike window
	SpikeWindow = TSet.SpikeWindow;
end

% calculate threshold - will depend on method
switch ThresholdMethod
	case 'RMS'
		Threshold = Threshold*mean_rms;
	case 'ABSOLUTE'
		Threshold = Threshold; %#ok<ASGSL>
end

%---------------------------------------------------------------------
% define filter for data
%---------------------------------------------------------------------
% local copy of sample rate
Fs = cInfo.Dinf.indev.Fs;
% build bandpass filter, store coefficients in BPfilt
if FilterData
	BPfilt.Fs = Fs;
	BPfilt.Fnyq = BPfilt.Fs / 2;
	BPfilt.Fc = [HPFreq LPFreq];
	BPfilt.cutoff = BPfilt.Fc / BPfilt.Fnyq;
	BPfilt.forder = Forder;
	[BPfilt.b, BPfilt.a] = butter(BPfilt.forder, BPfilt.cutoff, 'bandpass');
	sep_print('Filtering data...');
	[nr, nc] = size(tracesByStim);
	for r = 1:nr
		for c = 1:nc
			for n = 1:cInfo.nreps
					tracesByStim{r, c}(:, n) = filtfilt(BPfilt.b, BPfilt.a, ...
																tracesByStim{r, c}(:, n));
			end
		end
	end	
else
	BPfilt = [];
end

%---------------------------------------------------------------------
% find spikes!
%---------------------------------------------------------------------
% allocate spiketimes cell to match tracesByStim
spiketimes = cell(trRows, trCols);
% allocate snippets cell to match tracesByStim
snippets = cell(trRows, trCols);
% detect spikes
for r = 1:trRows
	for c = 1:trCols
		% use rms threshold to find spikes
		% need to pass the transpose of tracesByStim because, by default, the
		% samples go down rows and trials across cols, but spikeschmitt2
		% expects the opposite....
		spiketimes{r, c} = spikeschmitt2(tracesByStim{r, c}', ...
										Threshold, 1, Fs, 'ms');
									
		% locate spikes and get snippets of data
		snippets{r, c} = extract_snippets(spiketimes{r, c}, ...
														tracesByStim{r, c}, ...
														SpikeWindow, ...
														Fs);
	end
end

%---------------------------------------------------------------------
% outputs
%---------------------------------------------------------------------
if nargout
	% spiketimes struct that has spiketimes cell array,  snippets
	% cell array, tset (threshold settings struct)
	% **need to put field in curly braces to prevent struct() command from
	%   creating a struct array (which we don't want, in this case)
	varargout{1} = struct(	'ts', {spiketimes}, ...
									'snips', {snippets}, ...
									'tset', ...
										struct(	'netrmsvals', netrmsvals, ...
												'maxvals', maxvals, ...
												'mean_rms', mean_rms, ...
												'global_max', global_max, ...
												'ThresholdMethod', ThresholdMethod, ...
												'Threshold', Threshold, ...
												'SpikeWindow', SpikeWindow, ...
												'BPfilt', BPfilt		) ...
								);
	% return tracesByStim;
	if nargout > 1
		varargout{2} = tracesByStim;
	end
end

