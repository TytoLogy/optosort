function varargout = threshold_opto_data(Dinf, tracesByStim, varargin)
%------------------------------------------------------------------------
% [D, Dinf, S, T, P] = threshold_opto_data(see help)
%------------------------------------------------------------------------
% TytoLogy:Experiments:OptoAnalysis
%--------------------------------------------------------------------------
% Processes data collected by the opto program
%
%------------------------------------------------------------------------
%  Input Options:
%	 FILTER_DATA			if given, data will be bandpass filtered. default
%									is no filter applied
%   HPFREQ 					high pass filter cutoff frequency for neural data (Hz)
%   LPFREQ					low pass filter cutoff frequency for neural data (Hz)
%	 FORDER					filter order (typically 5 or lower)
%   THRESHOLD				RMS spike threshold (# RMS, default: 3)
%	 SPIKE_WINDOW			[pre_ts post_ts] window for grabbing spike
%								waveform snippets, in milliseconds
%	 FRAMETHOD				'OLD' uses original method for FRA data
%	 SHOW_DEFAULTS			show default values for options
% 
% 	Outputs:
% 
% 		spiketimes
% 		snippets
% 		CurveInfo		curve information object
% 		ThresholdInfo	thresholding information
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

% RMS spike threshold
Threshold = 3;
% Spike Window [preDetectTime postDetectTime] (ms)
SpikeWindow = [1 2];

FRAmethod = 'NEW';

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
					if tmp > 0
						Threshold = tmp;
					else
						error('%s: invalid threshold value: %.4f', mfilename, tmp)
					end
				else
					error('%s: invalid threshold value: %s', mfilename, tmp)
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

			case 'FRAMETHOD'
				tmp = varargin{argIndx + 1};
				if ischar(tmp)
					if strcmpi(tmp, 'OLD')
						fprintf('%s: using old FRA method threshold\n', mfilename);
						FRAmethod = 'OLD';
					elseif strcmpi(tmp, 'NEW')
						fprintf('%s: using new FRA method threshold\n', mfilename);
						FRAmethod = 'NEW';
					else
						error('%s: unknown FRAMETHOD option: %s', mfilename, tmp);
					end
				else
					error('%s: invalid FRAMETHOD option: %s', mfilename, tmp)
				end
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
				fprintf('\FRAMETHOD: %d\n', FRAmethod);
				return
			otherwise
				error('%s: unknown input arg %s', mfilename, varargin{argIndx});
		end
	end
end

%---------------------------------------------------------------------
% create CurveInfo from Dinf
%---------------------------------------------------------------------
cInfo = CurveInfo(Dinf);

%---------------------------------------------------------------------
% determine global RMS and max - used for thresholding
%---------------------------------------------------------------------
% first, get  # of stimuli (called ntrials by opto) as well as # of reps

%{ 
ORIG: uses test information from cInfo
if strcmpi(cInfo.testtype, 'WavFile')
	nstim = cInfo.Dinf.test.nCombinations;
	nreps = cInfo.Dinf.test.Reps;
else
	nstim = cInfo.ntrials;
	nreps = cInfo.nreps;
end
%}

% testing: use size of tracesByStim to get number of stimuli
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


%% allocate matrices
netrmsvals = zeros(nstim, nreps);
maxvals = zeros(nstim, nreps);
% find rms, max vals for each stim
for s = 1:nstim
	netrmsvals(s, :) = rms(tracesByStim{s});
	maxvals(s, :) = max(abs(tracesByStim{s}));
end
% compute overall mean rms for threshold
fprintf('Calculating mean and max RMS for data...\n');
mean_rms = mean(reshape(netrmsvals, numel(netrmsvals), 1));
fprintf('\tMean rms: %.4f\n', mean_rms);
% find global max value (will be used for plotting)
global_max = max(max(maxvals));
fprintf('\tGlobal max abs value: %.4f\n', global_max);

% get varlist
[varlist, nvars] = cInfo.varlist;

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

switch FRAmethod
	case 'OLD'
		% different approaches to storage depending on test type
		switch upper(cInfo.testtype)
			case 'FREQ+LEVEL'
				% for FRA data, nvars has values [nfreqs nlevels];
				spiketimes = cell(nvars(2), nvars(1));
				for v1 = 1:nvars(1)
					for v2 = 1:nvars(2)
						% use rms threshold to find spikes
						spiketimes{v2, v1} = ...
								spikeschmitt2(tracesByStim{v2, v1}', Threshold*mean_rms, ...
																					1, Fs, 'ms');
					end
				end
			otherwise
				% if test is not FREQ+LEVEL (FRA), nvars will be a single number
				spiketimes = cell(nvars, 1);
				for v = 1:nvars
					% use rms threshold to find spikes
					spiketimes{v} = spikeschmitt2(tracesByStim{v}', Threshold*mean_rms, ...
																						1, Fs, 'ms');
				end
		end
		
	case 'NEW'
		% allocate spiketimes cell to match tracesByStim
		spiketimes = cell(trRows, trCols);
		% detect spikes
		for r = 1:trRows
			for c = 1:trCols
				% use rms threshold to find spikes
				% need to pass the transpose of tracesByStim because, by default, the
				% samples go down rows and trials across cols, but spikeschmitt2
				% expects the opposite....
				spiketimes{r, c} = spikeschmitt2(tracesByStim{r, c}', ...
												Threshold*mean_rms, 1, Fs, 'ms');
			end
		end
end


%---------------------------------------------------------------------
% outputs
%---------------------------------------------------------------------
if nargout
	% raw traces
	varargout{1} = spiketimes;
	% data information
	varargout{2} = cInfo;
	% analysis output
	varargout{3} = struct(	'netrmsvals', netrmsvals, ...
									'maxvals', maxvals, ...
									'mean_rms', mean_rms, ...
									'global_max', global_max, ...
									'Threshold', Threshold, ...
									'BPfilt', BPfilt, ...
									'nvars', nvars, ...
									'varlist', varlist		);
	varargout{4} = tracesByStim;

end

% save snipindata.mat spiketimes cInfo tracesByStim 
