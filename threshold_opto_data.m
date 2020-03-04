function varargout = threshold_opto_data(D, Dinf, tracesByStim, varargin)
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
%   THRESHOLD				RMS spike threshold (# RMS, default: 3)
%	 SHOW_DEFAULTS			show default values for options
% 
% 	Outputs:
% 
% 		D			raw data, cell
% 		Dinf		data information struct
% 		S			spike data, struct
% 		T			data traces, sorted by stimulus, cell
% 		P			plot options, struct
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
SpikeWindow

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
				Forder = = varargin{argIndx + 1};
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

			case 'SHOW_DEFAULTS'
				% display default values for settings & options
				fprintf('%s: Default values:\n', mfilename)
				fprintf('\FILTER_DATA: %d\n', FilterData);
				fprintf('\tHPFREQ: %d\n', HPFreq);
				fprintf('\tLPFREQ: %d\n', LPFreq);
				fprintf('\tFORDER: %d\n', Forder);
				fprintf('\tTHRESHOLD: %d\n', Threshold);
				fprintf('\tCHANNEL: %d\n', channelNumber);
				fprintf('\tBINSIZE: %d\n', binSize);
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
if strcmpi(cInfo.testtype, 'WavFile')
	nstim = cInfo.Dinf.test.nCombinations;
	nreps = cInfo.Dinf.test.Reps;
else
	nstim = cInfo.ntrials;
	nreps = cInfo.nreps;
end
% allocate matrices
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

%---------------------------------------------------------------------
% Some test-specific things...
%---------------------------------------------------------------------
switch upper(cInfo.testtype)
	case 'FREQ'
		% list of frequencies, and # of freqs tested
		varlist = cInfo.varied_values;
		nvars = length(varlist);

	case 'LEVEL'
		% list of levels, and # of levels tested
		varlist = cInfo.varied_values;
		nvars = length(varlist);

	case 'FREQ+LEVEL'
		% list of freq, levels
		varlist = cell(2, 1);
		% # of freqs in nvars(1), # of levels in nvars(2)
		nvars = zeros(2, 1);
		tmprange = cInfo.varied_values;
		for v = 1:2
			varlist{v} = unique(tmprange(v, :), 'sorted');
			nvars(v) = length(varlist{v});
		end

	case 'OPTO'
		% not yet implemented
		
	case 'WAVFILE'
		% get list of stimuli (wav file names)
		varlist = cInfo.Dinf.test.wavlist;
		nvars = length(varlist);

	otherwise
		error('%s: unsupported test type %s', mfilename, cInfo.testtype);
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
else
	BPfilt = [];
end

sep_print('Filtering data...');
[nr, nc] = size(tracesByStim);
for r = 1:nr
	for c = 1:nc
		for 
		
	end
end

%---------------------------------------------------------------------
% find spikes!
%---------------------------------------------------------------------


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
									'varlist', {varlist}		);


end
