function varargout = threshold_opto_data(varargin)
%------------------------------------------------------------------------
% [D, Dinf, S, T, P] = optoproc(see help)
%------------------------------------------------------------------------
% TytoLogy:Experiments:OptoAnalysis
%--------------------------------------------------------------------------
% Processes data collected by the opto program
%
%------------------------------------------------------------------------
%  Input Options:
%   DATAFILE or FILE		path and name of .dat file
%   PLOTPATH				path (directory) for output of plot files
%   HPFREQ					high pass filter cutoff frequency for neural data (Hz)
%   LPFREQ					low pass filter cutoff frequency for neural data (Hz)
%   THRESHOLD				RMS spike threshold (# RMS, default: 3)
%   CHANNEL					Input data channel (1-16, default is 8)
%   BINSIZE					binsize in ms for PSTH (default: 5 ms)
%   PLOT_TRACES or		draw plots of all individual traces (1 = yes, 2 = no)
%		PLOTTRACES
%   PLOT_PSTH or			draw plots of PSTHs (1 = yes, 2 = no) as individual
%		PLOTPSTH				figures 
%   PLOT_PSTH_BY_LEVEL	draw plots of PSTHs (1 = yes, 2 = no) grouped by
%		or PLOTPSTHBYLEVEL			stimulus level
%   PLOT_PSTH_MATRIX		draw plots of PSTHs (1 = yes, 2 = no) as 1 figure
%		or PLOTPSTHMAT			 
%   PLOT_RLF				plot Rate-Level Function
%	 PLOT_FTC				plot Frequency Tuning Curve
%   PLOT_FRA				plot FrequencyResponseArea data
%   SAVEFIG					save .fig plots (1 = yes, 2 = no)
%   SAVEPNG					save plots as .png files (1 = yes, 2 = no)
%   SAVEPDF					save plots as .pdf files (1 = yes, 2 = no)
%   TIMELIMITS				limit time range to specified interval
%									(e.g., [0 1000], in ms);
%   YLIMITS					limit y axis to specified interval (e.g., [-1 1])
%   PLOT_FILENAME			user-specified plot file base name (e.g., '1142_unit1')
%   PLOTROWCOLS			# of rows and columns for plots	
%									([2 3] is 2 rows, 3 cols)
%   EXPLORE					open optexplore app (not yet working 23 Jul 2019)
%   EXPORT_CHANNEL		channel(s) to export	for spike sorting
%									[8 10 11] will export data for channels 8, 10 11
%								exported data will be in separate .mat files for
%								each channel
%	 EXPORT_PATH			path (directory) for export of files for spike sorting
%	 EXPORT_MODE			default: 'wave_clus'
%								options:
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
datafile = '';


% filter
HPFreq = 300;
LPFreq = 4000;
% RMS spike threshold
% Threshold = 4.5;
Threshold = 3;
% Channel Number (default: use 8 for single channel data)
channelNumber = 8;
% limits for channels (low high)
CHANNEL_LIMITS = [1 16];


%---------------------------------------------------------------------
% Parse inputs
%---------------------------------------------------------------------
if nargin
	argIndx = 1;
	while argIndx <= nargin
		fprintf('%s\n', upper(varargin{argIndx}))
		switch upper(varargin{argIndx})
			case {'DATAFILE', 'FILE'}
				datafile = varargin{argIndx + 1};
				[datapath, dfile, dext] = fileparts(datafile);
				datafile = [dfile dext];
				argIndx = argIndx + 2;
			case 'HPFREQ'
				HPFreq = varargin{argIndx + 1};
				argIndx = argIndx + 2;
			case 'LPFREQ'
				LPFreq = varargin{argIndx + 1};
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
			case 'CHANNEL'
				channelNumber = varargin{argIndx + 1};
				argIndx = argIndx + 2;
			case 'BINSIZE'
				binSize = varargin{argIndx + 1};
				argIndx = argIndx + 2;
			case 'SHOW_DEFAULTS'
				% display default values for settings & options
				fprintf('%s: Default values:\n', mfilename)
				fprintf('\tDATAFILE: %s\n', datafile);
				fprintf('\tHPFREQ: %d\n', HPFreq);
				fprintf('\tLPFREQ: %d\n', LPFreq);
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
% need to get information about system
%---------------------------------------------------------------------
[data_root_path, tytology_root_path] = optoanalysis_paths; %#ok<ASGLU>

%---------------------------------------------------------------------
% data file things
%---------------------------------------------------------------------
if isempty(datafile)
	% get data file from user
	[datafile, datapath] = uigetfile('*.dat', 'Select opto data file', ...
														data_root_path);
	% abort if cancelled
	if datafile == 0
		fprintf('Cancelled\n');
		return
	end	
end

%---------------------------------------------------------------------
% Read Data
%---------------------------------------------------------------------
[D, Dinf, tracesByStim] = getFilteredOptoData( ...
											fullfile(datapath, datafile), ...
											'Filter', [HPFreq LPFreq], ...
											'Channel', channelNumber);
if isempty(D)
	warning('%s: D is empty???!!!??!!', mfilename);
	return
end

% create CurveInfo from Dinf
cInfo = CurveInfo(Dinf);
cInfo.F = OptoFileName(fullfile(datapath, datafile));


%---------------------------------------------------------------------
% determine global RMS and max - used for thresholding
%---------------------------------------------------------------------
% first, get  # of stimuli (called ntrials by opto) as well as # of reps
if strcmpi(cInfo.testtype, 'WavFile')
	nstim = cInfo.test.nCombinations;
	nreps = cInfo.test.Reps;
else
	nstim = cInfo.test.stimcache.ntrials;
	nreps = cInfo.test.stimcache.nreps;
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
switch upper(cInfo.testtyoe)
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
		titleString = cell(nvars, 1);
		for v = 1:nvars
			if v == 1 
				titleString{v} = {fname, sprintf('wav name: %s', varlist{v})};
			else
				titleString{v} = sprintf('wav name: %s', varlist{v});
			end
		end
	otherwise
		error('%s: unsupported test type %s', mfilename, Dinf.test.Type);
end

%---------------------------------------------------------------------
% find spikes!
%---------------------------------------------------------------------
% local copy of sample rate
Fs = cInfo.Dinf.indev.Fs;

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
% Plot raw data
%---------------------------------------------------------------------
if plotTraces
	% new figure
	hF = figure;
	% name figure
	set(hF, 'Name', [fname '_sweeps']);
	% loop through variable
	for v = 1:nvars
		% time vector for plotting
		t = (1000/Fs)*((1:length(tracesByStim{v}(:, 1))) - 1);
		subplot(prows, pcols, v);
		% flip tracesByStim in order to have sweeps match the raster plots
		stackplot(t, fliplr(tracesByStim{v}), 'colormode', 'black', ...
												'ymax', global_max, ...
												'figure', hF, 'axes', gca);
		title(titleString{v}, 'Interpreter', 'none');
		xlabel('ms')
		ylabel('Trial')
	end
	% save plot
	if any([saveFIG savePDF savePNG])
		pname = fullfile(plotpath, [plotFileName '_traces']); 
		if saveFIG
			savefig(hF, pname, 'compact');
		end
		if savePDF
			print(hF, pname, '-dpdf');
		end
		if savePNG
			print(hF, pname, '-dpng', '-r300');
		end
	end
end

%---------------------------------------------------------------------
% raster, psths as matrix
%---------------------------------------------------------------------
if plotPSTHMAT
	% compute range of time for x axis
	if isempty(timeLimits)
		% time vector for plotting
		t = (1000/Fs)*((1:length(tracesByStim{1}(:, 1))) - 1);
		timeLimits = [0 ceil(max(t))];
	end
	
	hPR = plotPSTHMATRIX(spiketimes, Dinf, binSize, nvars, varlist, ...
									[prows pcols], timeLimits, yLimits, titleString);
	% set plot name
	set(hPR, 'Name', plotFileName)
	% save plot
	if any([saveFIG savePDF savePNG])
		pname = fullfile(plotpath, [plotFileName '_rp']);
		if saveFIG
			savefig(hPR, pname, 'compact');
		end
		if savePDF
			print(hPR, pname, '-dpdf');
		end
		if savePNG
			print(hPR, pname, '-dpng', '-r300');
		end
	end
end

%---------------------------------------------------------------------
% raster, psths individually ***** debugging 23 May 2019 ****
%---------------------------------------------------------------------
if plotPSTH
	
	% compute range of time for x axis
	if isempty(timeLimits)
		% time vector for plotting
		t = (1000/Fs)*((1:length(tracesByStim{1}(:, 1))) - 1);
		timeLimits = [0 ceil(max(t))];
	end
	
	if strcmpi(Dinf.test.Type, 'WAVFILE')
		if plotPSTH_BY_LEVEL
			hPR = optoproc_plotPSTH_WAVbyLevel(spiketimes, Dinf, binSize, ...
									[prows pcols], timeLimits, yLimits);			
		else
			hPR = optoproc_plotPSTH_byWAV(spiketimes, Dinf, binSize, ...
									timeLimits, yLimits);
		end
	end
	% save plot
	if any([saveFIG savePDF savePNG])
		for f = 1:length(hPR)
			if isempty(get(hPR{f}, 'FileName'))
				[pname, pfolder] = uiputfile(fullfile(plotpath, plotFileName), ...
														'Save plot');
				pname = fullfile(pfolder, pname);
			else
				pname = fullfile(plotpath, get(hPR{f}, 'FileName'));
			end
			if saveFIG
				savefig(hPR, pname, 'compact');
			end
			if savePDF
				print(hPR, pname, '-dpdf');
			end
			if savePNG
				print(hPR, pname, '-dpng', '-r300');
			end
			end
	end
end

%---------------------------------------------------------------------
% Rate-Level function plot
%---------------------------------------------------------------------
if plotRateLevelFun && strcmpi(Dinf.test.Type,'FREQ+LEVEL')
	% not FRA data
	warning('optoproc: Test type (%s) is FREQ+LEVEL (FRA)!', Dinf.test.Type);
elseif plotRateLevelFun && any(strcmpi(Dinf.test.Type, {'LEVEL', 'BBN'}))
	% time window for counting spikes - use Delay to Delay+Duration interval
	analysisWindow = [Dinf.audio.Delay (Dinf.audio.Delay + Dinf.audio.Duration)];
	RLF = computeRLF(spiketimes, Dinf.audio.Level, analysisWindow);
	hRLF = plotCurveAndCI(RLF, 'median');
	% build title string
	tStr = {sprintf('RLF [%d %d] dB SPL', min(RLF.xdata), max(RLF.xdata)), ...
					[datafile ', ' sprintf('Channel %d', channelNumber)]};
	title(tStr, 'Interpreter', 'none');
	% set plot name
	set(hRLF, 'Name', plotFileName);
    % save plot
	if any([saveFIG savePDF savePNG])
		pname = fullfile(plotpath, [plotFileName '_RLF']);
		if saveFIG
			savefig(hRLF, pname, 'compact');
		end
		if savePDF	
			print(hRLF, pname, '-dpdf');
		end
		if savePNG
			print(hRLF, pname, '-dpng', '-r300');
		end
	end	
	Dinf.RLF = RLF;
end


%---------------------------------------------------------------------
% Frequency tuning curve plot
%---------------------------------------------------------------------
if plotFreqTuningCrv && strcmpi(Dinf.test.Type,'FREQ+LEVEL')
	% not useful for FRA data (at moment - this could be used to plot a
	% "slice" of the FRA!!!!
	warning('optoproc: Test type (%s) is FREQ+LEVEL (FRA)!', Dinf.test.Type);
elseif plotFreqTuningCrv && ~strcmpi(Dinf.test.Type, 'FREQ')
	% other non freq test
	warning('optoproc: Test type (%s) is not compatible with FTC plot!', ...
					Dinf.test.Type);
elseif plotFreqTuningCrv && strcmpi(Dinf.test.Type, 'FREQ')
	% time window for counting spikes - use Delay to Delay+Duration interval
	analysisWindow = [Dinf.audio.Delay (Dinf.audio.Delay + Dinf.audio.Duration)];
	FTC = computeFTC(spiketimes, Dinf.audio.signal.Frequency, analysisWindow);
	hFTC = plotCurveAndCI(FTC, 'median');
	% build title string
	tStr = {sprintf('FTC [%d %d] kHz, %d dB SPL', ...
								min(FTC.xdata), max(FTC.xdata), Dinf.audio.Level), ...
					[datafile ', ' sprintf('Channel %d', channelNumber)]};
	title(tStr, 'Interpreter', 'none');
    % set plot name
	set(hFTC, 'Name', plotFileName);
	% save plot
	if any([saveFIG savePDF savePNG])
		pname = fullfile(plotpath, [plotFileName '_FTC']);
		if saveFIG
			savefig(hFTC, pname, 'compact');
		end
		if savePDF
			print(hFTC, pname, '-dpdf');
		end
		if savePNG
			print(hFTC, pname, '-dpng', '-r300');
		end
	end	
	% save FTC in Dinf... not ideal, but for now ok
	Dinf.FTC = FTC;
end

%---------------------------------------------------------------------
% FRA plot
%---------------------------------------------------------------------
if plotFreqRespArea && ~strcmpi(Dinf.test.Type,'FREQ+LEVEL')
	% not FRA data
	warning('optoproc: Test type (%s) is not FREQ+LEVEL (FRA)!', Dinf.test.Type);
elseif plotFreqRespArea
	% they are FRA data
	% window for spike count
	frawin = [Dinf.audio.Delay (Dinf.audio.Delay + Dinf.audio.Duration)];
	% calculate FRA stored in struct FRA
	FRA = computeFRA(spiketimes, varlist{1}, varlist{2}, frawin);
	% set fname to data file name
	FRA.fname = datafile;
	hFRA = plotFRA(FRA, 'dB');
    % set plot name
	set(hFRA, 'Name', plotFileName);
% 	set(hFRA, 'FileName', plotFileName);
	% save plot
	if any([saveFIG savePDF savePNG])
		pname = fullfile(plotpath, [plotFileName '_FRA']);
		if saveFIG
			savefig(hFRA, pname, 'compact');
		end
		if savePDF
			print(hFRA, pname, '-dpdf');
		end
		if savePNG
			print(hFRA, pname, '-dpng', '-r300');
		end
	end
	Dinf.FRA = FRA;
end


%---------------------------------------------------------------------
% outputs
%---------------------------------------------------------------------
if nargout
	% raw traces
	varargout{1} = D;
	% data information
	varargout{2} = Dinf;
	% analysis output
	varargout{3} = struct(	'spiketimes', {spiketimes}, ...
									'netrmsvals', netrmsvals, ...
									'maxvals', maxvals, ...
									'mean_rms', mean_rms, ...
									'global_max', global_max, ...
									'Threshold', Threshold, ...
									'nvars', nvars, ...
									'varlist', {varlist}		);
	% raw traces by stimulus
	varargout{4} = tracesByStim;
	% plot options
	if plotPSTH
		if exist('plotopts', 'var')
			varargout{5} = plotopts;
		else
			varargout{5} = [];
		end
	else
		varargout{5} = [];
	end
end
