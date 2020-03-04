%------------------------------------------------------------------------
% Definitions
%------------------------------------------------------------------------
% string used to separate text 
sepstr = '----------------------------------------------------';

% data locations
rawPath = '~/Work/Data/TestData/MT';
rawFile = '1382_20191212_02_02_3200_FREQ_TUNING.dat';
% Channels = [4, 5, 11, 14];

%------------------------------------------------------------------------
% get filtered data
%------------------------------------------------------------------------
% probably don't need full raw data in D....
[D, Dinf, tracesByStim] = getFilteredOptoData( ...
											fullfile(rawPath, rawFile), ...
											'Filter', [300 5000], ...
											'Channel', 4);

%------------------------------------------------------------------------
%% threshold data
%------------------------------------------------------------------------

sep_print('threshold_opto_data...');
[spiketimes, cI, threshSettings] = threshold_opto_data(D, Dinf, tracesByStim);


%% test snippetizing


%% load test data
load snipindata.mat
Fs = cI.Dinf.indev.Fs;
timestamps = spiketimes{1};
sweeps = tracesByStim{1};

[tsrows, tscols] = size(timestamps);
[sweeprows, sweepcols] = size(sweeps);

% check to make sure sweeps is in proper orientation
if sweeprows > sweepcols
	sweeps = sweeps';
	[sweeprows, sweepcols] = size(sweeps);
end

if sweeprows ~= tsrows
	error('%s: mismatch in number of reps between timestamps and sweeps', ...
																					mfilename);
end

% convert snippet window to bins
SpikeWindow = [1 2];
snippetwindow = ms2bin(SpikeWindow, Fs)
snippetbins = sum(snippetwindow)

% allocate snippet data
snippets = cell(sweeprows, 1);

for r = 1:sweeprows
	% check # of timestamps for this sweep
	nspikes = length(timestamps{r});
	if nspikes > 0
		% temporary storage for snippets - will convert to matrix later
		tmpsnip = cell(nspikes, 1);
		% loop through spiketimes
		for s = 1:nspikes
			spikebin = ms2bin(timestamps{r}(s), Fs);
			% need to make sure snippetbins are in bounds
			if spikebin < snippetbins(1)
				% if spike occurred too close to start of sweep, left-pad the
				% snippet with zeros
				lpad = zeros(1, snippetbins(1) - spikebin);
			else
				lpad = [];
			end
			if spikebin > (length(sweep(r, :)) + snippetbins(2))
				% if spike occurred too close to end of sweep, right-pad 
				% with zeros
				rpad = zeros(1, spikebin - snippetbins(2));
			else
				rpad = [];
			end
			% apply pads
			tmpsweep = [lpad sweeps(r, :) rpad];
			% grap the snippet for this timestamp
			tmpsnip{s} = tmpsweep(
		end
	else
		% no spikes, no snippets!
		snippets{r, 1} = [];
	end
end	