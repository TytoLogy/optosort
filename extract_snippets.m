%function snips = extract_snippets(timestamps, sweeps)

% timestamps is a nreps X 1 cell array with each row containing a list
% (vector) of timestamps, in milliseconds
%
% sweeps is an [nsamples X nreps] matrix of recorded data from electrode

%% load test data
load snipindata.mat
timestamps = spiketimes{1};
sweeps = tracesByStim{1};

[tsrows, tscols] = size(timestamps);
[sweeprows, sweepcols] = size(sweeps);

% check to make sure sweeps is in proper orientation
if sweeprows < sweepcols
	sweeps = sweeps'

%%
snippets = cell(nr, nc);
for r = 1:nr
	for c = 1:nc
		for n = 1:cInfo.nreps
				tracesByStim{r, c}(:, n) = filtfilt(BPfilt.b, BPfilt.a, ...
															tracesByStim{r, c}(:, n));

		end
	end
end	