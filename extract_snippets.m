%function snips = extract_snippets(snipwin, timestamps, sweeps, Fs)

% timestamps is a nreps X 1 cell array with each row containing a list
% (vector) of timestamps, in milliseconds
%
% sweeps is an [nsamples X nreps] matrix of recorded data from electrode

%% load test data
load snipindata.mat
timestamps = spiketimes{1};
sweeps = tracesByStim{1};
Fs = cInfo.Dinf.indev.Fs;

[tsrows, tscols] = size(timestamps);
[sweeprows, sweepcols] = size(sweeps);

% check to make sure sweeps is in proper orientation (samples in columns,
% trials/sweeps in rows)
if sweeprows > sweepcols
	sweeps = sweeps';
	[sweeprows, sweepcols] = size(sweeps);
end
% make sure # of rows match in the two arrays
if sweeprows ~= tsrows
	error('%s: mismatch in # of rows in sweeps and timestamps', mfilename);
end
% pre-ts post-ts window
snipwin = [1 2];
snipbin = ms2bin(snipwin, Fs)
lpad = zeros(1, snipbin(1));
rpad = zeros(1, snipbin(2));
ctrbin = snipbin(1) + 1;

%% test algorithm on one sweep
snippets = cell(sweeprows, 1);
for r = sweeprows
	ts = timestamps{sweeprows}(1)
	tsbin = ms2bin(ts, Fs)
	plot(1:sweepcols, sweeps(r, :));
	grid on
	hold on
		plot(tsbin, sweeps(r, tsbin), 'or')
	hold off
	drawnow
	
	% idea 1: by default, pad the left and right of the sweep
	padsweep = [lpad sweeps(r, :) rpad];
	hold on
		plot([1:length(padsweep)]-snipbin(1), padsweep, 'g.')
	hold off
	axis tight
end	