function snippets = extract_snippets(timestamps, sweeps, snipwin, Fs)

% timestamps is a nreps X 1 cell array with each row containing a list
% (vector) of timestamps, in milliseconds
%
% sweeps is an [nsamples X nreps] matrix of recorded data from electrode

%% load test data
load snipindata.mat
timestamps = spiketimes{1};
sweeps = tracesByStim{1};
Fs = cInfo.Dinf.indev.Fs;

%% check sizes of input arrays
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

% allocate storage for snippets
snippets = cell(sweeprows, 1);
% loop through traces
for r = sweeprows
	% plot sweep
	plot(1:sweepcols, sweeps(r, :), 'c.-');
	grid on
	
	% idea 1: by default, pad the left and right of the sweep
	padsweep = [lpad sweeps(r, :) rpad];
	% plot the padded data - need to offset the xaxis by length of lpad
	hold on
		plot( (1:length(padsweep))-snipbin(1), padsweep, 'k.')
	hold off
	axis tight

	% loop through timestamps
	for t = 1:length(timestamps{r})
		ts = timestamps{r}(t);
		tsbin = ms2bin(ts, Fs);

		% highlight timestamp
		hold on
			plot(tsbin, sweeps(r, tsbin), 'or')
		hold off
	
		% get snippet start and end bins - need to account for the shift due to
		% lpad!
		% this is the explicit routine: 
		% for the start bin, take original bin, apply the lpad offset
		% stored in snipbin(1) and then subtract the lpad bins to capture the
		% desired "pre-timestamp" waveform.
		% for the end bin, shift the timestamp bin by lpad length (snipbin(1))
		% and then add the post-timestamp bins (snipbin(2))
		% 	
		% 	snip_startbin = ((tsbin + snipbin(1)) - snipbin(1))
		% 	snip_endbin = ((tsbin + snipbin(1)) + snipbin(2))
		% 	
		% but we can simply use the original tsbin, since the padded sweep is
		% already shifted by snipbin(1)!!!
		snip_startbin = tsbin
		snip_endbin = ((tsbin + snipbin(1)) + snipbin(2))
		snip = padsweep( snip_startbin:snip_endbin);

		hold on
			plot((snip_startbin:snip_endbin) - snipbin(1), snip, 'g.:')
		hold off
		drawnow
		pause
	end
end	