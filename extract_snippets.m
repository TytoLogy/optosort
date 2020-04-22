function snippets = extract_snippets(timestamps, sweeps, snipwin, Fs)
%------------------------------------------------------------------------
% snippets = extract_snippets(timestamps, sweeps, snipwin, Fs)
%------------------------------------------------------------------------
% TytoLogy:Experiments:optoproc
%--------------------------------------------------------------------------
% given arrays of timestamps (in milliseconds) and sweeps, size of
% snippet window (in milliseconds) and sampling rate, extract_snippets will
% return a cell array of snippets from the sweep data
%
% timestamps is a nreps X 1 cell array with each row containing a list
% (vector) of timestamps, in milliseconds
%
% sweeps is an [nreps, nsamples] or [nsamples,  nreps]matrix of recorded
% data from electrode
%------------------------------------------------------------------------
%  Input Options:
% 		timestamps		{ntrials, 1} cell vector with each cell holding a vector
% 							of spike timestamps (milliseconds)
% 		sweeps			[ntrials, nsamples] matrix holding the recording traces 
% 							from each stimulus presentation.
% 		snipwin			[pre-timestamp time, post-timestamp time] 2 element vector
% 							providing the length of the snippet (in milliseconds) 
% 							before and after the timestamp
% 		Fs					Sampling rate (samples/second) used for the sweeps data
% 
% 	Outputs:
%		snippets		{ntrials} cell vector with each element holding a matrix of
%						snippets [ntimestamps for trial, nsamples]
%------------------------------------------------------------------------
% See Also: opto, plotFRA, plotRLF, plotFTC
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%  Sharad Shanbhag
%   sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 5 March 2020 (SJS)
%
% Revisions:
%--------------------------------------------------------------------------



%% check sizes of input arrays
[tsrows, ~] = size(timestamps);
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
snipbins = ms2bin(snipwin, Fs);
% create left, right padding vectors
lpad = zeros(1, snipbins(1));
rpad = zeros(1, snipbins(2));

% allocate storage for snippets
snippets = cell(sweeprows, 1);

% loop through traces (aka sweeps)
for r = 1:sweeprows
	% pad the left and right of the sweep - this is to deal with the
	% possibility of a timestamp being too close to the start or end of the
	% actual data sweep
	padsweep = [lpad sweeps(r, :) rpad];

	% cell array to hold snippets for this sweep prior to matrix conversion
	tmpS = cell(length(timestamps{r}));
	
	% loop through timestamps
	for t = 1:length(timestamps{r})
		ts = timestamps{r}(t);
		tsbin = ms2bin(ts, Fs);
	
		% get snippet start and end bins - need to account for the shift due to
		% lpad!
		% this is the explicit routine: 
		% for the start bin, take original bin, apply the lpad offset
		% stored in snipbins(1) and then subtract the lpad bins to capture the
		% desired "pre-timestamp" waveform.
		% for the end bin, shift the timestamp bin by lpad length (snipbins(1))
		% and then add the post-timestamp bins (snipbins(2))
		% 	
		% 	snip_startbin = ((tsbin + snipbins(1)) - snipbins(1))
		% 	snip_endbin = ((tsbin + snipbins(1)) + snipbins(2))
		% 	
		% but we can simply use the original tsbin, since the padded sweep is
		% already shifted by snipbins(1)!!!
		tmpS{t} = padsweep(tsbin:((tsbin + snipbins(1)) + snipbins(2)));

	end
	snippets{r} = cell2mat(tmpS);
end	