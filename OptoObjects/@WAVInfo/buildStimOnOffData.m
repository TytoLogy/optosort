% function [obj, varargout] = buildStimOnOffData(obj)
%------------------------------------------------------------------------
% [obj, onsetI, offsetI] = WAVInfo.buildStimOnOffData()
%------------------------------------------------------------------------
% TytoLogy:optosort:WAVInfo Object method
%------------------------------------------------------------------------
% 
% builds vector of stimulus onset and offset bins referenced to start
% of file
% 
%------------------------------------------------------------------------
% Input Arguments:
% 	Dinf			data information struct from opto program output file
%
% Output Arguments:
%	obj		updated copy of object
% 	onsetI	[# sweeps, 1] array with each element stimulus onset timestamp 
%					for each sweep
% 	offsetI	[# sweeps, 1] array with each element stimulus onset timestamp 
%					for each sweep
%------------------------------------------------------------------------
% See also: CurveInfo, WAVInfo, buildChannelData, export_for_plexon
%------------------------------------------------------------------------


%------------------------------------------------------------------------
% Sharad J. Shanbhag
% sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 15 April, 2020 from CurveInfo.buildStimOnOffData (SJS)
%
% Revisions:
%------------------------------------------------------------------------
% TO DO:
%	include option to have onset include/exclude stim delay???
%
%------------------------------------------------------------------------


%% read in test data
load('a.mat');

obj = a;



%------------------------------------------------------------------------
%% settings
%------------------------------------------------------------------------
% check to make sure startSweepBin is defined
if isempty(obj.startSweepBin)
	error('%s: startSweepBin is empty array', mfilename);
end

%------------------------------------------------------------------------
% process data: get onset and offset of each stimulus (will need to be 
% slightly modified for WAV stimuli - do this in WAVInfo subclass of
% CurveInfo
% assume consistent across channels!!!!!
%------------------------------------------------------------------------
%{
idea for algorithm:

	loop through trials

	get stimulus onset, offset info:
		onset = stimulus delay
		offset = stimulus delay + stimulus duration

	wavOnsetBin will be defined as for stimStartBin in CurveInfo


	use indev.Fs to convert to samples from milliseconds

	To reference to absolute "time"/samples, add CurveInfo.startSweepBin(n) to
	n'th onset and offset

some caveats:

	will be different synthesized (tone, bbn) vs. wav stimuli, so will need to
	handle those different cases - look at optoproc_plotPSTH_WAVbyLevel to get ideas!

10 Apr 2020 (SJS):
Roadmap:
 work here on more simple approach for BBN, freq tuning, etc. Then
 incorporate as method in CurveInfo object. Then override for use with WAV
 data in WAVInfo object
%}


% quick check on stimindices and nread
if obj.Dinf.nread ~= length(obj.Dinf.test.stimIndices)
	tmpstr = sprintf('Dinf.nread (%d) ~= length(Dinf.test.stimIndices (%d)', ...
							obj.Dinf.nread, obj.Dinf.test.stimIndices);
	error('%s: %s', mfilename, tmpstr);
end

% information about the audio stimulus for each sweep is in
% stimList.audio struct
% stimList(s).audio:
% 	ISI: 100
% 	Delay: 55
% 	Duration: 100
% 	Level: 60
% 	Ramp: 5
% 	Frozen: 0
% 	signal: [1x1 struct] - will be either 'null', 'noise' or 'wav'
% 				Type: 'wav'
% 				WavPath: 'C:\TytoLogy\Experiments\WAVs'
% 				WavFile: 'Flat_adj.wav'

% quick check on stimindices and nread
if obj.Dinf.nread ~= length(obj.Dinf.test.stimIndices)
	tmpstr = sprintf('Dinf.nread (%d) ~= length(Dinf.test.stimIndices (%d)', ...
							obj.Dinf.nread, obj.Dinf.test.stimIndices);
	error('%s: %s', mfilename, tempstr);
end

% allocate storage for times
onsetTime = zeros(1, obj.Dinf.nread);
offsetTime = zeros(1, obj.Dinf.nread);
% initialize onsetIs and offsetI to store stop/start locations
onsetI = zeros(1, obj.Dinf.nread);
offsetI = zeros(1, obj.Dinf.nread);
% wavonset will hold start of wav stimulus playback (i.e. same as stimonset
% in CurveInfo...
wavonsetI = zeros(1, obj.Dinf.nread);

% loop through each sweep
for s = 1:obj.Dinf.nread
	% get the index into stimList from stimIndices vector in Dinf.test
	stim_index = obj.Dinf.test.stimIndices(s);
	% then get the stimulus out of Dinf.stimList
	stim = obj.Dinf.stimList(stim_index);
	
	% what happens next will depend on stimulus type ('null', 'noise' or
	% 'wav')
	fprintf('stimulus %d: %s\n', s, stim.audio.signal.Type);
	switch upper(stim.audio.signal.Type)
		case {'NULL', 'NOISE'}
			% if NULL or NOISE stim, delay and duration are just plain old
			% values
			delay_ms = stim.audio.Delay;
			duration_ms = stim.audio.Duration;
			% onset Time is just the delay time (like in CurveInfo)
			onsetTime(s) = delay_ms;
			offsetTime(s) = delay_ms + duration_ms;
			
		case 'WAV'
			% for wav stim, delay will be stim.audio.Delay plus the
			% onset time in the wav file itself
			% to compute this, we need the data from 
			% wavInfo about the wav stimulus...
			wavdata = lookup_wav(stim.audio.signal.WavFile, obj.wavInfo.wavs);
			% convert onset and offset bins to time in ms
			wavonset_ms = round(bin2ms(wavdata.OnsetBin, wavdata.SampleRate));
			wavoffset_ms = bin2ms(wavdata.OffsetBin, wavdata.SampleRate);
			% onset time is overall stimulus delay + wavonset_ms
			onsetTime(s) = stim.audio.Delay + wavonset_ms;
			offsetTime(s) = stim.audio.Delay + wavoffset_ms;
	end
	
	% now, we can reference the onset and offset times to the sweep times
	% onset will be startSweepBin + onsetTime sample (in ADfs!)
	onsetI(s) = obj.startSweepBin{1}(s) + ms2bin(onsetTime(s), obj.ADFs);
	% stim offset will onset bin + duration
	offsetI(s) = onsetI(s) + ms2bin(offsetTime(s), obj.ADFs);
	% waveonset is simply startSweepBin + delay.
	wavonsetI(s) = obj.startSweepBin{1}(s) + ms2bin(stim.audio.Delay, obj.ADFs);
end



% assign to obj.stimStartBin and obj.stimEndBin
obj.stimStartBin = onsetI;
obj.stimEndBin = offsetI;


if nargout > 1
	varargout{1} = onsetI;
	varargout{2} = offsetI;
	varargout{3} = sweepLen;
end


%%









%------------------------------------------------------------------------
%% process data
%------------------------------------------------------------------------


% get stimulus delay - but do a check for consistency
if obj.Dinf.audio.Delay ~= obj.wavInfo.audio.Delay
	warning('%s: unequal Dinf.audio.Delay and obj.wavInfo.audio.Delay', ...
							mfilename);
end
delay_ms = obj.Dinf.audio.Delay;
% convert to samples - need to use spike sampling rate (ADFs)
delay = ms2bin(delay_ms, obj.ADFs);



% get wav onset times
wav_onsetTimes_ms = bin2ms([obj.wavInfo.wavs.OnsetBin], ...
														obj.wavInfo.wavs(1).SampleRate);
% and offset bins
wav_offsetTimes_ms = bin2ms([obj.wavInfo.wavs.OffsetBin], ...
														obj.wavInfo.wavs(1).SampleRate);
% loop through each sweep
for s = 1:obj.Dinf.nread
	
	% need to get stimulus index from stimIndices
	stim_index = obj.Dinf.test.stimIndices(s);
	
	
	% build list of stimulus onset and offset indices 
	% **** referenced to start of each sweep ********
	% wav onset will be startSweepBin + stimulus delay
	onsetI(s) = obj.startSweepBin{1}(s) + delay;
	% stim offset will onset bin + duration
	offsetI(s) = onsetI(s) + duration;
	
	% waveonset is simple startSweepBin + delay.
	wavonsetI(s) = obj.startSweepBin{1}(s) + delay;
end


%{

%% test

plot(obj.startSweepBin{1}, ones(obj.Dinf.nread, 1), 'g.')
hold on
plot(obj.endSweepBin{1}, ones(obj.Dinf.nread, 1), 'r.')
plot(obj.stimStartBin, ones(obj.Dinf.nread, 1), 'c*')
plot(obj.stimEndBin, ones(obj.Dinf.nread, 1), 'k*')

hold off
%}
