function [obj, varargout] = buildStimOnOffData(obj)
%------------------------------------------------------------------------
% [obj, onsetI, offsetI] = CurveInfo.buildStimOnOffData()
%------------------------------------------------------------------------
% TytoLogy:optosort:CurveInfo Object method
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
% Created: 9 April, 2020 using some code from buildChannelData (SJS)
%
% Revisions:
% 	10 Apr 2020 (SJS):
%		Roadmap:
%		 work here on more simple approach for BBN, freq tuning, etc. Then
%		 incorporate as method in CurveInfo object. Then override for use with WAV
%		 data in WAVInfo object
%------------------------------------------------------------------------
% TO DO:
%	include option to have onset include/exclude stim delay???
%
%------------------------------------------------------------------------

%{
%% read in test data
load('~/Work/Data/TestData/MT/1382_20191212_02_02_3200_FREQWAVMERGE_nexinfo.mat');

obj = nexInfo.FileData{2};
%}

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



%------------------------------------------------------------------------
%% process data
%------------------------------------------------------------------------

% initialize onsetIs and offsetI to store stop/start locations
onsetI = zeros(1, obj.Dinf.nread);
offsetI = zeros(1, obj.Dinf.nread);

% get stimulus delay and duration
delay_ms = obj.Dinf.audio.Delay;
duration_ms =  obj.Dinf.audio.Duration;

% convert to samples - need to use spike sampling rate (ADFs)
delay = ms2bin(delay_ms, obj.ADFs);
duration = ms2bin(duration_ms, obj.ADFs);

% loop through each sweep
for s = 1:obj.Dinf.nread
	% build list of stimulus onset and offset indices 
	% **** referenced to start of each sweep ********
	% stim onset will be startSweepBin + stimulus delay
	onsetI(s) = obj.startSweepBin{1}(s) + delay;
	% stim offset will onset bin + duration
	offsetI(s) = onsetI(s) + duration;
end
% assign to obj.stimStartBin and obj.stimEndBin
obj.stimStartBin = onsetI;
obj.stimEndBin = offsetI;


if nargout > 1
	varargout{1} = onsetI;
	varargout{2} = offsetI;
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
