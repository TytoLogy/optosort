function varargout = getSpikesByStim(obj, fIndx, chanNum, unitNum)
% [spike_struct, info_struct] = getSpikesByStim(fIndx, chanNum, unitNum)
%------------------------------------------------------------------------
% TytoLogy:optosort:OptoObjects
%--------------------------------------------------------------------------
%	method: SpikeData.getSpikesByStim
%--------------------------------------------------------------------------
%
%	basically a front end for CurveInfo(and subclasses WAVInfo and FRAInfo)
%	that takes care of some checks and then calls
%	SpikeData.spikesForAnalysis to get a list of spikes by sweep (aka trial)
%	and then, for given file,
%	CurveInfo/WAVInfo/FRAInfo.convertSpikeTableToSpikeTimes
%	
%------------------------------------------------------------------------
% Input Arguments:
%   fIndx    index to SpikeData.FileInfo{} cell for desired test
% 	 chanNum  A/D channel for requested data
% 	 unitNum  unit id number (from sorting program) 
%
% Output Arguments:
% if only 1 output requested, spike_struct will be returned:
%  spike_struct	output struct. organization of elements within the struct will
%                 vary depending on test type (particularly for FRA), but the 
%                 fields are (this is for FRA):
% 				spiketimes   {# levels, # frequencies} cell array, each element is
% 								 an {nreps, 1} cell array of spiketimes (in milliseconds, 
% 								 relative to start of sweep) for each rep of this stimulus
% 				stimindex    {# levels, # frequencies} cell array, each element is
% 								 an [nreps, 1] vector of indices into spiketable{} 
% 								 for each instance of this stimulus combination
% 				stimvar      { [# frequencies, 1], [# levels, 1]} array of stimulus 
% 								 frequency and level (dB SPL) values for each of the
% 								 elements in spiketimes and stimindex. This really should be
% 								 in order of { levels, freqs} ...
% 				unique_stim  {2, 1} cell, simply vectors of frequency and level
% 				nstim        [2, 1] array, with # of freqs, levels
% 				spiketable   {total # stimulus combinations * # reps, 1} cell array
% 								 with elements from the SpikeData.Spikes table for this
% 								 unit's spikes for each stimulus. Has waveforms for spikes
% 				fileIndex    index into SpikeData.FileInfo{} cell of file/test
% 								 information
% 				channel      A/D channel
% 				unit         unit ID #
% 				fileName	    source raw .dat file from opto program
%
% With two output arguements, both spike_struct and info_struct will be
% returned. Fields in the structs will be different:
% 
% 		spike_struct	output struct with spike information/data
% 			Fields (this is for FRA data):
% 				spiketimes   {# levels, # frequencies} cell array, each element is
% 								 an {nreps, 1} cell array of spiketimes (in milliseconds, 
% 								 relative to start of sweep) for each rep of this stimulus
% 				spiketable   {total # stimulus combinations * # reps, 1} cell array
% 								 with elements from the SpikeData.Spikes table for this
% 								 unit's spikes for each stimulus. Has waveforms for spikes
% 				channel      A/D channel
% 				unit         unit ID #
% 
% 		info_struct	output struct with information about test
% 			Fields (for FRA data):
% 				stimindex    {# levels, # frequencies} cell array, each element is
% 								 an [nreps, 1] vector of indices into spiketable{} 
% 								 for each instance of this stimulus combination
% 				stimvar      { [# frequencies, 1], [# levels, 1]} array of stimulus 
% 								 frequency and level (dB SPL) values for each of the
% 								 elements in spiketimes and stimindex. This really should be
% 								 in order of { levels, freqs} ...
% 				unique_stim  {2, 1} cell, simply vectors of frequency and level
% 				nstim        [2, 1] array, with # of freqs, levels
% 				fileIndex    index into SpikeData.FileInfo{} cell of file/test
% 								 information
% 				fileName	    source raw .dat file from opto program
%------------------------------------------------------------------------
% See Also: SpikeData, SpikeInfo
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%  Sharad Shanbhag
%   sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 15 June, 2020 (SJS)
%	- adapted from OptoAnalysis:getFilteredOptoData
%
% Revisions:
%	18 Jun 2020 (SJS):
%		- modifying to return values for a single file/test, channel and unit
%		- using test fake data for now
%		- moved into method for SpikeData
%	24 Jan 2020 (SJS): moved different methods for each curve type
%	(subclass) into its respective file
%  3 Oct 2020 (SJS): mods for APANalyze, future analysis
%   - if looking across units, there is redundant information in the str
%     struct. modified output behavior to deal with this
%------------------------------------------------------------------------
% TO DO:
%   *Documentation!
%--------------------------------------------------------------------------

%{
--------------------------------------------------------------------------
--------------------------------------------------------------------------
ALGORITHM

Table in obj.Spikes has SpikeTimes in sequential order, aligned to start of
entire file

obj.spikesForAnalysis will return spikes for a specific file (can also
limit to specific channel and/or unit). for now, get spiketimes organized
into sweeps and across all channels
** NOT ACROSS ALL CHANNELS!

for each test, will need:
	list of stimuli and indices
		these are available via the CurveInfo.getStimulusIndices method


output will be in format that can be passed to computeRLF, computeFRA, etc:
% assumes spikeTimes is in format (for standard 1-D curves)
% 		spikeTimes{nLevels, 1}
% 			spikeTimes{n} = {nTrials, 1}
% 				spikeTimes{n}{t} = [spike1_ms spike2_ms spike3ms ...
%
% and
% 		spikeTimes{nLevels, nFreqs}
% 			where spikeTimes{l, f} = {nTrials, 1}
% 			and spikeTimes{l, f}{t} = [spike1_ms spike2_ms spike3ms ...]
% 
% for 2-D data (FRA, eventually WAV)


--------------------------------------------------------------------------
--------------------------------------------------------------------------
%}

%------------------------------------------------------------------------
% definitions
%------------------------------------------------------------------------
% string used in error messages
mname = 'SpikeData.getSpikesByStim';

%------------------------------------------------------------------------
% parse inputs
%------------------------------------------------------------------------
% need to use nargin == 4 since obj counts as an input in 
% object methods!!!!!
if nargin < 4
	error('%s: need file index #, channel and unit', mname);
end

% file index
if (fIndx < 1) || (fIndx > obj.Info.nFiles)
	error('%s: fIndx out of bounds', mname);
end

% channels
if length(chanNum) > 1
	error('%s: only works for single channel at a time', mname);
elseif obj.check_channels(chanNum) == -1
	error('%s: invalid channel %d', mname, chanNum);
end

% unit
if length(unitNum)  > 1
	error('%s: only works for single unit at a time', mname);
elseif ~obj.check_units(chanNum, unitNum)
	error('%s: invalid unit %d', mname, unitNum);
end	

%------------------------------------------------------------------------
% get spikes for each stimulus sweep for desired file, channel and unit
%------------------------------------------------------------------------
spiketable = obj.spikesForAnalysis(fIndx, 'channel', chanNum, ...
											'Unit', unitNum, 'Align', 'Sweep');
%------------------------------------------------------------------------
% convert to spiketimes array using this class' method
%------------------------------------------------------------------------
str = obj.Info.FileInfo{fIndx}.convertSpikeTableToSpikeTimes(spiketable);

if nargout == 1
	% everything goes in spike_struct (varargout{1})
	% add fields to str
	str.fileIndex = fIndx;
	str.channel = chanNum;
	str.unit = unitNum;
	str.fileName = obj.Info.FileInfo{fIndx}.F.file;
	% assign output
	varargout{1} = str;
else
	%{
	varargout{1} = struct(	'spiketimes', str.spiketimes, ...
									'spiketable', str.spiketable, ...
									'channel', chanNum, ...
									'unit', unitNum );
	varargout{2} = struct(	'stimindex', str.stimindex, ...
									'stimvar', str.stimvar, ...
									'unique_stim', str.unique_stim, ...
									'nstim', str.nstim, ...
									'fileIndex', fIndx, ...
									'fileName', obj.Info.FileInfo{fIndx}.F.file);
	%}
	% assign spike_struct fields to varargout{1}
	varargout{1}.spiketimes = str.spiketimes;
	varargout{1}.spiketable = str.spiketable;
	varargout{1}.channel = chanNum;
	varargout{1}.unit = unitNum;
	% assign info_struct to varargout{2}
	varargout{2}.stimindex = str.stimindex;
	varargout{2}.stimvar = str.stimvar;
	varargout{2}.unique_stim = str.unique_stim;
	varargout{2}.nstim = str.nstim;
	varargout{2}.fileIndex = fIndx;
	varargout{2}.fileName = obj.Info.FileInfo{fIndx}.F.file;

end
									