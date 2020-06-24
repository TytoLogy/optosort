function varargout = getSpikesByStim(obj, fIndx, chanNum, unitNum)
%------------------------------------------------------------------------
% [data, datainfo, tracesByStim] = viewOptoData(varargin)
%------------------------------------------------------------------------
% TytoLogy:opto
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
% Output Arguments:
%
%
%------------------------------------------------------------------------
% See Also:
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

% add fields to str
str.fileIndex = fIndx;
str.channel = chanNum;
str.unit = unitNum;
str.fileName = obj.Info.FileInfo{fIndx}.F.file;
% assign output
varargout{1} = str;

	