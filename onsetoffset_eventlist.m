%{



need overall list of:
1) stimuli
		name (xxx.wav, bbn, tone, etc.)
		level
		other parameter that's varied
2)	time of onset of stimulus (relative to merged file 0)
3) offset of stimulus (relative to merged file 0)

then, for each stimulus type, create events
	xxx_onset
	xxx_offset

%}



%% load test data - bbn and freq resp

% this will load the nI (SpikeInfo) object
load MERGALLobj.mat

%%
%{
SpikeInfo already has stimStartBin and stimEndBin for each file in a
{nfiles, 1} cell array. each element of the cell will be a [1, nsweeps]
vector of timestamps (in samples) indicating start and end of stimuli.
	stimStartBin: {[1×160 double]  [1×140 double]}
	stimEndBin: {[1×160 double]  [1×140 double]}

So, what's needed is a list of: stimulus combinations (stimtype, level,
etc). along with the sweep associated with each combination

CurveInfo.getStimulusIndices returns (for FREQ and LEVEL) test types
	stimindex 
		{nstimulus types, 1} cell array of sweep indices, where each
		element is a (# trials, 1) vector of indices into arrays dimensioned
		like stimStartBin and stimEndBin
	stimlist
		(nstimulus types, 1) vector of stimulus variables
		(LEVEL or FREQ)


For WAV and FRA, these are a little different

For wav types
	stimindex 
		{nstimulus types, 1} cell array of sweep indices, where each
		element is a (# trials, 1) vector of indices into arrays dimensioned
		like stimStartBin and stimEndBin
	stimlist
		{nstimulus types, 1} cell array of wav file name (no ext)

for FRA:
	stimindex 
		{n levels, nfreqs} cell array of sweep indices, where each
		element is a (# trials, 1) vector of indices into arrays dimensioned
		like stimStartBin and stimEndBin
	stimlist
		{1, 2} cell array with elements holding stimulus frequencies and
		stimulus levels stimlist{1} = list of frequencies, stimlist{2} = list
		of levels


%}

%------------------------------------------------------------------------
%------------------------------------------------------------------------
%% BBN, leveldata
%------------------------------------------------------------------------
%------------------------------------------------------------------------

findx = 1;

obj = nI.FileInfo{findx};

% get stimulus indices
[stimindex, stimlist] = obj.getStimulusIndices;
% get varied values
varied_values = obj.varied_values;
% get stimulus onset offset bins
% to align to appended/merged file, will need to add
% SpikeInfo.fileStartBin(findx) - 1
onsetbins = obj.stimStartBin;
offsetbins = obj.stimEndBin;

%------------------------------------------------------------------------
% create eventList as struct array
%------------------------------------------------------------------------
% {# stimulus levels, 2}
% 1st col will hold event name, second col will hold event times
nevents = length(varied_values);
events = repmat(	struct(	'name', '', ...
									'samples', [], ...
									'timestamps', [] ), ...
						nevents, 1);
for n = 1:nevents
	events(n).name = sprintf('%s_%s_%ddB', obj.testtype, obj.testname, varied_values(n));
	events(n).samples = onsetbins(stimindex{n});
end

% to write to nex file, this would be in SpikeData...
% 1) loop through nevents
for n = 1:nevents
	% 2) event times  are ebins + (filestartbin -1) / sampling rate
	events(n).timestamps = ((nI.fileStartBin(findx) - 1) + events(n).samples) ...
							./ nI.Fs;
			
	% 3) then add to nex struct/file
% 	 nD = nexAddEvent(nD, force_col(etimes{n}), eventSamplesByStim{n, 1})
end



%{
%------------------------------------------------------------------------
% create eventList as cell array
%------------------------------------------------------------------------
% {# stimulus levels, 2}
% 1st col will hold event name, second col will hold event times
nevents = length(varied_values);
eventSamplesByStim = cell(nevents, 2);
for n = 1:nevents
	ename = sprintf('%s_%ddB', obj.testtype, varied_values(n));
	ebins = onsetbins(stimindex{n});
	eventSamplesByStim{n, 1} = ename;
	eventSamplesByStim{n, 2} = ebins;
end

% to write to nex file, this would be in SpikeData...
etimes = cell(nevents, 1);
% 1) loop through nevents
for n = 1:nevents
	% 2) event times  are ebins + (filestartbin -1) / sampling rate
	etimes{n} = ((nI.fileStartBin(findx) - 1) + eventSamplesByStim{n, 2}) ...
							./ nI.Fs;
			
	% 3) then add to nex struct/file
% 	 nD = nexAddEvent(nD, force_col(etimes{n}), eventSamplesByStim{n, 1})
end
%}



%{
%------------------------------------------------------------------------
%------------------------------------------------------------------------
%% FREQ, tuning data
%------------------------------------------------------------------------
%------------------------------------------------------------------------

% get stimulus indices
[stimindex, stimlist] = nI.FileInfo{2}.getStimulusIndices;
% get varied_values - freqs
varied_values = nI.FileInfo{2}.varied_values;
% get stimulus onset offset bins
onsetbins = nI.FileInfo{2}.stimStartBin;
offsetbins = nI.FileInfo{2}.stimEndBin;

% create eventList
% {# stimulus levels, 2}
% 1st col will hold event name, second col will hold event times
eventList = cell(length(varied_values), 2);


%------------------------------------------------------------------------
%------------------------------------------------------------------------
%% wav data
%------------------------------------------------------------------------
%------------------------------------------------------------------------

% get stimindices
[stimindex, stimlist] = nI.FileInfo{3}.getStimulusIndices;
% get matching levels (convert to matrix/vector from cell array)
levellist = cell2mat(nI.FileInfo{3}.getlevelList);


%------------------------------------------------------------------------
%------------------------------------------------------------------------
%% FRA data
%------------------------------------------------------------------------
%------------------------------------------------------------------------
% get stimindices
[stimindex, stimlist] = nI.FileInfo{4}.getStimulusIndices;
%}