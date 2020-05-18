function [cSweeps, sInfo] = read_data_for_export(F, Channels, BPfilt)
%------------------------------------------------------------------------
% [cSweeps, sInfo] = read_data_for_export(F, Channels, BPfilt)
%------------------------------------------------------------------------
% TytoLogy:Experiments:optosort
%------------------------------------------------------------------------
% 
% reads data from multiply raw .dat files generated by opto program
%
% called by export_xxx functions for exporting data to spike sorting
% programs
%% 
%------------------------------------------------------------------------
% Input Arguments:
%
%	F	file object or array of OptoFileName objects (see defineSampleData())
%		
%	Channels:
% 	 neural A/D channels to export from data file(s). 
% 	 note that channels must be consistent across all the data files to be
% 	 merged in the exported file!
% 	 If blank, all channels present in file will be included
% 		Channels = [11 9 14];
% 
% 	BPfilt:
% 	 specifies filter for processing output data.
% 	 if blank/unspecified, no filter will be applied to data
% 		[highpass lowpass] cutoff frequencies in Hz:
% 		 BPfilt.Fc = [300 4000];
% 		order of filter. note that the filtfilt() function in MATLAB is used,
% 		so the effective order is doubled. typically use 5:
% 		 BPfilt.forder = 5;
% 		ramp time (ms) to apply to each sweep in order to cutdown on onset/offset
% 		transients from filtering:
% 		 BPfilt.ramp = 1;
%		filter type is either 'bessel' or 'butter' (butterworth, default)
%		 BPfilet.type = 'bessel';
% 
% Output Arguments:
% 	cSweeps	{nfiles, 1} cell array of sweep data
%	sInfo	SpikeInfo object

%------------------------------------------------------------------------
% See also: export_for_plexon, export_raw
%				OptoFileName, SpikeInfo, CurveInfo, WAVInfo classes
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%  Sharad J. Shanbhag
%	sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 18 May, 2020 (extracted from export_for_plexon) (SJS)
%
% Revisions:
%------------------------------------------------------------------------
% TO DO:
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% some initial things
%------------------------------------------------------------------------
% determine # of files
nFiles = length(F);

%------------------------------------------------------------------------
% create sInfo object (SpikeInfo) to hold sweep/file bin and time data
%------------------------------------------------------------------------
sInfo = SpikeInfo();
%{
sInfo.FileName = fullfile(NexFilePath, NexFileName);
% create output _nexinfo.mat file name - base is same as .nex file
[~, baseName] = fileparts(sInfo.FileName);
sInfo.InfoFileName = fullfile(NexFilePath, [baseName '_nexinfo.mat']);
%}
% store channel information
sInfo.ADchannel = Channels;


%------------------------------------------------------------------------
% pre-allocate some things
%------------------------------------------------------------------------
% bins for start and end of each file's data
sInfo.fileStartBin = zeros(1, nFiles);
sInfo.fileEndBin = zeros(1, nFiles);
% bins for all sweep starts and ends
sInfo.sweepStartBin = cell(1, nFiles);
sInfo.sweepEndBin = cell(1, nFiles);
% bins for stim onset, offset
sInfo.stimStartBin = cell(1, nFiles);
sInfo.stimEndBin = cell(1, nFiles);

% each file's sampling rate for neural data
tmpFs = zeros(nFiles, 1);
% cell array to hold sweep data - this will be converted to a single
% "vector" of values per channel that will be added to the .nex file
% (algorithm will be outlined in the next section)
cSweeps = cell(nFiles, 1);
% allocate cell array to hold information for each file
% we'll assign/allocate/initialize objects (CurveInfo, WAVINfo) as we loop
% through the files in the next section
cInfo = cell(nFiles, 1);

%------------------------------------------------------------------------
% Read and process data from raw files
%------------------------------------------------------------------------

% loop through files
for f = 1:nFiles
	sendmsg(sprintf('Reading data for file %s', F(f).file));
	
	% get data for each file and channel and convert to row vector format
	% algorithm:
	%		(1) put each sweep for this channel in a {1, # sweeps} cell array
	%				cSweeps
	%		(2) make a note of the length of each sweep to use for
	%				markers/timestamps
	%		(3) after cSweeps is built, convert to a row vector using cell2mat

	% use readOptoData to read in raw data. 
	[D, tmpDinf] = readOptoData(fullfile(F(f).path, F(f).file));
	% Fix test info
	tmpDinf = correctTestType(tmpDinf);
	
	% CurveInfo (or WavInfo) to hold everything for each file
	switch(upper(tmpDinf.test.Type))
		case {'FREQ', 'LEVEL', 'FREQ+LEVEL', 'OPTO'}
			cInfo{f} = CurveInfo(tmpDinf);
			
		case {'WAV', 'WAVFILE'}
			% build wavinfo file and load it
			wavinfo_filename = [F(f).base '_wavinfo.mat'];
			if ~exist(fullfile(F(f).path,wavinfo_filename), 'file')
				warning('%s: wavinfo file %s not found', mfilename, ...
																wavinfo_filename);
				cInfo{f} = WAVInfo(tmpDinf);
			else
				tmpW = load(fullfile(F(f).path,wavinfo_filename));
				cInfo{f} = WAVInfo(tmpDinf, tmpW);
			end
	
		otherwise
			error('%s: unknown test.Type %s', mfilename, tmpDinf.test.Type);
	end
	
	% build filter for neural data
	if ~isempty(BPfilt)
		BPfilt.Fs = cInfo{f}.Dinf.indev.Fs;
		BPfilt.Fnyq = cInfo{f}.Dinf.indev.Fs / 2;
		BPfilt.cutoff = BPfilt.Fc / BPfilt.Fnyq;
		if strcmpi(BPfilt.type, 'bessel')
			[BPfilt.b, BPfilt.a] = besself(BPfilt.forder, BPfilt.cutoff, 'bandpass');
		elseif strcmpi(BPfilt.type, 'butter')
			[BPfilt.b, BPfilt.a] = butter(BPfilt.forder, BPfilt.cutoff, 'bandpass');
		else
			error('%s: unknown filter type %s', mfilename, BPfilt.type)
		end
	end
	
	% check to make sure consistent # of sweeps (aka trials)
	if cInfo{f}.Dinf.nread ~= length(D)
		error('%s: mismatch in Dinf.nread (%d) and length(D) (%d)', ...
					mfilename, cInfo{f}.Dinf.nread, length(D));
	end
	
	% build into sweeps by channel format
	fprintf('Test type: %s\n', cInfo{f}.testtype);
	[cInfo{f}, cSweeps{f}] = cInfo{f}.buildChannelData(Channels, BPfilt, D);
	
	% store sample for start of this file (should be 1); use channel 1 value
	 cInfo{f}.fileStartBin = cInfo{f}.startSweepBin{1}(1);
	% store sample for end of this file
	 cInfo{f}.fileEndBin = cInfo{f}.endSweepBin{1}(end);

	% to avoid any issues, should make sure sample rates are consistent
	% to do this, store list of sample rates and check once out of this loop
	tmpFs(f) = cInfo{f}.Dinf.indev.Fs;
	
	% calculate overall start and end bins for each file's data
	if f == 1
		sInfo.fileStartBin(f) = cInfo{f}.fileStartBin;
		sInfo.fileEndBin(f) = cInfo{f}.fileEndBin;
	else
		% add 1 to prior end bin for start
		sInfo.fileStartBin(f) = sInfo.fileEndBin(f-1) + 1;
		sInfo.fileEndBin(f) = sInfo.fileStartBin(f) + cInfo{f}.fileEndBin - 1;
	end
	
	% compute stimulus onset/offset bins for this file
	cInfo{f} = cInfo{f}.buildStimOnOffData;
		
end

% assign cInfo to sInfo.FileData
sInfo.FileInfo = cInfo;
% store filter info
sInfo.dataFilter = BPfilt;

sendmsg('Checking sample rates across files');
% check sampling rates
if ~all(tmpFs(1) == tmpFs)
	error('%s: Sample Rate Mismatch!!!', mfilename);
else
	% store overall sample rate
	Fs = tmpFs(1);
	sInfo.Fs = Fs;
	sendmsg(sprintf('Neural A/D Fs = %.4f', Fs));
end

%------------------------------------------------------------------------
% Create file, sweep, stimulus start/end bin and timestamp vectors
%------------------------------------------------------------------------
sendmsg('Building start and end sweep indices:');
% assign values for bins
for f = 1:nFiles
	% calculate start and end sweep bins for each file's data
	if f == 1
		sInfo.sweepStartBin{f} = cInfo{f}.startSweepBin{1};
		sInfo.sweepEndBin{f} = cInfo{f}.endSweepBin{1};
	else
		% add previous file's final endSweepBin value as offset
		sInfo.sweepStartBin{f} = cInfo{f}.startSweepBin{1} + ...
												sInfo.sweepEndBin{f-1}(end);
		sInfo.sweepEndBin{f} = cInfo{f}.endSweepBin{1} + ...
												sInfo.sweepEndBin{f-1}(end);
	end
	
	% for stim onset/offset, align to file start bin *** would this work for
	% the sweep start and end computation above???? need to test!!!!!
	sInfo.stimStartBin{f} = sInfo.fileStartBin(f) + cInfo{f}.stimStartBin;
	sInfo.stimEndBin{f} = sInfo.fileStartBin(f) + cInfo{f}.stimEndBin;
end
