%------------------------------------------------------------------------
% plot_curve_demo.m
%------------------------------------------------------------------------
% TytoLogy:optosort
%--------------------------------------------------------------------------
% example script for generating plots of data
%	plots shown:
%		
%------------------------------------------------------------------------
% See Also: optoproc, opto (TytoLogy:opto program)
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%  Sharad Shanbhag
%   sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 25 June 2020 (SJS)
%	 
% Revisions:
%
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%% add path to readPLXFileC if needed
%------------------------------------------------------------------------
% readPLXFileC is a function downloaded from MATLAB Central that allows
% direct reading of PLX file data in Matlab
if ~exist('readPLXFileC', 'file')
	fprintf('plot_curve_demo: adding readPLXFile to path\n');
	try
		addpath('readPLXFileC');
	catch
		error('Cannot find readPLXFileC function. Please add to path!')
	end
end

%------------------------------------------------------------------------
%------------------------------------------------------------------------
%% define paths and data files
%------------------------------------------------------------------------
%------------------------------------------------------------------------
% data locations - adjust this for your setup
%------------------------------------------------------------------------

% location of .plx file (from Plexon OfflineSorter)
plxFilePath = '/Users/sshanbhag/Work/Data/TestData/working';
% location of _nexInfo.mat file (generated by export_for_plexon)
nexInfoPath = plxFilePath;

% sorted data file name
plxFile = '1407_20200309_03_01_1350_MERGEEVENTS2.plx';
% nexinfo file
nexInfoFile = '1407_20200309_03_01_1350_MERGEEVENTS2_nexinfo.mat';

sendmsg(sprintf('Using data from file: %s', plxFile));

%------------------------------------------------------------------------
%------------------------------------------------------------------------
%------------------------------------------------------------------------
%% load sorted data
%------------------------------------------------------------------------
% How to use:
% import_from_plexon(<plx file name>, <nexinfo file name>, 
%								<'continuous'/'nocontinuous'>)
%
%	will return a SpikeData object containing data from plx file:
% 	- spike times/unit information if sorted
% 	- continuous data, if saved in plx and 'continuous' is 
%		specified (default)
% 	- file stimulus info
%------------------------------------------------------------------------
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% there are two options here - if you feel that you'll want to look at the
% continuous data (recordings from electroeds), you can specify the
% 'continuous' option to load them.  
%------------------------------------------------------------------------
S = import_from_plexon(fullfile(plxFilePath, plxFile), ...
							fullfile(nexInfoPath, nexInfoFile), 'continuous');
%------------------------------------------------------------------------
% or, specify 'nocontinuous' to tell import_from_plexon to ignore the
% continuous channel data. this should save on memory
% S = import_from_plexon(fullfile(plxFilePath, plxFile), ...
%------------------------------------------------------------------------
% 							fullfile(nexInfoPath, nexInfoFile), 'nocontinuous');

% The returned SpikeData object, S, contains:
%	Info
%		information about the experiment, stimuli, etc in the Info object
%	Spikes
%		results from spike sorting, stored as a MATLAB table object
% 		Typically, data in the Spikes table will be accessed through methods
%		that are part of the SpikeData class.
%	Continuous
%		
%	plxvar
%		holdover from attempted solution for missing A/D channel info
%		when exporting from Plexon OFS.
%	hasContinuousData
% 		1 if continuous data are loaded, 
%------------------------------------------------------------------------
%------------------------------------------------------------------------
%% information about file
%------------------------------------------------------------------------
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% show file, channel, unit
%------------------------------------------------------------------------
% display file, channel, unit info, store information about files, channels
% and units
[fileList, channelList, unitList] = S.printInfo;

%------------------------------------------------------------------------
% Show how to list files and curve types
%------------------------------------------------------------------------
% loop through the number of files merged into file to be sorted
sendmsg('File and test information:')
for f = 1:S.Info.nFiles
	% for each file number, display the test type and test name (i know why
	% two different things? an effect of an old kludge...)
	fprintf('File %d:\n', f);
	fprintf('\tTest Type: %s\n', S.Info.FileInfo{f}.testtype);
	fprintf('\tTest Name: %s\n', S.Info.FileInfo{f}.testname);
end

%{
alternative: use fileList
%}


%------------------------------------------------------------------------
%------------------------------------------------------------------------
%% test, channel and unit to process and plot
%------------------------------------------------------------------------
%------------------------------------------------------------------------

%-------------------------------------------------------
% What test data do you want to plot?
% Options:
% testToPlot = 'FREQ_TUNING';
% testToPlot = 'FRA';
% testToPlot = 'BBN';
% testToPlot = 'WAV';
%-------------------------------------------------------
testToPlot = 'BBN';

% specify channel to plot
channel = 4;
% specify unit:
unit = 1;
% binsize (in milliseconds) for psth
psth_bin_size = 5;

%-------------------------------------------------------
% figure out file index for this test. the indexForTestName method of
% SpikeData class allows an easy way to do this.
%-------------------------------------------------------
findx = S.indexForTestName(testToPlot);
if isempty(findx)
	fprintf('Test %s not found in %s\n', testToPlot, plxFile);
	error('%s: bad testToPlot', mfilename);
end

%------------------------------------------------------------------------
% get spikes times struct (store in st) for this test, channel and unit
% spiketimes will be aligned to start of each sweep
%------------------------------------------------------------------------
fprintf('Getting data for file %d (%s), channel, %d unit %d\n', ...
								findx, S.listFiles{findx}, channel, unit);
st = S.getSpikesByStim(findx, channel, unit);
% st  struct with fields:
% st = 
%      spiketimes: {53�1 cell}
%       stimindex: {53�1 cell}
%         stimvar: {53�1 cell}
%     unique_stim: {53�1 cell}
%           nstim: 53
%      spiketable: {1060�1 cell}
%       fileIndex: 3
%         channel: 5
%            unit: 1
%        fileName: '1407_20200309_03_01_1350_WAV.dat'
		 
% make a local copy of Dinf for this file to make things a little simpler
Dinf = S.Info.FileInfo{findx}.Dinf;


%------------------------------------------------------------------------
%------------------------------------------------------------------------
%% Two ways to plot waveforms:
%------------------------------------------------------------------------
%------------------------------------------------------------------------
%------------------------------------------------------------------------
%% to plot spike waveforms for this test
%------------------------------------------------------------------------
% (1) vertically concatenate data stored in cell array of sweeps in
% st.spiketable into a single table
tmpT = vertcat(st.spiketable{:});

% (2) extract just the wave field 
tmpwav = tmpT.Wave;

% (3) plot overlaid waveforms

% plot in new figure
figure

% need time vector for x-axis
[~, nBins] = size(tmpwav);
t_ms = (1000/S.Info.Fs) * (0:(nBins - 1));

% plot waveforms, with mean and legend
% need tmpwav to be in column format - time in rows, indv. units by column
% so send the function the transpose of the tmpwav matrix.
hWF = plot_spike_waveforms(t_ms, tmpwav', 'MEAN', true, 'LEGEND', true);

% add title to plot
% create title string with 2 rows:
%	filename (either from st struct or S.Info.FileInfo{findx}.F.file
%	channel and unit
tstr = {	st.fileName, ...
			sprintf('Channel %d Unit %d', channel, unit)};
title(tstr, 'Interpreter', 'none');	

% set figure filename - use the base from the FreqTuningInfo.F object
set(gcf, 'Name', sprintf('%s_Ch%d_Un%d', ...
					S.Info.FileInfo{findx}.F.base, st.channel, st.unit));


%------------------------------------------------------------------------
%% plot this unit for ALL tests in sorted file
%plot waveforms for this channel and unit - note that these are extracted
%from entire file
%------------------------------------------------------------------------
S.plotUnitWaveforms(channel, unit);

%------------------------------------------------------------------------
%------------------------------------------------------------------------
%% LEVEL/BBN tests:
%	plot rate-level curves (firing rate- level function)
%------------------------------------------------------------------------
if any(strcmpi(testToPlot, {'LEVEL', 'BBN'}))
	% set analysis window to [stimulus onset   stimulus offset]
	analysisWindow = [Dinf.audio.Delay ...
									(Dinf.audio.Delay + Dinf.audio.Duration)];
	% compute rate level function
	RLF = computeRLF(st.spiketimes, st.unique_stim, analysisWindow);
	% plot
	hRLF = plotCurveAndCI(RLF, 'mean');
	% create title string with 2 rows:
	%	filename;
	%	channel and unit
	tstr = {	st.fileName, ...
				sprintf('Channel %d Unit %d', channel, unit)};
	% add title to plot
	title(tstr, 'Interpreter', 'none');
	% set filename - use the base from the FreqTuningInfo.F object
	set(hRLF, 'Name', sprintf('%s_Ch%d_Un%d', ...
					S.Info.FileInfo{findx}.F.base, st.channel, st.unit));
end
%------------------------------------------------------------------------
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%------------------------------------------------------------------------
%% FREQ_TUNING tests:
%	plot frequency-tuning curves
%------------------------------------------------------------------------
%------------------------------------------------------------------------
if any(strcmpi(testToPlot, 'FREQ_TUNING'))
	% set analysis window to [stimulus onset   stimulus offset]
	analysisWindow = [Dinf.audio.Delay ...
									(Dinf.audio.Delay + Dinf.audio.Duration)];
	% compute rate level function
	FTC = computeFTC(st.spiketimes, st.unique_stim, analysisWindow);
	% plot
	hFTC = plotCurveAndCI(FTC, 'median');
	% create title string with 2 rows:
	%	filename
	%	channel and unit
	tstr = {	st.fileName, ...
				sprintf('Channel %d Unit %d', channel, unit)};
	% add title to plot
	title(tstr, 'Interpreter', 'none');
	% set filename - use the base from the FreqTuningInfo.F object
	set(hFTC, 'Name', ...
					sprintf('%s_Ch%d_Un%d', S.Info.FileInfo{findx}.F.base, ...
										st.channel, st.unit));
end


%------------------------------------------------------------------------
%------------------------------------------------------------------------
%% for WAV data, plot PSTH and rasters
% use different plots/pages for each different stimulus level
%------------------------------------------------------------------------
%------------------------------------------------------------------------
if strcmpi(testToPlot, 'WAV')
	% plotPSTH is a method in the WAVInfo class, stored at the findx
	% element within the FileInfo cell array. 'LEVEL' tells the method to
	% separate plots according to stimulust level
	hWAV = S.Info.FileInfo{findx}.plotPSTH(st, psth_bin_size, 'LEVEL');

	% rename plots with filename, level, channel and unit
	for h = 1:length(hWAV)
		% build on original plot name, since this has the level information
		% already in it
		origname = get(hWAV{h}, 'Name');
		set(hWAV{h}, 'Name', ...
					sprintf('%s_Ch%d_Un%d', origname, st.channel, st.unit));
	end
end


%------------------------------------------------------------------------
%------------------------------------------------------------------------
%% FRA (frequency-response area) tests:
%	plot heat map and "waterfall" plot 
%------------------------------------------------------------------------
%------------------------------------------------------------------------
if any(strcmpi(testToPlot, 'FRA'))
	% window for spike count
	frawin = [Dinf.audio.Delay (Dinf.audio.Delay + Dinf.audio.Duration)];
	% calculate FRA stored in struct FRA
	FRA = computeFRA(st.spiketimes, st.unique_stim{1}, ...
															st.unique_stim{2}, frawin);
	% set fname in FRA struct to data file name
	FRA.fname = st.fileName;
	hFRA = plotFRA(FRA, 'dB');	
	% set filename - use the base from the FreqTuningInfo.F object
	set(hFRA, 'Name', sprintf('%s_Ch%d_Un%d', ...
					S.Info.FileInfo{findx}.F.base, st.channel, st.unit));

end

