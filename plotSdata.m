% plotSdata
%
% script to plot/select sweeps


%------------------------------------------------------------------------
%% load data object
%------------------------------------------------------------------------
sortedPath = '/Users/sshanbhag/Work/Data/TestData/MT/1407';
rawPath = sortedPath;
nexPath = sortedPath;
nexInfoFile = '1407_20200309_03_01_1350_BBN_nexinfo.mat';
nexFile = '1407_20200309_03_01_1350_BBN.nex';
plxFile = '1407_20200309_03_01_1350_BBN-sorted.ch4,5,7,15.plx';
sfile = '1407_20200309_03_01_1350_BBN-sorted.ch4,5,7,15_Sobj.mat';
%------------------------------------------------------------------------
% get base from one of the file objects in S.Info
fprintf('\n%s\n', sepstr);
fprintf('loading Sobj from file\n\t%s\n', fullfile(sortedPath, sfile));
fprintf('%s\n', sepstr);
load(fullfile(nexPath, sfile), '-MAT', 'S');
% get plx data
Plx = PLXData(fullfile(sortedPath, plxFile), 'all', 'continuous');
%------------------------------------------------------------------------
%% get data by file
%------------------------------------------------------------------------
%% extract timestamps for use in analysis, raster plots, etc., separated by channel
% convert to timestamps
spikesForAnalysis = cell(S.Info.nFiles, S.Info.nChannels);

% loop through files
for f = 1:S.Info.nFiles
	% loop through channels
	for c = 1:S.Info.nChannels
		fprintf('Getting spikes for file %d, channel %d\n', f, S.Info.ADchannel(c));
		% extract timestamp data from each table, store in cell matrix
		spikesForAnalysis{f, c} = S.spikesForAnalysis(f, ...
									'Channel', S.Info.ADchannel(c), 'align', 'sweep');
	end
end

%------------------------------------------------------------------------
%% resample continuous data (demo for export... function)
%------------------------------------------------------------------------

cIndx = 1;
channel = S.Info.ADchannel(cIndx);
sweepStartBin = S.Info.sweepStartBin{1};
sweepEndBin = S.Info.sweepEndBin{1};



% convert i16 data to double & scale by continuous channel max magnitude
dData = double(S.Continuous(cIndx).Values) / double(Plx.P.ContMaxMagnitudeMV);




%%

%

% specify sample rates
Fs_old = S.Info.Fs;
Fs_new = 48820;
fprintf('Original Sample Rate: %.4f\n', Fs_old);
fprintf('Resampled SampleRate: %.4f\n', Fs_new);

% determine ratio via closest rational approximation
[p, q] = rat(Fs_new/Fs_old, 1e-9)
fprintf('New Ratio:\n');
fprintf('p = %.2f\nq = %.2f\n', p, q);
fprintf('Error abs( (p/q)*Fs_old - Fs_new): %.12f\n', abs( (p/q)*Fs_old - Fs_new));
cData = resample(dData, p, q);

%------------------------------------------------------------------------
%% plot
%------------------------------------------------------------------------
fH = figure(142);


for s = 1:10
	% for each sweep, will need:
	%	channel to display
	%	unit(s) to display
	%	continuous data snipped (sweepStart, sweepEnd)
	tvec = (sweepStartBin(s):sweepEndBin(s)) / Fs_new;
	plot(tvec, cData(sweepStartBin(s):sweepEndBin(s)))
	
	pause(1)
end




