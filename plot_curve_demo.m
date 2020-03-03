% function varargout = build_curve_from_spikes(spikesByStim, unit, curveInfo)

%% load data

% data locations
sortedPath = '~/Work/Data/TestData/MT';
rawPath = '~/Work/Data/TestData/MT';
nexPath = '~/Work/Data/TestData/MT';

% sorted data file
% sortedFile = '1323_20190722_03_02_632_MERGE.mat';
sortedFile = '1382_20191212_02_02_3200.mat';

% nexinfo file
% nexInfoFile = '1323_20190722_03_02_632_MERGE_nexinfo.mat';
nexInfoFile = '1382_20191212_02_02_3200_MERGE_nexinfo.mat';

% nex file
nexFile = '1382_20191212_02_02_3200_MERGE.nex';

% load the data - returned as a SpikeData object
S = load_plexon_data(fullfile(sortedPath, sortedFile), ...
							fullfile(sortedPath, nexInfoFile));
						
%% test things

% show unit id numbers
fprintf('Found spike data for units: ')
fprintf('%d  ', S.listUnits);
fprintf('\n\n');

% Show how to list files and curve types
% loop through the number of files merged into file to be sorted
for f = 1:S.Info.nFiles
	% for each file number, display the test type 
	fprintf('File %d:\n', f);
	fprintf('\tTest Type: %s\n', S.Info.FileData(f).testtype);
end

%% plot freq-tuning curves for non-zero units





