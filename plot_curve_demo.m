% function varargout = build_curve_from_spikes(spikesByStim, unit, curveInfo)

%% load data

% S object file
Sfile = '~/Work/Data/TestData/MT/1382_20191212_02_02_3200_Sobj.mat';

% load data
load(Sfile)

%% test things

% need to select a unit
allUnits = S.listUnits;
fprintf('Found spike data for units: ')
fprintf('%d  ', S.listUnits);
fprintf('\n\n');

%% Show how to list files and curve types

for f = 1:S.Info.nFiles
	fprintf('File %d:\n', f);
	fprintf('\tTest Type: %s\n', S.Info.FileData(f).testtype);
	
end

%%

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

S = load_plexon_data(fullfile(sortedPath, sortedFile), ...
							fullfile(sortedPath, nexInfoFile))