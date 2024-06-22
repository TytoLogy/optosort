%------------------------------------------------------------------------
% SI_binarydemo.m
%------------------------------------------------------------------------
% TytoLogy:Experiments:optosort
%------------------------------------------------------------------------
% sample script to demonstrate binary data for use with SpikeInterface
% 
%------------------------------------------------------------------------
% See also: 
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%  Sharad J. Shanbhag
%	sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 21 June, 2024 from code in exportTest_SpikeInterface and
%          exportTestRaw (SJS)
%
% Revisions:
%  14 May, 2020 exportTestRaw created from exportTest (SJS)
%   1 May, 2024 exportTest_SpikeInterface created from exportTestRaw (SJS)
%               uses OutputShape->'SamplesChannels' option to write output
%               data in format [samples, channels] needed for Kilosort
%------------------------------------------------------------------------
% TO DO:
%------------------------------------------------------------------------

% check to make sure tytoLogy (esp opto) things are on path
if ~exist('readOptoData', 'file')
	fprintf(['readOptoData (and possibly other TytoLogy/opto functions' ...
						'not found!\n'])
	fprintf('Please check Matlab path(s)\n');
	fprintf(['e.g.:\n' ...
				'addpath(''~/Work/Code/Matlab/dev/TytoLogy/Experiments/Opto'')\n']);
end

if ~exist('SpikeData', 'file')
	fprintf('SpikeData.m class definition file not found!\n')
	fprintf('This is usually found in the OptoObjects folder\n')
	fprintf('Please check Matlab path(s)\n');
end

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% PATHS TO DATA FILES
%------------------------------------------------------------------------
%------------------------------------------------------------------------


% Wav Only, BLA data
%{
exportOpts.DataPath = {'/media/Data/NeuroData/Mouse/Raw/BLA/1500/20230213'};
exportOpts.OutputPath = '/media/Data/NeuroData/TestData/exports/SI';
exportOpts.DataFile = {	'1500_20230213_01_0_3352_WAV.dat'};
exportOpts.TestFile = {	'1500_20230213_01_0_3352_WAV_testdata.mat'};
exportOpts.OutputFile = '1500_20230213_01_0_3352_WAV.bin';
exportOpts.Channels = 1:16;
%}

% IC data for testing
exportOpts.DataPath = {'/Users/sshanbhag/Work/Data/IC/working/Raw/1429/20200701'};
exportOpts.OutputPath = '/Volumes/SSData/Data/Mouse/Testing/SI';
exportOpts.DataFile = {	'1429_20200701_01_0_817_WAV.dat'};
exportOpts.TestFile = {	'1429_20200701_01_0_817_WAV_testdata.mat'};
exportOpts.OutputFile = '1429_20200701_01_0_817_WAV.bin';
exportOpts.Channels = 1:16;
%---------------------------------------
%---------------------------------------

%------------------------------------------------------------------------
% Specify output shape as [samples, channels] as needed for SpikeInterface
% i.e., samples in rows, channels in columns
%------------------------------------------------------------------------
exportOpts.OutputShape = 'SamplesChannels';

%------------------------------------------------------------------------
% SpikeInterface wants data in 'double' format
%------------------------------------------------------------------------
exportOpts.OutputFormat = 'double';


%---------------------------------------------------------------------	
%---------------------------------------------------------------------	
% this is code from export_raw
%---------------------------------------------------------------------	
%---------------------------------------------------------------------	

%---------------------------------------------------------------------	
% define path to data file and data file for testing
%---------------------------------------------------------------------
F = defineSampleData(exportOpts.DataPath, exportOpts.DataFile, exportOpts.TestFile);

% # channel(s) of data to obtain
nChannels = length(exportOpts.Channels);
nFiles = 1;

%------------------------------------------------------------------------
%% get the data and information from the raw files
%------------------------------------------------------------------------
% no filter or resample
BPfilt = [];
resampleData = [];
[cSweeps, expInfo] = read_data_for_export(F, exportOpts.Channels, ...
                                                BPfilt, resampleData);

%%

OutputFile = fullfile(exportOpts.OutputPath, exportOpts.OutputFile);

%{
from https://spikeinterface.readthedocs.io/en/latest/how_to/load_matlab_data.html

----------------------------------------------------
Writing data for SpikeInterface:
----------------------------------------------------
% Define the size of your data
numSamples = 1000;
numChannels = 384;

% Generate random data as an example
data = rand(numSamples, numChannels);

% Write the data to a binary file
fileID = fopen('your_data_as_a_binary.bin', 'wb');
fwrite(fileID, data, 'double');
fclose(fileID);
%}

%% select subset of data, assign dummy values
nReps = 2;
swpdata{1} = cSweeps{1}(:, 1:nReps);

dumbdata{1} = swpdata{1};
for c = 1:nChannels
   for s = 1:nReps
      dumbdata{1}{c, s} = (c+0.1*s)*ones(size(swpdata{1}{c, s}))++1e-6*(1:19531);
   end
end

%------------------------------------------------------------------------
%%
% convert (concatenate) cSweeps to [nchannels, nsamples] matrix and write
% to binary output file
%------------------------------------------------------------------------
sendmsg(sprintf('Exporting raw binary data to %s', expInfo.FileName));

% open raw file
fp = fopen(OutputFile, 'wb');

fprintf('Writing data for %s\n', F(f).file);
% write data for this file, specify shape based on outputShape setting
% 1 May 2024: this may not have any effect since binary files store
% everything in 1D format!!!!
switch exportOpts.OutputShape
   case 'ChannelsSamples'
      tmpmat = cell2mat([dumbdata{f}]);
      nw = fwrite(fp, tmpmat, exportOpts.OutputFormat);
   case 'SamplesChannels'
      tmpmat = cell2mat([dumbdata{f}])';
      nw = fwrite(fp, tmpmat, exportOpts.OutputFormat);
end

% could also:
%	concatenate: tmp = [cVector{:}];
%	fwrite(fp, cell2mat(tmp), 'float64');

fprintf('\t%d values written\n', nw);
fclose(fp);

%------------------------------------------------------------------------
%% as a check, read in binary data
%------------------------------------------------------------------------
nC = length(exportOpts.Channels);

fp = fopen(fullfile(exportOpts.OutputPath, exportOpts.OutputFile), ...
                  'rb');
M = fread(fp, exportOpts.OutputFormat);
fclose(fp);
p
sM = size(M);
fprintf('read %d elements\n', numel(M));
fprintf('Size of data read from file: [%d, %d]\n', sM(1), sM(2));
M2 = reshape(M, [numel(M)/nC, nC]);
size(M2)




%{
%% do a simple test
nC = length(exportOpts.Channels);

a = ones(nC, 20);
for n = 1:nC
   a(n, :) = n*a(n, :);
end

%%

% open raw file
fp = fopen('test.bin', 'wb');
nw = fwrite(fp, a, exportOpts.OutputFormat);
fclose(fp);

% read raw file
fp = fopen('test.bin', 'rb');
b  = fread(fp, exportOpts.OutputFormat);
fclose(fp);
%%
reshape(b, size(a))

reshape(b, [nC, numel(b)/nC])
%}

