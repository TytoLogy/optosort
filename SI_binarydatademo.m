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

%% create syndata cell with synthetic arrays of data - this is how data
% are read in by our matlab code, although with many more Reps

% use 2 reps, 19531 samples per sweep
nReps = 2;
swpPoints = 19531;

syndata{1} = cell(nChannels, nReps);
for c = 1:nChannels
   for s = 1:nReps
      syndata{1}{c, s} = (c+0.1*s)*ones(1, 19531)++1e-6*(1:19531);
   end
end

%------------------------------------------------------------------------
%%
% convert (concatenate) syndata to [nchannels, nsamples] matrix and write
% to binary output file
%------------------------------------------------------------------------
sendmsg(sprintf('Exporting raw binary data to %s', 'syndata.bin'));

% open raw file
fp = fopen('syndata.bin', 'wb');

fprintf('Writing data\n');
% write data for this file, specify shape based on outputShape setting
% 1 May 2024: this may not have any effect since binary files store
% everything in 1D format!!!!
switch exportOpts.OutputShape
   case 'ChannelsSamples'
      tmpmat = cell2mat([syndata{f}]);
      nw = fwrite(fp, tmpmat, 'double');
   case 'SamplesChannels'
      tmpmat = cell2mat([syndata{f}])';
      nw = fwrite(fp, tmpmat, 'double');
end

% could also:
%	concatenate: tmp = [cVector{:}];
%	fwrite(fp, cell2mat(tmp), 'float64');

fprintf('\t%d values written\n', nw);
fclose(fp);

%------------------------------------------------------------------------
%% as a check, read in binary data
%------------------------------------------------------------------------
fp = fopen('syndata.bin', 'rb');
M = fread(fp, 'double');
fclose(fp);
sM = size(M);
fprintf('read %d elements\n', numel(M));
fprintf('Size of data read from file: [%d, %d]\n', sM(1), sM(2));
M2 = reshape(M, [numel(M)/nC, nC]);
size(M2)

% plot the data, these should be 16 parallel lines
plot(M2)
