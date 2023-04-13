%------------------------------------------------------------------------
% testCommonReference
%------------------------------------------------------------------------
% TytoLogy:Experiments:optosort
%------------------------------------------------------------------------
% 
% script to test common avg/median reference for spike data
%
% 
%------------------------------------------------------------------------
% See also: exportTest (script)
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%  Sharad J. Shanbhag
%	sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 12 April, 2023(SJS)
%
% Revisions:
%------------------------------------------------------------------------

%{
Things to do:

1) specify and load data
2) apply referencing 
3) display and compare pre/post data across channels
questions:
- how to exclude channels?
%}

%------------------------------------------------------------------------
%% specify file to examine
% path, dat file, _testdata.mat file
%------------------------------------------------------------------------
dpath = '/media/sshanbhag/SSData/Data/Mouse/Raw/BLA/1500/20230213';
dname = '1500_20230213_01_0_3352_BBN.dat';
tname = '1500_20230213_01_0_3352_BBN_testdata.mat';

% define path to data file and 
F = defineSampleData({dpath}, {dname}, {tname});

%------------------------------------------------------------------------
%% load data
% get the data and information from the raw files
%	cSweeps is a {nfiles, 1} cell array with each element holding an
% 	{nChannels X ntrials} cell array of data for each file
%	nexInfo is an instance of a SpikeInfo object
%
% use empty values for BPfilt (arg 3) and resampleData (arg 4) so 
% that no filtering or resampling is done
%------------------------------------------------------------------------
Channels = 1:16;
nChannels = length(Channels);
BPfilt = [];
resampleData = [];
[cSweeps, nexInfo] = read_data_for_export(F, Channels, ...
                                                [], []);
%%
% set up raw data plot
[r.fH, ar.X, r.pH] = tracePlot([], cell2mat(cSweeps{1}(:, 1))', nexInfo.Fs);
r.fH.Name = 'RAW';
% set up avg data plot
[a.fH, a.aX, a.pH] = tracePlot([], cell2mat(cSweeps{1}(:, 1))', nexInfo.Fs);
a.fH.Name = 'AVG';
% set up med data plot
[m.fH, m.aX, m.pH] = tracePlot([], cell2mat(cSweeps{1}(:, 1))', nexInfo.Fs);
m.fH.Name = 'MED';

%% plot raw data



trialN = 13;

D = cSweeps{1}(:, trialN);

% define filter
BPfilt = buildFilter(nexInfo.Fs);
% filter the data
for c = 1:length(D)
   D{c} = filtfilt(BPfilt.b, BPfilt.a, ...
							sin2array(D{c}, nexInfo.Fs, BPfilt.ramp));
end

% get data for one trial, all channels (rows)
T = cell2mat(D)';

% plot data

tracePlot(r.pH, T, nexInfo.Fs)


Ta = common_avg_ref(T);
tracePlot(a.pH, Ta, nexInfo.Fs)
Tm = common_med_ref(T);
tracePlot(m.pH, Tm, nexInfo.Fs)





function Tout = common_avg_ref(Tin)
   A = mean(Tin, 2);
   Tout = Tin - A;
end


function Tout = common_med_ref(Tin)
   A = median(Tin, 2);
   Tout = Tin - A;
end


%------------------------------------------------------------------------
%------------------------------------------------------------------------
%------------------------------------------------------------------------
% Nested Functions
%------------------------------------------------------------------------
%------------------------------------------------------------------------
%------------------------------------------------------------------------


















function BPfilt = buildFilter(Fs)
   %---------------------------------------------------------------------
   % define filter
   %---------------------------------------------------------------------
   % [highpass lowpass] cutoff frequencies in Hz
   BPfilt.Fc = [250 4000];
   % order of filter. note that the filtfilt() function in MATLAB is used,
   % so the effective order is doubled. typically use 5
   BPfilt.forder = 5;
   % ramp time (ms) to apply to each sweep in order to cutdown on
   % onset/offset transients from filtering
   BPfilt.ramp = 1;
   % filter type. 'bessel' or 'butter' (for butterworth). typically use
   % butter as bessel is low-pass only and we need to remove low frequency
   % crap from the raw data
   BPfilt.type = 'butter';
   
   BPfilt.Fs = Fs;
   BPfilt.Fnyq = Fs / 2;
   BPfilt.cutoff = BPfilt.Fc / BPfilt.Fnyq;
   if strcmpi(BPfilt.type, 'bessel')
	   [BPfilt.b, BPfilt.a] = besself(BPfilt.forder, BPfilt.cutoff, ...
                                          'bandpass');
   elseif strcmpi(BPfilt.type, 'butter')
	   [BPfilt.b, BPfilt.a] = butter(BPfilt.forder, BPfilt.cutoff, ...
                                          'bandpass');
   else
	   error('%s: unknown filter type %s', mfilename, BPfilt.type)
   end

end