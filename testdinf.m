DataPath = '~/Work/Data/TestData/MT';


%% FREQ file info
DataFile = '1382_20191212_02_02_3200_FREQ_TUNING.dat';
% use readOptoData to read in raw data. 
[~, tmpDinf] = readOptoData(fullfile(DataPath, DataFile));
% Fix test info
freqDinf = correctTestType(tmpDinf);

%% FRA file info
DataFile = '1382_20191212_02_02_3200_FRA.dat';
% use readOptoData to read in raw data. 
[~, tmpDinf] = readOptoData(fullfile(DataPath, DataFile));
% Fix test info
fraDinf = correctTestType(tmpDinf);

%% BBN file info
DataFile = '1382_20191212_02_02_3200_BBN.dat';
% use readOptoData to read in raw data. 
[~, tmpDinf] = readOptoData(fullfile(DataPath, DataFile));
% Fix test info
bbnDinf = correctTestType(tmpDinf);

%% WAV file info
DataFile = '1382_20191212_02_02_3200_WAV.dat';
WAVInfoFile = '1382_20191212_02_02_3200_WAV_wavinfo.mat';
% use readOptoData to read in raw data. 
[~, tmpDinf] = readOptoData(fullfile(DataPath, DataFile));
% Fix test info
wavDinf = correctTestType(tmpDinf);
tmpWinfo = load(fullfile('~/Work/Data/TestData/MT', '1382_20191212_02_02_3200_WAV_wavinfo.mat'))

WI = WAVInfo(wavDinf, tmpWinfo) 

