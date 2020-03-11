
% data locations

%{
% Data for 1382_20191212_02_02_3200 has good recordings on 
% Channels 4, 5, 11, 14
rawPath = '~/Work/Data/TestData/MT';
rawFile = '1382_20191212_02_02_3200_FREQ_TUNING.dat';
Channels = [4 5];
%}

%
% Test WAV data
% 1382_20191212_02_02_3200_WAV.dat
% 1382_20191212_02_02_3200_WAV_PSTH.fig
% 1382_20191212_02_02_3200_WAV_wavinfo.mat
% Channels 4, 5, 11, 14
rawPath = '~/Work/Data/TestData/MT';
rawFile = '1382_20191212_02_02_3200_WAV.dat';
Channels = [4 5];
%{
% data for FRA - this is correct data to use for testing as it was
collected using updated FRA routine in opto program

rawPath = '/Volumes/Wenstrup Laboratory/By User/SJS/Data/SpikeSort/IC-probe';
rawFile = '1407_20200305_01_01_550_FRA.dat';
Channels = 8;
%}

% create opto file name object (helps to create new file names for these
% data)

fI = OptoFileName(fullfile(rawPath, rawFile));

%% load data

tmpF = fullfile('~/Work/Data/TestData/MT', ...
			fI.newname(  sprintf('Chan%d-%d_TracesByStim', ...
														Channels(1), Channels(end)), ...
 							'mat') );
load(tmpF);

%% plot data in explorable form

stimNames = unique(Dinf.test.wavlist, 'stable');

tracesByStim = Traces{1};
spikesByStim = spikes{1};

stimNum = 1;
sweepNum = 1;
stimTraces = tracesByStim{stimNum};
sweepData = stimTraces(:, sweepNum)';

[nPts, nTraces] = size(tracesByStim{currentStim});
fH = figure('Units','Normalized','Position',[0.25 0.25 0.5 0.5]);
aH =   axes('Units','Normalized','Position',[0.05 0.15, 0.75 0.75]);
stimCtrlH = uicontrol(fH, 'Style','listbox','Units','Normalized', ...
							'Position',[0.85 0.15, 0.1, 0.75],...
							'String', stimNames, ...
							 'Callback',{@changeStim, aH, sweepNum, stimTraces} );



%%
%{
function changeChannel(l,evtData,a,s,data)
	cla(a);
	chanNum = str2double(l.String{l.Value});
	%500Hz
	sR = 500;  
	%Reshape each epoch into a column
	tempData = reshape(data(chanNum,:,:),[],size(data,3)); 
	%Build time array
	tempTime = [0:1/sR:(size(data,2)-1)/sR]' + (0:1:size(data,3)-1)*2; 
	%plot all the lines
	plot(a,tempTime,tempData)
	%Rest Slider Position
	s.Value = 1; 
end
%}