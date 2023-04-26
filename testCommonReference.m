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
- idea: 
   select a channel
   export three "channels" of data
      raw
      c.a.r.
      c.m.r.
   then sort and compare
%}

%------------------------------------------------------------------------
%% specify file to examine
% path, dat file, _testdata.mat file
%------------------------------------------------------------------------
dpath = '/media/sshanbhag/SSData/Data/Mouse/Raw/BLA/1500/20230213';
dpath = '/Volumes/SSData/Data/Mouse/Raw/BLA/1500/20230213';
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
% get stimulus information
onset = nexInfo.FileInfo{1}.Dinf.audio.Delay;
offset = onset + nexInfo.FileInfo{1}.Dinf.audio.Duration;
%%
%{
% set up raw data plot
[r.fH, ar.X, r.pH] = tracePlot([], cell2mat(cSweeps{1}(:, 1))', nexInfo.Fs);
r.fH.Name = 'RAW';
% set up avg data plot
[a.fH, a.aX, a.pH] = tracePlot([], cell2mat(cSweeps{1}(:, 1))', nexInfo.Fs);
a.fH.Name = 'AVG';
% set up med data plot
[m.fH, m.aX, m.pH] = tracePlot([], cell2mat(cSweeps{1}(:, 1))', nexInfo.Fs);
m.fH.Name = 'MED';
%}

% define filter
BPfilt = buildFilter(nexInfo.Fs);

%% process data

% select trial
trialN = 100;
D = cSweeps{1}(:, trialN);

% filter the data
for c = 1:length(D)
   D{c} = filtfilt(BPfilt.b, BPfilt.a, ...
							sin2array(D{c}, nexInfo.Fs, BPfilt.ramp));
end

% transform cell array to matrix with time bins in rows, channels by column
T = cell2mat(D)';
% apply common average referencing to data
Ta = common_avg_ref(T);
% apply common median referencing to data
Tm = common_med_ref(T);

% Plot data
figure(1)
tl = tiledlayout(1, 3);
tl.TileSpacing = 'compact';
tl.Padding = 'compact';

% plot raw data
nexttile
[r.pH, r.aX, r.fH] = tracePlot(T, nexInfo.Fs);
title(r.aX, 'RAW');
% draw onset/offset lines
line([onset onset], r.aX.YLim, 'Color', 'g');
line([offset offset], r.aX.YLim, 'Color', 'r');

% plot avg data
nexttile
[a.pH, a.aX, a.fH] = tracePlot(Ta, nexInfo.Fs);
title(a.aX, 'AVG');
line([onset onset], r.aX.YLim, 'Color', 'g');
line([offset offset], r.aX.YLim, 'Color', 'r');
% plot med data
nexttile;
[m.pH, m.aX, m.fH] = tracePlot(Tm, nexInfo.Fs);
title(m.aX, 'MED');
line([onset onset], r.aX.YLim, 'Color', 'g');
line([offset offset], r.aX.YLim, 'Color', 'r');

r.fH.Name = sprintf('%s_Sweep%d', F.base, trialN);
r.fH.FileName = sprintf('%s_Sweep%d', F.base, trialN);
r.fH.PaperOrientation = 'landscape';
%{
   % assign data to plot
   for c = 1:nchan
      % update plot
      set(pH(c), 'YData', T(:, c) + c*yabsmax);
   end
   % update plots
   refreshdata(ancestor(pH(1), 'figure'));
end
%}



%------------------------------------------------------------------------
%------------------------------------------------------------------------
%------------------------------------------------------------------------
% Nested Functions
%------------------------------------------------------------------------
%------------------------------------------------------------------------
%------------------------------------------------------------------------

function varargout = tracePlot(T, Fs, varargin)
   % figure and axes background grey level ( 0 = black, 1 = white)
   bgcol = 1;
   % max range for each trace
   yabsmax = 5;
   
   % size of T
   [npts, nchan] = size(T);
   
   % time vector for xdata
   tv = timevec(npts, Fs, 'ms');
   % generate shifted channel data
   Ts = zeros(size(T));
   for c = 1:nchan
      Ts(:, c) = c*(yabsmax) + T(:, c);
   end
   % and plot it, storing handles to lines
   pH = plot(tv, Ts, 'k-');
   % get handles to axes and figure
   aX = gca;
   fH = gcf;
   
   % rename/scale ticks
   yticks_yvals = yabsmax*(1:nchan);
   yticks_txt = cell(nchan, 1);
   for n = 1:nchan
      yticks_txt{n} = num2str(n);
   end
   set(aX, 'YTick', yticks_yvals);
   set(aX, 'YTickLabel', yticks_txt);
   set(aX, 'TickDir', 'out');
   
   set(aX, 'Box', 'off');
   
   % adjust x, y limits
   xlim([0 ceil(max(pH(1).XData))]);
   ylim(yabsmax*[0 nchan+1]);
   
   grid on

   % labels for axes
   xlabel('Time (ms)')
   ylabel('Channel')
   % set bg color
   set(aX, 'Color', bgcol*[1 1 1]);
   set(fH, 'Color', bgcol*[1 1 1]);
   % turn off toolbar
   % set(fH, 'ToolBar', 'none');
   
   varargout{1} = pH;
   varargout{2} = aX;
   varargout{3} = fH;
end


% apply common average reference to matrix of channel data [samples,
% channels]
function Tout = common_avg_ref(Tin)
   A = mean(Tin, 2);
   Tout = Tin - A;
end

% apply common median reference to matrix of channel data [samples,
% channels]
function Tout = common_med_ref(Tin)
   A = median(Tin, 2);
   Tout = Tin - A;
end


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