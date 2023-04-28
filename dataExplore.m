function dataExplore(cSweeps, nexInfo)
% used to examine data from opto program.
%
%	cSweeps     {nfiles, 1} cell array with each element holding an
% 	            {nChannels X ntrials} cell array of data for each file
%	nexInfo     instance of a SpikeInfo object

% check input
if nargin ~= 2
   help dataExplore
   error('%s: needs cSweeps cell array and nexInfo object', mfilename);
end

% get number of files in cSweeps and set file index number
nFiles = length(cSweeps);
fileN = 1;

% get number of channels, number of trials for current file, set 
% current trial to 1;
[nchannels, ntrials] = size(cSweeps{fileN});
channels = 1:nchannels;
trialN = 1;

% get sample rate (samples/s)
Fs = nexInfo.Fs;


% to initialize plot/gui, get a single trial data in a matrix
trialdata = cell2mat(cSweeps{fileN}(:, trialN))';

% length of trial data
L = size(trialdata, 1);
% x (time) axis vector
T = (1:L) / Fs;

% create 
fig = uifigure;
fig.Name = 'dataExplore';

% use grid layout to set format of app window
gl = uigridlayout(fig);
% set row and column heights
%     '1x' specifies variable width, equal to size of first column
gl.RowHeight = {22, 22, '1x'};
gl.ColumnWidth = {150, '2x'};

% set up button controls

% previous trial button
trialUp = uibutton(

% set up axes
ax = uiaxes(gl);
ax.Layout.Row = [1, 3];
ax.Layout.Column = 2;

return

% need to generate y data based on min max and # of channels
yLim = [min(trialdata(:)) max(trialdata(:))];
yInterv = yLim(2)-yLim(1);
nChan = length(channels);
yTotal = yInterv * nChan;
yStart = (nChan:-1:1)*yInterv-(yInterv/2);
% scale data appropriately
data = bsxfun(@plus, trialdata, yStart);

% create axes
h_ax = axes(fig, 'Position',[0.1 0.125 0.8 0.75]);

% plot data
plot(h_ax, T, data);
% set yticks
set(h_ax,'ytick',yStart(end:-1:1),'yticklabel',channels(end:-1:1))
% adjust limits
ylim([0 yTotal]);
xlim([T(1) T(end)]);

h_tUp = uicontrol('Style','pushbutton','Parent',h_fig,...
      'Units','normalized','Position',[0.01 0.01 0.05 0.05],...
      'Value',0,'Callback',{@next_trial});


%% Nested Functions

function next_trial(hObject, eventdata, handles)
   if trialN + 1 > ntrials
      warning('at last trial');
   else
      trialN = trialN + 1;
   end
end


function prev_trial(hObject, eventdata, handles)
   if trialN - 1 < 1
      warning('at first trial');
   else
      trialN = trialN - 1;
   end
end






end