function testplot(cSweeps, Fs)
% plot and step through trials in data - uses programmatic GUI 
%
% 27 Apr 2023: creatd dataExplore.m to try different approach using new app
% building techniques in MATLAB

if ~iscell(cSweeps)
   error
end
fileN = 1;
nFiles = length(cSweeps);
[nchannels, ntrials] = size(cSweeps{fileN});
channels = 1:nchannels;
trialN = 1;

% get a single trial data in a matrix
data = cell2mat(cSweeps{fileN}(:, trialN))';

L = size(data, 1);

% x (time) axis vector
T = (1:L) / Fs;

h_fig = figure;

% need to generate y data based on min max  and # of channels
yLim = [min(data(:)) max(data(:))];
yInterv = yLim(2)-yLim(1);
nChan = length(channels);
yTotal = yInterv * nChan;

yStart = (nChan:-1:1)*yInterv-(yInterv/2);

dataWin = bsxfun(@plus, data, yStart);

% create axes
h_ax = axes(h_fig, 'Position',[0.1 0.125 0.8 0.75]);
plot(h_ax, T, dataWin);
set(h_ax,'ytick',yStart(end:-1:1),'yticklabel',channels(end:-1:1))
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