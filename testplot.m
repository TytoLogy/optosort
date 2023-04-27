function testplot(cSweeps, Fs)

if ~iscell(cSweeps)
   error
end
fileN = 1;
nFiles = length(cSweeps);
[nchannels, ntrials] = size(cSweeps{fileN});
trialN = 1;

% get a single trial data in a matrix
data = cell2mat(cSweeps{fNum}(:, trialN))';

L = size(data, 1);

% x (time) axis vector
T = (1:L) / srate;

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

%% Nested Functions

function next_trial
   if trialN + 1 > ntrials
      warning('at last trial');
   else
      trialN = trialN + 1;
   end
end


function prev_trial
   if trialN - 1 < 1
      warning('at first trial');
   else
      trialN = trialN - 1;
   end
end






end