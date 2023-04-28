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

% create figure
fig = figure;
fig.Name = 'dataExplore';
fig.Units = 'Normalized';
fig.Position = [0.2236    0.7151    0.1626    0.2229];


% create axes
ax = axes(fig);

% set up button controls

% place trial buttons in panel
trialPanel = uipanel(fig, 'Title', 'Trial', 'Tag', 'trialPanel');

% prev trial button
trialDown = uicontrol(trialPanel, 'Style', 'pushbutton', ...
                        'Tag', 'trialDown', 'String', '<<');
trialDown.Units = 'normalized';
trialDown.Callback = @change_trial;
% next trial button
trialUp = uicontrol(trialPanel, 'Style', 'pushbutton', ...
                        'Tag', 'trialUp', 'String', '>>');
trialUp.Units = 'normalized';
trialUp.Callback = @change_trial;
% trial indicator control
trialNum = uicontrol(trialPanel, 'Style', 'edit', ...
                        'Tag', 'trialNum', 'String', '1');
trialNum.Units = 'normalized';
trialNum.Callback = @change_trial;


dataExploreLayout

update_plot


%% Nested Functions

   function update_plot
      % need to generate y data based on min max and # of channels
      trialdata = cell2mat(cSweeps{fileN}(:, trialN))';
      yLim = [min(trialdata(:)) max(trialdata(:))];
      yInterv = yLim(2)-yLim(1);
      nChan = length(channels);
      yTotal = yInterv * nChan;
      yStart = (nChan:-1:1)*yInterv-(yInterv/2);
      % scale data appropriately
      data = bsxfun(@plus, trialdata, yStart);
      
      % plot data
      plot(ax, T, data);
      % set yticks
      set(ax,'ytick',yStart(end:-1:1),'yticklabel',channels(end:-1:1))
      % adjust limits
      ylim([0 yTotal]);
      xlim([T(1) T(end)]);      
   end


   function change_trial(hObject, eventdata, handles)
      switch(hObject.Tag)
         case 'trialDown'
            if trialN == 1
               beep
               sendmsg('At first trial');
            else
               trialN = trialN - 1;
               update_ui_str(trialNum, trialN);
               update_plot
            end
         case 'trialUp'
            if trialN == ntrials
               beep
               sendmsg('At last trial');
            else
               trialN = trialN + 1;
               update_ui_str(trialNum, trialN);
               update_plot
            end

         case 'trialNum'
            % get the value from string
            tmp = read_ui_str(hObject, 'n');
            if ~between(tmp, 1, ntrials)
               beep
               sendmsg('Trial # out of bounds');
               update_ui_str(trialNum, trialN);
            else
               trialN = tmp;
               update_plot
            end
      end
   end

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