function varargout = dataExplore(cSweeps, nexInfo)
%------------------------------------------------------------------------
%------------------------------------------------------------------------
% used to examine data from opto program.
%
%	cSweeps     {nfiles, 1} cell array with each element holding an
% 	            {nChannels X ntrials} cell array of data for each file
%	nexInfo     instance of a SpikeInfo object
%------------------------------------------------------------------------
%------------------------------------------------------------------------

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


% reference type
refmode = 'RAW';

% variables for stimulus settings
stimtype = '';
testname = '';
stimOnset = [];
stimOffset = [];
levels_by_trial = [];
get_stimulus_parameters;

% get sample rate (samples/s)
Fs = 0.001*nexInfo.Fs;


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


% filename text
filenameText = uicontrol(fig, 'Style', 'Text', ...
                           'String', nexInfo.FileInfo{1}.F.file, ...
                           'FontWeight', 'Bold', ...
                           'HorizontalAlignment', 'center');
% set up trial controls
% place trial buttons in panel
trialPanel = uipanel(fig, 'Title', sprintf('Trial (1 - %d)', ntrials), ...
                           'Tag', 'trialPanel');
% prev trial button
trialDownButton = uicontrol(trialPanel, 'Style', 'pushbutton', ...
                        'Tag', 'trialDownButton', 'String', '<<');
trialDownButton.Callback = @change_trial;
% next trial button
trialUpButton = uicontrol(trialPanel, 'Style', 'pushbutton', ...
                        'Tag', 'trialUpButton', 'String', '>>');
trialUpButton.Callback = @change_trial;
% trial indicator control
trialnEdit = uicontrol(trialPanel, 'Style', 'edit', ...
                        'Tag', 'trialnEdit', 'String', '1');
trialnEdit.Callback = @change_trial;


% set up reference radio buttons
% create ref button group
refButtonGroup = uibuttongroup(fig, 'SelectionChangedFcn', @change_ref, ...
                                 'Title', 'Common Reference');
rawRadioButton = uicontrol(refButtonGroup, 'Style', 'radiobutton', ...
                              'String', 'Raw', ...
                              'Tag', 'RAW', ...
                              'HandleVisibility', 'off');
avgRadioButton = uicontrol(refButtonGroup, 'Style', 'radiobutton', ...
                              'String', 'Average', ...
                              'Tag', 'AVG', ...
                              'HandleVisibility', 'off');
medRadioButton = uicontrol(refButtonGroup, 'Style', 'radiobutton', ...
                              'String', 'Median', ...
                              'Tag', 'MED', ...
                              'HandleVisibility', 'off');


% panel with stimulus information
stimInfoPanel = uipanel(fig, 'Title', 'Stimulus Information', ...
                           'Tag', 'stimInfoPanel');
stimulusText = uicontrol(stimInfoPanel, 'Style', 'Text', ...
                           'String', 'stim', ...
                           'HorizontalAlignment', 'center');
stimuluslevelText = uicontrol(stimInfoPanel, 'Style', 'Text', ...
                                    'String', 'level', ...
                                    'HorizontalAlignment', 'center');
onsetoffsetText = uicontrol(stimInfoPanel, 'Style', 'Text', ...
                                    'String', '[onset offset]', ...
                                    'HorizontalAlignment', 'center');



dataExploreLayout

update_plot
draw_onsetoffset

varargout{1} = fig;

%------------------------------------------------------------------------
%------------------------------------------------------------------------
%% Nested Functions
%------------------------------------------------------------------------
%------------------------------------------------------------------------

   %---------------------------------------------------------------------
   %---------------------------------------------------------------------
   function change_ref(hObject, eventdata, handles)
      refmode = hObject.SelectedObject.Tag;
      update_plot;
   end

   %---------------------------------------------------------------------
   %---------------------------------------------------------------------
   function apply_reference
      data = cell2mat(cSweeps{fileN}(:, trialN))';
      switch refmode
         case 'RAW'
            trialdata = data;
         case 'AVG'
            trialdata = common_avg_ref(data);
         case 'MED'
            trialdata = common_med_ref(data);
      end
   end

   %---------------------------------------------------------------------
   %---------------------------------------------------------------------
   function update_plot
      % need to generate y data based on min max and # of channels
      apply_reference;
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
      set(ax, 'ytick', yStart(end:-1:1), 'yticklabel', channels(end:-1:1))
      % adjust limits
      ylim([0 yTotal]);
      xlim([T(1) T(end)]);
      xlabel(ax, 'Time (ms)');
      ylabel(ax, 'Channel');
      draw_onsetoffset
      update_stimulus_info
      % change figure name - eases saving of figure as png, or jpg
      figname = sprintf('%s_%s_t%d_%ddB_%s', ...
                                    nexInfo.FileInfo{1}.F.base, ...
                                    stimtype, ...
                                    trialN, ...
                                    levels_by_trial(trialN), ...
                                    refmode);
      set(fig, 'Name', figname);
      set(fig, 'FileName', figname);
   end


   %---------------------------------------------------------------------
   %---------------------------------------------------------------------
   function change_trial(hObject, eventdata, handles)
      switch(hObject.Tag)
         case 'trialDownButton'
            if trialN == 1
               beep
               sendmsg('At first trial');
            else
               trialN = trialN - 1;
               update_ui_str(trialnEdit, trialN);
               update_plot
            end
         case 'trialUpButton'
            if trialN == ntrials
               beep
               sendmsg('At last trial');
            else
               trialN = trialN + 1;
               update_ui_str(trialnEdit, trialN);
               update_plot
            end

         case 'trialnEdit'
            % get the value from string
            tmp = read_ui_str(hObject, 'n');
            if ~between(tmp, 1, ntrials)
               beep
               sendmsg('Trial # out of bounds');
               update_ui_str(trialnEdit, trialN);
            else
               trialN = tmp;
               update_plot
            end
      end

   end

   function update_stimulus_info
      update_ui_str(stimulusText, sprintf('%s  %s', testname, stimtype));
      update_ui_str(stimuluslevelText, sprintf('%d', ...
                                                levels_by_trial(trialN)));
      update_ui_str(onsetoffsetText, sprintf('[%d %d]', ...
                                                   stimOnset, stimOffset));
   end

   %---------------------------------------------------------------------
   %---------------------------------------------------------------------
   function get_stimulus_parameters
      stimtype = nexInfo.FileInfo{1}.Dinf.test.stimcache.stimtype;

      % build stimulus information string
      testname = nexInfo.FileInfo{1}.Dinf.test.Name;

      switch testname
         case 'BBN'
            stimOnset = nexInfo.FileInfo{1}.Dinf.audio.Delay;
            stimOffset = stimOnset + ...
                                 nexInfo.FileInfo{1}.Dinf.audio.Duration;
            levels_by_trial = ...
               cell2mat(nexInfo.FileInfo{1}.Dinf.test.stimcache.stimvar);
         case WAV
            
            levels_by_trial = ...
               cell2mat(nexInfo.FileInfo{1}.Dinf.test.stimcache.stimvar);
            

      end
   end

   %---------------------------------------------------------------------
   %---------------------------------------------------------------------
   function draw_onsetoffset
      hold on
      plot(ax, stimOnset, ax.YLim(1), 'g.', 'MarkerSize', 12);
      plot(ax, stimOffset, ax.YLim(1), 'r.', 'MarkerSize', 12);
      hold off
   end



end