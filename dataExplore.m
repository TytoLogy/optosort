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
stimtype = {};
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

% set locations and size of controls
dataExploreLayout
% update plot
update_plot
% draw dots for stimulus onset and offset
draw_onsetoffset
% assign output
varargout{1} = fig;
%------------------------------------------------------------------------
%------------------------------------------------------------------------
%------------------------------------------------------------------------
%------------------------------------------------------------------------

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
                                    stimtype{trialN}, ...
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
      update_ui_str(stimulusText, ...
                           sprintf('%s  %s', testname, stimtype{trialN}));
      update_ui_str(stimuluslevelText, sprintf('%d', ...
                                                levels_by_trial(trialN)));
      update_ui_str(onsetoffsetText, sprintf('[%d %d]', ...
                                                   stimOnset, stimOffset));
   end

   %---------------------------------------------------------------------
   %---------------------------------------------------------------------
   function get_stimulus_parameters
      % build stimulus information string
      testname = nexInfo.FileInfo{1}.Dinf.test.Name;

      switch testname
         case 'BBN'
            stimtype = cell(ntrials, 1);
            for t = 1:ntrials
               stimtype{t} = 'BBN';
            end
            stimtype = nexInfo.FileInfo{1}.Dinf.test.stimcache.stimtype;
            stimOnset = nexInfo.FileInfo{1}.Dinf.audio.Delay;
            stimOffset = stimOnset + ...
                                 nexInfo.FileInfo{1}.Dinf.audio.Duration;
            levels_by_trial = ...
               cell2mat(nexInfo.FileInfo{1}.Dinf.test.stimcache.stimvar);
         case 'WAV'
            [stimtype, stimOnset, stimOffset, levels_by_trial] = ...
                                 get_stimulus_parameters_WAV;
      end
   end

   function [sName, sOn, sOff, sLevel] = get_stimulus_parameters_WAV
      % nexInfo.FileInfo{1}.Dinf.test.stimIndices(1:# of trials) 
      % holds indices into nexInfo.FileInfo{1}.Dinf.stimList struct that 
      % has stimulus information     
      stimIndices = nexInfo.FileInfo{1}.Dinf.test.stimIndices;
      nstim = length(stimIndices);
      
      sName = cell(nstim, 1);
      sLevel = zeros(nstim, 1);
      % for WAV, this will need to be obtained from stimulus struct
      % see APANalyze - the values in stimList are for TDT HW trigger and don't 
      % account for the offset within the .wav files 
      sOn = zeros(nstim, 2);
      sOff = zeros(nstim, 2);
      
      % init ANApath object
      AP = ANApath;
      %------------------------------------------------------------------------
      % stimulus onset/offset data from IC Paper Table 1 (with adjusted
      % onsets/durations)- saved by checkStimDurations script
      %------------------------------------------------------------------------
      table1Path = AP.ssPath;
      % table1FileName = 'ICPaper1_Table1.mat';
      table1FileName = 'ICPaper1_Table1_StimTiming_13Mar2023.mat';
      % load Table 1
      sendmsg('Loading stimulus timing information');
      load(fullfile(table1Path, table1FileName), 'Table1');
      %{
      Assume onset is 100 ms for everything
      For WAV stimuli:
      offset = onset + Table1.Durations_AdjStim(tI);
      %}
      
      % loop through the stimulus indices
      for s = 1:nstim
         % get audio portion of stimulus struct for this trial (don't need opto)
         S = nexInfo.FileInfo{1}.Dinf.stimList(stimIndices(s)).audio;
         % save stimulus name and stimulus onset/offset
         % need to account for 'null', 'noise' and 'wav'
         switch lower(S.signal.Type)
            case 'null'
               sName{s} = 'null';
               sOn(s) = 100;
               sOff(s) = 200;
            case 'noise'
               sName{s} = 'BBN';
               sOn(s) = S.Delay;
               sOff(s) = S.Delay + [0 S.Duration];
            case 'wav'
               % remove _adj
               sName{s} = deAdjWAVName(S.signal.WavFile);
               % find stimulus name in Table1.Syllable
               sindx = strcmp(sName{s}, Table1.Syllable);
               sOn(s) = 100;
               sOff(s) = 100 + Table1.Durations_AdjStim(sindx);
            otherwise
               error('Unknown stimulus type: %s', S.signal.Type);
         end
         % save stimulus level
         sLevel(s) = S.Level;
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