function testCommonReference(varargin)
%------------------------------------------------------------------------
% testCommonReference
%------------------------------------------------------------------------
% TytoLogy:Experiments:optosort
%------------------------------------------------------------------------
% script to test common avg/median reference for spike data
%------------------------------------------------------------------------
% See also: exportTest (script), export_for_plexon, common_avg_ref, 
%           common_med_ref
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%  Sharad J. Shanbhag
%	sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 12 April, 2023(SJS)
%
% Revisions:
%  3 May 2023 (SJS): implemented ref viewer, now need to create exported
%  test data where 1 channel is exported in 3 versions:
%     raw
%     common avg ref
%     common med ref
% 11 May 2023: functionalized
%------------------------------------------------------------------------

   %---------------------------------------------------------------------
   %---------------------------------------------------------------------
   %% specify file to examine
   % path, dat file, _testdata.mat file
   %---------------------------------------------------------------------
   %---------------------------------------------------------------------
   
   if ~isempty(varargin)
      if length(varargin) ~= 3
         error('need .dat file, ,nex output path, and channel as input');
      end
      % path and file for data file
      [dpath, dfile, dext] = fileparts(varargin{1});
      dname = [dfile dext];
      % _testdata.mat filename
      tname = [dfile '_testdata.mat'];
      % output (export) path 
      NexFilePath = varargin{2};
      % specify channel to export to test sorting file (.nex) for plexon
      exportChannel = varargin{3};
   else
      % get info from user
      % data file
      [dname, dpath] = uigetfile('*.dat', 'Select .dat file for analysis');
      if dname == 0
         sendmsg('Cancelling');
         return
      end
      [~, dfile, dext] = fileparts(dname);

      % _testdata.mat filename
      [~, dbase] = fileparts(dname);
      tname = [dbase '_testdata.mat'];

      % channel
      exportChannel = uiaskvalue('value', 1, ...
                        'valuetext', 'Channel (1-16)', ...
                        'valuetype', 'num', ...
                        'questiontext', 'Select Channel for Export', ...
                        'figurename', 'Export Channel'   );
      if isnan(exportChannel)
         sendmsg('Cancelling');
         return;
      elseif ~between(exportChannel, 1, 16)
         error('Invalid channel value - must be 1 through 16');
      else
      end

      % nexpath
      NexFilePath = uigetdir(dpath, 'path for exporting .nex file');
      if NexFilePath == 0
         sendmsg('Cancelling');
         return
      end
   
   end

   % define path to data file and return as OptoFileName object
   F = defineSampleData({dpath}, {dname}, {tname});

   % output (export) path and file
   NexFileName = sprintf('%s_CH%d_reftest.nex', F.base, exportChannel);

   sendmsg({'Reading from:', fullfile(dpath, dname)});
   sendmsg(sprintf('Exporting channel %d', exportChannel));
   sendmsg({'Writing to:', fullfile(NexFilePath, NexFileName)});
 
   
   %---------------------------------------------------------------------
   %% load data
   % get the data and information from the raw files
   %	cSweeps is a {nfiles, 1} cell array with each element holding an
   % 	{nChannels X ntrials} cell array of data for each file
   %	nexInfo is an instance of a SpikeInfo object
   %
   % use empty values for BPfilt (arg 3) and resampleData (arg 4) so 
   % that no filtering or resampling is done
   %---------------------------------------------------------------------
   % determine # of files
   nFiles = length(F);
   
   % specify channel to import
   Channels = 1:16;
   nChannels = length(Channels);
   
   BPfilt = buildFilter;
   resampleData = [];
   [cSweepsRaw, nexInfo] = read_data_for_export(F, Channels, ...
                                                   BPfilt, resampleData);
   
   %---------------------------------------------------------------------
   %% specify export file for plexon
   %---------------------------------------------------------------------   
   % some of this code is pulled from export_for_plexon
   nexInfo.FileName = fullfile(NexFilePath, NexFileName);
   % create output _nexinfo.mat file name - base is same as .nex file
   [~, baseName] = fileparts(nexInfo.FileName);
   nexInfo.InfoFileName = fullfile(NexFilePath, [baseName '_nexinfo.mat']);
   
   %---------------------------------------------------------------------
   %% apply common referencing before concatenating data - this should save on
   % memory and time... hopefully
   %---------------------------------------------------------------------
   cSweepsAvg = applyCommonReference(cSweepsRaw, @common_avg_ref);
   cSweepsMed = applyCommonReference(cSweepsRaw, @common_med_ref);
   
   %---------------------------------------------------------------------
   %---------------------------------------------------------------------
   %% export to nex file
   %-------------------------
   % using code from export_for_plexon()
   %---------------------------------------------------------------------
   %---------------------------------------------------------------------
   sendmsg('Adding continuous and event data to nex struct:');
   fprintf('Exporting data to %s\n', nexInfo.FileName);
   
   % write all events
   eventsToWrite = {'all'};
   
   % start new nex file data struct
   nD = nexCreateFileData(nexInfo.Fs);
   
   %---------------------------------------------------------------------
   % RAW channel data
   %---------------------------------------------------------------------
   % this will be a matrix of format
   % 	[# channels, (# sweeps) * (# samples per sweep)
   cVector = cell(1, nFiles);
   for f = 1:nFiles
	   cVector{1, f} = cSweepsRaw{f}(exportChannel, :);
   end
   % for each channels's data, concatenate cell array, convert to vector, 
   % add to nex struct
   % steps:
   %	concatenate: tmp = [cVector{:}];
   %	tmpVector = cell2mat(tmp);
   %  add to nex struct:
   %     [nexFile] = nexAddContinuous(nexFile, startTime, ...
   %                                    adFreq, values, name)
   nD = nexAddContinuous(nD, nexInfo.fileStartTime(1), nexInfo.Fs, ...
								   cell2mat([cVector{:}]), ...
								   sprintf('RAW_chan_%d', exportChannel));
   % clear cVector to save memory
   clear cVector
   
   %---------------------------------------------------------------------
   % AVG channel data
   %---------------------------------------------------------------------
   cVector = cell(1, nFiles);
   for f = 1:nFiles
	   cVector{1, f} = cSweepsAvg{f}(exportChannel, :);
   end
   nD = nexAddContinuous(nD, nexInfo.fileStartTime(1), nexInfo.Fs, ...
								   cell2mat([cVector{:}]), ...
								   sprintf('AVG_chan_%d', exportChannel));
   % clear cVector to save memory
   clear cVector
   
   %---------------------------------------------------------------------
   % MED channel data
   %---------------------------------------------------------------------
   cVector = cell(1, nFiles);
   for f = 1:nFiles
	   cVector{1, f} = cSweepsMed{f}(exportChannel, :);
   end
   nD = nexAddContinuous(nD, nexInfo.fileStartTime(1), nexInfo.Fs, ...
								   cell2mat([cVector{:}]), ...
								   sprintf('MED_chan_%d', exportChannel));
   % clear cVector to save memory
   clear cVector
   
   
   %---------------------------------------------------------------------
   % Add Events (aka timestamps, markers)
   %---------------------------------------------------------------------
   % add start sweep time stamps as event - assume consistent across
   % channels!
   if any(strcmpi('all', eventsToWrite) | ...
                  strcmpi('startsweep', eventsToWrite))
      %  [nexFile] = nexAddEvent( nexFile, timestamps, name )
      % events must be in column format...?
      nD = nexAddEvent(nD, ...
                     force_col(nexInfo.startTimeVector), 'startsweep');
   end
   
   % add end sweep time stamps as event - assume consistent across
   % channels! this is technically redundant, as startsweep event should be
   % 1 sample or time interval greater than endsweep. there is little
   % overhead involved in adding it so for now keep it here
   if any(strcmpi('all', eventsToWrite) | ...
                  strcmpi('endsweep', eventsToWrite))
      nD = nexAddEvent(nD, force_col(nexInfo.endTimeVector), 'endsweep');
   end
   
   % add file start times as single event type
   if any(strcmpi('all', eventsToWrite) | ...
                                    strcmpi('filestart', eventsToWrite))
      nD = nexAddEvent(nD, force_col(nexInfo.fileStartTime), 'filestart');
   end
   
   % add file end times as single event type
   if any(strcmpi('all', eventsToWrite) | ...
                              strcmpi('fileend', eventsToWrite))
      nD = nexAddEvent(nD, force_col(nexInfo.fileEndTime), 'fileend');
   end
   
   % Add individual event for each file with filename as event name
   if any(strcmpi('all', eventsToWrite) | ...
                                    strcmpi('filename', eventsToWrite))
      for f = 1:nFiles
         nD = nexAddEvent(nD, nexInfo.fileStartTime(f), ...
                           nexInfo.FileInfo{f}.F.base);
      end
   end
   
   % add stimulus onset ...
   if any(strcmpi('all', eventsToWrite) | ...
                                    strcmpi('stimstart', eventsToWrite))
      nD = nexAddEvent(nD, force_col(nexInfo.stimStartTimeVector), ...
                           'stimstart');
   end
   % ... and offset times
   if any(strcmpi('all', eventsToWrite) | ...
                                    strcmpi('stimend', eventsToWrite))
      nD = nexAddEvent(nD, ...
                        force_col(nexInfo.stimEndTimeVector), 'stimend');
   end
   
   % add stimulus-specific onset times
   if any(strcmpi('all', eventsToWrite) | ...
                                       strcmpi('stim_id', eventsToWrite))
      for f = 1:nFiles
         events = nexInfo.stimEventTimesForFile(f);
         fprintf('Adding events from file %s\n', ...
                           nexInfo.FileInfo{f}.F.file);
         for n = 1:length(events)
            fprintf('\t%s\n', events(n).name);
            nD = nexAddEvent(nD, force_col(events(n).timestamps), ...
                                 events(n).name);
         end
      end
   end
   
   %---------------------------------------------------------------------
   % write to nexfile
   %---------------------------------------------------------------------
   sendmsg(sprintf('Writing nex file %s:', nexInfo.FileName));
   writeNexFile(nD, nexInfo.FileName);
   
   %---------------------------------------------------------------------
   % write useful information to _nexinfo.mat file 
   %---------------------------------------------------------------------
   % save to matfile
   sendmsg(sprintf('Writing _nexinfo.mat file %s:', ...
											   nexInfo.InfoFileName));
   save(nexInfo.InfoFileName, 'nexInfo', '-MAT');

end

%------------------------------------------------------------------------
%------------------------------------------------------------------------
%------------------------------------------------------------------------
% Nested Functions
%------------------------------------------------------------------------
%------------------------------------------------------------------------
%------------------------------------------------------------------------

function varargout = plotRawAvgMedRef(figH, raw, avg, med, Fs)

   % transform cell array to matrix with time bins in columns, channels by
   % rows
   T = cell2mat(raw);
   % apply common average referencing to data
   Ta = cell2mat(avg);
   % apply common median referencing to data
   Tm = cell2mat(med);
   
   % Plot data
   if isempty(figH)
      figH = figure;
   else
      figure(figH);
   end

   tl = tiledlayout(1, 3);
   tl.TileSpacing = 'compact';
   tl.Padding = 'compact';
   
   % plot raw data
   nexttile
   [r.pH, r.aX, r.fH] = tracePlot(T', Fs);
   title(r.aX, 'RAW');
%    % draw onset/offset lines
%    line([onset onset], r.aX.YLim, 'Color', 'g');
%    line([offset offset], r.aX.YLim, 'Color', 'r');
   
   % plot avg data
   nexttile
   [a.pH, a.aX, a.fH] = tracePlot(Ta', Fs);
   title(a.aX, 'AVG');
%    line([onset onset], r.aX.YLim, 'Color', 'g');
%    line([offset offset], r.aX.YLim, 'Color', 'r');
   % plot med data
   nexttile;
   [m.pH, m.aX, m.fH] = tracePlot(Tm', Fs);
   title(m.aX, 'MED');
%    line([onset onset], r.aX.YLim, 'Color', 'g');
%    line([offset offset], r.aX.YLim, 'Color', 'r');
   
   varargout{1} = figH;
end

% plot traces of data
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
   A = mean(Tin);
   Tout = Tin - A;
end

% apply common median reference to matrix of channel data [samples,
% channels]
function Tout = common_med_ref(Tin)
   A = median(Tin);
   Tout = Tin - A;
end

function out = applyCommonReference(sweepCell, refFunctionHandle)
   % make a copy of input cell
   out = sweepCell;
   nFiles = length(sweepCell);
   % loop through files
   for f = 1:nFiles
      nSweeps = size(sweepCell{f}, 2);
      % loop through sweeps (aka trials)
      for s = 1:nSweeps
         % channels are in rows, trials in columns of the cell array 
         % convert to mat (channels, samples)
         m = cell2mat(sweepCell{f}(:, s));
         % apply function
%          a = refFunctionHandle(m);
         % convert to cell (channels, 1) and assign to cSweepsAvg
         out{f}(:, s) = mat2cell(refFunctionHandle(m), ...
                                       ones(1, size(m, 1)));
      end
   end
end

function BPfilt = buildFilter(varargin)
   %---------------------------------------------------------------------
   % define filter parameters
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
   
   %---------------------------------------------------------------------
   % if sample rate was provided, generate coefficients
   %---------------------------------------------------------------------
   if ~isempty(varargin)
      Fs = varargin{1};
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
	      error('%s: unsupported filter type %s', mfilename, BPfilt.type)
      end
   end
end