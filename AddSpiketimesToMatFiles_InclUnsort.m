% AddSpiketimesToMatFiles.m
% 
% Prior to this:
% 1) Use ExportTytologyDataToNex.m to concatenate all .mat files for each paradigm into one big .nex file for offline sorting.
% 2) Use Plexon Offline Sorter to create a concatenated .mat file containing sorted spikes.
% 3) Run this program separately for each cell/unit:
%    In the activepath folder, place:
%    - .mat file with Plexon OfflineSorter sorted spike times, CONCATENATED for all paradigms.   Format: animalID_cellID.mat
%    - The _cellinfo.mat file created during concatenation by ExportTytologyDataToNex.m.         Format: animalID_cellID_cellinfo.mat
%      NB: You can load the _cellinfo.mat file to see which files were concatenated, that is:
%    - All .mat and .dat files for each paradigm, collected via Tytology or PresentStimCurve.    (Listed in the struct cellinfo.curvetype.)
%
% HPSearch and PresentStimCurve both use Tytology scripts and save data similarly.
% .mat files contain structs curvedata and curveinfo
% .dat files contain raw data for each trial
%
% 1) Read in .mat files exported by Plexon's Offline Sorter
% Column 1: unit number (where 0 is unsorted)
% Column 2: timestamp where spike crosses threshold (in seconds)
% Columns 3-34 (assuming waveform window of 1311us / 32 samples):
%   waveform snippet, with or without prewindow as set in Offline Sorter
%   (prewindow default: 494us / 12 samples)
%   (window default: 1311us / 32 samples)
%   This is in units of samples/sec of raw data file 
%   (24414.063 Hz based on settings in data acquisition program 
%        HPSearch or PresentStimCurve in RosenLab)
% 2) Load matching raw data file. Error check to ensure correct match.
% 3) Map spike times to trials.
% 4) Scroll through file trial by trial. User can notate which trials need
%    to be excluded. Later script will exclude those trials (if we decide
%    that's necessary).
% 5) Append sorted spikes from Plexon's Offline Sorter to the raw .mat 
%    data file (e.g. M331_020_041013-16-22_RISEFALL.mat)
%

% TO DO:
% exclude trials (not essential at the moment)
% But at least exclude all spikes that are invalid (unitID = -1)
% process data (start by plotting each paradigm)


%% Load .mat file exported by Offline Sorter. 
% Loads matrix adc001.

clear all; close all;

% processpath = ['/Users/Merri/Documents/Matlab_work/Gerbil/_TytoLogy/'];
% activepath = ['/Users/Merri/Documents/Gerbil/Data&Analysis/headfixed_physiology/runme/'];
% activepath = ['/Users/merri/Dropbox/MjR_Ethologeek/EXPERIMENTS/Headfixed_extracell_Tytology/Data_FilesForMR/April2016/F415-004/'];
processpath = ['C:\Matlab_work\Tytology\'];
% activepath = ['C:\Experiments\Data\Headfixed_extracell_Tytology\Data_Files\CTR_Adult\M1019_CTR_Adult\'];
activepath = ['C:\Experiments\Data\Headfixed_extracell_Tytology\runme\'];
% activepath = ['C:\Ephys\Ctrl_Juv_Urethane\F938_CTR_Juv_urethane\004\'];
cd(activepath)

% folders = dir;
% flist = []; % this will become a list of valid folders to process (those that include animal_cell.mat files)
% for i = 3:length(folders)
%     cd(folders(i).name);
%     ci = dir('*cellinfo.mat');
%     u = strfind(ci(1).name,'_');
%     pos = u(2);
%     a = exist([ci(1).name(1:pos-1) '.mat'],'file');
%     if a == 2
%         flist = [flist folders(i).name];
%     end
% end
% Have it prompt user for a string of Y or N, whether to include unsorted
% spikes for each cell. Then batch process.

[filename,pathname] = uigetfile('*.mat','Choose CONCATENATED .mat file of exported spikes from Offline Sorter.');
cd(pathname)
disp(filename(1:end-4))

% rename pieces of OfflineSorter-created file to be meaningful
load([pathname filename])
if exist('adc001')
    unit = adc001(:,1);
    spiketime = adc001(:,2); % in Seconds
    snippet = adc001(:,3:end); % in Samples
    clear adc001
else
    error('The variable adc001 did not load as expected.')
end


%% Prep to load raw data files that were concatenated and exported as .nex files to Offline Sorter
cd(pathname);

% Load struct cellinfo that was created by ExportTytologyDataToNex.m
load([pathname, filename(1:end-4),'_cellinfo.mat']) 

% load(cellinfo_files.name);  % we've already changed to the correct directory
% disp(['Cell: ',cellinfo_files.name]) % NB: this assumes there is ONLY ONE _cellinfo.mat file per cell folder

disp(['Cell: ',filename]) % NB: this assumes there is ONLY ONE _cellinfo.mat file per cell folder

% User chooses whether or not to include unsorted spikes in future analyses:
cellinfo.includeUnsort = input('Do you want to include unsorted spikes (marked as zeros)? Y or N (default Y): ','s');
if isempty(cellinfo.includeUnsort)
    cellinfo.includeUnsort = 'Y';
end


dat_files = dir('*.dat');
if isempty(dat_files)
    warndlg('Put some .dat files - and matching .mat files - in the right folder!'); return
end

% Ensure that loaded .dat files match those originally sorted by Offline Sorter
check = [];
for i = 1:length(dat_files)
    if strcmp(dat_files(i).name(1:end-4), cellinfo.filename{i})
        check = [check 1];
    else
        check = [check 0];
    end
end

if mean(check)
else
    error('The .dat files do not match those used to create the .nex file exported to Offline Sorter.')
end
disp(['Loaded ', num2str(length(dat_files)), ' .dat files'])


%%

% setup for loop
signal=[]; timestampsSec=[]; timestampsN=[];  parameter=[];  count = 0;
% rawplots = figure(1); set(gcf,'position',[50,400,1200,500])

% Step through each raw data file
for FF=1:length(dat_files)
    
    cd(pathname);
    
    % Load raw .dat and .mat files
    disp(dat_files(FF).name(1:end-4))
    data = readHPData(dat_files(FF).name);
    % Load struct datainfo, which has fields .curvedata and .curvesettings
    datainfo = load([dat_files(FF).name(1:end-4), '.mat']);
    curvesettings = datainfo.curvesettings;
    curvedata = datainfo.curvedata;

    cd(processpath);
    
    % if data was collected prior to Fs being written to data file
    if datenum(curvesettings.time_start(1:end-9)) < datenum('13-Mar-2013')
        Fs = 2.4414e+004;
    else
        Fs = curvesettings.tdt.neuralFs;
    end
    
    if FF == 1
        % create filter to remove the annoying ~120Hz powered-speaker noise
        cutoff = 200/(Fs/2);
        numpole = 7;
        [b,a] = butter(numpole,cutoff,'high'); % set to 'high' for highpass
    end
    
%     rawplots = figure(1); set(gcf,'position',[50,400,1200,500])
    
    % Step through each trial
    for T = 1:length(data)
        
        count = count+1;
%         figure(rawplots);
%         set(gcf,'name',['RawTrace + SortedSpikes: ', filename(1:end-4), '  TRIAL ', num2str(T)]); set(gcf,'numbertitle','off');

        wave = data{T}.datatrace;  % waveform for each trial (single electrode)
        wavefilt = filter(b,a,wave);    % filter waveform to remove low freq noise
        signal = [signal; wavefilt];
        timestampsSec = [timestampsSec; length(signal)/Fs];
        timestampsN = [timestampsN; length(signal)];
        
        if T == 1 && FF == 1 % Process this separately because first trial in concatenated file starts at time 0
            sortedspikes{T} = spiketime(find(spiketime <= timestampsSec(T)));
            spikeidx = find(spiketime <= timestampsSec(T));
            offset = 0;
        else
            sortedspikes{T} = spiketime(find( spiketime <= timestampsSec(count) & spiketime > timestampsSec(count-1) ));
            spikeidx = (find( spiketime <= timestampsSec(count) & spiketime > timestampsSec(count-1) ));
            offset = timestampsSec(count-1);
        end
        
        unitID{T} = unit(spikeidx);
        
        % create time vector for waveform
        time = [0: (1/Fs) :length(wavefilt)/Fs]; 
        time = time(2:end);
        time = time+offset;
        
        % separate types of units
        invalididx{T}  = find(unitID{T} == -1);
        unsortedidx{T} = find(unitID{T} == 0);
        unit1idx{T}    = find(unitID{T} == 1);
        unit2idx{T}    = find(unitID{T} == 2);

        
        % PLOT
%         figure(rawplots); hold off; 
%         plot(time, wavefilt, 'k'); hold on;
%         
%         plot(sortedspikes{T}(invalididx{T}),  ones(length(invalididx{T}),1)*0.095, 'g.', 'markersize', 15)
%         plot(sortedspikes{T}(unsortedidx{T}), ones(length(unsortedidx{T}),1)*0.1,  'k.', 'markersize', 15)
%         plot(sortedspikes{T}(unit1idx{T}),    ones(length(unit1idx{T}),1)*0.105,   'r.', 'markersize', 15)
%         plot(sortedspikes{T}(unit2idx{T}),    ones(length(unit2idx{T}),1)*0.11,    'b.', 'markersize', 15)
%         xlim([time(1), time(end)])
%         ylim([-0.3, 0.3])
% pause

        % write useful info to struct (eventually to be saved over .mat data file)
        curvedata.sortedspikes{T} = sortedspikes{T};        % spike times in sec
        curvedata.offsetsec{T}    = offset;                 % offset indicating trial start (sec)
        curvedata.unitID{T}       = unitID{T};              % -1 = invalid, 0 = unsorted, 1 = unit1, 2 = unit2
        curvedata.invalididx{T}   = invalididx{T};          % indices of invalidated spikes
        curvedata.unsortedidx{T}  = unsortedidx{T};         % indices of unsorted spikes
        curvedata.unit1idx{T}     = unit1idx{T};            % indices of unit1 spikes
        curvedata.unit2idx{T}     = unit2idx{T};            % indices of unit2 spikes
       
%         parameter(T) = curvesettings.stimcache.stimvar{T};
        
    end % each trial T
    
    % CREATE NEW .mat file containing 2 structs: curvesettings, curvedata
    matfile = [dat_files(FF).name(1:end-4), 'new.mat'];
    save([pathname, matfile], 'curvesettings', 'curvedata');
    save([pathname,filename(1:end-4),'_cellinfo.mat'],'cellinfo')
    
    clear sortedspikes unitID unsortedidx unit1idx unit2idx invalididx
     
end % each data file FF



%% Plot waveform snippets
% figure
% for i=1:size(snippet,2)
%     subplot(10,4,i)
%     plot(snippet(i,:),'-')
% end

