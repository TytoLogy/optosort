% add plexon toolbox to path
addpath(['/Users/sshanbhag/Work/Code/Matlab/stable/Toolbox' ...
				'/Plexon/HowToReadAndWriteNexFilesInMatlab']);

%% Organise tytology data files into subfolders based on cell number

clear all; close all;

% startloc = 'C:\Ephys\runme\';
startloc = 'C:\Experiments\Data\Headfixed_extracell_Tytology\runme\';
cd(startloc);

% look for males...
throw1 = dir('M*');
% ... and females
throw2 = dir('F*');
% .. and combine into full list
critters = [throw1;  throw2];
clear throw1 throw2

% loop through all animals
for ii = 1:length(critters)
	% directory for source data for this animal
    critterloc = [startloc, critters(ii).name, '\'];
    cd(critterloc);
    throw1 = dir('M*');
    throw2 = dir('F*');
    datafiles = [throw1;  throw2];
    clear throw1 throw2
    for jj = 1:length(datafiles)
        tempstr = datafiles(jj).name(7:9);
        cellnumber{jj} = tempstr;
    end
    throw = cellnumber;
    uniquecells = unique(throw);
    for kk = 1:length(uniquecells)
       if exist(uniquecells{kk},'dir') ~= 7
        mkdir(uniquecells{kk}); 
       end
       copystr = [critterloc,'*','_',uniquecells{kk},'_','*'];
       deststr = [critterloc,uniquecells{kk}];
       movefile(copystr,deststr);
    end
    
    clear throw
    clear datafiles uniquecells cellnumber
end

%%
clearvars -except startloc

processpath = 'C:\Matlab_work\Tytology\';
cd(startloc);

critters = dir;
critters = critters(3:end);

for ii = 1:length(critters)
	critterpath = [startloc critters(ii).name '\'];
	cd(critterpath);

	neurons = dir;
	neurons = neurons(3:end);

	for jj = 1: length(neurons)
		neuronpath = [critterpath neurons(jj).name '\'];
		cd(neuronpath);

		dat_files = dir('*.dat');

		if isempty(dat_files)
		warndlg(['Put some .dat files - and matching .mat files - ' ...
		'in the right folder!']);
		return
		end

		disp(['Loaded ', num2str(length(dat_files)), ' .dat files'])

		signal = []; timestamps = []; parameter = [];

		% Step through each dat file and append them together
		try
			for FF=1:length(dat_files)
				% change to neuron data directory
				cd(neuronpath);

				disp(dat_files(FF).name(1:end-4))
				% read data file
				data = readHPData(dat_files(FF).name);
				% Load struct datainfo, which has fields .curvedata and .curvesettings
				datainfo = load([dat_files(FF).name(1:end-4), '.mat']);
				% change to exported/processed data directory
				cd(processpath);

				% if data was collected prior to Fs being written to data file
				if datenum(datainfo.curvesettings.time_start(1:end-9)) < datenum('13-Mar-2013')
					Fs = 2.4414e+004;
				else
					Fs = datainfo.curvesettings.tdt.neuralFs;
				end

				if FF == 1
					% start new nex file data
					nexFile = nexCreateFileData(Fs);

				%                 % create filter to remove the annoying ~120Hz powered-speaker noise
				%                 cutoff2 = 200/(Fs/2);
				%                 numpole = 7;
				%                 [B,A] = butter(numpole,cutoff2,'high'); % set to 'high' for highpass

						 % create bandpass filter, since raw traces are being saved as of 9/2016
						 % changed LPF to 5000 from 10,000 to remove noise from Kopf 01/2018

					cutoff = [200 5000]/(Fs/2);
					numpole = 7;
					[b,a] = butter(numpole,cutoff,'bandpass'); % set to 'high' for highpass, 'bandpass' for bandpass
				end

				cellinfo.filename{FF} = dat_files(FF).name(1:end-4);
				cellinfo.curvetype{FF} = datainfo.curvesettings.curve.curvetype;
				cellinfo.trialduration{FF} = datainfo.curvesettings.tdt.AcqDuration;
				cellinfo.sigstart{FF} = datainfo.curvesettings.tdt.StimDelay;
				cellinfo.sigend{FF} = datainfo.curvesettings.tdt.StimDelay + datainfo.curvesettings.tdt.StimDuration;

				% Step through each trial
				for T = 1:length(data)
					wave = data{T}.datatrace;  % waveform for each trial (single electrode)
					offset = wave(1);
					wave = wave - offset;
					wavefilt = filter(b,a,wave);

					%         wavename = [datainfo.curvesettings.stimcache.curvetype, ' ', num2str(datainfo.curvesettings.stimcache.stimvar{T})];
					timestamps = [timestamps; length(signal)/Fs];
					signal = [signal; wavefilt];

					parameter(T) = datainfo.curvesettings.stimcache.stimvar{T};
				end % each trial T

				cellinfo.parameter{FF} = parameter;

			end % each data file FF
		catch
			error('file error');
		end
        
		% use wavename to identify data for this file
		wavename = [dat_files(FF).name(1:9)];

		% SYNTAX: nexfile = nexAddContinuous(nexFile, starttime, Fs, wave, 'wavename');  HOWEVER: starttime must be 0 - otherwise it doesn't work
		nexFile = nexAddContinuous(nexFile, 0, Fs, signal, wavename);   % add continuous variable
		nexFile = nexAddEvent(nexFile, timestamps, 'signal_onsets');

		% save nex file
		writeNexFile(nexFile, [neuronpath, dat_files(FF).name(1:9), '.nex']);
		%  Copy .nex files into the main folder for processing by plexon
		%         copyfile([neuronpath, dat_files(FF).name(1:9), '.nex'],startloc );
		% save cellinfo file
		save([neuronpath, dat_files(FF).name(1:9), '_cellinfo.mat'],'cellinfo')

		clear neuronpath
	end % jj loop
	clear critterpath neurons
end % ii loop

