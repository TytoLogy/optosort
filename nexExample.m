% exampleSaveDataInNexAndNex5File.m
% file to test export of data to .nex files
%{
% based on code from 
% /Users/sshanbhag/Work/Code/Matlab/stable/Toolbox/NeuroExplorer/...
% HowToReadAndWriteNexAndNex5FilesInMatlab/README.txt


Writing .nex and .nex5 files
------------------

Use writeNexFile function to write the contents of nexFileData
structure to .nex file.

Use write5NexFile function to write the contents of 
nexFileData structure to .nex5 file.

You can read .nex file and then save nexFileData in .nex5 file and vice versa.

For example (see file exampleSaveDataInNexAndNex5File.m):

!!!!
note that nex5 file supports continuous values as float vs. nex file
which writes continuous values as int16. 
problem? Plexon OFS 3.3.5 does not support nex5 format...
¡¡¡¡
%}

addpath(['/Users/sshanbhag/Work/Code/Matlab/stable/Toolbox/NeuroExplorer' ...
				'/HowToReadAndWriteNexAndNex5FilesInMatlab']);

%% Initialize nex file data (nexFile struct)
%{
[nexFile] = nexCreateFileData(timestampFrequency) -- creates empty nex file data structure

INPUT:
  timestampFrequency - timestamp frequency in Hertz
%}
nexFile = nexCreateFileData(40000);

%% add continuous variable
%{
specify start time (t(1)), digitizing frequency (Fs), data (x2) and name
[nexFile] = nexAddContinuous( nexFile, startTime, adFreq, values, name ) 
        -- adds continuous variable to nexFile data structure

INPUT:
  nexFile - nex file data structure created in nexCreateFileData
  startTime - time of the first data point in seconds
  adFreq - A/D sampling rate of continuous variable in samples per second
  values - vector of continuous variable values in milliVolts
  name - continuous variable name  
%}
% digitizing frequency 1000 Hz
Fs = 1000;
% time interval from 1 to 5
t= 1:1/Fs:5;
% sin with frequency 2 Hz
contValues_2Hz = sin(2*pi*t*2);
% sin with frequency 2 Hz
contValues_4Hz = sin(2*pi*t*4);
nexFile = nexAddContinuous(nexFile, t(1), Fs, contValues_2Hz, 'sin2Hz');
nexFile = nexAddContinuous(nexFile, t(1), Fs, contValues_4Hz, 'sin4Hz');

%% add neuron spike train
%{
[nexFile] = nexAddNeuron( nexFile, timestamps, name ) -- adds a neuron 
            to nexFile data structure

INPUT:
  nexFile - nex file data structure created in nexCreateFileData
  timestamps - vector of neuron timestamps in seconds
  name - neuron name
  wireNumber - wire number (optional)
  unitNumber - unit number (optional). Zero unit number means unsorted.
%}
neuronTs = [0.5 0.9 2.1 2.3 2.5]';
nexFile = nexAddNeuron(nexFile, neuronTs, 'neuron1');

%% add event spike train
%{
[nexFile] = nexAddEvent( nexFile, timestamps, name ) -- adds an event 
            (series of timestamps) to nexFile data structure

INPUT:
  nexFile - nex file data structure created in nexCreateFileData
  timestamps - vector of event timestamps in seconds
  name - event name
%}
% events must be in column format...?
eventTs = [10 20 30 40]';
nexFile = nexAddEvent(nexFile, eventTs, 'event1');

%% add interval variable
% [nexFile] = nexAddInterval( nexFile, intervalStarts, intervalEnds, name ) 
%              -- adds an interval variable (series of intervals) to nexFile data structure
%
% INPUT:
%   nexFile - nex file data structure created in nexCreateFileData
%   intervalStarts - a vector of interval starts in seconds
%   intervalEnds - a vector of interval ends in seconds
%   name - interval variable name
intStarts = [5 10];
intEnds = [6 12];
nexFile = nexAddInterval(nexFile, intStarts, intEnds, 'interval1');

%% add  waveforms
%{
[nexFile] = nexAddWaveform( nexFile, WFreq, timestamps, waveforms, name, preThresholdTimeInSeconds, numberOfPointsInWaveform )
            -- adds waveform variable to nexFile data structure

INPUT:
  nexFile - nex file data structure created in nexCreateFileData
  startTime - time of the first data point in seconds
  WFreq - A/D sampling rate of waveform variable in samples per second
  timestamps - vector of wave timestamps (in seconds)
  waveforms - matrix of waveform variable values in milliVolts
              each waveform is a column (matrix can be empty if numberOfPointsInWaveform is specified)
  name - waveform variable name  
  preThresholdTimeInSeconds - (optional parameter) if waveform timestamp in seconds is t, 
        then the timestamp of the first point of waveform is t - prethresholdTimeInSeconds
  numberOfPointsInWaveform - (optional parameter) number of points in waveform; 
        if this parameter is omitted, the number of points in waveform is 
        the number of rows in waveforms matrix
  wireNumber - wire number (optional parameter)
  unitNumber - unit number (optional parameter). Zero unit number means unsorted.
%}
% waveform timestamps
waveTs = [1 2]';
% 2 waveforms (columns), 5 data points each
waves = [-10 0 10 20 30; -15 0 15 25 15]';
nexFile = nexAddWaveform(nexFile, 40000, waveTs, waves, 'wave1');

%% save nex file
%{
 [result] = writeNexFile(nexFile, fileName) -- write nexFile structure
to the specified .nex file. returns 1 if succeeded, 0 if failed.

INPUT:
  nexFile - a structure containing .nex file data

          SOME FIELDS OF THIS STRUCTURE (VERSIONS ETC.) ARE NOT DESCRIBED
          BELOW. IT IS RECOMMENDED THAT YOU READ A VALID .NEX FILE
          TO FILL THIS STRUCTURE, THEN MODIFY THE STRUCTURE AND SAVE IT.

          IF YOU WANT TO CREATE NEW .NEX FILE, USE nexCreateFileData.m,
          nexAddContinuous.m etc. See exampleSaveDataInNexFile.m.

  fileName - if empty string, will use File Save dialog

  nexFile - a structure containing .nex file data
  nexFile.version - file version
  nexFile.comment - file comment
  nexFile.tbeg - beginning of recording session (in seconds)
  nexFile.tend - end of recording session (in seconds)

  nexFile.neurons - array of neuron structures
          neurons{i}.name - name of a neuron variable
          neurons{i}.timestamps - array of neuron timestamps (in seconds)
              to access timestamps for neuron 2 use {n} notation:
              nexFile.neurons{2}.timestamps

  nexFile.events - array of event structures
          events{i}.name - name of event variable
          events{i}.timestamps - array of event timestamps (in seconds)

  nexFile.intervals - array of interval structures
          intervals{i}.name - name of interval variable
          intervals{i}.intStarts - array of interval starts (in seconds)
          intervals{i}.intEnds - array of interval ends (in seconds)

  nexFile.waves - array of wave structures
          waves{i}.name - name of waveform variable
          waves{i}.NPointsWave - number of data points in each wave
          waves{i}.WFrequency - A/D frequency for wave data points
          waves{i}.timestamps - array of wave timestamps (in seconds)
          waves{i}.waveforms - matrix of waveforms (in milliVolts), each
                            waveform is a column

  nexFile.contvars - array of continuous variable structures
          contvars{i}.name - name of continuous variable
          contvars{i}.ADFrequency - A/D frequency for data points

          Continuous (a/d) data for one channel is allowed to have gaps
          in the recording (for example, if recording was paused, etc.).
          Therefore, continuous data is stored in fragments.
          Each fragment has a timestamp and an index of the first data
          point of the fragment (data values for all fragments are stored
          in one array and the index indicates the start of the fragment
          data in this array).
          The timestamp corresponds to the time of recording of
          the first a/d value in this fragment.

          contvars{i}.timestamps - array of timestamps (fragments start times in seconds)
          contvars{i}.fragmentStarts - array of start indexes for fragments in contvar.data array
          contvars{i}.data - array of data points (in milliVolts)

  nexFile.popvectors - array of population vector structures
          popvectors{i}.name - name of population vector variable
          popvectors{i}.weights - array of population vector weights

  nexFile.markers - array of marker structures
          markers{i}.name - name of marker variable
          markers{i}.timestamps - array of marker timestamps (in seconds)
          markers{i}.values - array of marker value structures
              markers{i}.value.name - name of marker value
              markers{i}.value.strings - array of marker value strings

%}
nexFileDir = pwd;
nexFileName = 'nexExample.nex';
writeNexFile(nexFile, fullfile(nexFileDir, nexFileName));
%% Read nex file

nexData = readNexFile(fullfile(nexFileDir, nexFileName));
