load P
%%
%{
P = 
                 Version: 107
                 Comment: ''
                    Date: 7.3788e+05
        NumSpikeChannels: 4
        NumEventChannels: 4
         NumContChannels: 4
             ADFrequency: 48828
           NumPointsWave: 64
         NumPointsPreThr: 24
                FastRead: 0
            WaveformFreq: 48828
           LastTimestamp: 1953119
              Trodalness: 1
          DataTrodalness: 1
      BitsPerSpikeSample: 16
       BitsPerContSample: 16
     SpikeMaxMagnitudeMV: 3000
      ContMaxMagnitudeMV: 5000
         SpikePreAmpGain: 1000
       AcquiringSoftware: ''
      ProcessingSoftware: 'OFS 3.3.5'
    SpikeTimestampCounts: [27x4 double]
     SpikeWaveformCounts: [27x4 double]
             EventCounts: [160 160 1 1]
        ContSampleCounts: [0 0 0 0]
     ContSampleFragments: [0 0 0 0]
           SpikeChannels: [4x1 struct]
           EventChannels: [4x1 struct]
      ContinuousChannels: [4x1 struct]
                FullRead: 1
       DataStartLocation: 13952
                 PLXFile: '/Users/sshanbhag/Work/Data/TestData/MT/1407/140...'

P.SpikeChannels(1)
          Name: 'spikechan_4'
       Channel: 1
       SIGName: 'adc000'
           SIG: 0
      SourceID: 2
     ChannelID: 0
       Comment: ''
        NUnits: 3
           Ref: 0
        Filter: 0
          Gain: 0
     Threshold: -11508
        WFRate: 10
    SortMethod: 2
       SortBeg: 0
     SortWidth: 64
      Template: [5x64 double]
         Boxes: [5x2x4 double]
           Fit: [5x1 double]
    Timestamps: [78x1 uint32]
         Units: [78x1 uint8]
         Waves: [64x78 int16]
%}

c = 1;

% convert Timestamps to seconds (need to do conversion to double from
% uint32
ts_sec = double(P.SpikeChannels(c).Timestamps) / double(P.ADFrequency);
unit = double(P.SpikeChannels(c).Units);
channel = 4*ones(size(ts_sec));
PCA = zeros(length(ts_sec), 3);
snips = double(P.SpikeChannels(c).Waves') / double(P.SpikeMaxMagnitudeMV));
[channel unit ts_sec PCA]

