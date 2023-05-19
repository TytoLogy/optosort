% using this to develop ui for data explore (easier to deal with on its own
% vs in a function!!!)

% GUI Elements: Units
fig.Units = 'characters';
refButtonGroup.Units = 'characters';
stimInfoPanel.Units = 'characters';
filenameText.Units = 'characters';
trialPanel.Units = 'characters';
ax.Units = 'characters';
medRadioButton.Units = 'characters';
avgRadioButton.Units = 'characters';
rawRadioButton.Units = 'characters';
stimuluslevelText.Units = 'characters';
onsetoffsetText.Units = 'characters';
stimulusText.Units = 'characters';
trialnEdit.Units = 'characters';
trialUpButton.Units = 'characters';
trialDownButton.Units = 'characters';
% GUI Elements: Position
% fig.Position = [125.4286 59.8750 99 53.1875];
refButtonGroup.Position = [71.1429 0.6250 27.2857 6.3125 ];
stimInfoPanel.Position = [41.4286 1.6875 27.2857 5.2500 ];
filenameText.Position = [27.0000 51.3750 50.1429 1.3125 ];
trialPanel.Position = [1.2857 3.1250 38.7143 3.8125 ];
ax.Position = [8.2857 9.9375 88.7143 41.3125 ];
medRadioButton.Position = [1.2857 0.8125 11.5714 1.0625 ];
avgRadioButton.Position = [1.2857 2.3750 11.5714 1.0625 ];
rawRadioButton.Position = [1.2857 3.9375 11.5714 1.0625 ];
stimuluslevelText.Position = [1.2857 0.1250 24.5714 1.3125 ];
onsetoffsetText.Position = [1.1429 1.3750 24.7143 1.3125 ];
stimulusText.Position = [1.1429 2.6250 24.7143 1.3125 ];
trialnEdit.Position = [28.2857 0.6250 8.7143 1.9375 ];
trialUpButton.Position = [14.1429 0.5000 11.5714 2.0000 ];
trialDownButton.Position = [1.2857 0.5625 11.5714 2.0000 ];
%{

%% set positions of trial objects (easier to do together)
% these values are [x y width height]
trialPanel.Position = [0.02 0.01 0.35 0.15];
trialDown.Position = [0.02 0.12 0.3 0.75];
trialUp.Position = [0.35 0.12 0.3 0.75];
trialNum.Position = [0.7 0.12 0.26 0.75];


%% axes position and size
ax.Position = [0.1 0.245 0.85 0.68];
ax.TickDir = 'out';

%}