% loading exported .mat data from Plexon

% two ways:

% (1) use the load command or function:
%
% either
% 
% load(<string with path and filename>)
% 
% or
% 
% load <string with path and filename>
% 
% e.g.:

load('/Volumes/Lexar/Work/DATA/Plexon-ICdata/1407_20200309_03_01_1350_BBN-sorted.ch4,5,7,15.mat')

load '/Volumes/Lexar/Work/DATA/Plexon-ICdata/1407_20200309_03_01_1350_BBN-sorted.ch4,5,7,15.mat'


% To see what's variables were loaded, either check the workspace window, 
% or use the matlab who() command:
% 
% >> who
% 
% Your variables are:
% 
% adc000  adc001  adc002  adc003  ans     

% if you want to see what's  in the .mat file without loading it:
	
who('-file', '/Volumes/Lexar/Work/DATA/Plexon-ICdata/1407_20200309_03_01_1350_BBN-sorted.ch4,5,7,15.mat')
% 
% who('-file', '/Volumes/Lexar/Work/DATA/Plexon-ICdata/1407_20200309_03_01_1350_BBN-sorted.ch4,5,7,15.mat')
% 
% Your variables are:
% 
% adc000  adc001  adc002  adc003  


% and to load a specified variable:

load('/Volumes/Lexar/Work/DATA/Plexon-ICdata/1407_20200309_03_01_1350_BBN-sorted.ch4,5,7,15.mat', 'adc000')

% >> who
% 
% Your variables are:
% 
% adc000  

% okay so variable is loaded.  to look at the first few columns 
% the original files had:
%		Column 1: channel (?)
%		Column 2: unit #
%		Column 3: timestamp (in seconds)
%		Column 4: PCA1 weight
%		Column 5: PCA2 weight
%		Column 6: PCA3 weight
%		Column 7-end : waveform

% To see the values for the first five rows and the columns 1:3 for adc000:

adc000(1:5, 1:3)

% this should enable you to see whether the exported data have info about
% the channel. alternatively, just look at the variable in the workspace
% browser window...double clicking on the variable name in the workspace browser
% will bring up a spreadsheet-type view of the values stored in that
% variable

