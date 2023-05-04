function Tout = common_med_ref(Tin)
%---------------------------------------------------------------------
% apply common median reference to matrix of channel data [samples,
% channels]
%---------------------------------------------------------------------
% Created: 26 April 2023 (SJS)
%     Originally nested in testCommonReference script
% Revised:
%  4 May 2023 (SJS): making changes to account for row-oriented data
%---------------------------------------------------------------------
% To Do:
%---------------------------------------------------------------------
A = median(Tin);
Tout = Tin - A;
