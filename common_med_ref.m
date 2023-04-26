function Tout = common_med_ref(Tin)
%---------------------------------------------------------------------
% apply common median reference to matrix of channel data [samples,
% channels]
%---------------------------------------------------------------------
% Created: 26 April 2023 (SJS)
% Originally nested in testCommonReference script
%
%---------------------------------------------------------------------
% To Do:
%---------------------------------------------------------------------
A = median(Tin, 2);
Tout = Tin - A;
