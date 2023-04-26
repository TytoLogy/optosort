function Tout = common_avg_ref(Tin)
%---------------------------------------------------------------------
% apply common average reference to matrix of channel data [samples,
% channels]
%---------------------------------------------------------------------
% Created: 26 April 2023 (SJS)
% Originally nested in testCommonReference script
%
%---------------------------------------------------------------------
% To Do:
%---------------------------------------------------------------------
A = mean(Tin, 2);
Tout = Tin - A;
