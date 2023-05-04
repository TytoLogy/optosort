function Tout = common_avg_ref(Tin)
%---------------------------------------------------------------------
% apply common average reference to matrix of channel data [samples,
% channels]
%---------------------------------------------------------------------
% Created: 26 April 2023 (SJS)
%     Originally nested in testCommonReference script
% Revised:
%  4 May 2023 (SJS): making changes to account for row-oriented data
%
%---------------------------------------------------------------------
% To Do:
%---------------------------------------------------------------------
A = mean(Tin);
Tout = Tin - A;
