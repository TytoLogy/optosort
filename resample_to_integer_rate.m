function [reD, reDinf] = resample_to_integer_rate(D, Dinf)
%------------------------------------------------------------------------
% [cSweeps, sInfo] = read_data_for_export(F, Channels, BPfilt, resampleData)
%------------------------------------------------------------------------
% TytoLogy:Experiments:optosort
%------------------------------------------------------------------------
% Takes data in cell array D and resamples to nearest lower integer value
% Updates Fs field in Dinf to new value
%
% e.g., if original sample rate is 48828.125 samples/second, new rate will
% be 48818 samples/second
% 
%------------------------------------------------------------------------
% Input Arguments:
%	D			contains data in a cell structure array.
%	Dinf		has the file header information
% 
% Output Arguments:
%	reD			contains resampled data in a cell structure array.
%	reDinf		has the file header information with update Fs
%------------------------------------------------------------------------
% See also: read_data_for_export, readOptoData
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%  Sharad J. Shanbhag
%	sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 9 June, 2020 (SJS) will be called by read_data_for_export 
%
% Revisions:
%------------------------------------------------------------------------
% TO DO:
%------------------------------------------------------------------------







%

% specify sample rates
Fs_old = S.Info.Fs;
Fs_new = 48800;
fprintf('Original Sample Rate: %.4f\n', Fs_old);
fprintf('Resampled SampleRate: %.4f\n', Fs_new);

% determine ratio via closest rational approximation
[p, q] = rat(Fs_new/Fs_old, 1e-9)
fprintf('Error: %.12f\n', abs( (p/q)*Fs_old - Fs_new))
cData = resample(dData, p, q);







