function [reD, reDinf] = resample_data(D, Dinf, Fs_new)
%------------------------------------------------------------------------
% [reD, reDinf] = resample_data(D, Dinf, Fs_new)
%------------------------------------------------------------------------
% TytoLogy:Experiments:optosort
%------------------------------------------------------------------------
% Takes data in cell array D and resamples to rate specified in Fs_new
%
% 
%------------------------------------------------------------------------
% Input Arguments:
%	D			contains data in a cell structure array.
%	Dinf		has the file header information
%	Fs_new	new sampling rate in samples/second
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

% first, some checks

% if sample rates are equivalent, return input vars unchanged
if (Dinf.indev.Fs == Fs_new) || isempty(Fs_new)
	warning('resampling continuous data: Dinf.indev.Fs == Fs_new');
	reD = D;
	reDinf = Dinf;
	return
end



% specify sample rates
Fs_old = Dinf.indev.Fs;
fprintf('Original Sample Rate: %.4f\n', Fs_old);
fprintf('Resampled SampleRate: %.4f\n', Fs_new);

% determine ratio via closest rational approximation
[p, q] = rat(Fs_new/Fs_old, 1e-9);
fprintf('New Ratio:\n');
fprintf('p = %.2f\nq = %.2f\n', p, q);
fprintf('Error abs( (p/q)*Fs_old - Fs_new): %.12f\n', abs( (p/q)*Fs_old - Fs_new));

% beta value for kaiser window (default)
bta = 5;
% antialiasing filter will be order 2 X n X max(p, q)
% default for N is 10
N = 10;

% length (# rows) of datatrace matrix (assume all are same for this file)
Lx = size(D{1}.datatrace, 1);

% build filter
% get filter
[h, delay] = getMatrixFIR(p, q, Lx, N, bta);

% to get rid of trailing and leading data so input and output signals line up
% temporally:
Ly = ceil(Lx*p/q);  % output length

% allocate output
reD = D;
% assign new value for Fs
reDinf = Dinf;
reDinf.indev.Fs = Fs_new;


% loop through sweeps, resample
for c = 1:numel(D)
	y = upfirdn(D{c}.datatrace, h, p, q);
	y(1:delay,:) = [];
	y(Ly+1:end,:) = [];
	reD{c}.datatrace = stripZeros(y, p, q, Lx, delay);
end







