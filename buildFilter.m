function BPfilt = buildFilter(Fs, varargin)
%---------------------------------------------------------------------
% define filter
%---------------------------------------------------------------------
% Created: 26 April 2023 (SJS)
% Originally nested in testCommonReference script
%
%---------------------------------------------------------------------
% To Do:
%  add input args for Fc, forder, ramp time, filter type
%---------------------------------------------------------------------
% [highpass lowpass] cutoff frequencies in Hz
BPfilt.Fc = [250 4000];
% order of filter. note that the filtfilt() function in MATLAB is used,
% so the effective order is doubled. typically use 5
BPfilt.forder = 5;
% ramp time (ms) to apply to each sweep in order to cutdown on
% onset/offset transients from filtering
BPfilt.ramp = 1;
% filter type. 'bessel' or 'butter' (for butterworth). typically use
% butter as bessel is low-pass only and we need to remove low frequency
% crap from the raw data
BPfilt.type = 'butter';

BPfilt.Fs = Fs;
BPfilt.Fnyq = Fs / 2;
BPfilt.cutoff = BPfilt.Fc / BPfilt.Fnyq;
if strcmpi(BPfilt.type, 'bessel')
   [BPfilt.b, BPfilt.a] = besself(BPfilt.forder, BPfilt.cutoff, ...
                                       'bandpass');
elseif strcmpi(BPfilt.type, 'butter')
   [BPfilt.b, BPfilt.a] = butter(BPfilt.forder, BPfilt.cutoff, ...
                                       'bandpass');
else
   error('%s: unknown filter type %s', mfilename, BPfilt.type)
end
