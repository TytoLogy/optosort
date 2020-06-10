function [h, delay] = getMatrixFIR(p, q, Lx, N, bta)
% p, q	numerator, denominator of rational approximation to Fs_new/Fs_old
% Lx		num rows in matrix to be filtered
% N		order for antialiasing filter will be  2 X n X max(p, q)
%				default is 10
% bta		beta value for kaiser window
%				default is 5
%
% Uses code from RESAMPLE()

if nargin < 5,  bta = 5;  end   %--- design parameter for Kaiser window LPF
if nargin < 4,   N = 10;   end

% max(X,Y) returns an array the same size as X and Y with the largest
% elements taken from X or Y. Either one can be a scalar.
pqmax = max(p,q);

% design filter (assume N > 0)
fc = 1/2/pqmax;
L = 2*N*pqmax + 1;
h = firls( L-1, [0 2*fc 2*fc 1], [1 1 0 0]).*kaiser(L,bta)' ;
h = p*h/sum(h);
Lhalf = (L-1)/2;

% Need to delay output so that downsampling by q hits center tap of filter.
nz = floor(q-mod(Lhalf,q));
z = zeros(1,nz);
h = [z h(:).'];  % ensure that h is a row vector.
Lhalf = Lhalf + nz;

% Number of samples removed from beginning of output sequence 
% to compensate for delay of linear phase filter:
delay = floor(ceil(Lhalf)/q);

% Need to zero-pad so output length is exactly ceil(Lx*p/q).
nz1 = 0;
while ceil( ((Lx-1)*p+length(h)+nz1 )/q ) - delay < ceil(Lx*p/q)
    nz1 = nz1+1;
end
h = [h zeros(1,nz1)];


