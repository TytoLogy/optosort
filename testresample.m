
DataPath = '/Users/sshanbhag/Work/Data/TestData/MT/1407';
DataFile = '1407_20200309_03_01_1350_BBN.dat';
TestFile = '1407_20200309_03_01_1350_BBN_testdata.mat';
Channels = [4, 5, 7, 15];
Fs_new = 48820;

%% use readOptoData to read in raw data. 
[D, Dinf] = readOptoData(fullfile(DataPath, DataFile));
% Fix test info
Dinf = correctTestType(Dinf);

%% resample data if needed
% if ~isempty(resampleRate)
% 	[D, tmpDinf] = resample_data(D, tmpDinf, resampleRate);
% end
% 


%% specify sample rates
Fs_old = Dinf.indev.Fs;
fprintf('Original Sample Rate: %.4f\n', Fs_old);
fprintf('Resampled SampleRate: %.4f\n', Fs_new);

% determine ratio via closest rational approximation
[p, q] = rat(Fs_new/Fs_old, 1e-12);
fprintf('New Ratio:\n');
fprintf('p = %.2f\nq = %.2f\n', p, q);
fprintf('Error abs( (p/q)*Fs_old - Fs_new): %.12f\n', abs( (p/q)*Fs_old - Fs_new));

% allocate output
reD = D;
stime = zeros(size(D));


% loop through sweeps, resample
for c = 1:numel(D)
	if c == 1
		fprintf('first run, no B\n');
		tic
		[reD{c}.datatrace, B] = resample(D{c}.datatrace, p, q);
		stime(c) = toc; 
		fprintf('Run %d took\t\t%.2f\n', c, stime(c));
	else
		tic
		reD{c}.datatrace = resample(D{c}.datatrace, p, q, B);
		stime(c) = toc; 
		fprintf('Run %d took\t\t%.2f\n', c, stime(c));
	end
end

%% try using upfirdn
% code adapted from the uniformResample function in the matlab resample()
% function

Fs_old = Dinf.indev.Fs;
fprintf('Original Sample Rate: %.4f\n', Fs_old);
fprintf('Resampled SampleRate: %.4f\n', Fs_new);

% determine ratio via closest rational approximation
[p, q] = rat(Fs_new/Fs_old, 1e-9);
fprintf('New Ratio:\n');
fprintf('p = %.2f\nq = %.2f\n', p, q);
fprintf('Error abs( (p/q)*Fs_old - Fs_new): %.12f\n', abs( (p/q)*Fs_old - Fs_new));

x = D{2}.datatrace;

%% new functions
% beta value for kaiser window (default)
bta = 5;
% antialiasing filter will be order 2 X n X max(p, q)
% default for N is 10
N = 10;

% length (# rows) of matrix
Lx = size(x, 1);

% get filter
[h1, delay] = getMatrixFIR(p, q, Lx, N, bta);

tic
% ----  HERE'S THE CALL TO UPFIRDN  ----------------------------
y1 = upfirdn(x, h1, p, q);

y1out = stripZeros(y1, p, q, Lx, delay);
t1 = toc
%%

tic

% beta value for kaiser window (default)
bta = 5;
% antialiasing filter will be order 2 X n X max(p, q)
% default for N is 10
N = 10;

pqmax = max(p,q);
% design filter
if( N>0 )
	fc = 1/2/pqmax;
	L = 2*N*pqmax + 1;
	h = firls( L-1, [0 2*fc 2*fc 1], [1 1 0 0]).*kaiser(L,bta)' ;
	h = p*h/sum(h);
else
	L = p;
	h = ones(1,p);
end

% matrix?
Lhalf = (L-1)/2;
isvect = any(size(x)==1);
if isvect
    Lx = length(x);
else
    Lx = size(x, 1);
end

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
h2 = [h zeros(1,nz1)];

% ----  HERE'S THE CALL TO UPFIRDN  ----------------------------
y = upfirdn(x,h,p,q);

% Get rid of trailing and leading data so input and output signals line up
% temporally:
Ly = ceil(Lx*p/q);  % output length
% Ly = floor((Lx-1)*p/q+1);  <-- alternately, to prevent "running-off" the
%                                data (extrapolation)
if isvect
    y(1:delay) = [];
    y(Ly+1:end) = [];
else
    y(1:delay,:) = [];
    y(Ly+1:end,:) = [];
end

t2 = toc

y2out = y;
% h([1:nz (end-nz1+1):end]) = [];  % get rid of leading and trailing zeros 
                                 % in case filter is output
%% use resample

tic
[y3out, h3] = resample(x, p, q);
t3 = toc
%%
figure(4)
plot(y1out(:, 1), '.')
hold on
plot(y2out(:, 1), '.')
plot(y3out(:, 1), '.')
hold off