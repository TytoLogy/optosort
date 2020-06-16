% function val = genTestData(dsize, spike, opts, Fs)
% 
dsize = [1, 10000];
Fs = 48828.125;
Channels = [4 5];

% define spike
spike.waveform = zeros(1, ms2bin(1, Fs));
spike.center = ms2bin(0.5, Fs);
spike.length = ms2bin(1, Fs);
spike.waveform(2:(spike.center - 1)) = 1;
spike.waveform((spike.center+1):spike.length) = -1;
spike.waveform(end) = 0;

% define parameters
testOpts = define_testOptsStruct(Channels);


% index into channels
c = 2;
opts = testOpts(c)

% create tmp matrix to hold units
tmpdata = zeros(opts.nUnits, max(dsize));

%
for u = 1:opts.nUnits
	% get isi
	isi = 1000 / opts.FiringRate
	% convert to bins
	isibins = ms2bin(isi, Fs)
	% create list of spike times
	startbin = ms2bin(opts.BurstOnsetTime, Fs);
	if startbin == 0
		startbin = 1;
	end
	endbin = ms2bin(opts.BurstOffsetTime, Fs);
	% spike times (in bins)
	spikebins =  startbin:isibins:endbin;
	% jitter bins if  ~= last unit
	if u ~= opts.nUnits
		spikebins = spikebins + ceil(ms2bin(30, Fs)*rand(size(spikebins)))
	end
	% set bins in tmpdata to 1
	tmpbins = zeros(dsize);
	tmpbins(spikebins) = 1;
	% convolve with waveform
	tmp = conv(tmpbins, spike.waveform);
	% assign proper sized, normalized data & scale by unit amplitude
	tmpdata(u, :) = opts.UnitAmp(u) * normalize(tmp(1:length(tmpdata)));
end

% add unit traces
val = sum(tmpdata);
% plot waveform
figure(1)
subplot(211)
plot(spike.waveform)
grid on

subplot(212)
plot(val)
%% try genTestData

val2 = genTestData(dsize, spike, testOpts(1), Fs);


% plot waveform
figure(1)
subplot(211)
plot(spike.waveform)
grid on

subplot(212)
plot(val2)