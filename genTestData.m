% function val = genTestData(dsize, spikebins, spike)
% val = zeros(dsize);
% for s = 1:length(spikebins)
% 	val( (spikebins(s)-spike.center):(spikebins(s)+(spike.length-spike.center) = spike.waveform;
% end
% 	

spikebins = 500;

dsize = [1, 12207];

Fs = 48828.125;
spike.waveform = zeros(1, ms2bin(1, Fs));
spike.center = ms2bin(0.5, Fs);
spike.length = ms2bin(1, Fs);

spike.waveform(2:(spike.center - 1)) = 1;
spike.waveform((spike.center+1):spike.length) = -1;
spike.waveform(end) = 0;

plot(spike.waveform)
grid on

tmp = zeros(dsize);
tmp(spikebins - spike.center + 1) = 1;

val = conv(tmp, spike.waveform);
plot(val)