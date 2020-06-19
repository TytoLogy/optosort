
%------------------------------------------------------------------------
%% check events
%------------------------------------------------------------------------


events = nD.events;

sendmsg('event names')
fprintf('index\tevent name\n');
for n = 1:length(nD.events)
	fprintf('%d\t%s\n', n, nD.events{n}.name);
	ev_names{n} = nD.events{n}.name;
end
fprintf('\n');

% first ten
sendmsg('events 1:10')
fprintf('\tstartsweep\t\tendsweep\t\tdiff\n');
for n = 1:10
	fprintf('%d\t%.6f\t\t%.6f\t\t%.6f\n', ...
					n, ...
					nD.events{1}.timestamps(n), ...
					nD.events{2}.timestamps(n), ...
					nD.events{2}.timestamps(n) - nD.events{1}.timestamps(n) ...
				);
end


% last ten
nevents = length(nD.events{1}.timestamps);

sendmsg('events (end-10):end:')
fprintf('\tstartsweep\t\tendsweep\t\tdiff\n');
for n = (nevents-10):nevents
	fprintf('%d\t%.6f\t\t%.6f\t\t%.6f\n', ...
					n, ...
					nD.events{1}.timestamps(n), ...
					nD.events{2}.timestamps(n), ...
					nD.events{2}.timestamps(n) - nD.events{1}.timestamps(n) ...
				);
end


% difference between startsweep(2) and endsweep(1)
sendmsg('difference between startsweep(2) and endsweep(1)');
fprintf('startsweep(2):\t%.6f\n', nD.events{1}.timestamps(2))
fprintf('endsweep(1):\t%.6f\n', nD.events{2}.timestamps(1))
fprintf('difference:\t%.6f\n', ...
			nD.events{1}.timestamps(2) - nD.events{2}.timestamps(1));
