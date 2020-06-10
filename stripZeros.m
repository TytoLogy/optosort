function y = stripZeros(y, p, q, Lx, delay)
% Uses code from RESAMPLE()

% Get rid of trailing and leading data so input and output signals line up
% temporally:
Ly = ceil(Lx*p/q);  % output length

% Ly = floor((Lx-1)*p/q+1);  <-- alternately, to prevent "running-off" the
%                                data (extrapolation)

% yy = y( (delay+1):Ly, :);

y(1:delay,:) = [];
y(Ly+1:end,:) = [];

% yy = y;
