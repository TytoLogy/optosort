% Generate and plot data
Fs = 1000;
dt = 1/Fs;
t = 0:dt:5;
y = sin(2*pi*10*t);


%{
Syntax:
scrollHandles = scrollplot(plotHandles, propName,propValue,...)

Examples:
>> scrollplot; % add scroll sub-window to the current axes (gca)
>> scrollplot(plot(xdata,ydata), 'WindowSize',50); % plot with initial zoom
>> scrollplot('Min',20, 'windowsize',70); % add x-scroll to current axes
>> scrollplot([h1,h2], 'axis','xy'); % scroll both X&Y of 2 plot axes
>> scrollplot('axis','xy', 'minx',20, 'miny',10); % separate scroll minima
%}

plot(t, y)
Hx = scrollplot(gca, 'Axes', 'X', 'WindowSize', 5);
