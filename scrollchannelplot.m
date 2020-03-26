function varargout = scrollchannelplot( data )
% from https://stackoverflow.com/questions/47250746/matlab-scrolling-plot

%{

I have an EEG data base that I would like to plot. The database is a
19*1000*134 matrix, with:

19 being the number of channel. On a first approach, I'm working with only
one channel. 1000 the size of a sample (1000 points for a sampling rate of
500 Hz, i.e. 2 sec of data) 134 the number of epochs (number of different 2
second experience)


%}
% data = rand(19,1000,134);
f = figure('Units','Normalized','Position',[0.25 0.25 0.5 0.5]);
a =   axes('Units','Normalized','Position',[0.05 0.15, 0.75 0.75]);
s =   uicontrol(f, 'Style','Slider', 'Units','Normalized', ...
							'Position',[0.05 0.025, 0.75 0.05],...
                   'Min',1,'Max',size(data,3),'Value',1, ...
						 'Callback',{@sliderChange,a} );
l =   uicontrol(f, 'Style','listbox','Units','Normalized', ...
							'Position',[0.85 0.15, 0.1, 0.75],...
                   'String',cellstr(num2str([1:size(data,1)]')), ...
						 'Callback',{@changeChannel,a,s,data} );

stepSize = 1/(s.Max - s.Min);
s.SliderStep = [stepSize 2*stepSize];               
changeChannel(l,[],a,s,data)

if nargout
	varargout{1} = f;
	varargout{2} = a;
	varargout{3} = s;
	varargout{4} = l;
end

end

function changeChannel(l,evtData,a,s,data)
	cla(a);
	chanNum = str2double(l.String{l.Value});
	%500Hz
	sR = 500;  
	%Reshape each epoch into a column
	tempData = reshape(data(chanNum,:,:),[],size(data,3)); 
	%Build time array
	tempTime = [0:1/sR:(size(data,2)-1)/sR]' + (0:1:size(data,3)-1)*2; 
	%plot all the lines
	plot(a,tempTime,tempData)
	%Rest Slider Position
	s.Value = 1; 
end

function sliderChange(s,evtData,a)
	viewI = round(s.Value);
	disp(['Current Epoch: ' num2str(viewI)])
	xlim(a,[(viewI-1)*2 viewI*2] + [-.1 .1])
end