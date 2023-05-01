function varargout = LAYOUT_CommonReferenceDataViewer(varargin)
%------------------------------------------------------------------------
% LAYOUT_COMMONREFERENCEDATAVIEWER MATLAB code for LAYOUT_CommonReferenceDataViewer.fig
%      LAYOUT_COMMONREFERENCEDATAVIEWER, by itself, creates a new LAYOUT_COMMONREFERENCEDATAVIEWER or raises the existing
%      singleton*.
%
%      H = LAYOUT_COMMONREFERENCEDATAVIEWER returns the handle to a new LAYOUT_COMMONREFERENCEDATAVIEWER or the handle to
%      the existing singleton*.
%
%      LAYOUT_COMMONREFERENCEDATAVIEWER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LAYOUT_COMMONREFERENCEDATAVIEWER.M with the given input arguments.
%
%      LAYOUT_COMMONREFERENCEDATAVIEWER('Property','Value',...) creates a new LAYOUT_COMMONREFERENCEDATAVIEWER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before LAYOUT_CommonReferenceDataViewer_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to LAYOUT_CommonReferenceDataViewer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
%------------------------------------------------------------------------
% See also: GUIDE, GUIDATA, GUIHANDLES
%------------------------------------------------------------------------

% Edit the above text to modify the response to help LAYOUT_CommonReferenceDataViewer

% Last Modified by GUIDE v2.5 01-May-2023 13:20:29

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @LAYOUT_CommonReferenceDataViewer_OpeningFcn, ...
                   'gui_OutputFcn',  @LAYOUT_CommonReferenceDataViewer_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT
%------------------------------------------------------------------------


%------------------------------------------------------------------------
% --- Executes just before LAYOUT_CommonReferenceDataViewer is made visible.
%------------------------------------------------------------------------
function LAYOUT_CommonReferenceDataViewer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to LAYOUT_CommonReferenceDataViewer (see VARARGIN)


fnames = fieldnames(handles);
% print all object units
fprintf('%% GUI Elements: Units\n')
for f = 1:length(fnames)
   if isprop(handles.(fnames{f}), 'Position')
      fprintf('%s.Units = ''%s'';\n', fnames{f}, handles.(fnames{f}).Units);
   end
end

% print all object positions
fprintf('%% GUI Elements: Position\n')
for f = 1:length(fnames)
   if isprop(handles.(fnames{f}), 'Position')
      pos = handles.(fnames{f}).Position;
      pstr = sprintf('%.4f ', pos);
      fprintf('%s.Position = [%s];\n', fnames{f}, pstr);
   end
end
% Choose default command line output for LAYOUT_CommonReferenceDataViewer
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);




% UIWAIT makes LAYOUT_CommonReferenceDataViewer wait for user response (see UIRESUME)
% uiwait(handles.fig);
%------------------------------------------------------------------------
%------------------------------------------------------------------------


% --- Outputs from this function are returned to the command line.
function varargout = LAYOUT_CommonReferenceDataViewer_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in trialDownButton.
function trialDownButton_Callback(hObject, eventdata, handles)
% hObject    handle to trialDownButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in trialUpButton.
function trialUpButton_Callback(hObject, eventdata, handles)
% hObject    handle to trialUpButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function trialnEdit_Callback(hObject, eventdata, handles)
% hObject    handle to trialnEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of trialnEdit as text
%        str2double(get(hObject,'String')) returns contents of trialnEdit as a double


% --- Executes during object creation, after setting all properties.
function trialnEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trialnEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
