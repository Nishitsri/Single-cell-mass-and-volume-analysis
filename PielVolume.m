
function varargout = PielVolume(varargin) 
% PIELVOLUME MATLAB code for PielVolume.fig
%      PIELVOLUME, by itself, creates a new PIELVOLUME or raises the existing
%      singleton*.
%
%      H = PIELVOLUME returns the handle to a new PIELVOLUME or the handle to
%      the existing singleton*.
%
%      PIELVOLUME('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PIELVOLUME.M with the given input arguments.
%
%      PIELVOLUME('Property','Value',...) creates a new PIELVOLUME or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before PielVolume_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to PielVolume_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help PielVolume

% Last Modified by GUIDE v2.5 15-Feb-2016 12:52:18

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @PielVolume_OpeningFcn, ...
                   'gui_OutputFcn',  @PielVolume_OutputFcn, ...
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


% --- Executes just before PielVolume is made visible.
function PielVolume_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to PielVolume (see VARARGIN)

% Choose default command line output for PielVolume
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
update(hObject);

% UIWAIT makes PielVolume wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = PielVolume_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Normalization();


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.figure1, 'pointer', 'watch') 
drawnow;
manualSegmentation();
% your computation
set(handles.figure1, 'pointer', 'arrow')

% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
global currentParameters
currentParameters=ParametersClass;
currentParameters.load(@update, hObject);


function update(hObject)

global currentParameters
%handles=guidata(gcbo);
%handles=guidata(hObject);
handles=guihandles(hObject);
if currentParameters.isNormalizationDone()
   set(handles.pushbutton3, 'Enable', 'on');
else
    set(handles.pushbutton3, 'Enable', 'off');
end
if currentParameters.isSegmentationDone()
   set(handles.pushbutton4, 'Enable', 'on');
else
    set(handles.pushbutton4, 'Enable', 'off');
end


% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global currentParameters
currentParameters.save();



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[FileName,PathName,~] = uigetfile('*.tif','Select movie file');
global currentParameters
% if iscell(FileName)
%     currentParameters.setRawMovieFilename(fullfile(PathName, FileName{1}));
%     currentParameters.setAdditionalRawMovieFilename(FileName(2:end),PathName);
% else
    if numel(FileName)>1
        currentParameters.setRawMovieFilename(fullfile(PathName, FileName));
        handles.edit2.String=FileName;
        handles.listbox1.String=currentParameters.getLastFileName();
        update(hObject);
    end
% end

% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1
global currentParameters
 i=get(hObject, 'Value');
if i>=0, 
    filename=currentParameters.selectInLast(i);
    if exist(filename, 'file')
        % choose the selected file
        currentParameters.setRawMovieFilename(filename);
        [~,name,~] = fileparts(filename);
        handles.edit2.String=name;
        handles.listbox1.Value=1;
        handles.listbox1.String=currentParameters.getLastFileName();
        
        update(hObject);
    end
end


% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
global currentParameters
hObject.String=currentParameters.getLastFileName();


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global currentParameters
i=get(handles.listbox1, 'Value');
if i>=0, 
    filename=currentParameters.selectInLast(i);
    currentParameters.removeFilenameToLastFiles(filename);
    handles.listbox1.Value=1;
    handles.listbox1.String=currentParameters.getLastFileName();

    update(hObject);
end
