function varargout = massManagment(varargin)
% MASSMANAGMENT MATLAB code for massManagment.fig
%      MASSMANAGMENT, by itself, creates a new MASSMANAGMENT or raises the existing
%      singleton*.
%
%      H = MASSMANAGMENT returns the handle to a new MASSMANAGMENT or the handle to
%      the existing singleton*.
%
%      MASSMANAGMENT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MASSMANAGMENT.M with the given input arguments.
%
%      MASSMANAGMENT('Property','Value',...) creates a new MASSMANAGMENT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before massManagment_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to massManagment_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help massManagment

% Last Modified by GUIDE v2.5 19-Mar-2019 16:30:22

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @massManagment_OpeningFcn, ...
                   'gui_OutputFcn',  @massManagment_OutputFcn, ...
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


% --- Executes just before massManagment is made visible.
function massManagment_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to massManagment (see VARARGIN)

% Choose default command line output for massManagment
handles.output = hObject;
userdata.refImage='';
userdata.phasicsAPI='C:\Program Files\SID4_SDK_x86_64\MatLabInterface';
if ~exist(userdata.phasicsAPI, 'dir')
    userdata.phasicsAPI='';
end
handles.edit3.String=userdata.phasicsAPI;
userdata.movieList={};
userdata.userProfile='';

set(handles.figure1, 'userdata', userdata);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes massManagment wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = massManagment_OutputFcn(hObject, eventdata, handles) 
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

userdata=get(handles.figure1, 'userdata');
[FileName,PathName,~] = uigetfile({'*.txt'}, 'Select user profile', userdata.userProfile, 'MultiSelect', 'off');

if isnumeric(FileName)
    return;
end
userdata.userProfile=fullfile(PathName, FileName);
handles.edit1.String=userdata.userProfile;
set(handles.figure1,'userdata', userdata);

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.figure1, 'pointer', 'watch')
drawnow;

userdata=get(handles.figure1, 'userdata');
addpath(userdata.phasicsAPI);
[FileName,PathName,~] = uiputfile({'*.txt'}, 'Select user profile', userdata.userProfile);
if isnumeric(FileName)
    set(handles.figure1, 'pointer', 'arrow')
    return;
end
userProfile=fullfile(PathName, FileName);
ret=createProfile(userProfile);
if ret==0
    userdata.userProfile=userProfile;
    handles.edit1.String=userProfile;
    set(handles.figure1, 'userdata', userdata);
end

set(handles.figure1, 'pointer', 'arrow')


function libalias=loadDll()

warning('off','MATLAB:loadlibrary:parsewarnings');
libalias=SID4SDKx86_LoadDll;
warning('on','MATLAB:loadlibrary:parsewarnings');


function ret=createProfile(userProfile)
libalias=loadDll();
currentPath=pwd;
UserProfileInit='';
[ ErrorCode,SDKSession ] = SID4SDKx86_OpenSID4( libalias,UserProfileInit);
ret=checkError(ErrorCode, SDKSession, libalias, currentPath); if ret<0, return; end
cd(currentPath);
SNPhasics = 'SID4-322';
userName='Piel';

answer = inputdlg({'Serial number','User name'}, 'New user',1,{SNPhasics, userName});
try 
    SNPhasics=strtrim(answer{1});
    userName=strtrim(answer{2});
catch
    msgbox('Error parsing arguments');
    ret=-1;
    return;
end

[ ErrorCode ]=SID4SDKx86_NewUserProfile_circle( SDKSession, SNPhasics, userProfile,userName, '', 0, '', 240, 240, 200 );
ret=checkError(ErrorCode, SDKSession, libalias, currentPath); if ret<0, return; end
cd(currentPath);
SID4SDKx86_CloseSID4( SDKSession );
cd(currentPath);
ret=0;

function checkTIFFTAG_SampleFormat(filename)
t = Tiff(filename,'r+');
t.setTag('SampleFormat',1); %1=Uint, 2=Int, 3=IEEEFP(floatting)
t.rewriteDirectory();
t.close();


function closeSDK(SDKSession, libalias)
SID4SDKx86_CloseSID4(SDKSession);
unloadlibrary(libalias);


function ret=checkError(ErrorCode, SDKSession, libalias, currentPath)
ret=0;
if ErrorCode~=0
    msg=SID4SDKx86_perror(SDKSession, ErrorCode);
    closeSDK(SDKSession, libalias);
    msgbox(msg, 'Error message');
    cd(currentPath);
    ret=-1;
end


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


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

userdata=get(handles.figure1, 'userdata');
[FileName,PathName,~] = uigetfile({'*.tif'}, 'Select reference image', userdata.refImage, 'MultiSelect', 'off');

if isnumeric(FileName)
    return;
end
userdata.refImage=fullfile(PathName, FileName);
handles.edit2.String=userdata.refImage;
set(handles.figure1,'userdata', userdata);


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1


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


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

userdata=get(handles.figure1, 'userdata');
[FileName,PathName,~] = uigetfile({'*.tif'}, 'Select any image of a sequence', 'MultiSelect', 'off');

if isnumeric(FileName)
    return;
end

if numel(userdata.movieList)==0
    userdata.movieList={};
end
[~, file, ext]=fileparts(FileName);
underscorePos=strfind(file, '_t');
if numel(underscorePos)==0,
    msgbox(['Unable to load the image sequence "', FileName,  '". The file should respect "XXXXX_t[num].tif" where [num] is the frame number']);
    return;
end
filename=fullfile(PathName, [file(1:underscorePos(end)),'t*',ext]);

list=dir(filename);
indexes=[];% "zeros(numel(list), 1);
for i=1:numel(list),
    [~, file, ~]=fileparts(list(i).name);
    
    num=str2num(file(underscorePos(end)+2:end));
    if numel(num)>0
        indexes=[indexes,num];
    end
end
[~, firstIndex]=min(indexes);
filename0=fullfile(PathName, [list(firstIndex).name]);
userdata.movieList=cat(1, userdata.movieList, {filename0});

userdata.movieList=unique(userdata.movieList);
handles.listbox1.String=userdata.movieList;
set(handles.figure1,'userdata', userdata);



% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
userdata=get(handles.figure1, 'userdata');

if numel(handles.listbox1.Value)==0
    return;
end
if handles.listbox1.Value<1
    return;
end
userdata.movieList(handles.listbox1.Value)=[];
handles.listbox1.String=userdata.movieList;
set(handles.figure1,'userdata', userdata);


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1

function ret = testIfAdministratorMode(phasicsAPI)
ret=-1;
try
    filename=fullfile(phasicsAPI, 'testAdministratorMode');
    mkdir(filename);
    if exist(filename, 'dir')
        rmdir(filename);
        ret=1;
    end
catch
    disp('No admin rights');
end

% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.figure1, 'pointer', 'watch')
drawnow;

doCrop=(handles.checkbox2.Value==1);
userdata=get(handles.figure1, 'userdata');

if ~exist(userdata.refImage, 'file')
    msgbox('Define reference image first');
    set(handles.figure1, 'pointer', 'arrow')
    return;
end
if ~exist(userdata.userProfile, 'file')
    msgbox('Define user profile first');
    set(handles.figure1, 'pointer', 'arrow')
    return;
end
if ~exist(userdata.phasicsAPI, 'dir')
    msgbox('Define Phasics API path first');
    set(handles.figure1, 'pointer', 'arrow')
    return;
end
if testIfAdministratorMode(userdata.phasicsAPI)<0
    msgbox('Matlab must be launch with admin rights');
    set(handles.figure1, 'pointer', 'arrow')
    return;
end
addpath(userdata.phasicsAPI);
if handles.checkbox4.Value==0
    ret=changeRef(userdata.userProfile, userdata.refImage, userdata.movieList{1});
    if ret<0 
        set(handles.figure1, 'pointer', 'arrow')
        return;
    end
end
lambda=str2num(handles.edit5.String);
objectSize=str2num(handles.edit4.String);
if handles.checkbox3.Value==0
    objectSize=0;
end

if handles.checkbox4.Value==1
    refROI=[];
    for i=1:numel(userdata.movieList)
        [centerX, centerY, height, width]=getROIRef(userdata.movieList{i});
        refROI=cat(1, refROI, struct('centerX', centerX, 'centerY', centerY, 'height', height, 'width', width));

    end
end

for i=1:numel(userdata.movieList)
    if handles.checkbox4.Value==1
        ret=changeRef(userdata.userProfile, userdata.refImage, refROI(i));
        if ret<0 
            set(handles.figure1, 'pointer', 'arrow')
            return;
        end
    end
    ret=processMovie(userdata.movieList{i}, userdata.userProfile, doCrop, objectSize, lambda);
    if ret<0 
        msgbox(['Fail to process in movie ', userdata.movieList{i}]);
        set(handles.figure1, 'pointer', 'arrow')
        return;
    end
end

set(handles.figure1, 'pointer', 'arrow')

function ret=processMovie(movieFilename, pathUserProfile, doCrop, objectSize, lambda)

currentFolder=pwd;
[PathName, file, ext]=fileparts(movieFilename);
analysisFolder=fullfile(PathName, 'Analysis');
if ~exist(analysisFolder, 'dir'),
    mkdir(analysisFolder)
end
underscorePos=strfind(file, '_t');
if numel(underscorePos)==0,
    msgbox(['Unable to load the image sequence "', FileName,  '". The file should respect "XXXXX_t[num].tif" where [num] is the frame number']);
    ret=-1;
    return;
end
fileTag=file(1:underscorePos(end));
filename=fullfile(PathName, [file(1:underscorePos(end)),'t*',ext]);

list=dir(filename);


indexes=[];% "zeros(numel(list), 1);
for i=1:numel(list),
    [~, file, ~]=fileparts(list(i).name);
    
    num=str2num(file(underscorePos(end)+2:end));
    if numel(num)>0
        indexes=[indexes,num];
    end
end

[~, sortIndex]=sort(indexes);
filenames={};
for i=1:numel(sortIndex)
    filenames=cat(1, filenames, {fullfile(PathName, list(sortIndex(i)).name)});
end
libalias=loadDll();
[ErrorCode, SDKSession]=SID4SDKx86_OpenSID4(libalias,pathUserProfile);
%disp(SDKSessionPar);
ret=checkError(ErrorCode, SDKSession, libalias, currentFolder); 
if ret<0, return; end

% wavelength associated to the acquisition
%lambda= 550;

for j=1:numel(filenames)
    try
        
        tic;
        checkTIFFTAG_SampleFormat(filenames{j});
        [ ErrorCode,Phase,Intensity,~ ] = SID4SDKx86_FileAnalysis( SDKSession,filenames{j});
        ret=checkError(ErrorCode, SDKSession, libalias, currentFolder); 
        if ret<0, return; end
        
        cd(currentFolder);
        disp(['loop # ', num2str(i), ', frame # ', num2str(j), ' duration SID4SDKx86_FileAnalysis ', num2str(toc)]);
        CalcPhase = lambda.*(max(max(Phase)).*ones(size(Phase))-double(Phase));

        
        if j==1,
            mode='overwrite';
            if doCrop
                BW = CalcPhase == max(max(CalcPhase));
                [row, col] = find(~BW);
                rect = [min(col) min(row) max(col)-min(col) max(row)-min(row)];
            end
        else
            mode='append';
        end
        

        OutputPhase = fullfile(analysisFolder, [fileTag,'_Phase.tif']);
        imwrite(uint16(CalcPhase), OutputPhase, 'WriteMode', mode);
        
        if doCrop
            img = imcrop(CalcPhase, rect);
            OutputPhase_crop = fullfile(analysisFolder, [fileTag,'_Phase_Cropped.tif']);
            imwrite(uint16(img), OutputPhase_crop, 'WriteMode', mode);
            if objectSize>0
                % background calculation
                background = imopen(img,strel('disk',objectSize));
                img2 = double(background)-double(img);
                OutputPhase_clean = fullfile(analysisFolder, [fileTag,'_Phase_Cleaned.tif']);
                imwrite(uint16(img2+ones(size(img2))*1000), OutputPhase_clean, 'WriteMode', mode); 
            end
            
            OutputIntensity = fullfile(analysisFolder, [fileTag,'_Intensity.tif']);
            Intensity = imcrop(Intensity, rect);
            imwrite(uint16(Intensity/100), OutputIntensity, 'WriteMode',mode);
        
        end
        
        
        

    catch Ex
        disp(['Ex ',Ex.message]);
        ret=-1;
        closeSDK(SDKSession, libalias);
        cd(currentFolder);
        return;
    end
end
ret=0;
closeSDK(SDKSession, libalias);
cd(currentFolder);



function [centerX, centerY, height, width]=getROIRef(ImageToShowPath)


ImageToShow=imread(ImageToShowPath);
hf=figure;
imagesc(ImageToShow);
colormap gray;
h = imrect();
ROIcoord = getPosition(h);
close(hf);
centerX=uint16(ROIcoord(1)+ROIcoord(3)/2);
centerY=uint16(ROIcoord(2)+ROIcoord(4)/2);
height=uint16(ROIcoord(4));
width=uint16(ROIcoord(3));

function ret=changeRef(UserProfile, RefPath, ImageToShowPath)
currentPath=pwd;
libalias=loadDll();
tic;
[ErrorCode, SDKSession]=SID4SDKx86_OpenSID4(libalias,UserProfile);
cd(currentPath);
ret=checkError(ErrorCode, SDKSession, libalias, currentPath); 
if ret<0 
    ret=-1;
    closeSDK(SDKSession, libalias);
    return;
end

checkTIFFTAG_SampleFormat(RefPath);
ErrorCode=SID4SDKx86_ChangeReference(SDKSession,1, RefPath);
cd(currentPath);
ret=checkError(ErrorCode, SDKSession, libalias, currentPath); 
if ret<0 
    ret=-1;
    closeSDK(SDKSession, libalias);
    return;
end

t=toc;
disp(['SID4SDKx86_ChangeReference duration: ', num2str(t)]);
tic;
if ischar(ImageToShowPath)
    ImageToShow=imread(ImageToShowPath);
    hf=figure;
    imagesc(ImageToShow);
    colormap gray;
    h = imrect();
    ROIcoord = getPosition(h);
    close(hf);
    centerX=uint16(ROIcoord(1)+ROIcoord(3)/2);
    centerY=uint16(ROIcoord(2)+ROIcoord(4)/2);
    height=uint16(ROIcoord(4));
    width=uint16(ROIcoord(3));
else
    centerX=ImageToShowPath.centerX;
    centerY=ImageToShowPath.centerY;
    height=ImageToShowPath.height;
    width=ImageToShowPath.width;
end

[ ErrorCode ] = SID4SDKx86_ChangeAnalysisMask_FromROI_rectangle( SDKSession, centerY, centerX, height, width);
cd(currentPath);
ret=checkError(ErrorCode, SDKSession, libalias, currentPath); 
if ret<0 
    ret=-1;
    closeSDK(SDKSession, libalias);
    return;
end

t=toc;
disp(['SID4SDKx86_ChangeAnalysisMask_FromROI_rectangle duration: ', num2str(t)]);
closeSDK(SDKSession, libalias);
cd(currentPath);
ret=0;

        

% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in checkbox2.
function checkbox2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox2



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

userdata=get(handles.figure1, 'userdata');

pathName = uigetdir(userdata.phasicsAPI);

if isnumeric(pathName)
    return;
end
userdata.phasicsAPI=pathName;
handles.edit3.String=pathName;
set(handles.figure1,'userdata', userdata);


% --- Executes on button press in checkbox3.
function checkbox3_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox3



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox4.
function checkbox4_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox4
