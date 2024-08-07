function varargout = massFiltering(varargin)
% MASSFILTERING MATLAB code for massFiltering.fig
%      MASSFILTERING, by itself, creates a new MASSFILTERING or raises the existing
%      singleton*.
%
%      H = MASSFILTERING returns the handle to a new MASSFILTERING or the handle to
%      the existing singleton*.
%
%      MASSFILTERING('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MASSFILTERING.M with the given input arguments.
%
%      MASSFILTERING('Property','Value',...) creates a new MASSFILTERING or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before massFiltering_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to massFiltering_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help massFiltering

% Last Modified by GUIDE v2.5 12-Jun-2020 10:33:46

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @massFiltering_OpeningFcn, ...
                   'gui_OutputFcn',  @massFiltering_OutputFcn, ...
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


% --- Executes just before massFiltering is made visible.
function massFiltering_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to massFiltering (see VARARGIN)

% Choose default command line output for massFiltering
handles.output = hObject;

userdata.movieList={};
userdata.threshold=-90;
userdata.doRemoveSmallObjects=true;
userdata.removeSmallObjects=30;

set(handles.edit5, 'String', num2str(userdata.threshold));
set(handles.edit4, 'String', num2str(userdata.removeSmallObjects));
set(handles.checkbox3, 'Value', userdata.doRemoveSmallObjects);
set(handles.figure1, 'userdata', userdata);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes massFiltering wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = massFiltering_OutputFcn(hObject, eventdata, handles) 
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

function ret=createProfile(userProfile)
libalias=loadDll;
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
updateSlider(handles);


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
[FileName,PathName,~] = uigetfile({'*.tif'}, 'Select any image of a sequence', 'MultiSelect', 'on');

if isnumeric(FileName)
    return;
end

if numel(userdata.movieList)==0
    userdata.movieList={};
end

if ~iscell(FileName)
    FileName={FileName};
end

for i=1:numel(FileName)
    filename=fullfile(PathName, FileName{i});
    userdata.movieList=cat(1, userdata.movieList, {filename});
end

userdata.movieList=unique(userdata.movieList);

userdata.movieList=unique(userdata.movieList);
handles.listbox1.String=userdata.movieList;
set(handles.figure1,'userdata', userdata);
updateSlider(handles);


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
if handles.listbox1.Value>numel(userdata.movieList)
    handles.listbox1.Value=1;
end
handles.listbox1.String=userdata.movieList;
set(handles.figure1,'userdata', userdata);
updateSlider(handles)

% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1



% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.figure1, 'pointer', 'watch')
drawnow;

userdata=get(handles.figure1, 'userdata');

tag=get(handles.edit6, 'String');
if numel(tag)==0
    tag='_filt';
    set(handles.edit6, 'String', tag);
end
if userdata.doRemoveSmallObjects
    removeSmallObjects=userdata.removeSmallObjects;
else
    removeSmallObjects=0;
end

medfiltMask=str2num(get(handles.edit8, 'String'));
sigma=str2num(get(handles.edit7, 'String'));
watershedDescend=str2num(get(handles.edit9, 'String'));
doGridNorm=(handles.radiobutton4.Value>0);
doTopHat=(handles.radiobutton3.Value>0);
disksize = str2num(get(handles.edit12, 'String')); % diameter of the rolling ball, should be smaller than most cells
ithresh = str2num(get(handles.edit10, 'String')); % intensity threshold
gthresh = str2num(get(handles.edit11, 'String')); % gradient threshold

for selected=1:numel(userdata.movieList)
    set(handles.figure1, 'Name', ['Processing ', num2str(selected), '/', num2str(numel(userdata.movieList))]);
    drawnow;
    %im3Out=processMassBackground(filename, threshold, removeSmallObjects, sigma, medfiltMask, tag)
    processMassBackground(userdata.movieList{selected}, userdata.threshold, removeSmallObjects, sigma, medfiltMask, tag, watershedDescend, doGridNorm, doTopHat, disksize, ithresh, gthresh);
    
end

set(handles.figure1, 'Name', ['Finished ', num2str(numel(userdata.movieList)), ' movies']);

set(handles.figure1, 'pointer', 'arrow')


        

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


function im=autoCrop(im)

BW = im == max(max(im));
if mean(BW(:))>0.01
    [row, col] = find(~BW);
    rect = [min(col) min(row) max(col)-min(col) max(row)-min(row)];
    im = imcrop(im, rect);
end



% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

userdata=get(handles.figure1, 'userdata');
if numel(userdata.movieList)
    selected=get(handles.listbox1, 'Value');
    if selected>=1
        if numel(userdata.movieList{selected})>0
            [im, map]=imread(userdata.movieList{selected});
            
            if handles.radiobutton4.Value>0
                im=autoCrop(im);
                im2=-getSimpleNormalization(double(im), str2num(handles.edit7.String), true);
                im=im2;
            end
            
            med=double(median(im(:)));
            userdata.threshold=str2num(get(handles.edit5, 'String'));
            newThresh=thresh_tool(im, map, userdata.threshold+med, 'Set cell detection threhsold');
            if numel(newThresh)>0
                newThresh=newThresh-med;
                userdata.threshold=newThresh;
                set(handles.edit5, 'String', num2str(newThresh));
                
            end
            set(handles.figure1, 'userdata', userdata);
        end
    else
        msgbox('Select a movie first', 'modal')
    end
else
    msgbox('Load a movie first', 'modal')
end


% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

userdata=get(handles.figure1, 'userdata');

removeSmallObjects=str2num(get(handles.edit4, 'String'));
userdata.removeSmallObjects=removeSmallObjects;
doRemoveSmallObjects=get(handles.checkbox3, 'Value');
userdata.doRemoveSmallObjects=doRemoveSmallObjects;
if doRemoveSmallObjects==false
    removeSmallObjects=0;
end
medfiltMask=str2num(get(handles.edit8, 'String'));
sigma=str2num(get(handles.edit7, 'String'));

watershedParameter=str2num(get(handles.edit9, 'String'));
useGridNorm=(handles.radiobutton4.Value>0);
useTopHat=(handles.radiobutton3.Value>0);
threshold=str2num(get(handles.edit5, 'String'));
if numel(userdata.movieList)
    selected=get(handles.listbox1, 'Value');
    if selected>=1
        
        if numel(userdata.movieList{selected})>0
            val =getSliderVal(handles);
            [im, map]=imread(userdata.movieList{selected}, val);
            
            im=double(im);
            
            if useGridNorm
                % do crop

                im = autoCrop(im);
                im=getGridNormalization(im, sigma, threshold, removeSmallObjects);
            end
            
            if useTopHat
                
                im = autoCrop(im);
                disksize = str2num(get(handles.edit12, 'String')); % diameter of the rolling ball, should be smaller than most cells
                ithresh = str2num(get(handles.edit10, 'String')); % intensity threshold
                gthresh = str2num(get(handles.edit11, 'String')); % gradient threshold
                im=getTopHatNormalization(im,disksize,ithresh,gthresh, true);
            end
            
            med=double(median(im(:)));
            th=userdata.threshold+med;

            
            [mask,~, ~]=autotrackParam(im, [], medfiltMask, 0, 0, removeSmallObjects, sigma, th, watershedParameter);
            %mask=maskRefilterMask(im, mask, medfiltMask, watershedParameter, th);
%             wat=watershed(imf);
%             statsWat=regionprops(wat, imf, 'MaxIntensity', 'PixelValues', 'PixelIdxList');
%             toRemove=[];
%             for wi=1:numel(statsWat)
%                 maxi=statsWat(wi).MaxIntensity;
%                 v=find(statsWat(wi).PixelValues>=maxi-10);
%                 toRemove=cat(1, toRemove, statsWat(wi).PixelIdxList(v));
%                 %wat(statsWat(wi).PixelIdxList(v))=
%             end
%             wat(toRemove)=0;
%             
%             wat2=removeSmallIntensity(wat, mask, 0);
%             wat2=bwlabel(max(wat2, mask)>0);
%             mask=wat2;
            
            
%             if userdata.doRemoveSmallObjects
%                 mask=removeSmall(mask, userdata.removeSmallObjects);
%             else
%                 mask=bwlabel(mask);
%             end
%             mask=imdilate(mask, ones(3));
            f=figure;
            set(f, 'Name', ['Frame: ',num2str(val),  ' Nb objects: ', num2str(max(mask(:)))]);
            
            im(mask==0)=1000;
            imagesc(im);
            colorbar
            colormap gray;
            drawRegions(mask, f, 0);
            %drawRegions(wat2, f, 0);
            %drawRegions(wat, f, 0);
            daspect([1 1 1])
            
            set(handles.figure1, 'userdata', userdata);
            
            
            
            
            %imagesc(f, mask);
%             if numel(removeSmallObjects)>0
%                 userdata.removeSmallObjects=removeSmallObjects;
%                 userdata.doRemoveSmallObjects=true;
%                 set(handles.checkbox3, 'Value', true);
%                 set(handles.edit4, 'String', num2str(removeSmallObjects));
%             end
        end
    else
        msgbox('Select a movie first', 'modal')
    end
else
    msgbox('Load a movie first', 'modal')
end


function im2=getSimpleNormalization(im, sigma, isNegative)
im2=im;
if isNegative
    imBG=imopen(im2, strel('disk', sigma*5));
else
    imBG=imclose(im2, strel('disk', sigma*5));
end
im2=im2-imBG;


% function imNorm=getTopHatNormalization(im, sigma, threshold, minSize)
% 
% im2=-getSimpleNormalization(im, sigma, true);
% 
% med=(median(im2(:)));
% %mask=im<threshold+med;
% level=threshold+med;
% 
% maskCells=(im2<level);
% maskCellsInit=maskCells;
% %maskCells=removeSmall(maskCells, minSize);
% intesnities=im2;
% if sigma>0
%     intesnities=imfilter(double(intesnities), fspecial('gaussian', (sigma)*6+1, sigma),    'symmetric', 'same'); 
% end
% %intesnities(maskCells==0)=-Inf;
% maskCells2=uint16(watershed(intesnities));
% maskCells2(maskCells==0)=0;
% maskCells=maskCells2;
% 
% maskCells=removeSmallIntensity(maskCells, maskCellsInit, 0);
% maskCells=removeSmall(maskCells, minSize);
% maskCells=imdilate(maskCells>0, strel('disk', round(sigma)));
% background=fillInMaskExtrapolation(im, maskCells, 2);
% imNorm=background-im+1000;

function imNorm=getGridNormalization(im, sigma, threshold, minSize)

im2=-getSimpleNormalization(im, sigma, true);

med=(median(im2(:)));
%mask=im<threshold+med;
level=threshold+med;

maskCells=(im2<level);
maskCellsInit=maskCells;
%maskCells=removeSmall(maskCells, minSize);
intesnities=im2;
if sigma>0
    intesnities=imfilter(double(intesnities), fspecial('gaussian', (sigma)*6+1, sigma),    'symmetric', 'same'); 
end
%intesnities(maskCells==0)=-Inf;
maskCells2=uint16(watershed(intesnities));
maskCells2(maskCells==0)=0;
maskCells=maskCells2;

maskCells=removeSmallIntensity(maskCells, maskCellsInit, 0);
maskCells=removeSmall(maskCells, minSize);
maskCells=imdilate(maskCells>0, strel('disk', round(sigma)));
background=fillInMaskExtrapolation(im, maskCells, 2);
imNorm=background-im+1000;


function mask=maskRefilterMask(im, mask, medfiltMask, offsetIntensity, originalThreshold)
if medfiltMask>0
    imf=medfilt2(im, [medfiltMask,medfiltMask], 'symmetric');
else
    imf=imInit;
end

    
wat=watershed(imf);
statsWat=regionprops(wat, imf, 'PixelValues', 'PixelIdxList');
toRemove=[];
for wi=1:numel(statsWat)
    %maxi=statsWat(wi).MaxIntensity;
    %v=find(statsWat(wi).PixelValues>=maxi-offsetIntensity);
    v=find(statsWat(wi).PixelValues>=originalThreshold+offsetIntensity);
    toRemove=cat(1, toRemove, statsWat(wi).PixelIdxList(v));
    %wat(statsWat(wi).PixelIdxList(v))=
end
wat(toRemove)=0;

mask=removeSmallIntensity(wat, mask, 0);
mask=bwlabel(max(mask, mask)>0);



function im3Out=processMassBackground(filename, threshold, removeSmallObjects, sigma, medfiltMask, tag, watershedParameter, doGridNorm, doTopHat, topHatRadius, topHatThreshold, topHatGradianThreshold)

im3=imreadTIFF3D(filename);
im3OutArr={};
maskOutArr={};


areaArr={};
frameArr={};
intensityArr={};
cellIndexArr={};
imArr={};
parfor z=1:size(im3, 3)
    im=double(im3(:,:,z));
    
    if doGridNorm
        im = autoCrop(im);  
        im=getGridNormalization(im, sigma, threshold, removeSmallObjects);
        imArr{z}=im;
    end

    if doTopHat
        im = autoCrop(im);  
        im=getTopHatNormalization(im,topHatRadius, topHatThreshold, topHatGradianThreshold, false);
        imArr{z}=im;
    end
    
    med=(median(im(:)));
    %mask=im<threshold+med;

    th=threshold+med;
    %mask=im<userdata.threshold+med;

    [mask,~, ~]=autotrackParam(im, [], medfiltMask, 0, 0, removeSmallObjects, sigma, th, watershedParameter);

    
    
    
    stats=regionprops(mask, im, 'Area', 'MeanIntensity');
    areaLoc=[stats.Area]';
    
    intensityLoc=[stats.MeanIntensity]'-med;
    
    frameLoc=ones(numel(areaLoc), 1)*z;
    areaArr{z}=areaLoc;
    frameArr{z}=frameLoc;
    intensityArr{z}=intensityLoc;
    cellIndexArr{z}=(1:numel(areaLoc))';
    im2=im;
    im2(mask==0)=1000;
    %if doGridNorm
    %    imArr{z}(mask==0)=1000;
    %end
    im3OutArr{z}=im2;
    maskOutArr{z}=mask;
end
 
frame=[]; 
area=[];
intensity=[];
cellIndex=[];
maskOut=zeros(size(im3OutArr{1}, 1), size(im3OutArr{1}, 2), size(im3, 3));
im3Out=maskOut;
imOut=maskOut;
for z=1:size(im3, 3)
    im3Out(:,:,z)=im3OutArr{z};
    maskOut(:,:,z)=maskOutArr{z};
    if doGridNorm || doTopHat
        imOut(:,:,z)=imArr{z};
    end
    area=cat(1, area, areaArr{z});
    frame=cat(1, frame, frameArr{z});
    intensity=cat(1, intensity, intensityArr{z});
    cellIndex=cat(1, cellIndex, cellIndexArr{z});
    
end

[path, file, ext]=fileparts(filename);
filenameNew=fullfile(path, [file, tag, ext]);
imwriteTIFF3D(uint16(im3Out), filenameNew, false);

if doGridNorm
    filenameNew=fullfile(path, [file, '_GridCleaned', ext]);
    imwriteTIFF3D(uint16(imOut), filenameNew, false);
end
if doTopHat
    filenameNew=fullfile(path, [file, '_TopHat', ext]);
    imwriteTIFF3D(uint16(imOut), filenameNew, false);
end

filenameNew=fullfile(path, [file, '_mask', ext]);
imwriteTIFF3D(uint16(maskOut), filenameNew, false);
intensityMean=intensity;
intensitySum=intensityMean.*area;
t=table(frame, cellIndex, area, intensityMean, intensitySum);
filenameNew=fullfile(path, [file, '_mask.csv']);
writetable(t, filenameNew);



function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double


% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double


% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double


% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
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


% --- Executes on button press in checkbox5.
function checkbox5_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox5


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

updateSlider(handles)
handles.text13.String=num2str(getSliderVal(handles));



function val=getSliderVal(handles)
val =handles.slider1.Value;
val=max(1, round(val));
val=min(handles.slider1.Max, val);

function updateSlider(handles)
userdata=get(handles.figure1, 'userdata');
if numel(userdata.movieList)>0
    selected=get(handles.listbox1, 'Value');
    if selected>=1
        if numel(userdata.movieList{selected})>0
            info=imfinfo(userdata.movieList{selected});
            val=handles.slider1.Value;
            if numel(info)<val
                val=0;
                handles.slider1.Value=val;
                handles.text13.String=num2str(val);
            end
            handles.slider1.Max=numel(info);
            handles.text14.String=num2str(numel(info));
        end
    end
end


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit10_Callback(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit10 as text
%        str2double(get(hObject,'String')) returns contents of edit10 as a double


% --- Executes during object creation, after setting all properties.
function edit10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit11_Callback(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit11 as text
%        str2double(get(hObject,'String')) returns contents of edit11 as a double


% --- Executes during object creation, after setting all properties.
function edit11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit12_Callback(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit12 as text
%        str2double(get(hObject,'String')) returns contents of edit12 as a double


% --- Executes during object creation, after setting all properties.
function edit12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
