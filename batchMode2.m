function varargout = batchMode2(varargin)
% BATCHMODE2 MATLAB code for batchMode2.fig
%      BATCHMODE2, by itself, creates a new BATCHMODE2 or raises the existing
%      singleton*.
%
%      H = BATCHMODE2 returns the handle to a new BATCHMODE2 or the handle to
%      the existing singleton*.
%
%      BATCHMODE2('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BATCHMODE2.M with the given input arguments.
%
%      BATCHMODE2('Property','Value',...) creates a new BATCHMODE2 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before batchMode2_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to batchMode2_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help batchMode2

% Last Modified by GUIDE v2.5 14-Dec-2018 18:14:01

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @batchMode2_OpeningFcn, ...
                   'gui_OutputFcn',  @batchMode2_OutputFcn, ...
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


% --- Executes just before batchMode2 is made visible.
function batchMode2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to batchMode2 (see VARARGIN)

% Choose default command line output for batchMode2
handles.output = hObject;

% Update handles structure
set(handles.uitable1, 'CellSelectionCallback', @cellSelect);
guidata(hObject, handles);
global currentParameters
currentParameters=ParametersClass;

% UIWAIT makes batchMode2 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = batchMode2_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function ret=isProcessing(hObject)

obj=get(hObject, 'UserData');
ret=false;
if isfield(obj, 'isProcessing')
    ret=obj.isProcessing;
end

function setIsProcessing(hObject, val)
obj=get(hObject, 'UserData');
obj.isProcessing=val;
set(hObject, 'UserData', obj);


function parsave(fname, imageFlat, deadZoneMask)
save(fname, 'imageFlat', 'deadZoneMask');

function setGrid(handles, row, col, str)
data=get(handles.uitable1, 'data');
data{row, col}=str;
set(handles.uitable1, 'data', data);
pause(0.1);
drawnow;

function val=getGrid(handles, row, col)
data=get(handles.uitable1, 'data');
val=data{row, col};

function processFrame(fi, f, frameCount, filename, parameters, normalizationFolder)
%f=frames(fi);
disp(['Frame ', num2str(fi), '/', num2str(frameCount)]);

[~, ~, imageFlat, ~, deadZoneMask]=backgroundEstimation2Params(filename, f, ...
    parameters.cellDiameter, parameters.noiseFactor, parameters.method, ...
    parameters.zoomFactor, parameters.distanceOfInfluence, parameters.alpha, ...
    parameters.alphaCells, 3);

parsave([normalizationFolder, '\frame', num2str(f), '.mat'], imageFlat, deadZoneMask);

function runLoadFilliation(filename)
global currentParameters
%global manualSeg
%currentParameters.setRawMovieFilename(filename);
%load(fullfile(currentParameters.trackingFolder, 'trackBeforeLoadingFilliations.mat'));
[path, file, ext]=fileparts(filename);

loadCurrentParamAndManualSeg(fullfile(path, file, 'Tracking', 'trackBeforeLoadingFilliations.mat'));
pause(0.1);
h=figure;
loadFilliation(h, fullfile(currentParameters.trackingFolder, 'trackFilliation_'));
close(h);
%save(fullfile(currentParameters.trackingFolder, 'trackFilliation.mat'), 'currentParameters','manualSeg', '-v7.3');
saveCurrentParamAndManualSeg(fullfile(currentParameters.trackingFolder, 'trackFilliation.mat'));
buildMultiPageTif(fullfile(currentParameters.trackingFolder, 'trackFilliation_*.tif'));

function runVolumeMeasurement(filename)
global currentParameters
global manualSeg
%currentParameters.setRawMovieFilename(filename);
%load(fullfile(currentParameters.trackingFolder, 'trackFilliation.mat'));
[path, file, ext]=fileparts(filename);
loadCurrentParamAndManualSeg(fullfile(path, file, 'Tracking', 'trackFilliation.mat'));
pause(0.1);
manualSeg.needUpdate=true;
updatePlot();
saveCurrentParamAndManualSeg(fullfile(currentParameters.trackingFolder, 'volumeMeasurement.mat'));
%save(fullfile(currentParameters.trackingFolder, 'volumeMeasurement.mat'), 'currentParameters','manualSeg', '-v7.3');


function runExportToCSV(filename, parametersCells)
global currentParameters
%currentParameters.setRawMovieFilename(filename);
%load(fullfile(currentParameters.trackingFolder, 'volumeMeasurement.mat'));
[path, file, ext]=fileparts(filename);
loadCurrentParamAndManualSeg(fullfile(path, file, 'Tracking', 'volumeMeasurement.mat'));
pause(0.1);

[data, headers]=exportToTable([], parametersCells.pilarHeight, parametersCells.pixelSize);
tab=cell2table(data, 'VariableNames', headers);
writetable(tab,fullfile(currentParameters.trackingFolder, 'results.csv'),'Delimiter',',');



function runTracking(filename, handles, fileId, parametersTrack, rawAdditionalMovieFilename)
global currentParameters
global manualSeg
manualSeg.cells=AllCells();
manualSeg.cellsData=[];
manualSeg.needUpdate=true;
manualSeg.stopAutoTrack=true;
manualSeg.image3DFiltered=[];

manualSeg.autoTrackInfo={};
manualSeg.sizePoly1Click=[];
if numel(manualSeg)==0
    manualSeg.cells=AllCells();
    manualSeg.cellsData=[];
    manualSeg.needUpdate=true;
    manualSeg.stopAutoTrack=true;
    manualSeg.image3DFiltered=[];
    manualSeg.autoTrackInfo={};
    manualSeg.sizePoly1Click=[];
else
    if ~isfield(manualSeg, 'cells')
        manualSeg.cells=AllCells();
        manualSeg.cellsData=[];
        manualSeg.needUpdate=true;
        manualSeg.stopAutoTrack=true;
        manualSeg.image3DFiltered=[];
        manualSeg.autoTrackInfo={};
        manualSeg.sizePoly1Click=[];
    end
end
currentParameters.setRawMovieFilename(filename);

currentParameters.openImage();

currentParameters.rawAdditionalMovieFilename={};
for i=1:numel(rawAdditionalMovieFilename)
    if numel(rawAdditionalMovieFilename{i})>0
        currentParameters.rawAdditionalMovieFilename=cat(1, currentParameters.rawAdditionalMovieFilename, rawAdditionalMovieFilename{i});
    end
end
if numel(currentParameters.rawAdditionalMovieFilename)>0
    currentParameters.rawAdditionalMovieFilename=unique(currentParameters.rawAdditionalMovieFilename);
end

manualSeg.stopAutoTrack=false;
h=figure;
autoTracking(parametersTrack.medfiltMask, parametersTrack.minDistance, parametersTrack.distanceDilation, parametersTrack.minSize, parametersTrack.sigma,parametersTrack.level, h, fullfile(currentParameters.trackingFolder, 'autotracking_'));
currentParameters.hObjectInMainFigure=[];
close(h);

saveCurrentParamAndManualSeg(fullfile(currentParameters.trackingFolder, 'trackBeforeLoadingFilliations.mat'));
%save(fullfile(currentParameters.trackingFolder, 'trackBeforeLoadingFilliations.mat'), 'currentParameters','manualSeg', '-v7.3');
buildMultiPageTif(fullfile(currentParameters.trackingFolder, 'autotracking_*.tif'));

function buildMultiPageTif(filename)
list=dir(filename);
[path, file, ext]=fileparts(filename);
newFilename=fullfile(path, [file(1:end-2), ext]);
for i=1:numel(list),
    im=imread(fullfile(path, list(i).name));
    if i==1,
        imwrite(im, newFilename);
    else
        imwrite(im, newFilename, 'WriteMode', 'append');
    end
end
for i=1:numel(list),
    delete(fullfile(path, list(i).name));
end

function runNormalization(filename, handles, fileId, parameters)
global currentParameters
currentParameters.setRawMovieFilename(filename);
%currentParameters.rawMovieFilename=fileList{fileId};
setGrid(handles, fileId, 2, 'Processing 1');
set(handles.pushbutton1, 'String', 'Cancel');
drawnow
pause(0.1) % try to force the GUI to be updated
info=imfinfo(currentParameters.rawMovieFilename);
frames=1:numel(info);

tic
currentParameters.deleteFileNormalizationDone();
rawMovieFilename=currentParameters.rawMovieFilename;
normalizationFolder=currentParameters.normalizationFolder;
if handles.checkbox2.Value==1
    alreadyExist=[];
    for fi=1:numel(frames)
        f=frames(fi);
        if exist([normalizationFolder, '\frame', num2str(f), '.mat'], 'file')
            alreadyExist=[alreadyExist, f];
        end
    end
    frames(alreadyExist)=[];
end      

frameCount=numel(frames);
if get(handles.checkbox1, 'Value')==1, % parallel computing
    parfor fi=1:numel(frames)
        processFrame(fi,frames(fi), frameCount, rawMovieFilename, parameters, normalizationFolder);
    end
else
    for fi=1:numel(frames)
        processFrame(fi,frames(fi), frameCount, rawMovieFilename, parameters, normalizationFolder);
    end
end
currentParameters.createFileNormalizationDone();
            
            
% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global currentParameters
isProc=isProcessing(hObject);
if isProc,
    setIsProcessing(hObject, false);
    set(hObject, 'String', 'Processing ...');
    return;
end
userdata=get(handles.figure1, 'UserData');

%if handles.checkbox3.Value % batrchMode
fileList=userdata.files(:,1);
initFilename='';


tic
setIsProcessing(hObject, true);
pause(0.1)
drawnow
for fileId=1:numel(fileList)
    try
        [pathLog, fileLog, extLog]=fileparts(fileList{fileId});
        filenameLog=fullfile(pathLog, fileLog, 'crash.txt');
        if exist(filenameLog, 'file')
            delete(filenameLog);
        end
        currentProc=0;
        if numel(userdata.configs)>=fileId
            if numel(userdata.configs{fileId})>0
                if getGrid(handles, fileId, 3)==1
                    parameters=userdata.configs{fileId}{1};
                    currentProc=1;
                    setGrid(handles, fileId, 2, ['Processing ', num2str(currentProc)]);
                    runNormalization(fileList{fileId}, handles, fileId, parameters);
                    setGrid(handles, fileId, 2, [num2str(currentProc), ' done']);
                end
            end

            if numel(userdata.configs{fileId})>1
                if getGrid(handles, fileId, 4)==1
                    parametersTrack=userdata.configs{fileId}{2};
                    currentProc=2;
                    setGrid(handles, fileId, 2, ['Processing ', num2str(currentProc)]);
                    runTracking(fileList{fileId}, handles, fileId, parametersTrack, userdata.files(fileId,2:end));
                    setGrid(handles, fileId, 2, [num2str(currentProc), ' done']);
                end
            end
            
            if getGrid(handles, fileId, 5)==1
                currentProc=3;
                setGrid(handles, fileId, 2, ['Processing ', num2str(currentProc)]);
                runLoadFilliation(fileList{fileId}); %fileList{fileId}, handles, fileId);
                setGrid(handles, fileId, 2, [num2str(currentProc), ' done']);
            end
            
            if getGrid(handles, fileId, 6)==1
                currentProc=4;
                setGrid(handles, fileId, 2, ['Processing ', num2str(currentProc)]);
                runVolumeMeasurement(fileList{fileId});
                setGrid(handles, fileId, 2, [num2str(currentProc), ' done']);
            end
            
            if numel(userdata.configs{fileId})>1
                if getGrid(handles, fileId, 7)==1
                    currentProc=5;
                    setGrid(handles, fileId, 2, ['Processing ', num2str(currentProc)]);
                    parametersTrack=userdata.configs{fileId}{2};
                    runExportToCSV(fileList{fileId}, parametersTrack); %fileList{fileId}, handles, fileId);
                    setGrid(handles, fileId, 2, [num2str(currentProc), ' done']);
                end
            end
        end
    catch Ex
        str=[Ex.message, char(10)];
        disp(Ex.message);
        for i=1:numel(Ex.stack)
            str=[str, 'file: ', Ex.stack(i).file, char(10)];
            str=[str, 'name: ', Ex.stack(i).name, char(10)];
            str=[str, 'line: ', num2str(Ex.stack(i).line), char(10)];
            disp(Ex.stack(i));
        end
        [pathLog, fileLog, extLog]=fileparts(fileList{fileId});
        filenameLog=fullfile(pathLog, fileLog, 'crash.txt');
        fid=fopen(filenameLog, 'wt');
        fwrite(fid, str);
        fclose(fid);
        setGrid(handles, fileId, 2, [num2str(currentProc), ' fail!']);
    end
end

setIsProcessing(hObject, false);
set(hObject, 'String', 'Run batch');
disp(['Info: Done in ', num2str(toc), 's']);




% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname, filterindex] = uigetfile( ...
{  '*.xls;*.xlsx','Excel (*.xls, *.xlsx)'; ...
   '*.*',  'All Files (*.*)'}, ...
   'Pick a file', ...
   'MultiSelect', 'off');
if numel(filename)>1
    userdata=get(handles.figure1, 'UserData');
    userdata.filename=filename;
    [num, txt, raw]=xlsread(fullfile(pathname, filename));
    userdata.files=txt;
    data={};
    for i=1:size(txt, 1),
        [path, file, ext]=fileparts(txt{i,1});
        data{i,1}=file;
        data{i,2}='';
        data{i,3}=true;
        data{i,4}=true;
        data{i,5}=true;
        data{i,6}=true;
        data{i,7}=true;
    end
    
    set(handles.uitable1, 'data', data);
    set(handles.figure1, 'UserData', userdata);
end
updateStats(handles);


function updateStats(handles)
userdata=get(handles.figure1, 'UserData');

if isfield(userdata, 'files'),
    for i=1:size(userdata.files, 1),
        [path, file, ext]=fileparts(userdata.files{i,1});
        pathNew=fullfile(path, file);
        configExist=false;
        if exist(pathNew, 'dir')
            paramFile=fullfile(pathNew, 'configuration parameters.mat');
            if exist(paramFile, 'file')
                configExist=true;
                
                config=load(paramFile);
                userdata.configs{i}=config.config;
                if numel(userdata.configs{i})==0
                    configExist=false;
                end
            end
        end
        if ~configExist
            setGrid(handles, i, 2, 'No Config');
        else
            normFile=fullfile(pathNew, 'Normalization', 'NormalizationDone.txt');
            normExist=exist(normFile, 'file');
            if ~normExist
                setGrid(handles, i, 2, 'Config OK');
            else
                
                trackFile=fullfile(pathNew, 'Tracking', 'trackBeforeLoadingFilliations.mat');
                trackExist=exist(trackFile, 'file');
                if ~trackExist
                    setGrid(handles, i, 2, '1 Done');
                else
                    filliationFile=fullfile(pathNew, 'Tracking', 'trackFilliation.mat');
                    filliationExist=exist(filliationFile, 'file');
                    if ~filliationExist
                        setGrid(handles, i, 2, '2 Done');
                    else
                        volumeFile=fullfile(pathNew, 'Tracking', 'volumeMeasurement.mat');
                        volumeExist=exist(volumeFile, 'file');
                        if ~volumeExist
                            setGrid(handles, i, 2, '3 Done');
                        else
                            csvFile=fullfile(pathNew, 'Tracking', 'results.csv');
                            csvExist=exist(csvFile, 'file');
                            if ~csvExist
                                setGrid(handles, i, 2, '4 Done');
                            else
                                setGrid(handles, i, 2, '5 Done');
                            end
                        end
                    end
                end
            end
        end
        set(handles.figure1, 'userdata', userdata);
    end
end


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
updateStats(handles)

% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1

function cellSelect(hObject, eventdat)
set(hObject,'UserData',eventdat.Indices);

function index=getSelectedRow(handles)
index=get(handles.uitable1, 'userdata');
if numel(index)>0
    index=index(:, 1);
end

function index=getSelectedCol(handles)
index=get(handles.uitable1, 'userdata');
if numel(index)>1
    index=index(:, 2);
end


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global currentParameters
global manualSeg
selectedImage=getSelectedRow(handles);
selectedCol=getSelectedCol(handles);
userdata=get(handles.figure1, 'UserData');
if numel(userdata)==0
    return;
end
if numel(selectedImage)==0
    return;
end

if selectedCol==3
    
    imageFilename=userdata.files{selectedImage, 1};
    if exist(imageFilename, 'file')
        % choose the selected file
        currentParameters.setRawMovieFilename(imageFilename);
    end
    
    %trackFilename=fullfile(currentParameters.trackingFolder, 'trackBeforeLoadingFilliations.mat');
    set(handles.figure1, 'pointer', 'watch') 
    drawnow;
    %if exist(trackFilename, 'file')
    %    load(trackFilename);
    %    manualSeg.doReset=false;
    %    pause(0.1);
    %end
    manualSeg.doReset=true;
    manualSegmentation();
    set(handles.figure1, 'pointer', 'arrow')
    return;
end

if selectedCol==4
    
    imageFilename=userdata.files{selectedImage, 1};
    if exist(imageFilename, 'file')
        % choose the selected file
        currentParameters.setRawMovieFilename(imageFilename);
    end
    
    trackFilename=fullfile(currentParameters.trackingFolder, 'trackBeforeLoadingFilliations.mat');
    set(handles.figure1, 'pointer', 'watch') 
    drawnow;
    if exist(trackFilename, 'file')
        loadCurrentParamAndManualSeg(trackFilename);
        manualSeg.doReset=false;
        pause(0.1);
    end
    manualSegmentation();
    set(handles.figure1, 'pointer', 'arrow')
    return;
end

if selectedCol==5
    
    imageFilename=userdata.files{selectedImage, 1};
    if exist(imageFilename, 'file')
        % choose the selected file
        currentParameters.setRawMovieFilename(imageFilename);
    end
    
    trackFilename=fullfile(currentParameters.trackingFolder, 'trackFilliation.mat');
    set(handles.figure1, 'pointer', 'watch') 
    drawnow;
    if exist(trackFilename, 'file')
        loadCurrentParamAndManualSeg(trackFilename);
        manualSeg.doReset=false;
        pause(0.1);
    end
    handlesMS=manualSegmentation();
    
    set(handles.figure1, 'pointer', 'arrow')
    return;
end


imageFilename=userdata.files{selectedImage, 1};
additionalFiles={};
for i=2:size(userdata.files, 2), 
    if numel(userdata.files{selectedImage, i})>0
        additionalFiles{i-1}=userdata.files{selectedImage, i};
    end
end
config=backgroundConfiguration(0, imageFilename);
if ~isfield(userdata, 'configs')
    userdata.configs={};
end

if numel(config)==0
    return;
end

userdata.configs{selectedImage}=config;
[path, file, ext]=fileparts(imageFilename);
pathNew=fullfile(path, file);
if ~exist(pathNew, 'dir')
    mkdir(pathNew);
end
save(fullfile(pathNew, 'configuration parameters.mat'), 'config');
set(handles.figure1, 'UserData', userdata);
updateStats(handles);

% --- Executes on button press in checkbox2.
function checkbox2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox2


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

selectedImage=getSelectedRow(handles);
selectedCol=getSelectedCol(handles);
userdata=get(handles.figure1, 'UserData');
if numel(userdata)==0
    return;
end
if numel(selectedImage)==0
    return;
end

%if selectedCol==3
defaultName='configuration parameters.mat';
[filename, pathname, filterindex] = uigetfile( ...
{  '*.mat','Matlab Variables (*.mat)'; ...
   '*.*',  'All Files (*.*)'}, ...
   'Select config file', defaultName,...
   'MultiSelect', 'off');
if numel(filename)>1
    filename=fullfile(pathname, filename);
    for i=1:numel(selectedImage)
        imageFilename=userdata.files{selectedImage(i), 1};
        [path, file, ext]=fileparts(imageFilename);
        matFilename=fullfile(path, file, defaultName);
        pathLoc=fullfile(path, file);
        if ~exist(pathLoc, 'dir'),
            mkdir(pathLoc);
        end
        if exist(matFilename, 'file')
            choice = questdlg(['The param file ', file, ' already exists. Do you want to overwrite ?'], 'Overwrite ?', 'Overwrite', 'Skip', 'Overwrite');
            if strcmpi(choice, 'Overwrite')
                copyfile(filename, matFilename);
            end
        else
            
            copyfile(filename, matFilename);
        end
    end
end


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data=get(handles.uitable1, 'data');
for i=1:size(data,1)
    for j=3:size(data,2),
        data{i, j}=false;
    end
end
set(handles.uitable1, 'data', data);


% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

userdata=get(handles.figure1, 'UserData');

path=uigetdir();
if numel(path)>1
    list=dir(fullfile(path, '*.tif'));
    
    
    listFiles={};
    for i=1:numel(list),
        listFiles{i, 1}=fullfile(path, list(i).name);
    end
    [selection, ok]=listdlg('ListString', listFiles);
    if ok~=1
        return;
    end
    
    listFiles=listFiles(selection);
    
    userdata.files=listFiles;
    data={};
    for i=1:size(listFiles, 1),
        [path, file, ext]=fileparts(listFiles{i});
        data{i,1}=file;
        data{i,2}='';
        data{i,3}=true;
        data{i,4}=true;
        data{i,5}=true;
        data{i,6}=true;
        data{i,7}=true;
    end
    
    set(handles.uitable1, 'data', data);
    set(handles.figure1, 'UserData', userdata);
end
updateStats(handles);


% --- Executes when entered data in editable cell(s) in uitable1.
function uitable1_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to uitable1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
