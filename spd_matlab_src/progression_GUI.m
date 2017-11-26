function varargout = progression_GUI(varargin)
% PROGRESSION_GUI M-file for progression_GUI.fig
%      PROGRESSION_GUI, by itself, creates a new PROGRESSION_GUI or raises the existing
%      singleton*.
%
%      H = PROGRESSION_GUI returns the handle to a new PROGRESSION_GUI or the handle to
%      the existing singleton*.
%
%      PROGRESSION_GUI('CALLBACK',hObject,eventData,handles,...) calls the
%      local
%      function named CALLBACK in PROGRESSION_GUI.M with the given input arguments.
%
%      PROGRESSION_GUI('Property','Value',...) creates a new PROGRESSION_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before progression_GUI_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to progression_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help progression_GUI

% Last Modified by GUIDE v2.5 08-Jun-2009 10:40:49

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @progression_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @progression_GUI_OutputFcn, ...
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


% --- Executes just before progression_GUI is made visible.
function progression_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to progression_GUI (see VARARGIN)

% Choose default command line output for progression_GUI
handles.output = hObject;
set(handles.Botton_Fit_module_mst,'Enable','off');
set(handles.Button_gene_module_adj,'Enable','off');
set(handles.Button_Add_handpicked_progression,'Enable','off');
set(handles.Botton_remove_selected_progression,'Enable','off');
set(handles.Botton_view_handpicked,'Enable','off');
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes progression_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = progression_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;




function Edit_raw_data_filename_Callback(hObject, eventdata, handles)
% hObject    handle to Edit_raw_data_filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% handles.data_filename = get(handles.Edit_raw_data_filename,'String');
% guidata(hObject,handles); 



% --- Executes during object creation, after setting all properties.
function Edit_raw_data_filename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Edit_raw_data_filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on key press over Button_Browse_raw_data_file with no controls selected.
function Button_Browse_raw_data_file_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to Button_Browse_raw_data_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in Button_Browse_raw_data_file.
function Button_Browse_raw_data_file_Callback(hObject, eventdata, handles)
% hObject    handle to Button_Browse_raw_data_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile('*.mat', 'Pick a .mat file');
if isequal(filename,0) || isequal(pathname,0)
    1;
else
    set(handles.Edit_raw_data_filename,'String',fullfile(pathname, filename))
end


% --- Executes on button press in Button_Load_raw_data.
function Button_Load_raw_data_Callback(hObject, eventdata, handles)
% hObject    handle to Button_Load_raw_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.data_filename = get(handles.Edit_raw_data_filename,'String');
if exist(handles.data_filename)~=2
    return
end
if ~isempty(handles.data_filename) && isequal(handles.data_filename(end-3:end),'.mat')
    load(handles.data_filename)
    handles.data = data; 
    handles.probe_names = probe_names;
    handles.exp_names = exp_names;
    if exist('color_code_names') && exist('color_code_vectors') && size(color_code_names,1)==size(color_code_vectors,1)
        handles.color_code_names = color_code_names;
        handles.color_code_vectors = color_code_vectors;
    else
        handles.color_code_names = [];
        handles.color_code_vectors = [];
    end
    handles.filter_gene_ind = 1:size(data,1);  % index of genes kept in the analysis
    handles.filtered_data = [];
    handleds.filtered_probe_names = [];
    handles.filter_std = 0;  
    handles.filter_acceptable_nulls = NaN; 
    handles.filter_throw_x_at = 0;
    handles.x_at_probes_ind=[]; for i=1:size(probe_names,2), handles.x_at_probes_ind = union(handles.x_at_probes_ind, isInList(probe_names(:,i),'_x_at')); end
    handles.noname_probes_ind=[]; for i=1:size(probe_names,1), for j=1:size(probe_names,2), if length(probe_names{i,j})==0, handles.noname_probes_ind=[handles.noname_probes_ind,i]; break; end; end;  end
%     handles.x_at_probes_ind = union(handles.x_at_probes_ind,handles.noname_probes_ind);
    handles.kmeans_iter =  str2num(get(handles.Edit_clustering_num_iter,'string'));
    handles.kmeans_clusters = 100;
    handles.min_module_size = str2num(get(handles.Edit_clustering_num_cluster,'string'));
    handles.module_coherence_thres = str2num(get(handles.Edit_clustering_coherence_threshold,'string'));
    handles.adj=[];
    handles.idx=[];
    handles.permutation_iter = str2num(get(handles.Edit_permutation_iter,'string'));
    handles.iter_merge_pvalue = 0.001;
    handles.pvalue_module_tree_thres = str2num(get(handles.Edit_pvalue_module_tree_thres,'string'));
    handles.candidate_trees = [];
    handles.p_value = [];
    handles.progression_modules=[];
    handles.handpicked_progression=[];
    handles.view_tree_handle = []; % this is useless though
    set(handles.Edit_std_threshold,'String',num2str(handles.filter_std));
    set(handles.Edit_num_of_acceptable_nulls,'String',num2str(handles.filter_acceptable_nulls ));
    set(handles.Checkbox_throw_x_at,'Value',handles.filter_throw_x_at);
    set(handles.Edit_num_genes,'String',num2str(size(data,1)));
    set(handles.Edit_num_samples,'String',num2str(size(data,2)));
    set(handles.Edit_num_null_entries,'String',num2str(sum(sum(isnan(data)))));
    set(handles.Edit_save_filename,'String',[handles.data_filename(1:end-4),'_result_',num2str(length(handles.filter_gene_ind)),'_',num2str(handles.filter_std),'_',num2str(handles.filter_acceptable_nulls),'_',num2str(handles.kmeans_iter),'_',num2str(handles.min_module_size),'_',num2str(handles.module_coherence_thres),'_',num2str(handles.permutation_iter),'_',num2str(handles.iter_merge_pvalue),'.mat']);
    set(handles.Edit_clustering_num_iter,'String',num2str(handles.kmeans_iter));
    set(handles.Edit_clustering_num_cluster,'String', num2str(handles.min_module_size)); % this box was designed for the number of clusters in kmeans, but this was later fixed to 2, and this box was used to specify another parameter, which is the minimun module size allowed
    set(handles.Edit_clustering_coherence_threshold,'String', num2str(handles.module_coherence_thres));
    set(handles.Edit_permutation_iter,'String',num2str(handles.permutation_iter));
    set(handles.Edit_iter_merge_pvalue,'String',num2str(handles.iter_merge_pvalue));
    set(handles.Listbox_progression_modules,'String',[]);
    set(handles.Edit_pvalue_module_tree_thres, 'String', num2str(handles.pvalue_module_tree_thres));
    set(handles.Edit_input_handpicked,'String',[]);
    set(handles.Listbox_handpicked_progression,'String',[]);
    set(handles.text20,'visible','off');
    set(handles.Button_show_module_quality,'Enable','off');
    set(handles.Button_view_modules_expr,'Enable','off');
    set(handles.Botton_Fit_module_mst,'Enable','off');
    set(handles.Button_gene_module_adj,'Enable','off');
    set(handles.Button_Add_handpicked_progression,'Enable','off');
    set(handles.Botton_remove_selected_progression,'Enable','off');
    set(handles.Botton_view_handpicked,'Enable','off');
    set(handles.text17,'visible','off'); set(handles.text18,'visible','off');
    guidata(hObject,handles); 
end


function Edit_Result_filename_Callback(hObject, eventdata, handles)
% hObject    handle to Edit_Result_filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Edit_Result_filename as text
%        str2double(get(hObject,'String')) returns contents of Edit_Result_filename as a double
handles.result_filename = get(handles.Edit_Result_filename,'String');
guidata(hObject,handles); 


% --- Executes during object creation, after setting all properties.
function Edit_Result_filename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Edit_Result_filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in Button_browse_result_file.
function Button_browse_result_file_Callback(hObject, eventdata, handles)
% hObject    handle to Button_browse_result_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename, pathname] = uigetfile('*.mat', 'Pick a .mat file');
if isequal(filename,0) || isequal(pathname,0)
    1;
else
    set(handles.Edit_Result_filename,'String',fullfile(pathname, filename))
end


% --- Executes on button press in Button_Load_result_file.
function Button_Load_result_file_Callback(hObject, eventdata, handles)
% hObject    handle to Button_Load_result_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.result_filename = get(handles.Edit_Result_filename,'String');
if exist(handles.result_filename)~=2
    return
end
if ~isempty(handles.result_filename) && isequal(handles.result_filename(end-3:end),'.mat')
    load(handles.result_filename)
    if ~isequal(handles.result_filename, result_filename)  % this means that the whole folder of data/result may have been moved to another directory or computer / or something wrong happened
        if ~isequal(handles.result_filename(max([1,find(handles.result_filename=='\')+1]):end),result_filename(max([1,find(result_filename=='\')+1]):end)) % if result filename part is not equal, something is wrong
            error('something wrong happened'); % this should never happen
        end
        if exist([handles.result_filename(1:max(find(handles.result_filename=='\'))),data_filename(max([1,find(data_filename=='\')+1]):end)])~=2   % if at the same directory of the user input, we can not find the data file, something wrong
            error('something wrong happened'); % this should never happen
        end
        % if the above two tests passed, this means that the whole thing is moved to another path
        data_filename = [handles.result_filename(1:max(find(handles.result_filename=='\'))),data_filename(max(find(data_filename=='\'))+1:end)];
        result_filename = handles.result_filename;
        save(result_filename,'data_filename','result_filename','-append');
    end
    handles.data_filename = data_filename;
    handles.data = data;                        clear data
    handles.probe_names = probe_names;          clear probe_names
    handles.exp_names = exp_names;
    handles.color_code_names = color_code_names;
    handles.color_code_vectors = color_code_vectors;
    handles.filter_gene_ind = filter_gene_ind;
    handles.filtered_data = filtered_data;
    handles.filtered_probe_names = handles.probe_names(filter_gene_ind,:);
    handles.filter_std = filter_std;
    handles.filter_acceptable_nulls = filter_acceptable_nulls;
    handles.filter_throw_x_at = filter_throw_x_at;
    handles.x_at_probes_ind = x_at_probes_ind;
    handles.kmeans_iter = kmeans_iter;
    if ~exist('min_module_size')
        [N,X] = hist(idx,1:max(idx));
        handles.min_module_size = min(N);
    else
        handles.min_module_size = min_module_size;
    end
    handles.kmeans_clusters = kmeans_clusters;
    handles.module_coherence_thres = module_coherence_thres;
    handles.adj = adj;                          clear adj
    handles.idx = idx;                          clear idx
    handles.permutation_iter = permutation_iter;
    handles.iter_merge_pvalue = iter_merge_pvalue;
    handles.pvalue_module_tree_thres = pvalue_module_tree_thres;
    handles.candidate_trees = candidate_trees;
    handles.p_value = p_value;          
    handles.progression_modules = progression_modules;
    handles.handpicked_progression = handpicked_progression;
    set(handles.Edit_std_threshold,'String',num2str(handles.filter_std));
    set(handles.Edit_num_of_acceptable_nulls,'String',num2str(handles.filter_acceptable_nulls ));
    set(handles.Checkbox_throw_x_at,'Value',handles.filter_throw_x_at);
    set(handles.Edit_num_genes,'String',num2str(size(handles.filtered_data,1)));
    set(handles.Edit_num_samples,'String',num2str(size(handles.filtered_data,2)));
    set(handles.Edit_num_null_entries,'String',num2str(sum(sum(isnan(handles.data(handles.filter_gene_ind,:))))));
    set(handles.Edit_save_filename,'String',handles.result_filename);
%     if ~isequal(handles.result_filename, [handles.data_filename(1:end-4),'_result_',num2str(length(handles.filter_gene_ind)),'_',num2str(handles.filter_std),'_',num2str(handles.filter_acceptable_nulls),'_',num2str(handles.kmeans_iter),'_',num2str(handles.kmeans_clusters),'_',num2str(handles.module_coherence_thres),'_',num2str(handles.permutation_iter),'_',num2str(handles.iter_merge_pvalue),'.mat'])
%         error('something wrong happened') % this should never happen
%     end
    set(handles.Edit_clustering_num_iter,'String',num2str(handles.kmeans_iter));
    set(handles.Edit_clustering_num_cluster,'String', num2str(handles.min_module_size));
    set(handles.Edit_clustering_coherence_threshold,'String', num2str(handles.module_coherence_thres));
    set(handles.Edit_permutation_iter,'String',num2str(handles.permutation_iter));
    set(handles.Edit_iter_merge_pvalue,'String',num2str(handles.iter_merge_pvalue));
    set(handles.Edit_pvalue_module_tree_thres, 'String', num2str(handles.pvalue_module_tree_thres));
    set(handles.Edit_input_handpicked,'String',[]);
    set(handles.text17,'visible','off'); set(handles.text18,'visible','off');
    set(handles.Button_show_module_quality,'Enable','on');
    set(handles.Button_view_modules_expr,'Enable','on');
    set(handles.Botton_Fit_module_mst,'Enable','on');
    set(handles.Button_gene_module_adj,'Enable','on');
    set(handles.Button_Add_handpicked_progression,'Enable','on');
    set(handles.Botton_remove_selected_progression,'Enable','on');
    set(handles.Botton_view_handpicked,'Enable','on');

% the following 3 needs special care
    set(handles.text20,'visible','on'); set(handles.text20,'string',['Total number of modules is ',num2str(max(handles.idx))]);
    tmp_for_listbox=[];
    for i=1:length(handles.progression_modules)
        tmp_for_listbox{i,1} = [num2str(size(handles.progression_modules(i).gene_expr,1)),' genes, gene modules:',num2str(handles.progression_modules(i).gene_modules(:)')];
    end
    set(handles.Listbox_progression_modules,'string',tmp_for_listbox); set(handles.Listbox_progression_modules,'value',1);
    tmp_for_listbox=[];
    for i=1:length(handles.handpicked_progression)
        tmp_for_listbox{i,1} = [num2str(size(handles.handpicked_progression(i).gene_expr,1)),' genes, gene modules:',num2str(handles.handpicked_progression(i).gene_modules(:)')];
    end
    set(handles.Listbox_handpicked_progression,'string',tmp_for_listbox); set(handles.Listbox_handpicked_progression,'value',1);    
    guidata(hObject,handles); 
end



% --- Executes on button press in radiobutton_load_raw_data.
function radiobutton_load_raw_data_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_load_raw_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(hObject,'Value')==0
    set(hObject,'Value',1);
end
% set the other option off
set(handles.Edit_Result_filename,'String',' '); set(handles.Edit_Result_filename,'Enable','off'); set(handles.Edit_Result_filename,'BackgroundColor',[0.9,0.9,0.9]);
set(handles.Button_browse_result_file,'Enable','off');
set(handles.Button_Load_result_file,'Enable','off');
% set self on
set(handles.Edit_raw_data_filename,'Enable','on'); set(handles.Edit_raw_data_filename,'BackgroundColor',[1,1,1]);
set(handles.Button_Browse_raw_data_file,'Enable','on');
set(handles.Button_Load_raw_data,'Enable','on');



% --- Executes on button press in radiobutton_load_result.
function radiobutton_load_result_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_load_result (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(hObject,'Value')==0
    set(hObject,'Value',1);
end
% set the other option off
set(handles.Edit_raw_data_filename,'String',' '); set(handles.Edit_raw_data_filename,'Enable','off'); set(handles.Edit_raw_data_filename,'BackgroundColor',[0.9,0.9,0.9]);
set(handles.Button_Browse_raw_data_file,'Enable','off');
set(handles.Button_Load_raw_data,'Enable','off');
% set self on
set(handles.Edit_Result_filename,'Enable','on'); set(handles.Edit_Result_filename,'BackgroundColor',[1,1,1]);
set(handles.Button_browse_result_file,'Enable','on');
set(handles.Button_Load_result_file,'Enable','on');



%%

% --- Executes on button press in Checkbox_throw_x_at.
function Checkbox_throw_x_at_Callback(hObject, eventdata, handles)
% hObject    handle to Checkbox_throw_x_at (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isnan(handles.filter_acceptable_nulls)
    num_acceptable_nulls = size(handles.data,2)+1;
else
    num_acceptable_nulls = handles.filter_acceptable_nulls;
end
handles.filter_gene_ind = find(nanstd(handles.data')>=handles.filter_std & sum(isnan(handles.data)',1)<=num_acceptable_nulls);
if get(handles.Checkbox_throw_x_at,'value')==1
    handles.filter_gene_ind = setdiff(handles.filter_gene_ind,handles.x_at_probes_ind);
end
% set its direct effect
set(handles.Edit_num_genes,'String',num2str(length(handles.filter_gene_ind)));
set(handles.Edit_num_null_entries,'String',num2str(sum(sum(isnan(handles.data(handles.filter_gene_ind,:))))));
set(handles.Edit_save_filename,'String',[handles.data_filename(1:end-4),'_result_',num2str(length(handles.filter_gene_ind)),'_',num2str(handles.filter_std),'_',num2str(handles.filter_acceptable_nulls),'_',num2str(handles.kmeans_iter),'_',num2str(handles.min_module_size),'_',num2str(handles.module_coherence_thres),'_',num2str(handles.permutation_iter),'_',num2str(handles.iter_merge_pvalue),'.mat']);
% set subsequent analyze empty
    handles.adj=[];
    handles.idx=[];
    handles.candidate_trees = [];
    handles.p_value = [];
    handles.progression_modules=[];
    handles.handpicked_progression=[];
    set(handles.text20,'visible','off');
    set(handles.Button_show_module_quality,'Enable','off');
    set(handles.Button_view_modules_expr,'Enable','off');
    set(handles.Botton_Fit_module_mst,'Enable','off');
    set(handles.Button_gene_module_adj,'Enable','off');
    set(handles.Button_Add_handpicked_progression,'Enable','off');
    set(handles.Botton_remove_selected_progression,'Enable','off');
    set(handles.Botton_view_handpicked,'Enable','off');
    set(handles.text17,'visible','off'); set(handles.text18,'visible','off');
    set(handles.Listbox_progression_modules,'String',[]);
    set(handles.Edit_input_handpicked,'String',[]);
    set(handles.Listbox_handpicked_progression,'String',[]);
guidata(hObject,handles); 



function Edit_num_of_acceptable_nulls_Callback(hObject, eventdata, handles)
% hObject    handle to Edit_num_of_acceptable_nulls (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tmp = str2num(get(hObject,'String'));
if length(tmp)~=1 || tmp~=round(tmp)  % this thing has to be integer
    set(hObject,'String',num2str(handles.filter_acceptable_nulls)); return
else
    handles.filter_acceptable_nulls = tmp;
end
% filter the genes
if isnan(handles.filter_acceptable_nulls)
    num_acceptable_nulls = size(handles.data,2)+1;
else
    num_acceptable_nulls = handles.filter_acceptable_nulls;
end
handles.filter_gene_ind = find(nanstd(handles.data')>=handles.filter_std & sum(isnan(handles.data)',1)<=num_acceptable_nulls);
if get(handles.Checkbox_throw_x_at,'value')==1
    handles.filter_gene_ind = setdiff(handles.filter_gene_ind,handles.x_at_probes_ind);
end
set(handles.Edit_num_genes,'String',num2str(length(handles.filter_gene_ind)));
set(handles.Edit_num_null_entries,'String',num2str(sum(sum(isnan(handles.data(handles.filter_gene_ind,:))))));
set(handles.Edit_save_filename,'String',[handles.data_filename(1:end-4),'_result_',num2str(length(handles.filter_gene_ind)),'_',num2str(handles.filter_std),'_',num2str(handles.filter_acceptable_nulls),'_',num2str(handles.kmeans_iter),'_',num2str(handles.min_module_size),'_',num2str(handles.module_coherence_thres),'_',num2str(handles.permutation_iter),'_',num2str(handles.iter_merge_pvalue),'.mat']);
% set subsequent analyze empty
    handles.adj=[];
    handles.idx=[];
    handles.candidate_trees = [];
    handles.p_value = [];
    handles.progression_modules=[];
    handles.handpicked_progression=[];
    set(handles.text20,'visible','off');
    set(handles.Button_show_module_quality,'Enable','off');
    set(handles.Button_view_modules_expr,'Enable','off');
    set(handles.Botton_Fit_module_mst,'Enable','off');
    set(handles.Button_gene_module_adj,'Enable','off');
    set(handles.Button_Add_handpicked_progression,'Enable','off');
    set(handles.Botton_remove_selected_progression,'Enable','off');
    set(handles.Botton_view_handpicked,'Enable','off');
    set(handles.text17,'visible','off'); set(handles.text18,'visible','off');
    set(handles.Listbox_progression_modules,'String',[]);
    set(handles.Edit_input_handpicked,'String',[]);
    set(handles.Listbox_handpicked_progression,'String',[]);
guidata(hObject,handles); 





% --- Executes during object creation, after setting all properties.
function Edit_num_of_acceptable_nulls_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Edit_num_of_acceptable_nulls (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Edit_std_threshold_Callback(hObject, eventdata, handles)
% hObject    handle to Edit_std_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tmp = str2num(get(hObject,'String'));
if length(tmp)~=1
    set(hObject,'String',num2str(handles.filter_std)); return
else
    handles.filter_std = tmp;
end
% filter the genes
if isnan(handles.filter_acceptable_nulls)
    num_acceptable_nulls = size(handles.data,2)+1;
else
    num_acceptable_nulls = handles.filter_acceptable_nulls;
end
handles.filter_gene_ind = find(nanstd(handles.data')>=handles.filter_std & sum(isnan(handles.data)',1)<=num_acceptable_nulls);
if get(handles.Checkbox_throw_x_at,'value')==1
    handles.filter_gene_ind = setdiff(handles.filter_gene_ind,handles.x_at_probes_ind);
end
set(handles.Edit_num_genes,'String',num2str(length(handles.filter_gene_ind)));
set(handles.Edit_num_null_entries,'String',num2str(sum(sum(isnan(handles.data(handles.filter_gene_ind,:))))));
set(handles.Edit_save_filename,'String',[handles.data_filename(1:end-4),'_result_',num2str(length(handles.filter_gene_ind)),'_',num2str(handles.filter_std),'_',num2str(handles.filter_acceptable_nulls),'_',num2str(handles.kmeans_iter),'_',num2str(handles.min_module_size),'_',num2str(handles.module_coherence_thres),'_',num2str(handles.permutation_iter),'_',num2str(handles.iter_merge_pvalue),'.mat']);
% set subsequent analyze empty
    handles.adj=[];
    handles.idx=[];
    handles.candidate_trees = [];
    handles.p_value = [];
    handles.progression_modules=[];
    handles.handpicked_progression=[];
    set(handles.text20,'visible','off');
    set(handles.Button_show_module_quality,'Enable','off');
    set(handles.Button_view_modules_expr,'Enable','off');
    set(handles.Botton_Fit_module_mst,'Enable','off');
    set(handles.Button_gene_module_adj,'Enable','off');
    set(handles.Button_Add_handpicked_progression,'Enable','off');
    set(handles.Botton_remove_selected_progression,'Enable','off');
    set(handles.Botton_view_handpicked,'Enable','off');
    set(handles.text17,'visible','off'); set(handles.text18,'visible','off');
    set(handles.Listbox_progression_modules,'String',[]);
    set(handles.Edit_input_handpicked,'String',[]);
    set(handles.Listbox_handpicked_progression,'String',[]);
guidata(hObject,handles); 




% --- Executes during object creation, after setting all properties.
function Edit_std_threshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Edit_std_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





% --- Executes during object creation, after setting all properties.
function Edit_num_genes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Edit_num_genes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function Edit_num_samples_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Edit_num_samples (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function Edit_num_null_entries_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Edit_num_null_entries (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Edit_save_filename_Callback(hObject, eventdata, handles)
% hObject    handle to Edit_save_filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Edit_save_filename as text
%        str2double(get(hObject,'String')) returns contents of Edit_save_filename as a double


% --- Executes during object creation, after setting all properties.
function Edit_save_filename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Edit_save_filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end






function Edit_clustering_num_iter_Callback(hObject, eventdata, handles)
% hObject    handle to Edit_clustering_num_iter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tmp = str2num(get(hObject,'String'));
if length(tmp)~=1 || tmp~=round(tmp)  || tmp<=0 % this thing has to be integer
    set(hObject,'String',num2str(handles.kmeans_iter));
else
    handles.kmeans_iter = tmp;
    set(handles.Edit_save_filename,'String',[handles.data_filename(1:end-4),'_result_',num2str(length(handles.filter_gene_ind)),'_',num2str(handles.filter_std),'_',num2str(handles.filter_acceptable_nulls),'_',num2str(handles.kmeans_iter),'_',num2str(handles.min_module_size),'_',num2str(handles.module_coherence_thres),'_',num2str(handles.permutation_iter),'_',num2str(handles.iter_merge_pvalue),'.mat']);
    % set subsequent analyze empty
    handles.adj=[];
    handles.idx=[];
    handles.candidate_trees = [];
    handles.p_value = [];
    handles.progression_modules=[];
    handles.handpicked_progression=[];
    set(handles.text20,'visible','off');
    set(handles.Button_show_module_quality,'Enable','off');
    set(handles.Button_view_modules_expr,'Enable','off');
    set(handles.Botton_Fit_module_mst,'Enable','off');
    set(handles.Button_gene_module_adj,'Enable','off');
    set(handles.Button_Add_handpicked_progression,'Enable','off');
    set(handles.Botton_remove_selected_progression,'Enable','off');
    set(handles.Botton_view_handpicked,'Enable','off');
    set(handles.text17,'visible','off'); set(handles.text18,'visible','off');
    set(handles.Listbox_progression_modules,'String',[]);
    set(handles.Edit_input_handpicked,'String',[]);
    set(handles.Listbox_handpicked_progression,'String',[]);
    guidata(hObject,handles); 
end


% --- Executes during object creation, after setting all properties.
function Edit_clustering_num_iter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Edit_clustering_num_iter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Edit_clustering_num_cluster_Callback(hObject, eventdata, handles)
% hObject    handle to Edit_clustering_num_cluster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tmp = str2num(get(hObject,'String'));
if length(tmp)~=1 || tmp~=round(tmp)  || tmp<=0 || tmp>=round(size(handles.data,1)/1) % this thing has to be integer
    set(hObject,'String',num2str(handles.min_module_size));
else
    handles.min_module_size = tmp;
    set(handles.Edit_save_filename,'String',[handles.data_filename(1:end-4),'_result_',num2str(length(handles.filter_gene_ind)),'_',num2str(handles.filter_std),'_',num2str(handles.filter_acceptable_nulls),'_',num2str(handles.kmeans_iter),'_',num2str(handles.min_module_size),'_',num2str(handles.module_coherence_thres),'_',num2str(handles.permutation_iter),'_',num2str(handles.iter_merge_pvalue),'.mat']);
    % set subsequent analyze empty
    handles.adj=[];
    handles.idx=[];
    handles.candidate_trees = [];
    handles.p_value = [];
    handles.progression_modules=[];
    handles.handpicked_progression=[];
    set(handles.text20,'visible','off');
    set(handles.Button_show_module_quality,'Enable','off');
    set(handles.Button_view_modules_expr,'Enable','off');
    set(handles.Botton_Fit_module_mst,'Enable','off');
    set(handles.Button_gene_module_adj,'Enable','off');
    set(handles.Button_Add_handpicked_progression,'Enable','off');
    set(handles.Botton_remove_selected_progression,'Enable','off');
    set(handles.Botton_view_handpicked,'Enable','off');
    set(handles.text17,'visible','off'); set(handles.text18,'visible','off');
    set(handles.Listbox_progression_modules,'String',[]);
    set(handles.Edit_input_handpicked,'String',[]);
    set(handles.Listbox_handpicked_progression,'String',[]);
    guidata(hObject,handles); 
end




% --- Executes during object creation, after setting all properties.
function Edit_clustering_num_cluster_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Edit_clustering_num_cluster (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Edit_clustering_coherence_threshold_Callback(hObject, eventdata, handles)
% hObject    handle to Edit_clustering_coherence_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tmp = str2num(get(hObject,'String'));
if length(tmp)~=1 || tmp<0 || tmp>1 
    set(hObject,'String',num2str(handles.module_coherence_thres));
else
    handles.module_coherence_thres = tmp;
    set(handles.Edit_save_filename,'String',[handles.data_filename(1:end-4),'_result_',num2str(length(handles.filter_gene_ind)),'_',num2str(handles.filter_std),'_',num2str(handles.filter_acceptable_nulls),'_',num2str(handles.kmeans_iter),'_',num2str(handles.min_module_size),'_',num2str(handles.module_coherence_thres),'_',num2str(handles.permutation_iter),'_',num2str(handles.iter_merge_pvalue),'.mat']);
    % set subsequent analyze empty
    handles.adj=[];
    handles.idx=[];
    handles.candidate_trees = [];
    handles.p_value = [];
    handles.progression_modules=[];
    handles.handpicked_progression=[];
    set(handles.text20,'visible','off');
    set(handles.Button_show_module_quality,'Enable','off');
    set(handles.Button_view_modules_expr,'Enable','off');
    set(handles.Botton_Fit_module_mst,'Enable','off');
    set(handles.Button_gene_module_adj,'Enable','off');
    set(handles.Button_Add_handpicked_progression,'Enable','off');
    set(handles.Botton_remove_selected_progression,'Enable','off');
    set(handles.Botton_view_handpicked,'Enable','off');
    set(handles.text17,'visible','off'); set(handles.text18,'visible','off');
    set(handles.Listbox_progression_modules,'String',[]);
    set(handles.Edit_input_handpicked,'String',[]);
    set(handles.Listbox_handpicked_progression,'String',[]);
    guidata(hObject,handles); 
end


% --- Executes during object creation, after setting all properties.
function Edit_clustering_coherence_threshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Edit_clustering_coherence_threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Button_perform_gene_clustering.
function Button_perform_gene_clustering_Callback(hObject, eventdata, handles)
if isfield(handles,'data')==0 || isempty(handles.data)
    display('Please load data first, before performing clustering !!!')
    return
end
set(handles.text18,'visible','on');
set(handles.text20,'visible','off');
set(handles.Button_show_module_quality,'Enable','off');
set(handles.Button_view_modules_expr,'Enable','off');
set(handles.Botton_Fit_module_mst,'Enable','off');
set(handles.Button_gene_module_adj,'Enable','off');
set(handles.Button_Add_handpicked_progression,'Enable','off');
set(handles.Botton_remove_selected_progression,'Enable','off');
set(handles.Botton_view_handpicked,'Enable','off');
handles.filtered_data = handles.data(handles.filter_gene_ind,:);
if sum(sum(isnan(handles.filtered_data)))~=0
    handles.filtered_data = knnimpute(handles.filtered_data);
end
% handles.filtered_data = per_gene_normalization(handles.filtered_data);
handles.filtered_probe_names = handles.probe_names(handles.filter_gene_ind,:);
handles.idx = agglomerative_clustering(handles.filtered_data, handles.module_coherence_thres, 0.9);
for i=1:max(handles.idx)
    ind = find(handles.idx==i); 
    if length(ind)<handles.min_module_size
        handles.idx(ind)=0;
    end
end
handles.idx = standardize_idx(handles.idx);
set(handles.Button_show_module_quality,'Enable','on');
set(handles.Button_view_modules_expr,'Enable','on');
set(handles.Botton_Fit_module_mst,'Enable','on');
set(handles.text18,'visible','off');
set(handles.text20,'visible','on');
set(handles.text20,'string',['Total number of modules is ',num2str(max(handles.idx))]);
guidata(hObject,handles); 


% --- Executes on button press in Button_perform_gene_clustering_method2.
function Button_perform_gene_clustering_method2_Callback(hObject, eventdata, handles)
% This botton is going to do something slightly different from above
% we make the consensus kmeans and the partitioning process more glued
% We use consensus 2-means to cluster the data into 2 clusters N times
% construct an adj matrix, and partition that into 2 pieces,
% check each resulting piece, see whether they met the stopping criterion
if isfield(handles,'data')==0 || isempty(handles.data)
    display('Please load data first, before performing clustering !!!')
    return
end
set(handles.text18,'visible','on');
set(handles.text20,'visible','off');
set(handles.Button_show_module_quality,'Enable','off');
set(handles.Button_view_modules_expr,'Enable','off');
set(handles.Botton_Fit_module_mst,'Enable','off');
set(handles.Button_gene_module_adj,'Enable','off');
set(handles.Button_Add_handpicked_progression,'Enable','off');
set(handles.Botton_remove_selected_progression,'Enable','off');
set(handles.Botton_view_handpicked,'Enable','off');
handles.filtered_data = handles.data(handles.filter_gene_ind,:);
if sum(sum(isnan(handles.filtered_data)))~=0
    handles.filtered_data = knnimpute(handles.filtered_data);
end
% handles.filtered_data = per_gene_normalization(handles.filtered_data);
handles.filtered_probe_names = handles.probe_names(handles.filter_gene_ind,:);
handles.idx = iterative_consensus_kmeans_graph_partition(handles.filtered_data, handles.kmeans_iter, handles.module_coherence_thres, 0.9);
for i=1:max(handles.idx)
    ind = find(handles.idx==i); 
    if length(ind)<=handles.min_module_size
        handles.idx(ind)=0;
    end
end
handles.idx = standardize_idx(handles.idx);
set(handles.Button_show_module_quality,'Enable','on');
set(handles.Button_view_modules_expr,'Enable','on');
set(handles.Botton_Fit_module_mst,'Enable','on');
set(handles.text18,'visible','off');
set(handles.text20,'visible','on');
set(handles.text20,'string',['Total number of modules is ',num2str(max(handles.idx))]);
guidata(hObject,handles); 


 

% --- Executes on button press in Button_show_module_quality.
function Button_show_module_quality_Callback(hObject, eventdata, handles)
if isempty(handles.idx), return; end
ss = zeros(2,max(handles.idx));
for i=1:max(handles.idx), 
    ss(1,i) = avg_center_gene_corr(handles.filtered_data(handles.idx==i,:)); 
    ss(2,i) = sum(handles.idx==i);
end; 
figure; subplot(2,1,1);stem(ss(1,:));title('average center gene correlation'); axis tight
subplot(2,1,2);bar(ss(2,:));title('module size');axis tight






function Edit_permutation_iter_Callback(hObject, eventdata, handles)
% hObject    handle to Edit_permutation_iter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tmp = str2num(get(hObject,'String'));
if length(tmp)~=1 || tmp<0 || tmp~=round(tmp)
    set(hObject,'String',num2str(handles.permutation_iter)); return
else
    handles.permutation_iter = tmp;
end
set(handles.Edit_save_filename,'String',[handles.data_filename(1:end-4),'_result_',num2str(length(handles.filter_gene_ind)),'_',num2str(handles.filter_std),'_',num2str(handles.filter_acceptable_nulls),'_',num2str(handles.kmeans_iter),'_',num2str(handles.min_module_size),'_',num2str(handles.module_coherence_thres),'_',num2str(handles.permutation_iter),'_',num2str(handles.iter_merge_pvalue),'.mat']);
% set subsequent analyze empty
    handles.candidate_trees = [];
    handles.p_value = [];
    handles.progression_modules=[];
    handles.handpicked_progression=[];
    set(handles.text17,'visible','off'); 
    set(handles.Listbox_progression_modules,'String',[]);
    set(handles.Edit_input_handpicked,'String',[]);
    set(handles.Listbox_handpicked_progression,'String',[]);

    set(handles.Button_gene_module_adj,'Enable','off');
    set(handles.Button_Add_handpicked_progression,'Enable','off');
    set(handles.Botton_remove_selected_progression,'Enable','off');
    set(handles.Botton_view_handpicked,'Enable','off');
guidata(hObject,handles); 


% --- Executes during object creation, after setting all properties.
function Edit_permutation_iter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Edit_permutation_iter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Edit_iter_merge_pvalue_Callback(hObject, eventdata, handles)
% hObject    handle to Edit_iter_merge_pvalue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tmp = str2num(get(hObject,'String'));
if length(tmp)~=1 || tmp<0 || tmp>1
    set(hObject,'String',num2str(handles.iter_merge_pvalue)); return
else
    handles.iter_merge_pvalue = tmp;
end
set(handles.Edit_save_filename,'String',[handles.data_filename(1:end-4),'_result_',num2str(length(handles.filter_gene_ind)),'_',num2str(handles.filter_std),'_',num2str(handles.filter_acceptable_nulls),'_',num2str(handles.kmeans_iter),'_',num2str(handles.min_module_size),'_',num2str(handles.module_coherence_thres),'_',num2str(handles.permutation_iter),'_',num2str(handles.iter_merge_pvalue),'.mat']);
% set subsequent analyze empty
    handles.candidate_trees = [];
    handles.p_value = [];
    handles.progression_modules=[];
    handles.handpicked_progression=[];
    set(handles.text17,'visible','off'); 
    set(handles.Listbox_progression_modules,'String',[]);
    set(handles.Edit_input_handpicked,'String',[]);
    set(handles.Listbox_handpicked_progression,'String',[]);
guidata(hObject,handles); 



% --- Executes during object creation, after setting all properties.
function Edit_iter_merge_pvalue_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Edit_iter_merge_pvalue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in Botton_Fit_module_mst.
function Botton_Fit_module_mst_Callback(hObject, eventdata, handles)
% hObject    handle to Botton_Fit_module_mst (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% GO_button
if isempty(handles.idx), return; end
set(handles.text17,'visible','on');
% % % task 1, generate progression modules
handles.progression_modules = [];% generate_progression_modules(handles.filtered_data,handles.idx, handles.iter_merge_pvalue, handles.permutation_iter);
% % % task 2, compute pairwise modules - candidate_trees fit
[handles.candidate_trees, handles.p_value] = build_candidate_trees_fit_modules_v3(handles.filtered_data,handles.idx, handles.permutation_iter);
% handles.candidate_trees=zeros(1,1,max(handles.idx));handles.p_value=zeros(max(handles.idx),max(handles.idx));
% % tmp_for_listbox=[];
% % for i=1:length(handles.progression_modules)
% %     tmp_for_listbox{i,1} = [num2str(size(handles.progression_modules(i).gene_expr,1)),' genes, gene modules:',num2str(handles.progression_modules(i).gene_modules(:)')];
% % end
% % set(handles.Listbox_progression_modules,'string',tmp_for_listbox); 
set(handles.Listbox_progression_modules,'value',1);
set(handles.text17,'visible','off');
set(handles.Button_gene_module_adj,'Enable','on');
set(handles.Button_Add_handpicked_progression,'Enable','on');
set(handles.Botton_remove_selected_progression,'Enable','on');
set(handles.Botton_view_handpicked,'Enable','on');
guidata(hObject,handles); 




% --- Executes on selection change in Listbox_progression_modules.
function Listbox_progression_modules_Callback(hObject, eventdata, handles)
% hObject    handle to Listbox_progression_modules (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns Listbox_progression_modules contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Listbox_progression_modules


% --- Executes during object creation, after setting all properties.
function Listbox_progression_modules_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Listbox_progression_modules (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in Botton_view_progression_module.
function Botton_view_progression_module_Callback(hObject, eventdata, handles)
% hObject    handle to Botton_view_progression_module (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ind = get(handles.Listbox_progression_modules,'value');
if ind<1 || ind>length(handles.progression_modules)
    ind=1;
    set(handles.Listbox_progression_modules,'value',1); return
end
handles.view_tree_handle = view_tree(handles.progression_modules(ind), handles.exp_names, handles.color_code_names, handles.color_code_vectors, handles);
guidata(hObject,handles)

 

% --- Executes on selection change in Listbox_handpicked_progression.
function Listbox_handpicked_progression_Callback(hObject, eventdata, handles)
% hObject    handle to Listbox_handpicked_progression (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns Listbox_handpicked_progression contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Listbox_handpicked_progression


% --- Executes during object creation, after setting all properties.
function Listbox_handpicked_progression_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Listbox_handpicked_progression (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Edit_input_handpicked_Callback(hObject, eventdata, handles)
% hObject    handle to Edit_input_handpicked (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tmp = str2num(get(hObject,'string'));
if length(tmp)==0
    set(hObject,'string',[]);
end
if min(tmp)<=0 || max(tmp)>max(handles.idx)
    set(hObject,'string',[]);
end
if sum(round(tmp)~=tmp)~=0
    set(hObject,'string',[]);
end



% --- Executes during object creation, after setting all properties.
function Edit_input_handpicked_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Edit_input_handpicked (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Button_Add_handpicked_progression.
function Button_Add_handpicked_progression_Callback(hObject, eventdata, handles)
% hObject    handle to Button_Add_handpicked_progression (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tmp = str2num(get(handles.Edit_input_handpicked,'string'));
if length(tmp)==0 || min(tmp)<=0 || max(tmp)>max(handles.idx) || sum(round(tmp)~=tmp)~=0
    set(handles.Edit_input_handpicked,'string',[]);
else
    if length(tmp)~=length(unique(tmp)), tmp = unique(tmp); end
%     gene_module_mean = get_module_mean(handles.filtered_data,handles.idx);   
    gene_module_mean = get_module_mean(per_gene_normalization(handles.filtered_data),handles.idx);   
    picked_progression_module.gene_modules = tmp; 
    picked_progression_module.gene_modules_mean = gene_module_mean(tmp,:);
    picked_progression_module.gene_expr = [];
    for i=1:length(tmp)
%         picked_progression_module.gene_expr = [picked_progression_module.gene_expr; handles.filtered_data(handles.idx==tmp(i),:)]; 
        picked_progression_module.gene_expr = [picked_progression_module.gene_expr; per_gene_normalization(handles.filtered_data(handles.idx==tmp(i),:))]; 
    end
    picked_progression_module_size = size(picked_progression_module.gene_expr,1);
    [adj1,adj2,cost] = mst(picked_progression_module.gene_expr');
%     [adj1,adj2,cost] = mst(picked_progression_module.gene_expr','corr');
    picked_progression_module.adj1 = adj1;
    picked_progression_module.adj2 = adj2;
    picked_progression_module.node_position = [];
    handles.handpicked_progression=[handles.handpicked_progression,picked_progression_module];

    tmp_for_listbox = get(handles.Listbox_handpicked_progression,'string');
    tmp_for_listbox{end+1,1} = [num2str(size(handles.handpicked_progression(end).gene_expr,1)),' genes, gene modules:',num2str(handles.handpicked_progression(end).gene_modules(:)')];
    set(handles.Listbox_handpicked_progression,'string',tmp_for_listbox)
    set(handles.Listbox_handpicked_progression,'value',length(handles.handpicked_progression));    
end
guidata(hObject, handles);


% --- Executes on button press in Botton_remove_selected_progression.
function Botton_remove_selected_progression_Callback(hObject, eventdata, handles)
% hObject    handle to Botton_remove_selected_progression (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ind = get(handles.Listbox_handpicked_progression,'value');
if length(handles.handpicked_progression)>=ind
    handles.handpicked_progression(ind)=[];
    tmp = get(handles.Listbox_handpicked_progression,'string');
    tmp(ind,:) = [];
    set(handles.Listbox_handpicked_progression,'string', tmp)
    set(handles.Listbox_handpicked_progression,'value',min(ind,size(tmp,1)));
    guidata(hObject, handles);
end

function Edit_pvalue_module_tree_thres_Callback(hObject, eventdata, handles)
% hObject    handle to Edit_pvalue_module_tree_thres (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tmp = str2num(get(hObject,'String'));
if length(tmp)~=1 || tmp<0 || tmp>1
    set(hObject,'String',num2str(handles.pvalue_module_tree_thres));
else
    handles.pvalue_module_tree_thres = tmp;
end
guidata(hObject,handles)



% --- Executes during object creation, after setting all properties.
function Edit_pvalue_module_tree_thres_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Edit_pvalue_module_tree_thres (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Button_gene_module_adj.
function Button_gene_module_adj_Callback(hObject, eventdata, handles)
% hObject    handle to Button_gene_module_adj (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles, 'p_value') && ~isempty(handles.p_value) && isfield(handles, 'pvalue_module_tree_thres') && ~isempty(handles.pvalue_module_tree_thres)
    p_threshold = handles.pvalue_module_tree_thres;
    adj_modules = zeros(max(handles.idx));
    p_value = handles.p_value;
    for i=1:size(p_value,1)
        ind_tmp = find(p_value(i,:)<=p_threshold);
        adj_modules(ind_tmp,ind_tmp)=adj_modules(ind_tmp,ind_tmp) + 1;
    end
    figure; perm_modules = HCC_heatmap(adj_modules,'a'); axis off
    for i=1:length(perm_modules)
        text(i,-0,num2str(perm_modules(i)));
        text(i,length(perm_modules)+1,num2str(perm_modules(i)));
    end
end



% --- Executes on button press in Botton_view_handpicked.
function Botton_view_handpicked_Callback(hObject, eventdata, handles)
% hObject    handle to Botton_view_handpicked (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ind = get(handles.Listbox_handpicked_progression,'value');
if ind<1 || ind>length(handles.handpicked_progression)
    ind=1;
    set(handles.Listbox_handpicked_progression,'value',1); return
end
handles.view_tree_handle = view_tree(handles.handpicked_progression(ind), handles.exp_names, handles.color_code_names, handles.color_code_vectors, handles);
% get(handles.view_tree_handle,'userdata')
guidata(hObject,handles)




% --- Executes on button press in Button_Save_results.
function Button_Save_results_Callback(hObject, eventdata, handles)
% hObject    handle to Button_Save_results (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


if sum(~isfield(handles,{'data','probe_names','exp_names','color_code_names','color_code_vectors','filter_gene_ind','filtered_data','filter_std','filter_acceptable_nulls','filter_throw_x_at','x_at_probes_ind','kmeans_iter','kmeans_clusters','min_module_size','module_coherence_thres','adj','idx','permutation_iter','iter_merge_pvalue','pvalue_module_tree_thres','p_value','progression_modules','handpicked_progression'})) ...
        || isempty(handles.data_filename) ...
        || isempty(handles.data) ...
        || isempty(handles.probe_names) ...
        || isempty(handles.exp_names) ...
        || isempty(handles.filter_gene_ind) ...
        || isempty(handles.filtered_data) ...
        || isempty(handles.filter_std) ...
        || isempty(handles.filter_acceptable_nulls) ...
        || isempty(handles.filter_throw_x_at) ...
        || isempty(handles.kmeans_iter) ...
        || isempty(handles.kmeans_clusters) ...
        || isempty(handles.module_coherence_thres) ...
        || isempty(handles.idx) ...
        || isempty(handles.permutation_iter) ...
        || isempty(handles.iter_merge_pvalue) ...
        || isempty(handles.pvalue_module_tree_thres) ...
        || isempty(handles.candidate_trees) ...
        || isempty(handles.p_value)
   
    display('Results not saved due to one of the following:');
    display('   (1) User did not load data/result file');
    display('   (2) User did not perform clustering');
    display('   (3) User did not click the GO button to fit modules and trees');
    display('   (4) User did input and add any modules in step 4');
    return
end
data_filename = handles.data_filename;
data = handles.data;
probe_names = handles.probe_names;
exp_names = handles.exp_names;
color_code_names = handles.color_code_names;
color_code_vectors = handles.color_code_vectors;
filter_gene_ind = handles.filter_gene_ind;
filtered_data = handles.filtered_data;
% filtered_probe_names = handles.filtered_probe_names;
filter_std = handles.filter_std;
filter_acceptable_nulls = handles.filter_acceptable_nulls;
filter_throw_x_at = handles.filter_throw_x_at;
x_at_probes_ind = handles.x_at_probes_ind;
kmeans_iter = handles.kmeans_iter;
kmeans_clusters = handles.kmeans_clusters; % this variable was used in concensus clustering, but later was discarded, this para determined how many clusters do we want in each run of kmean during the 200 runs, this number was later fixed to 2, so this is no longer meaningful
min_module_size = handles.min_module_size;
module_coherence_thres = handles.module_coherence_thres;
adj = handles.adj;   % this variable was used in concensus clustering, but later was discarded, each element is how many times gene i and j got assigned into the same cluster
idx = handles.idx;
permutation_iter = handles.permutation_iter;
iter_merge_pvalue = handles.iter_merge_pvalue;
pvalue_module_tree_thres = handles.pvalue_module_tree_thres;
candidate_trees = handles.candidate_trees;
p_value = handles.p_value;          
progression_modules = handles.progression_modules;
handpicked_progression = handles.handpicked_progression;
result_filename = get(handles.Edit_save_filename,'String');
save(result_filename,'data_filename','result_filename',...
    'data','probe_names','exp_names','color_code_names','color_code_vectors',...
    'filter_gene_ind','filtered_data','filter_std','filter_acceptable_nulls','filter_throw_x_at','x_at_probes_ind',...
    'kmeans_iter','kmeans_clusters','module_coherence_thres','adj','idx',...
    'permutation_iter','iter_merge_pvalue','pvalue_module_tree_thres',...
    'candidate_trees','p_value','progression_modules','handpicked_progression','min_module_size');



% --- Executes on button press in Button_view_modules_expr.
function Button_view_modules_expr_Callback(hObject, eventdata, handles)
% hObject    handle to Button_view_modules_expr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
view_module_expr_window(handles.filtered_data,handles.filtered_probe_names, handles.idx, handles.exp_names, handles.color_code_names, handles.color_code_vectors, handles);






% --- Executes on button press in Button_corr_structure.
function Button_corr_structure_Callback(hObject, eventdata, handles)
% hObject    handle to Button_corr_structure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles,'data')==0 || isempty(handles.data)
    display('Please load data first, before trying to view correlation structure of the data !!!')
    return
end
data = handles.data(handles.filter_gene_ind,:);
exclude_ind = find(sum(isnan(data'))>=(size(data,2)/10));
data(exclude_ind,:)=[];
if sum(isnan(data(:)))~=0
    data = knnimpute(data);
end
K = min(2000,size(data,1));
[Y,I] = sort(std(data'),'descend');
data = per_gene_normalization(data(1:K,:));
corr_matrix = data*data'/(norm(data(1,:))^2); 
corr_matrix = corr_matrix - diag(diag(corr_matrix));
for i=1:size(corr_matrix,1), corr_matrix(i,i)=0; end
figure;hist(squareform(corr_matrix),200);




