function path_objects = loadPathsFromExcel(file_name)

    % file_name='PathData.xlsx';

    % Loading general data about the paths
    if isfile('PathGeneralData.mat')
        load('PathGeneralData.mat','path_list');
    else
        error('⚠️ File (PathGeneralData.mat) not found!');
    end


    % Reading file data
    data = readtable(file_name, 'ReadVariableNames', false);

    % Extracting columns
    path_number_col = table2array(data(:,1));
    t_col = table2array(data(:,2));
    i_col = table2array(data(:,3));
    rho_col = table2array(data(:,4));

    % Identifying paths
    unique_paths = unique(path_number_col);
    num_paths = numel(unique_paths);

    % Initiallization 
    path_objects = cell(1, num_paths);  % Output in cell arrays

    % idx =2;
    for idx = 1:num_paths

        p_id = unique_paths(idx);
        
        general_info = path_list([path_list.ID] == p_id);

        % Choising the related paths
        mask = (path_number_col == p_id);
        i_vector = i_col(mask);
        rho_vector = rho_col(mask);

        % t_vector = t_col(mask);
        
        p = Path();
        p.ID = p_id;
        p.i = i_vector;
        p.rho = rho_vector;

        p.Length = general_info.Length;
        p.Duration = general_info.Duration;
        p.V_max = general_info.V_max;
        p.stop_time = general_info.stop_time;
        p.number_of_data = sum(mask); 
        path_objects{idx} = p;

        var_name = sprintf('path%d', p_id);
        assignin('base', var_name, p);
    end

    disp('✅ Path objects successfully created and assigned to workspace.');
end