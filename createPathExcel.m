function createPathExcel(station_number,t_s)
    try
        % Reset path ID counter
        getNextPathID(true);

        path_numbers = station_number - 1;
        path_list = Path.empty(path_numbers, 0);
        total_duration = 0;

        for j = 1:path_numbers
            fprintf('\n--- Enter general info for Path %d ---\n', j);
            len = input('Enter Path Length [m]: ');
            dur = input('Enter Experimental Duration [s]: ');
            v_max = input('Enter Allowed Speed [km/h]: ');
            stop_time = input('Enter Stop time [s]: ');

            total_duration = total_duration + dur + stop_time;
            p = Path(len, dur, v_max, stop_time);
            p.number_of_data = round((dur + stop_time) / t_s);  % total data points for this path
            path_list(j) = p;
        end
        
        % Always save path general data
        save('PathGeneralData.mat', 'path_list', 'station_number', 't_s');
        disp('✅ Path general data saved to "PathGeneralData.mat"');
        
        % Create Excel only if it doesn't exist
        file_name = 'PathData.xlsx';
        if ~isfile(file_name)
            t_col = (0:t_s:total_duration - t_s)';
            t_col = t_col + 1;
            number_of_data = length(t_col);

            path_number_col = zeros(number_of_data, 1);
            i_col = NaN(number_of_data, 1);
            rho_col = NaN(number_of_data, 1);

            start_idx = 1;
            for k = 1:length(path_list)
                dur = path_list(k).Duration + path_list(k).stop_time;
                data_count = round(dur / t_s);
                end_idx = start_idx + data_count - 1;

                if end_idx > number_of_data
                    end_idx = number_of_data;
                end

                path_number_col(start_idx:end_idx) = path_list(k).ID;
                start_idx = end_idx + 1;
            end

            T = table(path_number_col, t_col, i_col, rho_col, ...
                'VariableNames', {'Path number', 't[s]', 'i [m]', 'rho [m]'});

            writetable(T, file_name);

            fprintf('\x1b[32m✅ Excel file "%s" created successfully.\x1b[0m\n', file_name);
            fprintf('\x1b[34m📄 Please open the Excel file, fill in the "i [m]" and "rho [m]" columns,\nand then run the program again.\x1b[0m\n');
        else
            disp('ℹ️ Excel file already exists, skipped creation.');
        end

    catch
        fprintf('\x1b[31m❌ Something went wrong. The program will now exit safely.\x1b[0m\n');
        return;
    end
end