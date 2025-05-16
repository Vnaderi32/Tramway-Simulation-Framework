close all;
clear; 
clc;

% Constants
const = Constants();

%% Vehicle data

% Vehicle type
Vehicle_type = 'Tramva';

% Vehicle Name 
Vehicle_name = 'ATM Series 7000 (Eurotram)';

% Maximum motor power [kW]
P_max = 330;            % ( 4 * 106 kW )

% Maximum speed [km/h]
V_max_vehicle = 70;

% Supply system voltage [V_DC]
V_DC = 600;             % Overhead

% Mass [t]
m = 37.9;

% Gravitational mass [t]
G = m * const.g;

% Experimental value of beta
beta = 0.08;

% Equivalent mass [t]
m_e = m*(1+beta);

% Motorized wheelsets 
Nw = 6;

% Overall number of wheelsets
Mw = 8;

% Constant alpha
alpha = Nw/Mw;

% Wheel arrangement
Wheel_arr = 'Bo-2-Bo-Bo';

% Gauge [mm]
Gauge = 1445;

% Length [m]
len = 34.105; 

% Number of cars
n_cars = 7;

% Power absorber by auxiliaries [kW]
P_aux = 100;

%% Path data

if isfile('PathGeneralData.mat')
    data = load('PathGeneralData.mat');
    if isfield(data, 'station_number') && isfield(data, 't_s')
        station_number = data.station_number;
        t_s = data.t_s;
        fprintf('✅ Loaded station_number = %d and t_s = %.3f from saved data.\n', station_number, t_s);
    else
        t_s = input('-- Please enter the sampling time --\n');
        station_number = input('Also, please enter number of Stations: ');
        save('PathGeneralData.mat', 't_s', 'station_number');
        createPathExcel(station_number, t_s);
    end
else
    t_s = input('-- Please enter the sampling time --\n');
    station_number = input('Also, please enter number of Stations: ');
    save('PathGeneralData.mat', 't_s', 'station_number');
    createPathExcel(station_number, t_s);
end

file_name = 'PathData.xlsx';
path_objects = loadPathsFromExcel(file_name);

%% Path Calculations

% This is for summing data for all pathes ( Sum definition )        
v_km_h_fp =[];          % fp: Full Path
v_m_s_fp = [];
a_fp = [];
P_fp = [];
I_fp = [];

duration_full_path = 0;
all_data = 0;
for path_calculation_counter = 1:length(path_objects)
    current_path = path_objects{path_calculation_counter};
    
    V_max = min(current_path.V_max,V_max_vehicle);
    
    % Starting/Maximum tracktive effort [kN]
    F_tr_start = m_e * const.a_max;
    
    % Maximum adhesion force [kN]
    F_ad_max = const.f_ad * m * const.g * alpha;
    
    % Base speed
    V_base = P_max/F_tr_start;
    
    % Total distance 
    D_tot = current_path.Length;
    
    % Time Vector
    t = (0:t_s:current_path.Duration + current_path.stop_time - t_s)';  
    t = t + 1;
    number_of_data = length(t);
    all_data = all_data + number_of_data;
    %% Initial Values
    
    % Distence [m]
    d = zeros(number_of_data,1);
    
    % Velocity [m/s]
    v_m_s = zeros(number_of_data,1);
    
    % Velocity [km/h]
    v_km_h = v_m_s * 3.6;
    
    % Tractive force [kN]
    F_tr = zeros(number_of_data,1);
    
    %% Resistances
    
    % Relative resistance 
    r_0 = 2.7 + 0.04*(v_km_h/10).^2;
    
    % Resistance to forward motion [N]
    R_0 = r_0 * G;
    
    % Grade resistances [N]
    R_i = (current_path.i /1000)*G*1000;
    
    % Curviture resistances [N]
    R_c = zeros(number_of_data,1);
    valid = current_path.rho ~= 0;
    R_c(valid) = (const.c ./ current_path.rho(valid)) * G * 1000;
    
    % Total resistances [N]
    R = R_c + R_i + R_0;
    
    %% Tractive Force
    
    start_mode = v_m_s <= V_base;       %If V<V_base
    F_tr(start_mode) = F_tr_start;      % F_tr = F_tr_start
    
    F_tr(~start_mode) = P_max ./ v_m_s(~start_mode); % Else: F_tr = P_max / V
    
    %% Acceleration [m/s^2]
    a = (F_tr*1000 - R)/(m_e*1000);
    
    %% Velocity
    
    dt = diff(t);                        
    dv = a(2:end) .* dt;                 
    v_m_s(2:end) = v_m_s(1) + cumsum(dv);
    
    % Velocity [km/h]
    v_km_h = v_m_s * 3.6;
    
    %% Tractive Force Vector Correction
    start_mode = v_m_s <= V_base;       %If V<V_base
    F_tr(start_mode) = F_tr_start;      % F_tr = F_tr_start
    
    F_tr(~start_mode) = P_max ./ v_m_s(~start_mode); % Else: F_tr = P_max / V
    
    %% Distance
    for j=2:number_of_data
        dt = t(j)-t(j-1);
        d(j) = ((v_m_s(j)+v_m_s(j-1))/2)*dt+d(j-1);
    end
    
    %% Coasting Mode
    
    % Tractive force vector correction
     for j=1:number_of_data
         if v_km_h(j) > V_max
             F_tr(j)=0;
         end
     end
    
     % Acceleration Vector Correction
     a = (F_tr*1000 - R)/(m_e*1000);
    
    % Velocity vector correction
    dt = diff(t);                        
    dv = a(2:end) .* dt;                 
    v_m_s(2:end) = v_m_s(1) + cumsum(dv);
    
    % Distance Vector correction
    for j=2:number_of_data
        dt = t(j)-t(j-1);
        d(j) = ((v_m_s(j)+v_m_s(j-1))/2)*dt+d(j-1);
    end
    
    %% Breaking Mode
    d_br = 0.5*v_m_s.^2/const.a_br;
    
    break_moment = D_tot - d - d_br < 0;
    a(break_moment)= -const.a_br;
    
    %% Update
    % Velocity vector correction
    dt = diff(t);                        
    dv = a(2:end) .* dt;                 
    v_m_s(2:end) = v_m_s(1) + cumsum(dv);
    
    % Velocity [km/h]
    v_km_h = v_m_s * 3.6;
    
    
    % Distance Vector correction
    for j=2:number_of_data
        dt = t(j)-t(j-1);
        d(j) = ((v_m_s(j)+v_m_s(j-1))/2)*dt+d(j-1);
    end
    
    %% Break Mode Tractive Force
    
     F_tr(break_moment)= a(break_moment)*m_e;
    
     %% Adhesion Effect
     F_ad = (const.f_ad./(1+0.011.*v_km_h))*G*alpha;
     slip = abs(F_tr)>F_ad; 

     %% Full stop
     full_stop = v_m_s < 0;
     v_m_s(full_stop) = 0;
     a(full_stop)=0;
     F_tr(full_stop)=0;
     d_br(full_stop)=0;
    
     %% Update
    
     % Velocity [km/h]
     v_km_h = v_m_s * 3.6;
    
     % Distance Vector correction
    for j=2:number_of_data
        dt = t(j)-t(j-1);
        d(j) = ((v_m_s(j)+v_m_s(j-1))/2)*dt+d(j-1);
    end
    
    %%
    %% Electrical specification
    
    % Power [kW]
    P = zeros(number_of_data,1);
    tractive_mode = F_tr < 0;
    P(tractive_mode) = (F_tr(tractive_mode).*v_m_s(tractive_mode))./const.eta_tr + P_aux;
    P(~tractive_mode) = (F_tr(~tractive_mode).*v_m_s(~tractive_mode))*const.eta_br + P_aux;
    
    reg_break = v_km_h < 15;
    P(reg_break & break_moment)= P_aux;
    
    % Current [A]
    I = (P * 1000)/V_DC;
    
    % Energy [kW/h]
    E=zeros(number_of_data,1);
    
    E(1)= P(1)/3600;
    
    for j=2:number_of_data
        dt = t(j)-t(j-1);
        E(j)=(P(j)/3600)*dt + E(j-1);
    end
    

    duration_full_path = duration_full_path + length(t);
    %% 
    %% FIGURES

    %% Plot V
    figure
    plot(t, v_km_h)
    xlabel('Time [s]')
    ylabel('Velocity [km/h]')
    title(sprintf('Velocity [km/h] vs. Time - Path %d', current_path.ID));
    grid on

    %% Plot Acceleration
    figure
    plot(t, a)
    xlabel('Time [s]')
    ylabel('Acceleration [m/s^2]')
    title(sprintf('Acceleration [m/s^2] vs. Time - Path %d', current_path.ID));
    grid on

    %% Plot Power
    figure
    plot(t, P)
    xlabel('Time [s]')
    ylabel('Power [kW]')
    title(sprintf('Power [kW] vs. Time - Path %d', current_path.ID));
    grid on

    %% Plot Current
    figure
    plot(t, I)
    xlabel('Time [s]')
    ylabel('Current [A]')
    title(sprintf('Current [A] vs. Time - Path %d', current_path.ID));
    grid on

    %% Plot Energy
    figure
    plot(t, E)
    xlabel('Time [s]')
    ylabel('Energy [kW/h]')
    title(sprintf('Energy [kW/h] vs. Time - Path %d', current_path.ID));
    grid on

    %% Plot Distance
    figure
    plot(t, d)
    xlabel('Time [s]')
    ylabel('Distance [m]')
    title(sprintf('Distance [m] vs. Time - Path %d', current_path.ID));
    grid on

    %% Plot Tractive Force vs Velocity

    % Velocity used for plotting F_tr
    V_axis_km = (1:1:round(max(v_km_h))+20)';
    V_axis_m = V_axis_km/3.6;

    F_tr_axis = zeros(size(V_axis_m));

    % F_tr for plotting
    plot_index = V_axis_m <= V_base;       %If V<V_base
    F_tr_axis(plot_index) = F_tr_start;      % F_tr = F_tr_start

    F_tr_axis(~plot_index) = P_max ./ V_axis_m(~plot_index); % Else: F_tr = P_max / V
    F_tr_axis = F_tr_axis';
    figure
    plot(V_axis_km, F_tr_axis)
    xlabel('Velocity [km/h]')
    ylabel('Tractive Effort [kN]')
    title(sprintf('Tractive Effort [kN] vs. Velocity [km/h] - Path %d', current_path.ID));
    grid on
    hold on

    % F_ad for plotting
    F_ad_axis = (const.f_ad./(1+0.011.*V_axis_km))*G*alpha;
    F_ad_axis = F_ad_axis';
    plot(V_axis_km, F_ad_axis)

    %% Finding Full Path Parameters
    t_fp = (0:t_s:duration_full_path*t_s - t_s)';              % fp: Full Path
    v_km_h_fp =[v_km_h_fp;v_km_h];
    v_m_s_fp = [v_m_s_fp;v_m_s];
    a_fp = [a_fp;a];
    P_fp = [P_fp;P];
    I_fp = [I_fp;I];

end

%% Total Path 
    
    %% Plot Acceleration
    figure
    plot(t_fp, a_fp)
    xlabel('Time [s]')
    ylabel('Acceleration [m/s^2]')
    title('Acceleration [m/s^2] vs. Time - Full Path');
    grid on
    
    %% Plot Power
    figure
    plot(t_fp, P_fp)
    xlabel('Time [s]')
    ylabel('Power [kW]')
    title('Power [kW] vs. Time - Full Path');
    grid on
    
    %% Plot Current
    figure
    plot(t_fp, I_fp)
    xlabel('Time [s]')
    ylabel('Current [A]')
    title('Current [A] vs. Time - Full Path');
    grid on
    
    %% Plot Energy

     % Energy [kW/h]
    E_fp=zeros(all_data,1);
    
    E_fp(1)= P_fp(1)/3600;
    
    for j=2:all_data
        dt = t_fp(j)-t_fp(j-1);
        E_fp(j)=(P_fp(j)/3600)*dt + E_fp(j-1);
    end
    
    figure
    plot(t_fp, E_fp)
    xlabel('Time [s]')
    ylabel('Energy [kW/h]')
    title('Energy [kW/h] vs. Time - Full Path');
    grid on
    
    %% Plot Distance

    d_fp=zeros(all_data,1);

    % Full path distance 
    for j=2:all_data
        dt = t_fp(j)-t_fp(j-1);
        d_fp(j) = ((v_m_s_fp(j)+v_m_s_fp(j-1))/2)*dt+d_fp(j-1);
    end
    figure
    plot(t_fp, d_fp)
    xlabel('Time [s]')
    ylabel('Distance [m]')
    title('Distance [m] vs. Time - Full Path');
    grid on

    %% Plot Tractive Force vs Velocity
    
    % Velocity used for plotting F_tr
    V_axis_km = (1:1:round(max(v_km_h_fp)))';
    V_axis_m = V_axis_km/3.6;
    
    F_tr_axis = zeros(size(V_axis_m));

    % F_tr for plotting
    plot_index = V_axis_m <= V_base;       %If V<V_base
    F_tr_axis(plot_index) = F_tr_start;      % F_tr = F_tr_start
    
    F_tr_axis(~plot_index) = P_max ./ V_axis_m(~plot_index); % Else: F_tr = P_max / V
    F_tr_axis = F_tr_axis';

    % F_ad for plotting
    F_ad_axis = (const.f_ad./(1+0.011.*V_axis_km))*G*alpha;
    F_ad_axis = F_ad_axis';
    
    figure
    plot(V_axis_km, F_tr_axis)
    xlabel('Velocity [km/h]')
    ylabel('Tractive Effort [kN]')
    title('Tractive Effort [kN] vs. Velocity [km/h] - Full Path');
    grid on
    hold on
    plot(V_axis_km, F_ad_axis);
    
    %% Error Calculator
    error = 0;
    for i=1:length(path_objects)
        error = error + path_objects{i}.Length;
    end
    error = error - d_fp(end);
    fprintf('\n\n\n✅✅✅Tram stays with %.4fm distance from the station!✅✅✅\n',abs(error));

    %% Commercial Speed

    Commercial_Speed = (d_fp(end) / t_fp(end))*3.6;

    %% Plot V
    figure
    plot(t_fp, v_km_h_fp);
    hold on;
    xlabel('Time [s]')
    ylabel('Velocity [km/h]')
    title('Velocity [km/h] vs. Time - Full Path');
    grid on
    yline(Commercial_Speed,'red');
    legend('Velocity [km/h]','Commercial Speed [km/h]');
