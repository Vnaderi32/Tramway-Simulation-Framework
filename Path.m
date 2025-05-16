classdef Path
    properties
        ID              % Path ID(auto-generated, numeric)
        Length          % Path Length
        Duration        % Experimental Duration
        V_max           % Allowed speed during the path
        i               % Slope
        rho             % curviture
        number_of_data  % Number of samples
        stop_time       % Time for the stop
    end

    properties (Access = private)
        PathCount = 0   % Static counter for auto-generated numeric IDs
    end
    
    methods
        % Constructors
        function obj = Path(distance, duration, speed, t_stop)
            if nargin > 0
                id = getNextPathID();   % Using getNextPathID for generating ID
                obj.ID = id;

                obj.Length = distance;
                obj.Duration = duration;
                obj.V_max = speed;
                obj.stop_time = t_stop;
            end
        end

        % Path Info
        function displayInfo(obj)
            fprintf('Path: %s\n', obj.ID);
            fprintf('Length: %.2f [km]\n', obj.Length);
            fprintf('Experimental Duration: %.2f [s]\n', obj.Duration);
            fprintf('Allowed speed: %2.f [km/h]\n', obj.V_max);
            fprintf('Stop time: %2.f [s]\n', obj.stop_time);
        end
    end
end