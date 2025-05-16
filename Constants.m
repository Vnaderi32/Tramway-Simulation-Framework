classdef Constants
    properties
        a_max = 1;            % Max acceleration [m/s^2]
        a_br = 1.1;           % Max deceleration [m/s^2]
        eta_tr = 0.8;         % Traction efficiency
        eta_br = 0.8;         % Braking efficiency
        beta_br = 0.7;        % Regenerative braking portion
        f_ad = 0.25;          % Adhesion coefficient
        g = 9.80665;             % Gravity [m/s^2]
        c = 0.8;              % Curvature resistance constant
    end
end