% MTOW Estimation

function output = MTOW_estimation(mission_type)

% Check if input is 1 or 2
if ~(mission_type == 1 || mission_type == 2)
    error('Input must be 1 or 2');
end

% Define variables

% W_TO = Maximum Take-off weight [kg]
% W_fuel = Fuel weight [kg]
% W_empty = Empty weight [kg]
W_payload = 27*2.205; % Payload weight [lb]

n_prop = 0.67; % Propulsive efficiency
R_max = 306*1000*4 *3.281; % Maximum range [ft] (based on travelling at the cruise speed of 245 km/h for 4 hours)
SFC_cruise = 0.25; % Specific fuel consumption in cruise [lb/h/lb]
SFC_loiter = 0.25; % Specific fuel consumption in loiter [lb/h/lb]

C_cruise = SFC_cruise/(3600*550); % Range coefficient in cruise [1/ft]
C_loiter = SFC_loiter/(3600*550); % Range coefficient in cruise [1/ft]

E_max = 4*3600; % Maximum endurance [s]
LD_max = 10; % Maximum lift-to-drag ratio
V_s = 12*3.281; % Stall speed [ft/s]

a = 0.0015; % Coefficient for empty weight estimation - Estimate based on historical data
b = 0.3363; % Coefficient for empty weight estimation - Estimate based on historical data
Comp = 0.9; % Composite correction factor


%% Mission profiles

% Define fuel ratios

f1 = 0.98; % Taxi and take-off fuel ratio (W2/W1)
f2 = 0.97; % Climb fuel ratio (W3/W2)
f3 = 0.99; % Descent fuel ratio (W4/W3)
f4 = 0.997; % Approach and landing fuel ratio (W5/W4)

f_safety = 1.05; % Safety factor





if mission_type == 1
    % Mission profile 1 - surveillance mission (half cruise, half loiter)
    % Consists of 7 stages: 1. Taxi, 2. Take-off, 3. Climb, 4. Cruise, 5. Loiter,  6. Descent, 7. Landing

    % Cruise weight fraction

    f5 = exp(-R_max*0.5*C_cruise/(n_prop*LD_max)); % Cruise weight fraction

    % Loiter weight fraction

    f6 = exp(-E_max*0.5*C_loiter*V_s*1.3/(0.866*n_prop*LD_max)); % Loiter weight fraction

    % Calculate fuel weight fraction

    f = (f1)*(f2)*(f3)*(f4)*(f5)*(f6);
    F = f_safety*(1-f); % Fuel ratio - W_fuel/W_TO
elseif mission_type == 2
    % Mission profile 2 - payload delivery (main section is cruise)
    % Consists of 9 stages: 1. Taxi, 2. Take-off, 3. Climb, 4. Cruise, 5. Descent,  6. Ascent, 7. Cruise, 8. Descent, 9. Landing

    % Cruise weight fraction

    f5 = exp(-R_max*0.5*C_cruise/(n_prop*LD_max)); % Cruise weight fraction

    % Calculate fuel weight fraction

    f = (f1)*(f2)*(f2)*(f3)*(f3)*(f4)*(f5); % Fuel ratio - W_fuel/W_TO
    F = f_safety*(1-f); % Fuel ratio - W_fuel/W_TO
end


% Simultaneous equations:
% x is W_empty/W_TO, y is W_TO
syms x y

% Empty weight estimation

eqn1 = x == Comp*a*y + b; % Empty weight ratio - W_empty/W_TO - estimation [kg]

% MTOW equation

eqn2 = y == W_payload/(1 - x - F); % Maximum take-off weight [kg]

% Solve equations

sol = solve([eqn1, eqn2], [x, y]);

% Convert the solutions to double
sol.x = double(sol.x);
sol.y = double(sol.y);

% Check if the solutions are valid

for n = 1:length(sol.x)
    if sol.x(n)>0 && sol.x(n)<1
        sol_x = sol.x(n);
    else
    end
end

for n = 1:length(sol.y)
    if sol.y(n)>0
        sol_y = sol.y(n);
    else
        n=n+1;
    end
end

% Convert to metric

sol_y = sol_y*0.453592; % Convert to kg

% Display results

fprintf('The empty weight ratio is %.2f\n', sol_x);

fprintf('The maximum take-off weight is %.2f kg\n', sol_y);

fprintf('The fuel ratio is %.2f\n', F);

output = [sol_x, sol_y];

end




