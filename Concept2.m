function [tau_s_h, tau_s_v, q_s_v, q_s_h, axial_stress_h, axial_stress_v, V_z_h, V_y_v] = Concept2(Airfoil_w, y_interp_h, L_h_interp, y_interp, L_v_interp, t_v, t_h, span_v, span_h, C_r_v, C_r_h)

% Stress Analysis
[V_z_h, M_z_h, V_y_v, M_y_v, y_h, z_v] = Force_dist(y_interp_h, L_h_interp, y_interp, L_v_interp);

%% Calculate Ixx of tail skin

% Import airfoil coordinates

if Airfoil_w == 'NACA 4415'
    naca0012 = csvread('n0012-il.csv');
elseif Airfoil_w == 'NACA 4418'
    naca0012 = csvread('n0012-il(w4418).csv');
end
% 
% Remove trailing section of airfoil profile - this is taken up by the control surface
naca0012 = naca0012(naca0012(:,1) < max(naca0012(:,1)) * 0.7, :);

% Repeat first coordinate pair at end to make closed section
naca0012 = [naca0012; naca0012(1,:)];

% convert from mm to m
naca0012 = naca0012 /1000;

% Scale the airfoil data to the size of the tails
naca0012_h = naca0012 * (C_r_h/max(naca0012(:,1)));
naca0012_v = naca0012 * (C_r_v/max(naca0012(:,1)));

% Taper adjustment
% Compute the scaled chord length
C = C_r_h - (C_r_h * ((2 * y_h / span_h) * (1 - t_h)));
C_v = C_r_v - (C_r_v * ((2 * z_v / span_v) * (1 - t_v)));
c_taper_h = double(C/C_r_h);
c_taper_v = (C_v/C_r_v);

t_s_h = 0.001; % Thickness of the horizontal tail section [m]
t_s_v = 0.001; % Thickness of the vertical tail section [m]

% Scale the cross-sections
x_scaled_h = double(naca0012_h(:,1) * C / C_r_h);
z_scaled_h = double(naca0012_h(:,2) * C / C_r_h);
x_scaled_v = double(naca0012_v(:,1) * C_v / C_r_v);
y_scaled_v = double(naca0012_v(:,2) * C_v / C_r_v);

% As airfoil is symmetric, neutral x-axis is given by the chord line
naca_chord_h = naca0012_h(:,3);
naca_chord_v = naca0012_v(:,3);
naca_chord_h = naca_chord_h .* c_taper_h;
naca_chord_v = naca_chord_v .* c_taper_v;

% Calculate the moment of inertia of the cross-sections

% Horizontal tail
% Compute the length of the current section
L_j = sqrt((diff(x_scaled_h)).^2 + (diff(z_scaled_h).^2));
% Compute the angle of the current section to its local x-axis
theta_i = atan(abs((diff(z_scaled_h)) ./ (diff(x_scaled_h))));
% Compute the moment of inertia of the sections about their local axes
I_xx_s = ((L_j .* t_s_h) / 12) .* ((L_j.^2 .* cos(theta_i).^2) + (t_s_h.^2 .* sin(theta_i).^2));
% Compute the perpendicular distance between the straight line local x-axis and the neutral x-axis
y_bar = abs(((z_scaled_h(1:end-1, :) + z_scaled_h(2:end, :)) ./ 2) - ((naca_chord_h(1:end-1, :) + naca_chord_h(2:end, :)) ./ 2));
% Compute the moment of inertia of the sections about the cross-section centroid
I_xx_h_s = I_xx_s + (L_j .* t_s_h .* y_bar.^2);

% Compute the overall moments of inertia about the centroid of the cross sections
I_xx_h = sum(I_xx_h_s, 1);

% Vertical tail
% Compute the length of the current section
L_j = sqrt((diff(x_scaled_v)).^2 + (diff(y_scaled_v).^2));
% Compute the angle of the current section to its local x-axis
theta_i = atan(abs((diff(y_scaled_v)) ./ (diff(x_scaled_v))));
% Compute the moment of inertia of the sections about their local axes
I_xx_s = ((L_j .* t_s_v) / 12) .* ((L_j.^2 .* cos(theta_i).^2) + (t_s_v.^2 .* sin(theta_i).^2));
% Compute the perpendicular distance between the straight line local x-axis and the neutral x-axis
y_bar = abs(((y_scaled_v(1:end-1, :) + y_scaled_v(2:end, :)) ./ 2) - ((naca_chord_v(1:end-1, :) + naca_chord_v(2:end, :)) ./ 2));
% Compute the moment of inertia of the sections about the cross-section centroid
I_xx_v_s = I_xx_s + (L_j .* t_s_v .* y_bar.^2);

% Compute the overall moments of inertia about the centroid of the cross sections
I_xx_v = sum(I_xx_v_s, 1);


%% Adjust Ixx of section to account for reinforcing components - ignore effects of ribs as they are at the extreme ends of the wing



I_xx_h_bending = I_xx_h;
I_xx_v_bending = I_xx_v;
I_xx_h_shear = I_xx_h;
I_xx_v_shear = I_xx_v;

% Horizontal tail
y_h = y_interp_h; % Spanwise locations of the horizontal tail
% Vertical tail
z_v = y_interp; % Spanwise locations of the vertical tail

% Additional Ixx due to C beam spar - for simplicity, curvature of spar caps is ignored

t_spar_web_h = zeros(size(y_h)); % Thickness of the spar web [m]
t_spar_cap_h = zeros(size(y_h)); % Thickness of the spar cap [m]
w_spar_h = zeros(size(y_h)); % Total width of the spar [m]
t_spar_web_h(:) = 0.002; % Thickness of the spar web [m]
t_spar_cap_h(:) = 0.002; % Thickness of the spar cap [m]
target = 0.7 * C_r_h;
differences = abs(naca0012_h(:,1) - target);
% Sort the differences and get their indices
[sorted_differences, indices] = sort(double(differences));
% The indices of the two closest values are the first two elements of indices
index1 = indices(1);
index2 = indices(2); 
h_spar_web_h = abs(naca0012_h(index2,2) - naca0012_h(index1,2)) - 2 * t_s_h - 2 * t_spar_cap_h(1); % Height of the spar web [m]
w_spar_h(:) = h_spar_web_h / 2; % Total width of the spar [m]

t_spar_web_v = zeros(size(z_v)); % Thickness of the spar web [m]
t_spar_cap_v = zeros(size(z_v)); % Thickness of the spar cap [m]
w_spar_v = zeros(size(z_v)); % Total width of the spar [m]
t_spar_web_v(:) = 0.003; % Thickness of the spar web [m]
t_spar_cap_v(:) = 0.003; % Thickness of the spar cap [m]
target = 0.7 * C_r_v;
differences = abs(naca0012_v(:,1) - target);
% Sort the differences and get their indices
[sorted_differences, indices] = sort(double(differences));
% The indices of the two closest values are the first two elements of indices
index1 = indices(1);
index2 = indices(2);
h_spar_web_v = abs(naca0012_v(index2,2) - naca0012_v(index1,2)) - 2 * t_s_v - 2 * t_spar_cap_v(1); % Height of the spar web [m]
w_spar_v(:) = h_spar_web_v / 2; % Total width of the spar [m]
    
% Scale the cross-section
h_spar_web_h = h_spar_web_h .* c_taper_h; % Height of the spar web [m]
w_spar_h = w_spar_h .* c_taper_h; % Total width of the spar [m]
h_spar_web_h = h_spar_web_h .* c_taper_h; % Height of the spar web [m]
w_spar_h = w_spar_h .* c_taper_h; % Total width of the spar [m]

h_spar_web_v = h_spar_web_v .* c_taper_v; % Height of the spar web [m]
w_spar_v = w_spar_v .* c_taper_v; % Total width of the spar [m]
h_spar_web_v = h_spar_web_v .* c_taper_v; % Height of the spar web [m]
w_spar_v = w_spar_v .* c_taper_v; % Total width of the spar [m]

% Compute the additional Ixx due to the spar
I_xx_spar_h = (1/12) .* (t_spar_web_h .* h_spar_web_h.^3) + 2*(((w_spar_h .* t_spar_cap_h.^3)/12) + ((t_spar_cap_h .* w_spar_h .* (t_spar_cap_h + h_spar_web_h).^2)/4));
I_xx_h_bending = I_xx_h_bending + I_xx_spar_h;
I_xx_h_shear = I_xx_h_shear + I_xx_spar_h;

I_xx_spar_v = (1/12) .* (t_spar_web_v .* h_spar_web_v.^3) + 2*(((w_spar_v .* t_spar_cap_v.^3)/12) + ((t_spar_cap_v .* w_spar_v .* (t_spar_cap_v + h_spar_web_v).^2)/4));
I_xx_v_bending = I_xx_v_bending + I_xx_spar_v;
I_xx_v_shear = I_xx_v_shear + I_xx_spar_v;

% Additonal Ixx due to stringers
% Input the number of stringers per surface
n_stringers_h = 2*input('Enter the number of stringers per surface for the horizontal tail: ');
n_stringers_v = 2*input('Enter the number of stringers per surface for the vertical tail: ');

% Stringers assumed to have square cross-section with side length equal to the thickness of the skin
a_stringer_h = t_s_h;
a_stringer_v = t_s_v;

% Compute the indices at which to place the stringers
indices_h = round(linspace(1, length(x_scaled_h), n_stringers_h));
indices_v = round(linspace(1, length(y_scaled_v), n_stringers_v));

% Find the distance from the neutral axis of the centre of each stringer
y_bar_stringers_h = abs(z_scaled_h(indices_h) - (a_stringer_h/2) - naca_chord_h(indices_h));
y_bar_stringers_v = abs(y_scaled_v(indices_v) - (a_stringer_v/2) - naca_chord_v(indices_v));

% Compute the additional Ixx due to the stringers
I_xx_stringers_h = sum(a_stringer_h.^2 .* y_bar_stringers_h) .* ((a_stringer_h.^2)/12);
I_xx_stringers_v = sum(a_stringer_v.^2 .* y_bar_stringers_v) .* ((a_stringer_v.^2)/12);

% Add the additional Ixx to the existing Ixx
I_xx_h_bending = I_xx_h_bending + I_xx_stringers_h;
I_xx_v_bending = I_xx_v_bending + I_xx_stringers_v;
I_xx_h_shear = I_xx_h_shear + I_xx_stringers_h;
I_xx_v_shear = I_xx_v_shear + I_xx_stringers_v;

%% Calculate maximum axial stress distribution

% Initialize arrays to store the axial stress
axial_stress_h = zeros(size(y_h));
axial_stress_v = zeros(size(z_v));
axial_stress_h_s = zeros(length(naca0012(:,1)), length(y_h));
axial_stress_v_s = zeros(length(naca0012(:,1)), length(z_v));

% M_y_h and M_y_v are the bending moments for horizontal and vertical tails respectively

% Horizontal tail
% Compute the perpendicular distance between the straight line local x-axis and the neutral x-axis
[max_y_h, index] = max(z_scaled_h);
y_bar_tens = max_y_h - naca_chord_h(index(1),:);
[min_y_h, index] = min(z_scaled_h);
y_bar_comp = min_y_h - naca_chord_h(index(1),:);
y_bar_tens = y_bar_tens(:);
y_bar_comp = y_bar_comp(:);
M_z_h = M_z_h(:);
I_xx_h_bending = I_xx_h_bending(:);
% Compute the max axial stress
max_axial_stress_h_tens = (M_z_h .* y_bar_tens) ./ I_xx_h_bending;
max_axial_stress_h_comp = (M_z_h .* y_bar_comp) ./ I_xx_h_bending;

% Vertical tail

% Compute the perpendicular distance between the straight line local x-axis and the neutral x-axis
[max_z_v, index] = max(y_scaled_v);
z_bar_tens = max_z_v - naca_chord_v(index(1),:);
[min_z_v, index] = min(y_scaled_v);
z_bar_comp = min_z_v - naca_chord_v(index(1),:);
z_bar_tens = z_bar_tens(:);
z_bar_comp = z_bar_comp(:);
M_y_v = M_y_v(:);
I_xx_v_bending = I_xx_v_bending(:);
% Compute the max axial stress
max_axial_stress_v_tens = (M_y_v .* z_bar_tens) ./ I_xx_v_bending;
max_axial_stress_v_comp = (M_y_v .* z_bar_comp) ./ I_xx_v_bending;

ult_sigma_h = max(max_axial_stress_h_tens);
ult_sigma_v = max(max_axial_stress_v_tens);
ult_sigma_h_comp = max(max_axial_stress_h_comp);
ult_sigma_v_comp = max(max_axial_stress_v_comp);

fprintf('The ultimate tensile stress in the horizontal tail skin is %.2f MPa\n', ult_sigma_h/1e6);
fprintf('The ultimate tensile stress in the vertical tail skin is %.2f MPa\n', ult_sigma_v/1e6);
fprintf('The ultimate compressive stress in the horizontal tail skin is %.2f MPa\n', ult_sigma_h_comp/1e6);
fprintf('The ultimate compressive stress in the vertical tail skin is %.2f MPa\n', ult_sigma_v_comp/1e6);

% Plot the maximum compressive and tensile stress in each cross-section against the spanwise location
figure;
plot(y_h, max_axial_stress_h_comp, 'r', y_h, max_axial_stress_h_tens, 'b');
hold on;
plot(z_v, max_axial_stress_v_comp, 'm', z_v, max_axial_stress_v_tens, 'c');
title('Max Compressive and Tensile Stress vs. Spanwise Location');
xlabel('Spanwise Location');
ylabel('Stress');
legend('Max Compressive Stress - Horizontal Tail', 'Max Tensile Stress - Horizontal Tail', 'Max Compressive Stress - Vertical Tail', 'Max Tensile Stress - Vertical Tail');
hold off;

%% Shear distribution in skin

% Initialize the shear flow distributions
q_s_h = zeros(length(naca0012(:,1)), length(y_h));
q_s_v = zeros(length(naca0012(:,1)), length(y_h));
L_h = zeros(length(naca0012(:,1)), length(y_h));
L_v = zeros(length(naca0012(:,1)), length(y_h));
P_h = zeros(length(naca0012(:,1)), length(y_h));
P_v = zeros(length(naca0012(:,1)), length(y_h));
theta_h = zeros(length(naca0012(:,1)), length(y_h));
theta_v = zeros(length(naca0012(:,1)), length(y_h));
y_bar_h = zeros(length(naca0012(:,1)), length(y_h));
y_bar_v_s = zeros(length(naca0012(:,1)), length(y_h));

% Calculate the shear flow distribution for the horizontal tail
for i = 1:length(y_h)
    for j = 1:length(naca0012(:,1))-1
        % Compute the length of the current section
        L_h(j,i) = sqrt((x_scaled_h(j+1,i) - x_scaled_h(j,i))^2 + (z_scaled_h(j+1,i) - z_scaled_h(j,i))^2);
        % Compute the angle of the current section to its local x-axis
        theta_h(j,i) = atan((z_scaled_h(j+1,i) - z_scaled_h(j,i)) / (x_scaled_h(j+1,i) - x_scaled_h(j,i)));
        % Compute the moment arm of the current section
        P_h(j,i) = sqrt(((((z_scaled_h(j,i)+z_scaled_h(j+1,i))/2) - z_scaled_h(1,i))^2) + ((((x_scaled_h(j,i)+x_scaled_h(j+1,i))/2) - x_scaled_h(1,i))^2)) * cos((pi/2) - theta_h(j,i) - atan((((z_scaled_h(j,i)+z_scaled_h(j+1,i))/2) - z_scaled_h(1,i)) / (((x_scaled_h(j,i)+x_scaled_h(j+1,i))/2) - x_scaled_h(1,i))));
        % Compute the perpendicular distance between the straight line local x-axis and the neutral x-axis
        y_bar_h(j,i) = ((z_scaled_h(j+1,i) + z_scaled_h(j,i)) / 2) - naca_chord_h(j,i);
        % Compute the shear flow distribution
        q_s_h(j,i) = ((-V_z_h(i) / I_xx_h_shear(i)) * sum(t_s_h .* y_bar_h(1:j,i) .* L_h(1:j,i)));
    end
end

% Initialize the shear force at the cut point
q_s_0_h = zeros(size(y_h));
q_s_0_v = zeros(size(z_v));
eps_h = zeros(size(y_h));
eps_v = zeros(size(z_v));
A_h = zeros(size(y_h));
A_v = zeros(size(z_v));

% Calculate the shear flow at the cut point for the horizontal tail
for i = 1:length(y_h)
    % Compute the scaled chord length
    C = C_r_h - (C_r_h * ((2 * y_h(i) / span_v) * (1 - t_h)));
    % Compute the moment arm of the external load
    eps_h(i) = abs((C * 0.25) - x_scaled_h(1,i));
    % Compute the area using the trapezoidal rule
    A_h(i) = trapz(double(x_scaled_h(:,i)), double(z_scaled_h(:,i)));
    % Compute the shear force at the cut point
    q_s_0_h(i) = (V_z_h(i)*eps_h(i) - sum(P_h(:,i).*q_s_h(:,i).*L_h(:,i))) / (2*A_h(i));
end

% Calculate the shear flow distribution for the horizontal tail
for i = 1:length(y_h)
    q_s_h(:,i) = q_s_h(:,i) + q_s_0_h(i);
end

% Calculate the shear flow distribution for the vertical tail
for i = 1:length(z_v)
    for j = 1:length(naca0012(:,1))-1
        % Compute the length of the current section
        L_v(j,i) = sqrt((x_scaled_v(j+1,i) - x_scaled_v(j,i))^2 + (y_scaled_v(j+1,i) - y_scaled_v(j,i))^2);
        % Compute the angle of the current section to its local x-axis
        theta_v(j,i) = atan((y_scaled_v(j+1,i) - y_scaled_v(j,i)) / (x_scaled_v(j+1,i) - x_scaled_v(j,i)));
        % Compute the moment arm of the current section
        P_v(j,i) = sqrt(((((y_scaled_v(j,i)+y_scaled_v(j+1,i))/2) - y_scaled_v(1,i))^2) + ((((x_scaled_v(j,i)+x_scaled_v(j+1,i))/2) - x_scaled_v(1,i))^2)) * cos((pi/2) - theta_v(j,i) - atan((((y_scaled_v(j,i)+y_scaled_v(j+1,i))/2) - y_scaled_v(1,i)) / (((x_scaled_v(j,i)+x_scaled_v(j+1,i))/2) - x_scaled_v(1,i))));
        % Compute the perpendicular distance between the straight line local x-axis and the neutral x-axis
        y_bar_v_s(j,i) = ((y_scaled_v(j+1,i) + y_scaled_v(j,i)) / 2) - naca_chord_v(j,i);
        % Compute the shear flow distribution
        q_s_v(j,i) = ((-V_y_v(i) / I_xx_v_shear(i)) * sum(t_s_v .* y_bar_v_s(1:j,i) .* L_v(1:j,i)));
    end
end

% Calculate the shear force at the cut point for the vertical tail
for i = 1:length(z_v)
    % Compute the scaled chord length
    C = C_r_v - (C_r_v * ((2 * z_v(i) / span_v) * (1 - t_v)));
    % Compute the moment arm of the external load
    eps_v(i) = abs((C * 0.25) - x_scaled_v(1,i));
    % Compute the area using the trapezoidal rule
    A_v(i) = trapz(double(x_scaled_v(:,i)), double(y_scaled_v(:,i)));
    % Compute the shear force at the cut point
    q_s_0_v(i) = (V_y_v(i)*eps_v(i) - sum(P_v(:,i).*q_s_v(:,i).*L_v(:,i))) / (2*A_v(i));
end

% Calculate the shear flow distribution for the vertical tail
for i = 1:length(z_v)
    q_s_v(:,i) = q_s_v(:,i) + q_s_0_v(i);
end

tau_s_h = zeros(length(naca0012(:,1)), length(z_v));
tau_s_v = zeros(length(naca0012(:,1)), length(z_v));
tau_spar_h = zeros(length(z_v));
tau_spar_v = zeros(length(z_v));

% Calculate the shear stress distribution for the horizontal tail skin
for i = 1:length(y_h)
    tau_s_h(:,i) = q_s_h(:,i) ./ t_s_h;
end

% Calculate the shear stress in the horizontal tail spar web
for i = 1:length(y_h)
    tau_spar_h(i) = q_s_h(end,i) / t_spar_web_h(i);
end

% Calculate the shear stress distribution for the vertical tail skin
for i = 1:length(z_v)
    tau_s_v(:,i) = q_s_v(:,i) ./ t_s_v;
end

% Calculate the shear stress in the vertical tail spar web
for i = 1:length(z_v)
    tau_spar_v(i) = q_s_v(end,i) / t_spar_web_v(i);
end


% Calculate the maximum shear stress for the horizontal tail
max_shear_stress_h = zeros(1, length(y_h));
for i = 1:length(y_h)
    max_shear_stress_h(i) = max(abs(tau_s_h(:,i)));
end

% Plot the maximum shear stress for the horizontal tail
figure;
plot(y_h, max_shear_stress_h);
title('Maximum Shear Stress vs Spanwise Location for Horizontal Tail');
xlabel('Spanwise Location (y_h)');
ylabel('Maximum Shear Stress');

% Calculate the maximum shear stress for the vertical tail
max_shear_stress_v = zeros(1, length(z_v));
for i = 1:length(z_v)
    max_shear_stress_v(i) = max(abs(tau_s_v(:,i)));
end

% Plot the maximum shear stress for the vertical tail
figure;
plot(z_v, max_shear_stress_v);
title('Maximum Shear Stress vs Spanwise Location for Vertical Tail');
xlabel('Spanwise Location (z_v)');
ylabel('Maximum Shear Stress');

% Plot the maximum shear stress in the horizontal tail spar web
figure;
plot(y_h, tau_spar_h);
title('Shear Stress in Horizontal Tail Spar Web vs Spanwise Location');
xlabel('Spanwise Location (y_h)');
ylabel('Shear Stress');

% Plot the maximum shear stress in the vertical tail spar web
figure;
plot(z_v, tau_spar_v);
title('Shear Stress in Vertical Tail Spar Web vs Spanwise Location');
xlabel('Spanwise Location (z_v)');
ylabel('Shear Stress');

ult_tau_h = max(max_shear_stress_h);
ult_tau_v = max(max_shear_stress_v);
ult_tau_h_spar = max(tau_spar_h);
ult_tau_v_spar = max(tau_spar_v);

fprintf('The ultimate shear stress in the horizontal tail skin is %.2f MPa\n', ult_tau_h/1e6);
fprintf('The ultimate shear stress in the vertical tail skin is %.2f MPa\n', ult_tau_v/1e6);

% Compute the max axial stress
h_spar_web_v = h_spar_web_v(:);
h_spar_web_h = h_spar_web_h(:);

max_axial_stress_v_tens_spar = (M_y_v .* (h_spar_web_v/2)) ./ I_xx_v_bending;
max_axial_stress_v_comp_spar = (M_y_v .* (h_spar_web_v/2)) ./ I_xx_v_bending;

max_axial_stress_h_tens_spar = (M_z_h .* (h_spar_web_h/2)) ./ I_xx_h_bending;
max_axial_stress_h_comp_spar = (M_z_h .* (h_spar_web_h/2)) ./ I_xx_h_bending;

ult_ax_str_comp_spar_h = max(max_axial_stress_h_comp_spar);
ult_ax_str_tens_spar_h = max(max_axial_stress_h_tens_spar);
ult_ax_str_comp_spar_v = max(max_axial_stress_v_comp_spar);
ult_ax_str_tens_spar_v = max(max_axial_stress_v_tens_spar);

% Plot the maximum axial stress in the horizontal tail spar 
figure;
plot(y_h, max_axial_stress_h_comp_spar, 'r', y_h, max_axial_stress_h_tens_spar, 'b');
title('Max Axial Stress in Horizontal Tail Spar vs Spanwise Location');
xlabel('Spanwise Location (y_h)');
ylabel('Max Axial Stress');
legend('Max Compressive Stress - Horizontal Tail Spar', 'Max Tensile Stress - Horizontal Tail Spar');

% Plot the maximum axial stress in the vertical tail spar
figure;
plot(z_v, max_axial_stress_v_comp_spar, 'm', z_v, max_axial_stress_v_tens_spar, 'c');
title('Max Axial Stress in Vertical Tail Spar vs Spanwise Location');
xlabel('Spanwise Location (z_v)');
ylabel('Max Axial Stress');
legend('Max Compressive Stress - Vertical Tail Spar', 'Max Tensile Stress - Vertical Tail Spar');

%% Buckling analysis

% Compute the critical buckling stress for the horizontal tail skin between the spar and the leading edge
% Assume that the skin is clamped at the leading edge and clamped at the spar
% Assume loaded edges are simply supported by ribs
% Assume that the skin is rectangular with length a equal to the wing span and width b equal to the chord length at the mid point between the spar and leading edge

a_h = span_h;
pos = ceil(length(y_h)/2);
target = 0.7*C_r_h;

% Compute the absolute differences
differences = abs(naca0012_h(:,1) - target);
pos_x_1 = differences == min(differences);
b_h = x_scaled_h(pos_x_1,pos) - min(x_scaled_h(:,pos));
a_b_h = a_h/b_h;
fprintf('The buckling formula aspect ratio of the horizontal tail skin panel is %.2f\n', a_b_h);

a_v = span_v;
pos = ceil(length(z_v)/2);
target = 0.7*C_r_v;

% Compute the absolute differences
differences = abs(naca0012_v(:,1) - target);
pos_x_1 = differences == min(differences);
b_v = x_scaled_v(pos_x_1,pos) - min(x_scaled_v(:,pos));
a_b_v = a_v/b_v;
fprintf('The buckling formula aspect ratio of the vertical tail skin panel is %.2f\n', a_b_v);

% Initialize the loop control variable
continue_loop = true;

while continue_loop
    % Get the input from the user
    k_c_h = input('Enter the critical bending buckling stress coefficient for the horizontal tail skin panel: ');
    k_s_h = input('Enter the critical shear buckling stress coefficient for the horizontal tail skin panel: ');
    k_c_v = input('Enter the critical bending buckling stress coefficient for the vertical tail skin panel: ');
    k_s_v = input('Enter the critical shear buckling stress coefficient for the vertical tail skin panel: ');

    t_s_h = input('Enter the thickness of the horizontal tail skin panel: ');
    t_s_v = input('Enter the thickness of the vertical tail skin panel: ');

    E_skin = input('Enter the Young''s modulus of the skin material: ');
    nu_skin = input('Enter the Poisson''s ratio of the skin material: ');

    % Compute the critical buckling stress for the horizontal tail skin
    sigma_crit_buckling_h = ((k_c_h * (pi^2) * E_skin) / (12 * (1 - (nu_skin^2)))) * ((t_s_h/b_h)^2);
    tau_crit_buckling_h = ((k_s_h * (pi^2) * E_skin) / (12 * (1 - (nu_skin^2)))) * ((t_s_h/b_h)^2);

    % Compute the critical buckling stress for the vertical tail skin
    sigma_crit_buckling_v = ((k_c_v * (pi^2) * E_skin) / (12 * (1 - (nu_skin^2)))) * ((t_s_v/b_v)^2);
    tau_crit_buckling_v = ((k_s_v * (pi^2) * E_skin) / (12 * (1 - (nu_skin^2)))) * ((t_s_v/b_v)^2);

    % Compare the critical buckling stress to the ultimate stress
    if ult_sigma_h < sigma_crit_buckling_h
        fprintf('The horizontal tail skin is safe against buckling in bending\n');
    else
        fprintf('The horizontal tail skin is not safe against buckling in bending\n');
    end

    if ult_sigma_v < sigma_crit_buckling_v
        fprintf('The vertical tail skin is safe against buckling in bending\n');
    else
        fprintf('The vertical tail skin is not safe against buckling in bending\n');
    end

    if ult_tau_h < tau_crit_buckling_h
        fprintf('The horizontal tail skin is safe against buckling in shear\n');
    else
        fprintf('The horizontal tail skin is not safe against buckling in shear\n');
    end

    if ult_tau_v < tau_crit_buckling_v
        fprintf('The vertical tail skin is safe against buckling in shear\n');
    else
        fprintf('The vertical tail skin is not safe against buckling in shear\n');
    end

    % Ask the user if they want to continue
    user_input = input('Do you want to continue? (yes/no): ', 's');
    if strcmpi(user_input, 'no')
        continue_loop = false;
    end
end

% Estimate total volume of concept

% Horizontal tail
% Estimate the volume of the horizontal tail skin
V_skin_h = trapz(y_h, t_s_h * trapz(x_scaled_h, z_scaled_h));

% Estimate the volume of the horizontal tail spar
V_spar_h = trapz(y_h, h_spar_web_h.* t_spar_web_h) + trapz(y_h, 2 * w_spar_h .* t_spar_cap_h);

% Vertical tail
% Estimate the volume of the vertical tail skin
V_skin_v = trapz(z_v, t_s_v * trapz(x_scaled_v, y_scaled_v));

% Estimate the volume of the vertical tail spar
V_spar_v = trapz(z_v, h_spar_web_v.* t_spar_web_v) + trapz(z_v, 2 * w_spar_v .* t_spar_cap_v);

Vol_h = V_skin_h + V_spar_h;
Vol_v = V_skin_v + V_spar_v;
Vol_total = Vol_h + Vol_v;

fprintf('The total volume of the horizontal tail is %.2f m^3\n', Vol_h);
fprintf('The total volume of the vertical tail is %.2f m^3\n', Vol_v);
fprintf('The total volume of the concept is %.2f m^3\n', Vol_total);