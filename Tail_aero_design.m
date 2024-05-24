% Tail aerodynamic design

%% Horizontal tail design

% Tail volume coefficient

V_barH = 0.7; % Horizontal tail volume coefficient

fprintf('The horizontal tail volume coefficient is %f\n', V_barH);

% Optimum tail arm length

L = 3; % Total length of the aircraft (m)
t = 0.4; % Taper ratio
C = 0.32457; % mean aerodynamic chord
S = 0.7374; % Wing planform area (m^2)
D_fus = 0.4; % Fuselage diameter (m)
K_c = 1.1; % Correction factor for fuselage diameter

l_opt = K_c*sqrt((4*C*S*V_barH)/(pi*D_fus)); % Optimum tail arm length

fprintf('The optimum tail arm length is %f\n', l_opt);

l_opt = 0.6*L; % actual tail arm length due to design requirements - based on average l/L ratio of 0.6 for selected configuration

fprintf('The actual tail arm length is %f\n', l_opt);

% Calculate horizontal tail area

S_h = (V_barH*S*C)/l_opt; % Tail area (m^2)

fprintf('The horizontal tail area is %f\n', S_h);

% Calculate wing fuselage pitching moment coefficient

% Choose wing airfoil type - user input

Airfoil_w = input('Enter the wing airfoil type: ', 's');

C_l_c = 0.42317; % Aircraft cruise lift coefficient
V_c = 85; % Cruise speed (m/s)
rho_c = 0.7779; % Air density at cruise - 15,000ft (kg/m^3)
mu_c = 1.642e-5; % Dynamic viscosity at cruise - 15,000ft (Ns/m^2)
Re_c = (rho_c*V_c*C)/mu_c; % Reynolds number at cruise

fprintf('The Reynolds number at cruise is %f\n', Re_c);

if Airfoil_w == 'NACA 4415'
    C_m_w = -0.085; % Wing pitching moment coefficient at cruise - alpha = 0 degrees, Ncrit = 9, Reynolds number = Re_c, airfoil = NACA 4415
elseif Airfoil_w == 'NACA 4418'
    C_m_w = -0.097; % Wing pitching moment coefficient at cruise - alpha = -0.08 degrees, Ncrit = 9, Reynolds number = Re_c, airfoil = NACA 4418
else
    fprintf('Airfoil type not recognised\n');
end

AR_w = 7; % Wing aspect ratio 
span_w = sqrt(AR_w*S); % Wing span (m) 
Delta_w = 0; % wing sweep angle (typical value - (%%unknown at present%%))
alpha_t_w = -3; % Wing twist angle 

C_m_wf = C_m_w*((AR_w*cos(Delta_w)^2)/(AR_w + 2*cos(Delta_w))); % Wing fuselage pitching moment coefficient

fprintf('The wing fuselage pitching moment coefficient is %f\n', C_m_wf);

% Calculate desired horizontal tail lift coefficient

n_h = 0.9; % Tail efficiency factor

F_cg_ac = 0.05; % Distance factor between center of gravity and aerodynamic center
cg__ac = F_cg_ac*C; % Distance between wing aerodynamic center and centre of gravity
h__h_0 = -cg__ac/C; % Distance between wing aerodynamic center and horizontal tail aerodynamic center in terms of mean aerodynamic chord

pos_ac_wf = L-l_opt; % Position of wing aerodynamic center from fusealge nose (m)
pos_cg = L-l_opt+h__h_0*C; % Position of center of gravity from fuselage nose (m)

fprintf('The position of the wing aerodynamic center from the fuselage nose is %f\n', pos_ac_wf);
fprintf('The position of the center of gravity from the fuselage nose is %f\n', pos_cg);

C_l_h = (C_m_wf+(C_l_c*(h__h_0)))/(n_h*V_barH); % Desired horizontal tail lift coefficient

fprintf('The desired horizontal tail lift coefficient is %f\n', C_l_h);

% Select horizontal tail airfoil type 

C_l_h_i = 0.9*C_l_h; % Ideal horizontal tail lift coefficient
C_m_wf_max = -0.1; % Maximum wing fuselage pitching moment coefficient 
C_l_max = 1.45125; % Maximum lift coefficient


C_l_h_max = (C_m_wf_max+(C_l_c*(h__h_0)))/(n_h*V_barH); % Maximum horizontal tail lift coefficient

% Select horizontal tail airfoil type

fprintf('The selected horizontal tail airfoil type is NACA 0012\n'); % Chosen using qualitative criteria

% Select horizontal tail sweep angle

Delta_h = Delta_w; % Horizontal tail sweep angle (degrees) - chosen using qualitative criteria (%%unknown at present%%)

fprintf('The horizontal tail sweep angle is %f\n', Delta_h);

% Select horizontal tail dihedreal angle

Gamma_h = 7; % Horizontal tail dihedral angle (degrees) - chosen using qualitative criteria (%%unknown at present%%)

fprintf('The horizontal tail dihedral angle is %f\n', Gamma_h);

% Select horizontal tail aspect ratio

AR_h = (2/3)*AR_w; % Horizontal tail aspect ratio

fprintf('The horizontal tail aspect ratio is %f\n', AR_h);

% Select horizontal tail taper ratio

t_h = 0.3; % Horizontal tail taper ratio - chosen using qualitative criteria (%%unknown at present%%)

fprintf('The horizontal tail taper ratio is %f\n', t_h);

% Calculate the horizontal tail lift curve slope

C_l_alpha_h = 0.1067; % NACA 0012 lift curve slope (1/degrees)

% convert to per radian

C_l_alpha_h = C_l_alpha_h*(180/pi); % NACA 0012 lift curve slope (1/radian)

C_L_alpha_h = C_l_alpha_h/(1-(C_l_alpha_h/(pi*AR_h))); % Horizontal tail lift curve slope (1/radian)

fprintf('The horizontal tail lift curve slope is %f\n', C_L_alpha_h);

% Calculate horizontal tail angle of attack at cruise

alpha_h = (C_l_h_i/C_L_alpha_h); % Horizontal tail angle of attack at cruise (radians)

fprintf('The theoretical horizontal tail angle of attack at cruise is %f\n', alpha_h);

% Calculate downwash angle at horizontal tail

if Airfoil_w == 'NACA 4415'
    C_l_alpha_w = 0.1067; % NACA 4415 lift curve slope (1/radians)
elseif Airfoil_w == 'NACA 4418'
    C_l_alpha_w = 0.1059; % NACA 4418 lift curve slope (1/radians)
else
    fprintf('Airfoil type not recognised\n');
end

C_L_alpha_w = C_l_alpha_w/(1-(C_l_alpha_w/(pi*AR_w))); % Wing lift curve slope (1/radian)

alpha_w_c = 1.68; % Wing angle of attack at cruise (radians)
C_l_c_w = 0.44544; % Wing lift coefficient at cruise

epsilon = ((2*C_l_c_w)/(pi*AR_w))+(((2*C_L_alpha_w)/(pi*AR_w))*alpha_w_c); % Downwash angle at horizontal tail (radians)

fprintf('The downwash angle at the horizontal tail is %f\n', epsilon);

% Calculate the horizontal tail setting angle

i_t_h = alpha_h+epsilon; % Horizontal tail setting angle (radians)

fprintf('The initial value of the horizontal tail setting angle is %f\n', i_t_h);

% Calculate horizontal tail span, root chord, tip chord and mean aerodynamic chord

% Define the symbolic variables
syms span_h C_t_h C_r_h C_h

% Define the equations
eqn1 = AR_h == span_h/C_h;
eqn2 = t_h == C_t_h/C_r_h;
eqn3 = C_h == (2/3)*C_r_h*((1+t_h+t_h^2)/(1+t_h));
eqn4 = S_h == span_h*C_h;

% Solve the equations
sol = solve([eqn1, eqn2, eqn3, eqn4], [span_h, C_t_h, C_r_h, C_h]);

% Extract the solutions
span_h_sol = sol.span_h;
C_t_h_sol = sol.C_t_h;
C_r_h_sol = sol.C_r_h;
C_h_sol = sol.C_h;

% Filter out negative solutions
span_h = span_h_sol(span_h_sol > 0);
C_t_h = C_t_h_sol(C_t_h_sol > 0);
C_r_h = C_r_h_sol(C_r_h_sol > 0);
C_h = C_h_sol(C_h_sol > 0);

% Display the solutions
if isempty(span_h)
    fprintf('No positive solution for span_h\n');
else
    fprintf('The horizontal tail span is %f\n', span_h);
end

if isempty(C_t_h)
    fprintf('No positive solution for C_t_h\n');
else
    fprintf('The horizontal tail tip chord is %f\n', C_t_h);
end

if isempty(C_r_h)
    fprintf('No positive solution for C_r_h\n');
else
    fprintf('The horizontal tail root chord is %f\n', C_r_h);
end

if isempty(C_h)
    fprintf('No positive solution for C_h\n');
else
    fprintf('The horizontal tail mean aerodynamic chord is %f\n', C_h);
end

%% Lifting line theory

N = 9; % (number of segments - 1)
% Horizontal tail twist angle is zero (deg)
alpha_0_h = 0; % Horizontal tail zero lift angle of attack (deg)


theta = pi/(2*N):pi/(2*N):pi/2;
alpha = alpha_h*(180/pi); % segments angle of attack (deg)
% segment’s angle of attack
z = (span_h/2)*cos(theta);
c = C_r_h * (1 - (1-t_h)*cos(theta)); % Mean Aerodynamics Chord at each segment (m)
mu = c * C_L_alpha_h / (4 * span_h);
LHS = mu .* (alpha-alpha_0_h)/57.3; % Left Hand Side
% Solving N equations to find coefficients A(i):
for i=1:N
    for j=1:N
        B(i,j) = sin((2*j-1) * theta(i)) * (1 + (mu(i) * (2*j-1)) /sin(theta(i)));
    end
end
A=B\transpose(LHS);
for i = 1:N
    sum1(i) = 0;
    sum2(i) = 0;
    for j=1:N
        sum1(i) = sum1(i) + (2*j-1) * A(j)*sin((2*j-1)*theta(i));
        sum2(i) = sum2(i) + A(j)*sin((2*j-1)*theta(i));
    end
end
CL = 4*span_h*sum2 ./ c;
CL1=[0 CL(1) CL(2) CL(3) CL(4) CL(5) CL(6) CL(7) CL(8) CL(9)];
y_s=[span_h/2 z(1) z(2) z(3) z(4) z(5) z(6) z(7) z(8) z(9)];
plot(y_s,CL1,'-o', 'LineWidth', 2)
grid
set(gca, 'FontSize', 30, 'LineWidth', 2);
xlabel('Semi-span location (m)', 'FontSize', 30)
ylabel ('Lift coefficient', 'FontSize', 30)
CL_horizontal_tail = pi * AR_h * A(1);

fprintf('The lifting theory horizontal tail lift coefficient is %f\n', CL_horizontal_tail);
fprintf('The required horizontal tail lift coefficient is %f\n', C_l_h);

if CL_horizontal_tail ~= C_l_h
    fprintf('Therefore, the tail setting angle must be adjusted to achieve the required lift coefficient ...')
end

% Iteratively adjust horizontal tail effective angle of attack to match required lift coefficient

% Define the initial guesses for alpha_h
alpha_h_0 = rad2deg(alpha_h);
alpha_h_1 = 10;

% Define the tolerance for the difference between C_l_h and CL_horizontal_tail
tolerance = abs(0.1 * C_l_h);

% Define the initial values for C_l_h and CL_horizontal_tail
C_l_h_0 = C_l_h;
CL_horizontal_tail_0 = CL_horizontal_tail;

% Calculate the new values for C_l_h and CL_horizontal_tail
N = 9; % (number of segments - 1)
% Horizontal tail twist angle is zero (deg)
alpha_0_h = 0; % Horizontal tail zero lift angle of attack (deg)

theta = pi/(2*N):pi/(2*N):pi/2;
alpha = alpha_h_1;
% segment’s angle of attack
z = (span_h/2)*cos(theta);
c = C_r_h * (1 - (1-t_h)*cos(theta)); % Mean Aerodynamics Chord at each segment (m)
mu = c * C_L_alpha_h / (4 * span_h);
LHS = mu .* (alpha-alpha_0_h)/57.3; % Left Hand Side
% Solving N equations to find coefficients A(i):
for i=1:N
    for j=1:N
        B(i,j) = sin((2*j-1) * theta(i)) * (1 + (mu(i) * (2*j-1)) /sin(theta(i)));
    end
end
A=B\transpose(LHS);

CL_horizontal_tail_1 = pi * AR_h * A(1); % initial value of CL_horizontal_tail_1

% Iteratively adjust alpha_h
while true
    % Calculate i_t_h using the secant method
    alpha_h = alpha_h_1 - (CL_horizontal_tail - C_l_h) * (alpha_h_1 - alpha_h_0) / (CL_horizontal_tail_1 - CL_horizontal_tail_0);

    % Calculate the new values for C_l_h and CL_horizontal_tail
    N = 9; % (number of segments - 1)
    % Horizontal tail twist angle is zero (deg)
    alpha_0_h = 0; % Horizontal tail zero lift angle of attack (deg)
    
    theta = pi/(2*N):pi/(2*N):pi/2;
    alpha = alpha_h;
    % segment’s angle of attack
    z = (span_h/2)*cos(theta);
    c = C_r_h * (1 - (1-t_h)*cos(theta)); % Mean Aerodynamics Chord at each segment (m)
    mu = c * C_L_alpha_h / (4 * span_h);
    LHS = mu .* (alpha-alpha_0_h)/57.3; % Left Hand Side
    % Solving N equations to find coefficients A(i):
    for i=1:N
        for j=1:N
            B(i,j) = sin((2*j-1) * theta(i)) * (1 + (mu(i) * (2*j-1)) /sin(theta(i)));
        end
    end
    A=B\transpose(LHS);

    CL_horizontal_tail = pi * AR_h * A(1);

    % Check if the difference between C_l_h and CL_horizontal_tail is within the tolerance
    if abs(C_l_h - CL_horizontal_tail) <= tolerance
        break;
    end

    % Update alpha_h_0, alpha_h_1, CL_horizontal_tail_0, and CL_horizontal_tail_1 for the next iteration
    alpha_h_0 = alpha_h_1;
    alpha_h_1 = alpha_h;
    CL_horizontal_tail_0 = CL_horizontal_tail_1;
    CL_horizontal_tail_1 = CL_horizontal_tail;
end

% Calculate the horizontal tail setting angle, accounting for downwash

i_t_h = alpha_h+rad2deg(epsilon); % Horizontal tail setting angle (degrees)

% Display the final values of i_t_h, and C_l_h
fprintf('The final value of the horizontal tail setting angle is %f degrees \n', i_t_h);
fprintf('The final value of the horizontal tail lift coefficient is %f\n', C_l_h);

%% Horizontal tail - final check

% Longitudinal stability derivative

C_m_alpha = C_L_alpha_w*(h__h_0) - C_L_alpha_h*n_h*(S_h/S)*((1/C)-pos_cg/C)*(1-((2*C_L_alpha_w)/(pi*AR_w))) ; % Longitudinal stability derivative

if C_m_alpha < 0
    fprintf('The aircraft is longitudinally stable\n');
else
    fprintf('The aircraft is longitudinally unstable\n');
end


%% Vertical tail design

fprintf('... \n');
fprintf('Vertical tail design \n');
fprintf('... \n');

% Tail volume coefficient

V_barV = 0.05; % Vertical tail volume coefficient

fprintf('The vertical tail volume coefficient is %f\n', V_barV);

% Vertical tail arm length

l_v = 0.99; % tail arm length (m) based on manufacturing requirements

% Calculate vertical tail area

S_v = (V_barV*S*span_w)/l_v; % Tail area (m^2)

fprintf('The vertical tail area is %f\n', S_v);

% Select vertical tail airfoil type

fprintf('The selected vertical tail airfoil type is NACA 0012\n'); % Chosen using qualitative criteria

% Select vertical tail aspect ratio

AR_v = 1.5; % Vertical tail aspect ratio - assummed for now (%%unknown at present%%)

fprintf('The vertical tail aspect ratio is %f\n', AR_v);

% Select vertical tail taper ratio

t_v = 0.8; % Vertical tail taper ratio - chosen using qualitative criteria (%%unknown at present%%)

fprintf('The vertical tail taper ratio is %f\n', t_v);

% Select vertical tail incidence angle

i_v = 1.5; % Vertical tail incidence angle (degrees) - chosen using qualitative criteria (%%unknown at present%%)

fprintf('The vertical tail incidence angle is %f\n', i_v);

% Select vertical tail sweep angle

Delta_v = Delta_w; % Vertical tail sweep angle (degrees) - chosen using qualitative criteria (%%unknown at present%%)

% Select vertical tail dihedral angle

Gamma_v = 0; % Vertical tail dihedral angle (degrees) 

% Calculate vertical tail span, root chord, tip chord and mean aerodynamic chord

% Define the symbolic variables
syms span_v C_t_v C_r_v C_v

% Define the equations
eqn1 = AR_v == span_v/C_v;
eqn2 = t_v == C_t_v/C_r_v;
eqn3 = C_v == (2/3)*C_r_v*((1+t_v+t_v^2)/(1+t_v));
eqn4 = S_v == span_v*C_v;

% Solve the equations
sol = solve([eqn1, eqn2, eqn3, eqn4], [span_v, C_t_v, C_r_v, C_v]);

% Extract the solutions
span_v_sol = sol.span_v;
C_t_v_sol = sol.C_t_v;
C_r_v_sol = sol.C_r_v;
C_v_sol = sol.C_v;

% Filter out negative solutions
span_v = span_v_sol(span_v_sol > 0);
C_t_v = C_t_v_sol(C_t_v_sol > 0);
C_r_v = C_r_v_sol(C_r_v_sol > 0);
C_v = C_v_sol(C_v_sol > 0);

% Display the solutions
if isempty(span_v)
    fprintf('No positive solution for span_v\n');
else
    fprintf('The vertical tail span is %f\n', span_v);
end

if isempty(C_t_v)
    fprintf('No positive solution for C_t_v\n');
else
    fprintf('The vertical tail tip chord is %f\n', C_t_v);
end

if isempty(C_r_v)
    fprintf('No positive solution for C_r_v\n');
else
    fprintf('The vertical tail root chord is %f\n', C_r_v);
end

if isempty(C_v)
    fprintf('No positive solution for C_v\n');
else
    fprintf('The vertical tail mean aerodynamic chord is %f\n', C_v);
end

% Check spin recovery

% User inputs - position of horizontal tail in z-direction and vertical tail arm length

pos_h_z = input('Enter the position of the horizontal tail in the z-direction (metres from neutral axis): ');
bin = input('Is the vertical tail arm length the same as the horizontal tail arm length? (y/n): ', 's');
if bin == 'y'
    tail_arm_v = l_opt;
else
    tail_arm_v = input('Enter the vertical tail arm length (metres): ');
end

Spin_analysis(Airfoil_w, l_opt, C_r_v, C_r_h, C_t_v, span_v, S_v, pos_h_z, tail_arm_v)