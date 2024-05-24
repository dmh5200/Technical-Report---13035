% Spin_analysis(Airfoil_w,l_opt, C_r_v, C_r_h, C_t_v, span_v, S_v pos_h_z, tail_arm_v) performs spin recovery analysis.
% Airfoil_w is the wing airfoil type.
% l_opt is the length of the tail arm for the horizontal tail.
% C_r_v is the root chord of the vertical tail.
% C_r_h is the root chord of the horizontal tail.
% C_t_v is the tip chord of the vertical tail.
% Span_v is the span of the vertical tail.
% S_v is the total area of the vertical tail.
% pos_h_z is the z-coordinate of the root of the horizontal tail.
% tail_arm_v is the length of the tail arm for the vertical tail.
% If pos_h_z is not provided, it defaults to 0.
% If tail_arm_v is not provided, it defaults to the value of l_opt.

function Spin_analysis(Airfoil_w, l_opt, C_r_v, C_r_h, C_t_v, span_v, S_v, pos_h_z, tail_arm_v)

% Check if the inputs are provided, if not, use default values
if nargin < 1 || isempty(l_opt)
    error('Input argument "l_opt" is required.');
end
if nargin < 2 || isempty(C_r_v)
    error('Input argument "C_r_v" is required.');
end
if nargin < 3 || isempty(C_r_h)
    error('Input argument "C_r_h" is required.');
end
if nargin < 4 || isempty(C_t_v)
    error('Input argument "C_t_v" is required.');
end
if nargin < 5 || isempty(span_v)
    error('Input argument "span_v" is required.');
end
if nargin < 6 || isempty(S_v)
    error('Input argument "S_v" is required.');
end
if nargin < 7 || isempty(pos_h_z)
    pos_h_z = 0; % Default value for pos_h_z
end
if nargin < 8 || isempty(tail_arm_v)
    tail_arm_v = l_opt; % Default value for tail_arm_v
end

% Check spin recovery

% Calculate the position of the aerodynamic center
AC_v = 0.25 * C_r_v;  % For the vertical tail
AC_h = 0.25 * C_r_h;  % For the horizontal tail

% Adjust the position of the tails
x_v = tail_arm_v - AC_v;  % root x-coordinate of the vertical tail
x_h = l_opt - AC_h;  % root x-coordinate of the horizontal tail

% Load the NACA 0012 airfoil data from a CSV file
if Airfoil_w == 'NACA 4415'
    naca0012 = csvread('n0012-il.csv');
elseif Airfoil_w == 'NACA 4418'
    naca0012 = csvread('n0012-il(w4418).csv');
end

% Convert from mm to m
naca0012 = naca0012 / 1000;
% Scale the airfoil data to the size of the horizontal tail
naca0012_h = double(naca0012 .* (C_r_h/max(naca0012(:,1))));
% Adjust position to match ac position
naca0012_h(:,1) = naca0012_h(:,1) + x_h;
% Adjust height of horizontal tail
naca0012_h(:,2) = naca0012_h(:,2) + pos_h_z;

% Calculate the coordinates of the vertical tail
x_v_coords = [x_v, x_v + C_r_v - C_t_v, x_v + C_r_v, x_v + C_r_v, x_v];
z_v_coords = [0, span_v, span_v, 0, 0];

% Plot the tail planforms
figure;
hold on;
p1 = plot(naca0012_h(:,1), naca0012_h(:,2), 'b');  % Plot the horizontal tail
p2 = plot(x_v_coords, z_v_coords, 'r');  % Plot the vertical tail

% maximum and minimum x and z coordinates of the horizontal tail
[naca0012_h_max_x, i_max] = max(naca0012_h(:,1));
[naca0012_h_min_x, i_min] = min(naca0012_h(:,1));

% Calculate the coordinates of the dashed lines
x_le = min(naca0012_h(:,1));  % x-coordinate of the leading edge
x_le_line = linspace(x_le, x_v+C_r_v+0.2, 100);  % x-coordinates from the leading edge to the max x-coordinate of the horizontal tail
y_le = (tan(deg2rad(60)) * (x_le_line - x_le)) + naca0012_h(i_min,2);  % y-coordinate of the leading edge line

x_te = max(naca0012_h(:,1));  % x-coordinate of the trailing edge

% Calculate the x and y coordinates for the trailing edge line
x_te_line = linspace(x_te-0.0001, max(naca0012_h(:,1))+0.2, 100);  % x-coordinates from the trailing edge to the max x-coordinate of the horizontal tail
y_te = (tan(deg2rad(30)) * (x_te_line - x_te)) + naca0012_h(i_max,2);  % y-coordinate of the trailing edge line

% Add the dashed lines to the plot
p3 = plot(x_le_line, y_le, '--k');  % Leading edge line
p4 = plot(x_te_line, y_te, '--k');  % Trailing edge line

hold off;

% Label the axes
xlabel('x (distance from wing aerodynamic centre, m)');
ylabel('z (distance from neutral line, m)');
title('Planform view of the tail in the x-z plane');

% Add legend
legend([p1, p2, p3], {'Vertical Airfoil Planform', 'Horizontal Tail Root Section (NACA 0012)', 'Horizontal Tail Wake'});

% Convert inputs to type double

x_le_line = double(x_le_line);
y_le = double(y_le);
x_te_line = double(x_te_line);
y_te = double(y_te);
x_v_coords = double(x_v_coords);
z_v_coords = double(z_v_coords);

% Calculate the slopes and y-intercepts of the lines
slope_le = (y_le(end) - y_le(1)) / (x_le_line(end) - x_le_line(1));
slope_te = (y_te(end) - y_te(1)) / (x_te_line(end) - x_te_line(1));
intercept_le = y_le(1) - slope_le * x_le_line(1);
intercept_te = y_te(1) - slope_te * x_te_line(1);

% Define the extended lines
x_le_line_1 = linspace(min(x_le_line), max(x_le_line) + 20, 100); % Extend the x values by 20 units on the positive side
x_te_line_1 = linspace(min(x_te_line), max(x_te_line) + 20, 100); % Extend the x values by 20 units on the positive side
y_le_1 = slope_le * x_le_line_1 + intercept_le; % Calculate the corresponding y values
y_te_1 = slope_te * x_te_line_1 + intercept_te; % Calculate the corresponding y values

% Find intersection points of dashed lines with the vertical tail
[x_intersect_le_v, z_intersect_le_v] = polyxpoly(x_le_line, y_le, x_v_coords, z_v_coords);
[x_intersect_te_v, z_intersect_te_v] = polyxpoly(x_te_line, y_te, x_v_coords, z_v_coords);

% Find coordinates of the vertical tail that are enclosed by the dashed lines
in = inpolygon(x_v_coords, z_v_coords, [x_le_line_1(1); x_le_line_1(end); x_te_line_1(end); x_te_line_1(1)], [y_le_1(1); y_le_1(end); y_te_1(end); y_te_1(1)]);

% Find intersection points of dashed lines with the horizontal airfoil profile
[x_intersect_le_h, z_intersect_le_h] = polyxpoly(x_le_line, y_le, naca0012_h(:,1), naca0012_h(:,2));
[x_intersect_te_h, z_intersect_te_h] = polyxpoly(x_te_line, y_te, naca0012_h(:,1), naca0012_h(:,2));
x_intersect_h = [x_intersect_le_h; x_intersect_te_h];
z_intersect_h = [z_intersect_le_h; z_intersect_te_h];

% 1. Find the intersection points of the horizontal tail profile with the vertical tail
[x_intersect_h_v, z_intersect_h_v] = polyxpoly(naca0012_h(:,1), naca0012_h(:,2), x_v_coords, z_v_coords);

% 2. Check if any of the intersection points of the dashed lines and horizontal tail are outside of the vertical tail profile
in_h_v = inpolygon(x_intersect_h, z_intersect_h, x_v_coords, z_v_coords);

% 3. If they are, use the intersection points of the horizontal tail and vertical tail with the highest y-values to replace them
if any(~in_h_v)
    if ~isempty(z_intersect_h_v)
        % Find the intersection points with the highest y-values
        [~, idx_max] = max(z_intersect_h_v);
        x_intersect_h_max = x_intersect_h_v(idx_max);
        z_intersect_h_max = z_intersect_h_v(idx_max);

        % Replace the intersection points that are outside of the vertical tail profile
        x_intersect_h(~in_h_v) = x_intersect_h_max;
        z_intersect_h(~in_h_v) = z_intersect_h_max;
    else
        % Replace the intersection points that are outside of the vertical tail profile
        x_intersect_h = x_intersect_h(in_h_v);
        z_intersect_h = z_intersect_h(in_h_v);
    end
end

% Find intersection points that are enclosed by the dashed lines
in_h = inpolygon(x_intersect_h, z_intersect_h, [x_le_line_1(1); x_le_line_1(end); x_te_line_1(end); x_te_line_1(1)], [y_le_1(1); y_le_1(end); y_te_1(end); y_te_1(1)]);

% Convert to column vectors
x_v_coords = x_v_coords(:);
z_v_coords = z_v_coords(:);
% Concatenate the intersection points and the enclosed coordinates
x_vertices = [x_intersect_le_v(:); x_v_coords(in); x_intersect_h(in_h); x_intersect_te_v(:)];
z_vertices = [z_intersect_le_v(:); z_v_coords(in); z_intersect_h(in_h); z_intersect_te_v(:)];


% Sort the vertices to represent the perimeter enclosing the largest area
K = convhull(x_vertices, z_vertices);
x_vertices = x_vertices(K);
z_vertices = z_vertices(K);

% Plot the enclosed area

% Create polygon objects for the enclosed area and the horizontal airfoil section
enclosed_area = polyshape(x_vertices, z_vertices);
horizontal_airfoil = polyshape(naca0012_h(:,1), naca0012_h(:,2));

% Subtract the overlapping area
resulting_area = subtract(enclosed_area, horizontal_airfoil);

% Calculate the area of the resulting polygon
area_resulting = area(resulting_area);

% Calculate the area as a percentage of the vertical tail area
percentage = (area_resulting / S_v) * 100;

% Display the percentage
disp(['The area of the resulting polygon is ' num2str(percentage) '% of the vertical tail area.']);
if percentage < 50
    disp('The spin recovery characterisitics are acceptable.');
else
    disp('The spin recovery characterisitics are not acceptable.');
end

% Plot the resulting area
figure; % Create a new figure
hold on; % Hold the current axes so that new plots do not delete existing plots

% Plot the vertical tail
p1 = plot(x_v_coords, z_v_coords, 'b', 'LineWidth', 2);

% Plot the horizontal airfoil profile
p2 = plot(naca0012_h(:,1), naca0012_h(:,2), 'r', 'LineWidth', 2);

% Plot the dashed lines
p3 = plot(x_le_line, y_le, 'k--', 'LineWidth', 2);
p4 = plot(x_te_line, y_te, 'k--', 'LineWidth', 2);

% Plot the resulting area
p5 = plot(resulting_area, 'FaceColor', 'g', 'FaceAlpha', 0.3);

% Add labels
set(gca, 'FontSize', 30, 'LineWidth', 2);
xlabel('x (distance from wing aerodynamic centre, m)', 'FontSize', 30);
ylabel('z (distance from neutral line, m)', 'FontSize', 30);

% Add legend
h_legend = legend([p1, p2, p3, p5], {'Vertical Airfoil Planform', 'Horizontal Tail Root Section (NACA 0012)', 'Horizontal Tail Wake', sprintf('Percentage of Vertical Tail Plan Form Area Within Wake, %.2f%%', percentage)});
set(h_legend, 'FontSize', 30);

hold off; % Release the current axes