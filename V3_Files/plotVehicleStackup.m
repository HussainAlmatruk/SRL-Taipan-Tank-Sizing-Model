function plotVehicleStackup( vehicle_name, ...
    d_ox_tank_m, l_cyl_ox_m, l_ox_tank_total_m, ...
    d_fuel_tank_m, l_cyl_fuel_m, l_fuel_tank_total_m, ...
    d_pressurant_tank_outer_m, l_cyl_pressurant_m, l_pressurant_tank_total_m, ...
    gap_lox_fuel_m, gap_fuel_pres_m)
%{
---------------------------------------------------------------------------
plotVehicleStackup
  Generates a 2D and 3D visualization of the vehicle tank stackup.
  The tanks are stacked from bottom to top in the following order:
    1. LOX Tank (Bottom, Blue)
    2. Fuel Tank (Middle, Green)
    3. Pressurant Tank (Top, Yellow)

  All tanks share a common centerline. This centerline is shifted 
  positively so the left-most edge of the largest tank is at X = 0.
---------------------------------------------------------------------------
%}

% --- 1.0 Setup Colors and Geometry ---
% RGB triplets for specific tank colors
color_lox_blue    = [0.00, 0.45, 0.74]; % Blue for Liquid Oxygen
color_fuel_green  = [0.47, 0.67, 0.19]; % Green for Fuel (Jet-A)
color_pres_yellow = [0.93, 0.69, 0.13]; % Yellow for Pressurant (N2)

% Conversion factor for the secondary axis (m -> in)
m_to_in = 39.3701; 

% Calculate the max diameter to determine the common central axis
d_max_m = max([d_ox_tank_m, d_fuel_tank_m, d_pressurant_tank_outer_m]);
x_center_m = d_max_m / 2;
y_center_m = d_max_m / 2; % For 3D centering

% Calculate the starting Z-coordinate (height) for the bottom of each tank
z_bottom_lox_m  = 0;
z_bottom_fuel_m = z_bottom_lox_m + l_ox_tank_total_m + gap_lox_fuel_m;
z_bottom_pres_m = z_bottom_fuel_m + l_fuel_tank_total_m + gap_fuel_pres_m;

% Calculate the total height of the stackup
z_max_m = z_bottom_pres_m + l_pressurant_tank_total_m;

% --- 2.0 Initialize Figure ---
figure('Name', sprintf('%s Vehicle Tank Stackup', vehicle_name), 'NumberTitle', 'off', 'Position', [150, 150, 1000, 800]);

% --- 3.0 Generate 2D Stackup Plot ---
ax1 = subplot(1, 2, 1);
hold(ax1, 'on'); grid(ax1, 'on');

% Draw the 2D cross-sections
drawTank2D(ax1, x_center_m, z_bottom_lox_m,  d_ox_tank_m,  l_cyl_ox_m,  l_ox_tank_total_m,  color_lox_blue);
drawTank2D(ax1, x_center_m, z_bottom_fuel_m, d_fuel_tank_m, l_cyl_fuel_m, l_fuel_tank_total_m, color_fuel_green);
drawTank2D(ax1, x_center_m, z_bottom_pres_m, d_pressurant_tank_outer_m, l_cyl_pressurant_m, l_pressurant_tank_total_m, color_pres_yellow);

% Define Y-limits with a small margin top and bottom
y_min = -0.1;
y_max = z_max_m + 0.2;
y_range = y_max - y_min;

% Make the X-axis range a proportional width of the Y-axis range (approx 45%)
% This forces the white background grid to be a "normal" looking rectangle
% instead of shrinking tightly around the skinny rocket.
x_range = max(y_range * 0.45, d_max_m * 3); 

xlim(ax1, [x_center_m - x_range/2, x_center_m + x_range/2]);
ylim(ax1, [y_min, y_max]);

% Force true 1:1 data aspect ratio (circles will look like circles)
axis(ax1, 'equal'); 

xlabel(ax1, 'Position [m]');
ylabel(ax1, 'Height [m]');
title(ax1, '2D Vehicle Stackup');

% Create dummy patches strictly for a clean legend
patch(ax1, NaN, NaN, color_lox_blue, 'DisplayName', 'LOX Tank');
patch(ax1, NaN, NaN, color_fuel_green, 'DisplayName', 'Fuel Tank');
patch(ax1, NaN, NaN, color_pres_yellow, 'DisplayName', 'Pressurant Tank');
legend(ax1, 'Location', 'best');

% Force a draw to establish the limits of the meters axis before making the inches axis
drawnow; 

% Create a clean secondary axis on top/right for inches
ax2 = axes('Position', ax1.Position, 'XAxisLocation', 'top', ...
           'YAxisLocation', 'right', 'Color', 'none');
ax2.XLim = xlim(ax1) * m_to_in;
ax2.YLim = ylim(ax1) * m_to_in;

% CRITICAL: Setting axis equal on the secondary axis forces its box 
% to snap perfectly onto the primary axis box so it doesn't float away.
axis(ax2, 'equal');

xlabel(ax2, 'Position [in]');
ylabel(ax2, 'Height [in]');

% --- 4.0 Generate 3D Stackup Plot ---
ax3 = subplot(1, 2, 2);
hold(ax3, 'on'); grid(ax3, 'on');

% Draw the 3D surfaces
drawTank3D(ax3, x_center_m, y_center_m, z_bottom_lox_m,  d_ox_tank_m,  l_cyl_ox_m,  l_ox_tank_total_m,  color_lox_blue);
drawTank3D(ax3, x_center_m, y_center_m, z_bottom_fuel_m, d_fuel_tank_m, l_cyl_fuel_m, l_fuel_tank_total_m, color_fuel_green);
drawTank3D(ax3, x_center_m, y_center_m, z_bottom_pres_m, d_pressurant_tank_outer_m, l_cyl_pressurant_m, l_pressurant_tank_total_m, color_pres_yellow);

% Force true 1:1:1 scale and tight limits
axis(ax3, 'tight');
axis(ax3, 'equal');

% Add margin for clarity
zlim(ax3, [-0.1, z_max_m + 0.2]);

view(ax3, 3);
xlabel(ax3, 'X [m]');
ylabel(ax3, 'Y [m]');
zlabel(ax3, 'Height [m]');
title(ax3, '3D Vehicle Stackup');
camlight(ax3, 'headlight');
lighting(ax3, 'gouraud');
material(ax3, 'dull');

% --- 5.0 Save Figure ---
if ~exist('Outputs', 'dir')
    mkdir('Outputs');
end
saveas(gcf, fullfile('Outputs', sprintf('%s_Vehicle_Stackup.png', vehicle_name)));

end

% =========================================================================
% Local Helper Functions
% =========================================================================

function drawTank2D(ax, x_center_m, z_bottom_m, d_tank_m, l_cyl_m, l_total_m, tank_color)
    r_tank_m = d_tank_m / 2;
    h_cap_m  = (l_total_m - l_cyl_m) / 2;

    % Create the curve for the bottom cap
    theta_bot = linspace(pi, 2*pi, 30);
    x_bot = x_center_m + r_tank_m * cos(theta_bot);
    z_bot = z_bottom_m + h_cap_m + h_cap_m * sin(theta_bot);

    % Create the curve for the top cap
    theta_top = linspace(0, pi, 30);
    x_top = x_center_m + r_tank_m * cos(theta_top);
    z_top = z_bottom_m + h_cap_m + l_cyl_m + h_cap_m * sin(theta_top);

    % Combine points into a closed polygon
    x_poly = [x_bot, x_top];
    z_poly = [z_bot, z_top];

    % Draw and fill the 2D shape with 'HandleVisibility', 'off' to clean the legend
    fill(ax, x_poly, z_poly, tank_color, 'EdgeColor', 'k', 'LineWidth', 1.2, 'HandleVisibility', 'off');
end

function drawTank3D(ax, x_center_m, y_center_m, z_bottom_m, d_tank_m, l_cyl_m, l_total_m, tank_color)
    r_tank_m = d_tank_m / 2;
    h_cap_m  = (l_total_m - l_cyl_m) / 2;

    % 1. Bottom Cap
    [theta_bot, phi_bot] = meshgrid(linspace(0, 2*pi, 40), linspace(pi/2, pi, 20));
    x_bot = x_center_m + r_tank_m * cos(theta_bot) .* sin(phi_bot);
    y_bot = y_center_m + r_tank_m * sin(theta_bot) .* sin(phi_bot);
    z_bot = z_bottom_m + h_cap_m + h_cap_m * cos(phi_bot);
    surf(ax, x_bot, y_bot, z_bot, 'FaceColor', tank_color, 'EdgeColor', 'none', 'HandleVisibility', 'off');

    % 2. Cylindrical Body
    if l_cyl_m > 0
        [theta_cyl, z_cyl_grid] = meshgrid(linspace(0, 2*pi, 40), linspace(0, l_cyl_m, 2));
        x_cyl = x_center_m + r_tank_m * cos(theta_cyl);
        y_cyl = y_center_m + r_tank_m * sin(theta_cyl);
        z_cyl = z_bottom_m + h_cap_m + z_cyl_grid;
        surf(ax, x_cyl, y_cyl, z_cyl, 'FaceColor', tank_color, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    end

    % 3. Top Cap
    [theta_top, phi_top] = meshgrid(linspace(0, 2*pi, 40), linspace(0, pi/2, 20));
    x_top = x_center_m + r_tank_m * cos(theta_top) .* sin(phi_top);
    y_top = y_center_m + r_tank_m * sin(theta_top) .* sin(phi_top);
    z_top = z_bottom_m + h_cap_m + l_cyl_m + h_cap_m * cos(phi_top);
    surf(ax, x_top, y_top, z_top, 'FaceColor', tank_color, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end