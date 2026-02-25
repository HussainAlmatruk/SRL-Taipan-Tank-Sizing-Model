%{
---------------------------------------------------------------------------
CU Sounding Rocket Lab - Liquid Propulsion Team
Pressurant Tank MEOP Optimization Trade Study

Description:
This script iterates through various Nitrogen (N2) storage pressures 
to find the optimal Maximum Expected Operating Pressure (MEOP) that 
minimizes the overall mass of the pressurant system (Gas + COPV).

The final Result is the minimum 

---------------------------------------------------------------------------
%}

clear; clc; close all;

%% 1.0 - CONSTANTS & FIXED PARAMETERS
PSI_TO_PA = 6894.76;
PA_TO_PSI = 0.000145038;
IN_TO_M   = 0.0254;
r_universal_jmolk = 8.3145;
pressurant_molar_mass_kgmol = 0.0280134;

% --- Engine & Tank Requirements (from main model) ---
p_op_ox_tank_pa   = 1000 * PSI_TO_PA;
p_op_fuel_tank_pa = 1000 * PSI_TO_PA;
v_total_ox_tank_m3   = 0.0077; % [m^3] LOX Tank Volume (with ullage)
v_total_fuel_tank_m3 = 0.0042; % [m^3] Fuel Tank Volume (with ullage)

% --- Temperatures ---
ox_temp_k         = 90;
fuel_temp_k       = 294;
pressurant_temp_k = 294;

% --- Pressurant Tank (COPV) Properties ---
safety_factor                           = 1.5;
joint_efficiency_pressurant_tank        = 0.8;
corrosion_allowance_m                   = 0.0002;
material_density_pressurant_kgm3        = 2700; % Aluminum/Composite proxy
material_allowable_stress_pressurant_pa = 2.76e8;
material_density_liner_kgm3             = 2700;
t_liner_m                               = 0.002; % 2mm liner thickness proxy

% Fixed Outer Diameter Constraint
d_pressurant_tank_allowed_m = 4.5 * IN_TO_M; 
r_pressurant_outer_m        = d_pressurant_tank_allowed_m / 2;

%% 2.0 - CALCULATE REQUIRED GAS MASS (Isothermal Baseline)
% Calculate moles required to displace propellants
n_ox_mol = (p_op_ox_tank_pa * v_total_ox_tank_m3) / (r_universal_jmolk * (ox_temp_k + pressurant_temp_k)/2);
n_fuel_mol = (p_op_fuel_tank_pa * v_total_fuel_tank_m3) / (r_universal_jmolk * fuel_temp_k);

n_total_mol = n_ox_mol + n_fuel_mol;
m_pressurant_gas_kg = n_total_mol * pressurant_molar_mass_kgmol;

%% 3.0 - SWEEP STORAGE PRESSURES
% Create an array of test pressures from 1500 psi to 6000 psi
p_test_psi = 1500:50:6000;
p_test_pa  = p_test_psi * PSI_TO_PA;

% Pre-allocate arrays for storing results
m_copv_empty_kg  = zeros(size(p_test_psi));
m_total_sys_kg   = zeros(size(p_test_psi));
l_copv_total_m   = zeros(size(p_test_psi));

for i = 1:length(p_test_pa)
    p_storage = p_test_pa(i);
    p_design  = p_storage * safety_factor;
    
    % 1. Required Internal Volume at this pressure
    v_internal_m3 = (n_total_mol * r_universal_jmolk * pressurant_temp_k) / p_storage;
    
    % 2. Wall Thickness (Hoop Stress)
    t_wall_m = (p_design * r_pressurant_outer_m) / (material_allowable_stress_pressurant_pa * joint_efficiency_pressurant_tank) + corrosion_allowance_m;
    
    % 3. Internal Geometry
    r_internal_m = r_pressurant_outer_m - t_wall_m - t_liner_m;
    r_liner_outer_m = r_internal_m + t_liner_m;
    
    % Check if wall thickness exceeds radius (impossible geometry)
    if r_internal_m <= 0
        m_total_sys_kg(i) = NaN; l_copv_total_m(i) = NaN;
        continue;
    end
    
    % 4. Volumes and Lengths
    v_caps_m3 = (2/3) * pi * r_internal_m^3; % Two 2:1 ellipsoidal caps
    v_cyl_m3  = v_internal_m3 - v_caps_m3;
    
    if v_cyl_m3 < 0
        % Pressure is so high that the end caps alone are too big
        l_cyl_m = 0;
    else
        l_cyl_m = v_cyl_m3 / (pi * r_internal_m^2);
    end
    l_copv_total_m(i) = l_cyl_m + (2 * r_pressurant_outer_m);
    
    % 5. Mass Calculations
    m_shell_cyl  = pi * (r_pressurant_outer_m^2 - r_liner_outer_m^2) * l_cyl_m;
    m_shell_caps = (2/3) * pi * (r_pressurant_outer_m^3 - r_liner_outer_m^3);
    m_liner_cyl  = pi * (r_liner_outer_m^2 - r_internal_m^2) * l_cyl_m;
    m_liner_caps = (2/3) * pi * (r_liner_outer_m^3 - r_internal_m^3);
    
    m_copv_empty_kg(i) = material_density_pressurant_kgm3 * (m_shell_cyl + m_shell_caps) + ...
                         material_density_liner_kgm3 * (m_liner_cyl + m_liner_caps);
                     
    m_total_sys_kg(i) = m_copv_empty_kg(i) + m_pressurant_gas_kg;
end

%% 4.0 - FIND OPTIMUM AND PLOT
% Find the minimum total mass and corresponding pressure
[min_mass, idx_opt] = min(m_total_sys_kg);
opt_p_psi = p_test_psi(idx_opt);
opt_length_m = l_copv_total_m(idx_opt);

fprintf('====================================================\n');
fprintf('         N2 MEOP Optimization Results               \n');
fprintf('====================================================\n');
fprintf('Optimal Storage Pressure:  %8.0f psi\n', opt_p_psi);
fprintf('Minimum Total Mass:        %8.2f kg\n', min_mass);
fprintf('Resulting COPV Length:     %8.2f m (%.2f in)\n', opt_length_m, opt_length_m / IN_TO_M);
fprintf('====================================================\n');

% Plotting
figure('Name', 'Pressurant Tank MEOP Trade Study', 'Position', [100, 100, 800, 400]);

% Subplot 1: Mass vs Pressure
subplot(1,2,1);
plot(p_test_psi, m_total_sys_kg, 'b-', 'LineWidth', 2); hold on;
plot(p_test_psi, m_copv_empty_kg, 'k--', 'LineWidth', 1.5);
plot(opt_p_psi, min_mass, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
grid on;
title('System Mass vs. Storage Pressure');
xlabel('Storage Pressure (psi)');
ylabel('Mass (kg)');
legend('Total System Mass', 'Empty COPV Mass', 'Optimal MEOP', 'Location', 'Best');

% Subplot 2: Length vs Pressure
subplot(1,2,2);
plot(p_test_psi, l_copv_total_m, 'g-', 'LineWidth', 2); hold on;
plot(opt_p_psi, opt_length_m, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
grid on;
title('Tank Length vs. Storage Pressure');
xlabel('Storage Pressure (psi)');
ylabel('Total Length (m)');
legend('COPV Length', 'Optimal MEOP', 'Location', 'Best');