function C = getConstants()
%{
---------------------------------------------------------------------------
getConstants
  Returns a struct containing all fixed physical constants and unit
  conversion factors used throughout the Taipan model. Nothing in this
  file should ever need to change unless a constant definition is being
  corrected.

  To use a value from this struct in a calculation:
      C.g_earth_ms2, C.LBF_TO_N, C.PSI_TO_PA, etc.

Output:
  C - [struct] Struct containing all constants and conversion factors
---------------------------------------------------------------------------
%}

% --- Physical Constants ---
C.g_earth_ms2       = 9.80665;     % [m/s^2]      Standard gravitational acceleration
C.r_universal_jmolk = 8.3145;      % [J/(mol*K)]  Universal gas constant

% --- Unit Conversion Factors ---

% Length
C.FT_TO_M = 0.3048;
C.M_TO_FT = 3.28084;
C.IN_TO_M = 0.0254;
C.M_TO_IN = 39.3701;
C.M_TO_MM = 1000;
C.MM_TO_M = 0.001;

% Mass
C.LBM_TO_KG = 0.453592;
C.KG_TO_LBM = 2.20462;

% Force
C.LBF_TO_N = 4.44822;
C.N_TO_LBF = 0.224809;

% Pressure
C.PSI_TO_PA = 6894.76;
C.PA_TO_PSI = 0.000145038;
C.ATM_TO_PA = 101325;
C.PA_TO_ATM = 9.86923e-6;
C.BAR_TO_PA = 100000;
C.PA_TO_BAR = 1e-5;

% Temperature
% Note: These are for differences, not absolute temperatures
% For absolute: T_K = T_C + 273.15, T_R = T_F + 459.67
C.C_TO_K_OFFSET = 273.15;
C.F_TO_R_OFFSET = 459.67;

% Angle
C.DEG_TO_RAD = pi / 180;
C.RAD_TO_DEG = 180 / pi;

% Velocity
C.KPH_TO_MS = 1/3.6;
C.MS_TO_KPH = 3.6;
C.MPH_TO_MS = 0.44704;
C.MS_TO_MPH = 2.23694;

end
