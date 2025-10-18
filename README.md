# SRL Liquid Rocket - Tank Sizing and Mass Model

**Authors:** Hussain Almatruk

**Last Updated:** October 17, 2025

This repository contains the official MATLAB model for the design and analysis of the propellant and pressurant tank system for the Taipan rocket engine.

## Project Purpose

The purpose of this model is to determine the optimal dimensions, materials, and mass of our LOX, Jet-A, and Nitrogen pressurant tanks. It takes fundamental engine and vehicle parameters as inputs and calculates a full vehicle mass budget and key performance metrics.

## Guiding Philosophy

Our approach is **"Calculate, Don't Assume."** The script is designed to derive as much information as possible from a small set of core inputs. For example, we input the engine's thrust and efficiency (`Isp`), and the code calculates the required propellant mass not the other way around. This ensures that our design choices are directly tied to performance requirements.

# Rocket Sizing Script

A modular MATLAB script for rocket propulsion system sizing and performance calculations.

## Script Structure

The script is organized into sequential blocks for clear organization and maintainability:

- **`%% 0.0 - SETUP & CONSTANTS`**: Clears the workspace and defines universal physical constants
- **`%% 1.0 - INPUTS`**: All primary design choices, material properties, and assumptions are defined here. This is the only section that should be modified for new trade studies
- **`%% 2.0 - CALCULATIONS`**: (WORK IN PROGRESS) This section will contain all the engineering equations. It is currently a placeholder and needs to be filled in
- **`%% 3.0 - OUTPUTS`**: (WORK IN PROGRESS) This section will display the final results in a clean, formatted summary

## How to Use This Repository

1. **Pull**: Before starting work, open GitHub Desktop and Pull the latest version to ensure you are up-to-date
2. **Edit Inputs**: Open the `.m` file in MATLAB. Modify the variables in the `%% INPUTS` section as needed. currently this file is "tank sizing script.m"
3. **Add Calculations**: Implement the necessary formulas in the `%% CALCULATIONS` section
4. **Run & Analyze**: Execute the script and review the results
5. **Commit & Push**: Save your changes, go to GitHub Desktop, Commit them with a clear message, and Push them to the repository

## ✅ Current Tasks & To-Do List

Our first priority is to get the skeleton fully functional.

This might require having more inputs or outputs than the ones in the trade study documents. Make sure you: discuss beefore, and make it clear when you add any equation or a new input that is not already in the trade study document.

### 1. Define Input TODOs:
- [ ] Research and define the material properties for the LOX, Fuel, and Pressurant tanks (`material_density_*`, `material_allowable_stress_*`)
- [ ] Determine appropriate `joint_efficiency_*` values for our chosen materials and welding processes
- [ ] Create initial estimates for `m_misc_kg` and `m_plumbing_kg`

### 2. Implement Core Calculations:
- [ ] Add formulas to section 2.1 to calculate propellant mass flow rates and total propellant mass
- [ ] Add formulas to section 2.2 to calculate tank volumes, dimensions, wall thicknesses, and empty masses
- [ ] Add formulas to section 2.3 to sum all component masses into a total vehicle mass
- [ ] Add formulas to section 2.4 to calculate TWR and Delta-V

### 3. Create Formatted Outputs:
- [ ] Add `fprintf` statements to the OUTPUTS section to display the results clearly

---------------------------------------------------------


graph TD
%% Define Node Styles
classDef inputs fill:#E6F3FF,stroke:#B3D9FF,stroke-width:2px;
classDef outputs fill:#E6FFED,stroke:#B3FFC6,stroke-width:2px;
classDef calcs fill:#FFF3E6,stroke:#FFD9B3,stroke-width:2px;

%% === INPUTS BLOCK ===
subgraph "A. Core Inputs"
    direction TB
    A1["Engine Performance<br/>(Thrust, Isp, Burn Time, O/F Ratio)"]
    A2["Propellant Properties<br/>(LOX & Fuel Densities)"]
    A3["Pressurant Properties<br/>(Gas Temp, Molar Mass, Storage Pressure)"]
    A4["Vehicle Geometry<br/>(Tank Diameters)"]
    A5["Material Properties<br/>(Density, Allowable Stress, Joint Efficiency)"]
    A6["Design Margins & Estimates<br/>(Safety Factor, Ullage, Plumbing/Misc Mass)"]
    A7["Engine Requirements<br/>(Injector Feed Pressure, Plumbing Pressure Drop)"]
end
class A1,A2,A3,A4,A5,A6,A7 inputs;

%% === CALCULATIONS BLOCK 1: PROPULSION ===
subgraph "B. Step 1: Propulsion System Analysis"
    direction TB
    B1["Calculate Total Mass Flow Rate<br/>(ṁ_total)"]
    B2["Calculate Propellant Flow Rates<br/>(ṁ_ox, ṁ_fuel)"]
    B3["Calculate Total Propellant Masses<br/>(m_ox, m_fuel)"]
end
class B1,B2,B3 calcs;

%% Connect Inputs to Step 1
A1 --> B1;
B1 --> B2;
B2 --> B3;

%% === PARALLEL CALCULATIONS: TANK SIZING ===
subgraph "C. Steps 2 & 3: Parallel Tank Sizing"
    direction LR
    %% === LOX TANK SUBGRAPH ===
    subgraph "Step 2: LOX Tank"
        direction TB
        C1_1["Calculate LOX Volume & Total Tank Volume"]
        C1_2["Calculate LOX Tank Dimensions<br/>(Cylinder Length)"]
        C1_3["Calculate LOX Tank Pressures<br/>(Operating & Design)"]
        C1_4["<b>Output: LOX Tank Wall Thicknesses</b><br/>(t_cyl_ox, t_caps_ox)"]
        C1_5["<b>Output: Empty LOX Tank Mass</b><br/>(m_empty_tank_ox)"]
    end

    %% === FUEL TANK SUBGRAPH ===
    subgraph "Step 3: Fuel Tank"
        direction TB
        C2_1["Calculate Fuel Volume & Total Tank Volume"]
        C2_2["Calculate Fuel Tank Dimensions<br/>(Cylinder Length)"]
        C2_3["Calculate Fuel Tank Pressures<br/>(Operating & Design)"]
        C2_4["<b>Output: Fuel Tank Wall Thicknesses</b><br/>(t_cyl_fuel, t_caps_fuel)"]
        C2_5["<b>Output: Empty Fuel Tank Mass</b><br/>(m_empty_tank_fuel)"]
    end
end
class C1_1,C1_2,C1_3,C2_1,C2_2,C2_3 calcs;
class C1_4,C1_5,C2_4,C2_5 outputs;


%% Connect Inputs to Parallel Blocks
B3 -- m_ox --> C1_1;
B3 -- m_fuel --> C2_1;
A2 --> C1_1;
A2 --> C2_1;
A4 --> C1_2;
A4 --> C2_2;
A7 --> C1_3;
A7 --> C2_3;
A5 --> C1_4;
A5 --> C2_4;
A5 --> C1_5;
A5 --> C2_5;
A6 --> C1_1;
A6 --> C2_1;
A6 --> C1_3;
A6 --> C2_3;
C1_1 --> C1_2 --> C1_3 --> C1_4 --> C1_5;
C2_1 --> C2_2 --> C2_3 --> C2_4 --> C2_5;

%% === CALCULATIONS BLOCK 4: PRESSURANT ===
subgraph "D. Step 4: Pressurant System Sizing"
    direction TB
    D1["Calculate Required Pressurant Moles<br/>(n_total)"]
    D2["Calculate Pressurant Gas Mass"]
    D3["Calculate Pressurant Tank Internal Volume & Radius"]
    D4["<b>Output: Pressurant Tank Wall Thickness</b>"]
    D5["<b>Output: Empty Pressurant Tank Mass</b>"]
end
class D1,D2,D3 calcs;
class D4,D5 outputs;

%% Connect Parallel Blocks to Pressurant Block
C1_1 -- Tank Volume & Pressure --> D1;
C2_1 -- Tank Volume & Pressure --> D1;
C1_3 -- Tank Volume & Pressure --> D1;
C2_3 -- Tank Volume & Pressure --> D1;
A3 --> D1;
A3 --> D2;
A3 --> D3;
A5 --> D4;
A5 --> D5;
D1 --> D2 --> D3 --> D4 --> D5;


%% === FINAL CALCULATIONS & OUTPUTS ===
subgraph "E. Step 5: Vehicle Mass Buildup & Performance"
    direction TB
    E1["<b>Output: Total Tank System Mass</b>"]
    E2["<b>Output: Vehicle Liftoff Mass</b>"]
    E3["<b>Output: Thrust-to-Weight Ratio (TWR)</b>"]
    E4["Calculate Vehicle Dry Mass"]
    E5["<b>Output: Ideal Delta-V</b>"]
    E6["<b>Output: Estimated Apogee</b>"]
end
class E4 calcs;
class E1,E2,E3,E5,E6 outputs;

%% Connect all masses to final block
B3 -- Propellant Masses --> E1;
C1_5 -- Empty Tank Masses --> E1;
C2_5 -- Empty Tank Masses --> E1;
D2 -- Pressurant Mass --> E1;
D5 -- Empty Tank Masses --> E1;
A6 -- Plumbing/Misc Mass --> E2;
E1 --> E2 --> E3;
E2 -- Total Mass --> E4;
B3 -- Propellant Mass --> E4;
E4 -- Dry Mass --> E5;
E5 --> E6;
