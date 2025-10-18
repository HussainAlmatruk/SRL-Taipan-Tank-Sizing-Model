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
2. **Edit Inputs**: Open the `.m` file in MATLAB. Modify the variables in the `%% INPUTS` section as needed
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
