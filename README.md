# SRL Liquid Rocket - Tank Sizing and Mass Model

**Authors:** Hussain Almatruk, Jonathan Forte

**Last Updated:** November 2, 2025

This repository contains the official MATLAB model for the design and analysis of the propellant and pressurant tank system for the Taipan rocket engine.

Official Project's Document For Phase 1 (Overleaf): https://www.overleaf.com/read/jwzdsjwwgnff#6aaf7d

## Project Purpose

The purpose of this model is to determine the optimal dimensions, materials, and mass of our LOX, Jet-A, and Nitrogen pressurant tanks. It takes fundamental engine and vehicle parameters as inputs and calculates a full vehicle mass budget and key performance metrics.

## Guiding Philosophy

Our approach is **"Calculate, Don't Assume."** The script is designed to derive as much information as possible from a small set of core inputs. For example, we input the engine's thrust and efficiency (`Isp`), and the code calculates the required propellant mass not the other way around. This ensures that our design choices are directly tied to performance requirements.

## Rocket Sizing Script

A modular MATLAB script for rocket propulsion system sizing and performance calculations.

## Script Structure

The script is organized into sequential blocks for clear organization and maintainability:

- **`%% 0.0 - SETUP & CONSTANTS`**: Clears the workspace and defines universal physical constants.
- **`%% 1.0 - INPUTS`**: All primary design choices, material properties, and assumptions are defined here. This is the only section that should be modified for new trade studies
- **`%% 2.0 - CALCULATIONS`**: This section will contain all the engineering equations.
- **`%% 3.0 - MODEL PHASES OF FLIGHT`**: Contains time-dependant simulations, AKA the flight sim.
- **`%% 4.0 - VALIDATION AND CHECKS`**:  This section implements the checks defined in the technical plan.
- **`%% 5.0 - OUTPUTS`**:  displays the final calculated values in a clean format.

## How to Use This Repository

1. **Pull**: Before starting work, open GitHub Desktop and Pull the latest version to ensure you are up-to-date
2. **Edit Inputs**: Open the `.m` file in MATLAB. Modify the variables in the `%% INPUTS` section as needed. currently this file is "tank sizing script.m"
3. **Add Calculations**: Implement the necessary formulas in the `%% CALCULATIONS` section
4. **Run & Analyze**: Execute the script and review the results
5. **Commit & Push**: Save your changes, go to GitHub Desktop, Commit them with a clear message, and Push them to the repository

## ✅ Current Tasks & To-Do List




