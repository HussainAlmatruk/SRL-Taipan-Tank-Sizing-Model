# SRL Liquid Rocket - Tank Sizing and Mass Model

**Authors:** Hussain Almatruk
**Last Updated:** October 17, 2025

This repository contains the official MATLAB model for the design and analysis of the propellant and pressurant tank system for the Taipan rocket engine.

## Project Purpose

The purpose of this model is to determine the optimal dimensions, materials, and mass of our LOX, Jet-A, and Nitrogen pressurant tanks. It takes fundamental engine and vehicle parameters as inputs and calculates a full vehicle mass budget and key performance metrics.

## Guiding Philosophy

Our approach is **"Calculate, Don't Assume."** The script is designed to derive as much information as possible from a small set of core inputs. For example, we input the engine's thrust and efficiency (`Isp`), and the code calculates the required propellant mass not the other way around. This ensures that our design choices are directly tied to performance requirements.

## Calculation Workflow

The script executes in a clear, logical sequence:

1.  **INPUTS:** All primary design choices and constants are defined in this section. This is the section that should mainly be modified for new trade studies.
2.  **Propulsion System Calculations:** Derives mass flow rates and total propellant mass from engine thrust, `Isp`, burn time, and O/F ratio.
3.  **Propellant Tank Sizing:** Calculates the required volume, dimensions, wall thickness, and empty mass for both the LOX and Fuel tanks.
4.  **Pressurant System Sizing:** Uses the Ideal Gas Law to determine the required mass of pressurant gas and then sizes a tank to hold it at high pressure.
5.  **Vehicle Mass Buildup:** Sums the mass of all components (full tanks, payload, plumbing, etc.) to find the total liftoff mass.
6.  **Performance Analysis & OUTPUTS:** Calculates the final Thrust-to-Weight Ratio (TWR), ideal ΔV, and an estimated apogee. All key results are printed clearly to the command window.

## How to Use

1.  **Pull:** Before starting work, always pull the latest version from GitHub to ensure you are up-to-date.
2.  **Edit Inputs:** Open the `.m` file in MATLAB. Modify the variables in the `%% INPUTS` section only.
3.  **Run:** Execute the script.
4.  **Analyze:** Review the summary printed in the command window.
5.  **Commit & Push:** Save your changes, commit them with a clear message in GitHub Desktop, and push them to the repository.

## Future Work & To-Do

- [ ] Add a more detailed model for `m_misc` and `m_plumbing` instead of using a single estimate.
- [ ] Integrate a basic atmospheric model for a more accurate apogee prediction.

