# SRL Liquid Flight Vehicle Tank Sizing Model

**Authors:** Hussain Almatruk, Jonathan Forte

**Current model:** Version 3 MATLAB model

This repository contains the MATLAB sizing and analysis model for the SRL
Liquid Flight Vehicle tank stack. The model is currently set up to support
multiple related vehicle/engine configurations, including Taipan and
Scorpion, using the same core calculation flow.

Phase 1 reference document:
https://www.overleaf.com/read/jwzdsjwwgnff#6aaf7d

## Purpose

The purpose of this model is to give the team a repeatable way to estimate
propellant loads, tank volumes, tank geometry, pressurant requirements,
vehicle mass properties, and first-pass flight performance from a controlled
set of vehicle inputs.

The model is not intended to replace detailed CAD, test data, structural
review, operations planning, or safety review. It is a working engineering
tool for early design trades and for keeping the team aligned on the current
assumptions.

## Design Philosophy

The model should calculate as much as practical from a small set of clear
inputs. For example, thrust, specific impulse, burn time, and mixture ratio
drive propellant mass rather than entering propellant mass by hand.

When possible, keep assumptions in the vehicle input files and keep equations
inside the shared calculation functions. This makes it easier to compare
vehicles without changing the math.

## Repository Structure

```text
.
|-- TaipanModelVersion3_0.m      # Main script / run order
|-- Vehicles/                    # Vehicle-specific inputs
|   |-- Taipan.m
|   `-- Scorpion.m
|-- V3_Files/                    # Shared calculation, plotting, and output functions
|-- Extra_Function/              # Trade studies and helper analyses
|-- Outputs/                     # Generated Excel sheets and figures
|-- Archive/                     # Older model versions and templates
`-- 4in_RASAero_CD_data.CSV      # Current aerodynamic lookup data
```

## Main Files

- `TaipanModelVersion3_0.m` is the main script. It sets up paths, selects the
  vehicle configuration, runs the calculation functions, displays results,
  writes output files, and generates plots.
- `Vehicles/Taipan.m` and `Vehicles/Scorpion.m` hold vehicle-specific design
  inputs. Start here when changing assumptions for a vehicle.
- `V3_Files/getConstants.m` holds unit conversions and physical constants.
- `V3_Files/calcMassAndSizing.m` contains the main propellant, tank, pressurant,
  mass, and geometry calculations.
- `V3_Files/runFlightSimulation.m` contains the current 1D flight simulation.
- `V3_Files/displayResults.m`, `writeResults.m`, and `plotVehicleStackup.m`
  handle outputs and visualization.

## How To Run

1. Pull the latest version of the repository.
2. Open MATLAB from the repository root.
3. Open `TaipanModelVersion3_0.m`.
4. Select the vehicle near the top of the script:

   ```matlab
   vehicle_name = 'Scorpion';
   ```

   Valid current options are:

   ```matlab
   vehicle_name = 'Taipan';
   vehicle_name = 'Scorpion';
   ```

5. Edit vehicle assumptions in the matching file under `Vehicles/`.
6. Run `TaipanModelVersion3_0.m`.
7. Review the command-window output, generated Excel file, and generated plots.

## Generated Outputs

When output writing is enabled, the model writes vehicle-specific outputs to
the `Outputs/` folder. Typical generated files include:

- `<Vehicle>_Outputs.xlsx`
- `<Vehicle>_Vehicle_Stackup.png`
- `<Vehicle>_Pressurant_Thermal_Analysis.png`

These files are useful for review, but they are generated artifacts. Before
committing them, make sure the outputs are intentionally part of the update.

## Development Workflow

1. Pull before starting work.
2. Make focused changes.
3. Run the model for the affected vehicle configuration.
4. Check for warnings in the command-window output.
5. Review generated plots and spreadsheets if your change affects outputs.
6. Commit with a clear message describing what changed and why.
7. Push only the files that belong in the update.

## Current Work Areas

The model is under active development. Current work generally falls into these
areas:

- Keeping Taipan and Scorpion support clean and consistent.
- Improving tank sizing, pressurant sizing, and vehicle stackup calculations.
- Making output tables and plots easier to read and review.
- Replacing temporary assumptions with reviewed design values as the team
  makes decisions.
- Improving aerodynamic and flight-performance modeling as better data becomes
  available.
- Keeping documentation current enough that new team members can understand
  the model without needing the full project history.

## Notes For Contributors

- Put vehicle-specific assumptions in `Vehicles/` whenever possible.
- Put shared equations and reusable logic in `V3_Files/`.
- Use `Extra_Function/` for trade studies, experiments, and helper analyses
  that should not be part of the main run path yet.
- Avoid hard-coding values inside calculation functions if they belong in a
  vehicle input file.
- If a number is temporary, label it clearly and include enough context for
  someone else to improve it later.
- Keep generated output changes separate from source-code changes when
  practical.

## Status

This is a working design model, not a frozen release. Treat the results as
engineering estimates that depend on the current assumptions in the vehicle
files. Any result used for design decisions should be checked against the
latest team requirements, P&ID, vendor data, and test results.
