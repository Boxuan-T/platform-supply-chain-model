# Data Accumulation and Revenue Sharing in Platform-Mediated Supply Chains

MATLAB implementation of **Data Accumulation and Revenue Sharing in Platform-Mediated Supply Chains: A Coupled Evolutionary Game Model**.

The model couples platform investment, manufacturer and retailer data sharing, and an endogenous data stock. Numerical experiments evaluate cooperation stability, attraction basins, revenue-sharing allocations, and transition payoffs.

## Requirements

MATLAB. Tested with R2024a.

## Quick start

Open the repository folder in MATLAB and run:

```matlab
result = minimal_example();
```

Run an experiment:

```matlab
result = run_experiment('matched_data');
```

Run all experiments and numerical validation:

```matlab
run_all();
```

Run tests:

```matlab
run_tests();
```

## Experiments

| Name | Analysis |
| --- | --- |
| `reference` | Reference configuration |
| `baseline` | Baseline strategy trajectories |
| `initial_data` | Initial data availability |
| `basin_sections` | Sections of attraction basins |
| `accumulation` | Data integration and depreciation |
| `allocation` | Revenue-sharing grid |
| `mechanisms` | Data dynamics and revenue functions |
| `sensitivity` | Parameter sensitivity |
| `adaptation` | Strategy adjustment rates |
| `matched_data` | Finite and instantaneous adjustment with matched initial data |
| `asymmetry` | Manufacturer costs and retailer data contributions |

## Outputs

Experiment results are saved to `results/` as MAT files and CSV tables. Checkpoints support resuming completed configurations.

```matlab
run_experiment('allocation', 'results/allocation');
```

## Authors

Boxuan Tian and Jing Gao

Contact: [13994012222@163.com](mailto:13994012222@163.com)
