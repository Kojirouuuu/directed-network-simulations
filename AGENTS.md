# Repository Guidelines

## Project Structure & Module Organization

This repository is a Java 21 Gradle project for directed-network SAR and percolation simulations, with supporting C and Python analysis code.

Paths in this document are relative to the repository root unless stated otherwise. Run the listed commands from the repository root. Use repository-relative paths in documentation and instructions rather than machine-specific absolute paths.

- `app/src/main/java/sim/`: Java application code. Main entry points are `App.java`, `SAR.java`, `GraphGen.java`, and `BondPercolation.java`.
- `app/src/main/java/sim/network/`: CSR-style graph representation and topology generators.
- `app/src/main/java/sim/simulation/`: SAR and bond-percolation simulators and result types.
- `app/src/main/java/sim/utils/`: shared utilities and topology dispatch logic.
- `app/src/test/java/`: JUnit tests, mirroring Java package names.
- `c-lang/ebcm/`: C EBCM numerical programs.
- `notebooks/`: Jupyter notebooks and Python EBCM helpers.
- `test-notebooks/`: local analysis notebooks and helpers, when present; ignored by Git. Do not force-add these files unless explicitly requested.
- `scripts/`: convenience wrappers for common Java and C runs.

Generated outputs go under `build/`, `app/out/`, or `out/`; keep them out of commits unless explicitly required.

### Java Responsibilities

- `SAR` coordinates batches, simulation parameters, engine selection, and run summaries. Keep experiment settings in `SAR.SimulationConfig`.
- `SARGraphPreparation` handles graph generation/loading, edge-list discovery, randomization dispatch, and SAR network output paths.
- `DirectedGraph` retains the CSR representation, construction, validation, and public graph API. Internal metrics and transformations live in `GraphOperations` and `FflRewiringState`.
- `SARSimulator` and `SARGillespieSimulator` implement their respective simulation algorithms. Preserve random draw order and event ordering when refactoring.
- `SARResult` holds results and exposes CSV methods; `SARResultCsvWriter` implements serialization.
- `InitialAdopters` handles seeded initial-node selection, and `ProgressBar` handles console rendering.

## Build, Test, and Development Commands

- `./gradlew build`: compile, test, and assemble the project.
- `./gradlew test`: run all JUnit tests.
- `./gradlew :app:runSAR`: run the SAR simulation using `SAR.SimulationConfig`.
- `./gradlew :app:runGraphGen`: generate edge lists for later simulations.
- `./gradlew :app:runBondPercolation`: run bond-percolation simulations.
- `./gradlew run`: run the `sim.App` demo entry point.
- `bash scripts/run-c.sh`: compile and run the selected C EBCM program.

Use the Gradle wrapper rather than a system Gradle installation.

## Simulation Workflow

Follow these rules whenever the user requests a simulation experiment:

- Reuse the existing simulation code and change its parameters. Do not create new Java files or experiment entry points for a simulation request. This restriction does not apply when code changes or refactoring are the main task.
- Run simulations through the existing scripts in `scripts/`, such as `bash scripts/run-sar.sh`, `bash scripts/run-graphgen.sh`, or `bash scripts/run-bondpercolation.sh`.
- Use the same simulation implementation for every condition in a comparison. For example, to study the effect of a scale-free power-law exponent, run that implementation for each confirmed `gamma`, keeping other comparison conditions consistent.
- Ask the user before execution if any simulation parameter is uncertain, including parameters for the small-scale test. Do not treat current source defaults as confirmed choices. Choose an appropriate `optionPath` yourself so the experiment is identifiable and its outputs do not mix with existing results.
- Use `test-notebooks/` as the reference for visualization. If useful, copy an existing notebook into `test-notebooks/temp/<experiment-name>/`, change its parameters and input paths, and execute its cells. Preserve the original notebook.
- First test the entire workflow on a small configuration: simulation, output loading, aggregation, and plotting. Verify that the configuration, data, and figures serve the intended comparison before starting the large run.
- If results are clearly abnormal, show the user the configuration, observations, and concerns, and wait for clarification before proceeding to the large run.
- If the small-scale test has no apparent problems, delete only the test data, figures, and other result artifacts created by that test, then run the confirmed large-scale configuration without asking for another execution confirmation. Keep test outputs separate so cleanup cannot remove existing research data or original notebooks.

## Coding Style & Naming Conventions

Java code uses 4-space indentation, package-private tests, and descriptive `UpperCamelCase` class names. Methods, fields, and local variables use `lowerCamelCase`; constants use `UPPER_SNAKE_CASE`. Prefer deterministic seeds for simulation code and tests. Add new topology generators under `sim.network.topology` and wire them through `SwitchUtils.generateGraph` and `SwitchUtils.buildNetworkPath`.

Use `SwitchUtils.GraphGeneratorParams.builder(networkType, n)` for named network parameters. Use `withN(actualN)` when constructing output paths from the generated or loaded graph's actual vertex count.

When changing network types or path formats, compare Java's `SwitchUtils` with `notebooks/utils/switch.py` and, when present, `test-notebooks/utils/switch.py`. Check parameter names and order, decimal formatting, and randomization suffixes. These implementations are not fully synchronized: the `notebooks/` helper lacks `DirectedER`, and Java/Python decimal rounding can differ. Java path formatting also currently depends on the default locale. Do not assume parity without checking representative inputs and rounding boundaries.

## Testing Guidelines

Tests use JUnit Jupiter. Name test classes `*Test` and place them in the matching package under `app/src/test/java`. Focus tests on graph invariants, argument validation, reproducibility, and simulator edge cases. Run `./gradlew test` before submitting changes; use targeted tests during development, for example:

```bash
./gradlew :app:test --tests sim.network.topology.SchwartzDirectedSFTest
```

## Commit & Pull Request Guidelines

Git history follows concise conventional commits such as `feat: add bond percolation simulation`, `docs: update README`, and `chore: update SAR configs`. Use English, imperative subjects, and keep them under about 72 characters.

The repository uses `main` as its integration branch. Open pull requests against `main` unless instructed otherwise. Include a clear description, linked issue or experiment context, commands run, and notes about generated outputs or configuration changes. Add plots, notebook output, or screenshots only when they clarify scientific or user-facing results.
