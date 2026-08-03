# Repository Guidelines

## Project Structure & Module Organization

This repository is a Java 21 Gradle project for directed-network SAR and percolation simulations, with supporting C and Python analysis code.

- `app/src/main/java/sim/`: Java application code. Main entry points are `App.java`, `SAR.java`, `GraphGen.java`, and `BondPercolation.java`.
- `app/src/main/java/sim/network/`: CSR-style graph representation and topology generators.
- `app/src/main/java/sim/simulation/`: SAR and bond-percolation simulators and result types.
- `app/src/main/java/sim/utils/`: shared utilities and topology dispatch logic.
- `app/src/test/java/`: JUnit tests, mirroring Java package names.
- `c-lang/ebcm/`: C EBCM numerical programs.
- `notebooks/`: Jupyter notebooks and Python EBCM helpers.
- `scripts/`: convenience wrappers for common Java and C runs.

Generated outputs go under `build/`, `app/out/`, or `out/`; keep them out of commits unless explicitly required.

## Build, Test, and Development Commands

- `./gradlew build`: compile, test, and assemble the project.
- `./gradlew test`: run all JUnit tests.
- `./gradlew :app:runSAR`: run the SAR simulation using `SAR.SimulationConfig`.
- `./gradlew :app:runGraphGen`: generate edge lists for later simulations.
- `./gradlew :app:runBondPercolation`: run bond-percolation simulations.
- `./gradlew run`: run the `sim.App` demo entry point.
- `bash scripts/run-c.sh`: compile and run the selected C EBCM program.

Use the Gradle wrapper rather than a system Gradle installation.

## Coding Style & Naming Conventions

Java code uses 4-space indentation, package-private tests, and descriptive `UpperCamelCase` class names. Methods, fields, and local variables use `lowerCamelCase`; constants use `UPPER_SNAKE_CASE`. Prefer deterministic seeds for simulation code and tests. Add new topology generators under `sim.network.topology` and wire them through `SwitchUtils.generateGraph` and `SwitchUtils.buildNetworkPath`.

## Testing Guidelines

Tests use JUnit Jupiter. Name test classes `*Test` and place them in the matching package under `app/src/test/java`. Focus tests on graph invariants, argument validation, reproducibility, and simulator edge cases. Run `./gradlew test` before submitting changes; use targeted tests during development, for example:

```bash
./gradlew :app:test --tests sim.network.topology.SchwartzDirectedSFTest
```

## Commit & Pull Request Guidelines

Git history follows concise conventional commits such as `feat: add bond percolation simulation`, `docs: update README`, and `chore: update SAR configs`. Use English, imperative subjects, and keep them under about 72 characters.

Open pull requests against `develop` unless instructed otherwise. Include a clear description, linked issue or experiment context, commands run, and notes about generated outputs or configuration changes. Add plots, notebook output, or screenshots only when they clarify scientific or user-facing results.
