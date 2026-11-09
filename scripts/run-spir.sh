#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
if [[ $# -ne 0 ]]; then
    printf 'Usage: bash scripts/run-spir.sh\n' >&2
    exit 2
fi
./gradlew :app:build
./gradlew :app:runSPIR
