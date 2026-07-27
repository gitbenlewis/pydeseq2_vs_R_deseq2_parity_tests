#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
LOCAL_ENV_FILE="${LOCAL_ENV_FILE:-${REPO_ROOT}/config/local_env.sh}"

if [[ ! -f "${LOCAL_ENV_FILE}" ]]; then
  echo "[ERROR] Missing LOCAL_ENV_FILE: ${LOCAL_ENV_FILE}" >&2
  exit 1
fi

# shellcheck disable=SC1090
source "${LOCAL_ENV_FILE}"

: "${PYDESEQ2_REPO:?Set PYDESEQ2_REPO in ${LOCAL_ENV_FILE}}"
: "${PARITY_CONDA_ENV:?Set PARITY_CONDA_ENV in ${LOCAL_ENV_FILE}}"
: "${CONDA_BASE:?Set CONDA_BASE in ${LOCAL_ENV_FILE}}"

CONDA_SH="${CONDA_BASE}/etc/profile.d/conda.sh"
if [[ ! -f "${CONDA_SH}" ]]; then
  echo "[ERROR] Conda initialization script not found: ${CONDA_SH}" >&2
  exit 1
fi
if [[ ! -d "${PYDESEQ2_REPO}/pydeseq2" ]]; then
  echo "[ERROR] PyDESeq2 checkout not found: ${PYDESEQ2_REPO}" >&2
  exit 1
fi

# shellcheck disable=SC1090
source "${CONDA_SH}"
conda activate "${PARITY_CONDA_ENV}"

export PYTHONPATH="${PYDESEQ2_REPO}${PYTHONPATH:+:${PYTHONPATH}}"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export MPLCONFIGDIR="${REPO_ROOT}/.cache/matplotlib"

PYTHON_BIN="$(command -v python || true)"
if [[ -z "${PYTHON_BIN}" ]]; then
  echo "[ERROR] No Python interpreter found after activating ${PARITY_CONDA_ENV}." >&2
  exit 1
fi

mkdir -p "${SCRIPT_DIR}/logs" "${MPLCONFIGDIR}"
LOG_PATH="${SCRIPT_DIR}/logs/000_run_parity_$(date -u +%Y%m%dT%H%M%SZ).log"
exec > >(tee -a "${LOG_PATH}") 2>&1

cd "${REPO_ROOT}"

export PYTHONDONTWRITEBYTECODE=1

"${PYTHON_BIN}" -m pytest -q -p no:cacheprovider tests
"${PYTHON_BIN}" -m pytest -q -p no:cacheprovider \
  "${PYDESEQ2_REPO}/tests/test_transcript_length_normalization.py"
"${PYTHON_BIN}" "${SCRIPT_DIR}/run_parity.py"
