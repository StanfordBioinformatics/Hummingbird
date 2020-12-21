#!/bin/bash

set -e

PATH="/bin:/usr/bin:/sbin:/usr/sbin:/usr/local/bin:/usr/local/sbin"
BASENAME="${0##*/}"

usage () {
  if [ "${#@}" -ne 0 ]; then
    echo "* ${*}"
    echo
  fi
  exit 2
}

# Standard function to print an error and exit with a failing return code
error_exit () {
  echo "${BASENAME} - ${1}" >&2
  exit 1
}

which az >/dev/null 2>&1 || error_exit "Unable to find Azure CLI executable."

# Check what environment variables are set
if [ -z "${BLOB_NAME}" ]; then
  usage "BLOB_NAME not set, unable to determine target blob name."
fi

if [ -z "${AZURE_STORAGE_ACCOUNT}" ]; then
  usage "AZURE_STORAGE_ACCOUNT not set. No storage account to download from."
fi

if [ -z "${AZURE_STORAGE_CONTAINER}" ]; then
  usage "AZURE_STORAGE_CONTAINER not set. No container to download from."
fi

if [ -z "${AZURE_STORAGE_CONNECTION_STRING}" ]; then
  usage "AZURE_STORAGE_CONNECTION_STRING not set. No connection string to use for download."
fi

script="$(pwd)/run-script.sh"

# Fetch and run a script
fetch_and_run_script () {
  az storage blob download \
    --container-name "${AZURE_STORAGE_CONTAINER}" \
    --file "${script}" \
    --name "${BLOB_NAME}" \
    --account-name "${AZURE_STORAGE_ACCOUNT}" \
    --connection-string "${AZURE_STORAGE_CONNECTION_STRING}" \
    || error_exit "Failed to download file from ${AZURE_STORAGE_CONTAINER}"

  chmod u+x "${script}" || error_exit "Failed to chmod script."
  ( exec "${script}" "${@}" ) || error_exit " Failed to execute script."
}

fetch_and_run_script "${@}"
exit 0
