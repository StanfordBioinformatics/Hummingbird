#!/bin/bash

if [[ -n "$AWS_BATCH_JOB_ID" ]]; then
  source aws_fetch_and_run.sh "${@}"
  exit 0
fi

if [[ -n "$AZ_BATCH_ACCOUNT_NAME" ]]; then
  source azure_fetch_and_run.sh "${@}"
  exit 0
fi

echo "Invalid cloud provider (missing 'AZ_BATCH_ACCOUNT_NAME' and 'AWS_BATCH_JOB_ID' env vars." >&2
exit 1