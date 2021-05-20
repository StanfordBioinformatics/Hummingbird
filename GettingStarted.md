## Section 1: Getting Started

### Getting started on Google Cloud

Have [Google Cloud SDK](https://cloud.google.com/sdk/docs/quickstarts) installed and run:
```
gcloud init
```
This will set up your default project and grant credentials to the Google Cloud SDK. Also, provide credentials so that dsub can call Google APIs:
```
gcloud auth application-default login
```

### Getting started on AWS Batch
Install [AWS CLI](https://docs.aws.amazon.com/cli/latest/userguide/cli-chap-install.html) and configure:
```
aws configure
```
It will ask for `Access key ID` and `Secret access key`. This credential will be used for all resources on AWS. See more instructions [here](https://docs.aws.amazon.com/cli/latest/userguide/cli-configure-quickstart.html).

Follow instructions [here](https://docs.aws.amazon.com/batch/latest/userguide/service_IAM_role.html) to create the
AWS Batch Service Role.

Follow instructions [here](https://docs.aws.amazon.com/batch/latest/userguide/instance_IAM_role.html) to create the
AWS ECS Instance Role. Additionally, make sure that the instance has read/write accesss to the Input/Output buckets.

### Getting started on Azure Batch
Install [Azure CLI](https://docs.microsoft.com/en-us/cli/azure/install-azure-cli) and login:
```bash
az login
```

Set the Subscription ID:
```bash
# Show all accounts
az account list --output table

# Set the Subscription ID or Name gathered from the table above (replace Example Subscription)
az account set --subscription "Example Subscription"
```
