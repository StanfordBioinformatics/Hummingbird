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

1. Create the Hummingbird Cloudformation Stack in the target AWS Account and Region.

   [![Launch Stack](https://cdn.rawgit.com/buildkite/cloudformation-launch-stack-button-svg/master/launch-stack.svg)](
   https://console.aws.amazon.com/cloudformation/home?#/stacks/new?stackName=hummingbird&templateURL=https://cf-templates-gvgta4w56y1c-us-west-2.s3.us-west-2.amazonaws.com/hummingbird-cloudformation.template)

2. Install [AWS CLI](https://docs.aws.amazon.com/cli/latest/userguide/cli-chap-install.html) and configure:
    ```
    aws configure
    ```

   It will ask for `Access key ID` and `Secret access key`. This credential will be used for all resources on AWS. See
   more instructions [here](https://docs.aws.amazon.com/cli/latest/userguide/cli-configure-quickstart.html).

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
