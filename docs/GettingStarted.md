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
1. Install [AWS CLI](https://docs.aws.amazon.com/cli/latest/userguide/cli-chap-install.html) and configure:
    ```
    aws configure
    ```

    It will ask for `Access key ID` and `Secret access key`. This credential will be used for all resources on AWS.
    See more instructions [here](https://docs.aws.amazon.com/cli/latest/userguide/cli-configure-quickstart.html).

2. Follow instructions [here](https://docs.aws.amazon.com/batch/latest/userguide/service_IAM_role.html) to create the
AWS Batch Service Role.

3. Follow instructions [here](https://docs.aws.amazon.com/batch/latest/userguide/instance_IAM_role.html) to create the
AWS ECS Instance Role. Additionally, make sure that the instance has read/write access to the Input/Output buckets.

4. Follow instructions [here](https://docs.aws.amazon.com/AmazonECS/latest/developerguide/task_execution_IAM_role.html)
to create the AWS ECS Task Execution Role.

5. Edit [compute_environment.json](./Hummingbird/AWS/compute_environment.json) to update the following sections:
    - `subnets` to a list of subnets in the target VPC that the AWS Batch instance can be provisioned in.
    - `securityGroupIds` to one or more Security Groups that define the ingress/egress rules. Hummingbird, will not require any ingress rules but may require TCP `0.0.0.0/0` on egress rules.
    - `instanceRole` to be replaced by the value from step #3 (IAM Role Name).
    - `<AWSServiceRoleForBatch>` to be replaced by value from step #2 (full IAM Role ARN).
6. Edit [job-definition.json](./Hummingbird/AWS/job-definition.json) to update the following sections:
    - `jobRoleArn` to be replaced by the value from step #4 (full IAM Role ARN).

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
