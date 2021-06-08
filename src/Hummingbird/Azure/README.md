# Hummingbird on Azure Batch

## Prerequisites
1. You will need an [Azure Subscription](https://portal.azure.com/) to deploy Cromwell on Azure.
2. You must have the proper [Azure role assignments](https://docs.microsoft.com/en-us/azure/role-based-access-control/overview) to deploy Cromwell on Azure.  To check your current role assignments, please follow [these instructions](https://docs.microsoft.com/en-us/azure/role-based-access-control/check-access).  You must have one of the following combinations of [role assignments](https://docs.microsoft.com/en-us/azure/role-based-access-control/built-in-roles):
   1. `Owner` of the subscription<br/>
   2. `Contributor` and `User Access Administrator` of the subscription
   3. `Owner` of the resource group. *Note: this level of access will result in a warning during deployment, and will not use the latest VM pricing data.</i>  [Learn more](/docs/troubleshooting-guide.md/#How-are-Batch-VMs-selected-to-run-tasks-in-a-workflow?).  Also, you must specify the resource group name during deployment with this level of access (see below).*
   4.  Note: if you only have `Service Administrator` as a role assignment, please assign yourself as `Owner` of the subscription.
3. Install the [Azure Command Line Interface (az cli)](https://docs.microsoft.com/en-us/cli/azure/?view=azure-cli-latest), a command line experience for managing Azure resources.
4. Run `az login` to authenticate with Azure.

For jobs specified using Azure, Hummingbird will launch jobs using Azure Batch.  
You will need to prepare for a customized Docker Image including the Azure CLI for your applications. This tutorial will include the instructions and preparation steps for your jobs.

## Use published Docker Image
See [Hummingbird GitHub Packages](https://hub.docker.com/u/cloudhummingbird) for latest Docker Images.

If you need to use your private container registry or require custom packages, see below.


## Build your Docker image
This repository provides a sample `Dockerfile` and a `fetch_and_run.sh` script to help you create the customized Docker image. We will go through line by line for each step.

```
FROM broadinstitute/gatk:4.1.3.0
```
This specified the base image of your applications.

```
RUN apt-get update && apt-get install -y \
    bwa \
    git \
    curl

RUN curl -sL https://aka.ms/InstallAzureCLIDeb | bash
```
These commands install Azure CLI. Depend on your base image, you may need to install differently than this. Refer to [AWS CLI installing user guide](https://docs.aws.amazon.com/cli/latest/userguide/cli-chap-install.html) for more details.

```
RUN git clone https://github.com/lh3/seqtk.git
WORKDIR seqtk
RUN make install
```
Install any other applications you may need for your job. This step is optional and up to your need.

```
ADD fetch_and_run.sh /usr/local/bin/fetch_and_run.sh
WORKDIR /tmp

ENTRYPOINT ["/usr/local/bin/fetch_and_run.sh"]
```
Copy `fetch_and_run.sh` to the image and set it as the entrypoint for your docker image. `fetch_and_run.sh` script is essentially a helper script that fetch user's application script and execute it in the container.
