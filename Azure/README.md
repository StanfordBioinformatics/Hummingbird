# Hummingbird on Azure Batch

For jobs specified using Azure, Hummingbird will launch jobs using Azure Batch.  
You will need to prepare for a customized Docker Image including the Azure CLI for your applications. This tutorial will include the instructions and preparation steps for your jobs.

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
Copy `fetch_and_run.sh` to the image and set it as the entrypoint for your docker image. `fetch_and_run.sh` script is essentially a helper script that fetch user's application script and execute it in the container. The original repo is available as https://github.com/awslabs/aws-batch-helpers/blob/master/fetch-and-run/fetch_and_run.sh.
