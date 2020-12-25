# Hummingbird on AWS Batch

For jobs specified using AWS, Hummingbird will launch jobs using AWS Batch. You will need to prepare for a customized docker image including the AWS CLI for your applications. This tutorial will include the instructions and preparation steps for your jobs.

## Build your Docker image
This repository provides a sample `Dockerfile` and a `fetch_and_run.sh` script to help you create the customized Docker image. We will go through line by line for each step.

```
FROM broadinstitute/gatk:4.1.3.0
```
This specified the base image of your applications.

```
RUN curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"
RUN unzip awscliv2.zip
RUN ./aws/install
```
These commands install AWS CLI. Depend on your base image, you may need to install differently than this. Refer to [AWS CLI installing user guide](https://docs.aws.amazon.com/cli/latest/userguide/cli-chap-install.html) for more details.

```
RUN apt-get update && apt-get install -y \
    bwa \
    git

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
