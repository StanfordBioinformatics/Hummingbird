# Hummingbird: Efficient Performance Prediction for Executing Genomic Applications in the Cloud

## Overview

Hummingbird is a Python framework that gives you a variety of optimum instance configurations to run your favorite genomics pipeline on cloud platforms.

The framework takes as input the necessary information required to run a cloud job and outputs different instance configurations that the user can use to run his/her pipeline on the cloud. The user can choose from a variety of instance configurations, such as the fastest, the cheapest, and the most efficient. What these configurations exactly mean will be explained later on in this README.

The unique thing about Hummingbird is that it takes the input files, downsamples them and then runs the whole pipeline on the downwsampled files. This results in users getting the resulting configurations in a short amount of time, compared to the time it would have taken to run the entire pipeline on the whole input file for different configurations and then give the user different optimum configurations.

As of now Hummingbird supports Google Cloud (GCP) and Amazon Web Service (AWS) and we hope to add other cloud providers in the future.

## Getting Started

### Installation Instructions

Based on which cloud service provider you want to use, Hummingbird can be installed using
```
pip install Hummingbird[gcp]
```
for Google cloud, and
```
pip install Hummingbird[aws]
```
for AWS

Hummingbird requires pip and python 2.7 or python 3 to be installed before installation.

It is highly recommended to use a [virtual environment](https://packaging.python.org/guides/installing-using-pip-and-virtual-environments/) to isolate the execution environment. Follow the instructions from the previous link to create a virtual environment, and then activate it:
```
source <virtual-environment-name>/bin/activate
```

#### Getting started on Google Cloud
Have [Google Cloud SDK](https://cloud.google.com/sdk/docs/quickstarts) installed and run:
```
gcloud init
```
This will set up your default project and grant credentials to the Google Cloud SDK. Also provide credentials so dsub can call Google APIs:
```
gcloud auth application-default login
```

#### Getting started on AWS Batch
Install [AWS CLI](https://docs.aws.amazon.com/cli/latest/userguide/cli-chap-install.html) and configure:
```
aws configure
```
It will ask for `Access key ID` and `Secret access key`. This credential will be used for all resources on AWS. See more instructions [here](https://docs.aws.amazon.com/cli/latest/userguide/cli-configure-quickstart.html).

### Editing the configuration file

Hummingbird has a `conf` folder which contains configuration files for all tested pipelines. The configuration file has a naming convention `<pipeline-name>.conf.json`, and contains all the information needed to launch jobs on cloud. The format of the configuration file is very important and any unnecessary lines or spaces will throw an error. We will go through each required field in the configuration file below.

1. `Platform` Specifies information about the cloud computing platform.
    - `service` The cloud computing service. Specify`gcp`for Google cloud, and `aws` for AWS Batch.
    - `project` The cloud project id. Make sure the project has access to all needed functionalities and APIs.
    - `regions` The region where the computing resource is hosted.
    - `bucket` The name of the cloud storage bucket where all the log and output files generated by Hummingbird will be stored.

2. `Downsample` Options for Hummingbird's downsampling processes.
    - `input` The full gs/s3 path of files to be used as input to Hummingbird. Specify in key-value pairs format, the keys will be used as interpolation of values later.
    - `target` The number of reads in the original input files. This number will be used for prediction purpose.
    - `output` Path to a directory in your bucket to store output files for downsample. Do not include the bucket name.
    - `logging` (GCP only) Path to a directory in your bucket to store log files for downsample. Do not include the bucket name.
    - `size` (optional) A list of decimals representing the downsample size. The default list is `[0.0001, 0.01, 0.1]` which means the sample will be downsized to 0.1%, 1%, 10% of its original size.
    - `fullrun` (optional) Default to `false`. Set to `true` to run the whole input without downsampling.
    - `index` (optional) Default to `false`. Use `samtools` to generate index files for the downsampled input files.

3. `Profiling` Options for Hummingbird's Profiling processes.
    - `image` The Docker image on which your pipeline will be executed. For AWS backend, it requires you to build a customized image. See documentation [here](AWS/README.md).
    - `logging` (GCP only) Path to a directory in your bucket to store log files for profiling. Do not include the bucket name.
    - `result` Path to a directory in your bucket to store result files for profiling. Do not include the bucket name.
    - `thread` A list of numbers representing number of virtual CPUs on a machine. Default is [8]. It will cause Hummingbird to test the downsampled inputs on a 8 virtual CPUs machine.
    - WDL/Cromwell
        - `wdl_file` Path to the workflow wdl_file in your bucket to be submitted by Hummingbird to Cromwell. Do not include the bucket name in the path.
        - `backend_conf` Path to the workflow backend configuration file in your bucket to be submitted by Hummingbird to Cromwell. Do not include the bucket name in the path.
        - `json_input` Two dimensional array containing the input for each Cromwell call. For each value in the `thread` option, Hummingbird requires a list referencing each downsampled file specified in the `size` option. These input files vary based on pipelines.
    - Command line tool
        - `command` Command directly executed in the image.
        - `input` and/or `input-recursive` Add any additional input resource in key-value pairs format.
    - `output` and/or `output-recursive` Path in your bucket to where Hummingbird will output the memory and time profiling results. Specify in key-value pairs format, the keys will be used as interpolation of values later. Do not include the bucket name.
    - `force` (optional) Default to `false`. Set to `true` to force to re-execute pipeline even the result exists.
    - `tries` (optional) Default to 1. Specify the number of tries for each task, the result will be reported as average of multiple tries.
    - `disk` (optional) Default to 500. The size in GB for your data disk size of your instance.

### Executing Hummingbird

Once the pre-requisites have been installed and the chosen configuration file has been modified properly, you can execute Hummingbird.

To execute Hummingbird run the following command:
```
python hummingbird.py [options] <path to your configuration file>
```
For example:
```
python hummingbird.py conf/bwa.conf.json
```
Hummingbird has two options:

1. `--fa_downsample` (optional) specifies the tool used to downsample the input files. Choose between seqtk and zless, default is seqtk.

1. `-p` or `--profiler` (optional) specifies the profiling tool used to monitor memory and runtime information. Default is `time` which uses /usr/bin/time on local backend.

During execution, Hummingbird will first downsample the input file and place the outputs in the bucket. At this point, Hummingbird will ask for your input to continue. Enter 'N' to stop Hummingbird so you can configure the input json files for your pipeline using the newly downsampled files. You will need to write a separate input file for each thread (change cpu count to match each thread) and downsample size. Then upload these files to the Google cloud path specified in the `json_input` section of your configuration file.

Now start Hummingbird again using the same command. The previously downsampled files will be saved, so it should be quick. Type 'y' when the same prompt appears, and Hummingbird will begin profiling. No further user input is required. The expected runtime will depend on your pipeline and downsample sizes chosen.

Please keep in mind Hummingbird will download your input files, downsample them, and then upload the downsampled files to the Google bucket folder mentioned in the `output` field of `Downsample` in the configuration file, so make sure that there is enough space locally in the machine started by dsub to download the input files. In addition, make sure the boot disk is large enough to support your docker image. The default boot disk size is 50GB disk size is 1000GB.

### Result

For each stage of the pipeline, Hummingbird will print 4 different configurations on the terminal: the fastest, the cheapest, and the most efficient. The meaning of each term will be more clear in the following lines.

1. The fastest: This configuration indicates the fastest amongst all the different configurations that Hummingbird executed the pertaining stage of the pipeline on. Here fastest means the configuration with the least execution time.

2. The cheapest: The instance which costs the least is chosen as the instance with the cheapest configuration.

3. The most efficient: The inctance which has the most cost-efficient value, maximizing the computing power of a unit of spending.
