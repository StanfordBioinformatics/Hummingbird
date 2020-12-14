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

It is recommended to use the ```--install-option="--prefix=$PREFIX_PATH"``` along with pip while installing Hummingbird. This would give users easy access to the sample configuration files present under conf/examples which the users might need to refer to while writing their own configuration file for their own pipeline. Otherwise you will be able to find the configuration files under ```<virtualenv_name>/lib/<python_ver>/site-packages/Hummingbird/conf/examples```

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
##### Sample run on Google Cloud
In this section we will walk you through how to run Hummingbird on Google cloud for a sample pipeline(BWA). The assumption is that you have already created a project(https://cloud.google.com/resource-manager/docs/creating-managing-projects#creating_a_project), installed the Google cloud SDK and granted credentials to the SDK. Along with the Google cloud essentials, you also need to have installed Hummingbird. We will be modifying the bwa.conf.json file which is under conf/examples.

1. Get a list of all projects by executing ```gcloud projects list```. Make a note of the project name of the project in which you want to execute Hummingbird. Add that to the ```project``` field under the ```Platform``` section
2. Identify which region you want all of the computing resources to be launched in. This would ideally be the same region you provided to the gcloud sdk while setting it up. https://cloud.google.com/compute/docs/regions-zones has more information about regions and zones
3. Create a new storage bucket with the instructions provided in https://cloud.google.com/storage/docs/creating-buckets#storage-create-bucket-console. You can either create a bucket using the cloud storage browser in the Google Cloud Console, or execute ```gsutil mb gs://<BUCKET_NAME>``` from the command line. If creating the bucket from the command line, provide the ```-p```(project name), ```-c```(storage class) and ```-l```(location) flags to have greater control over the creation of your bucket. Once the bucket is created add it to the  ```bucket``` field under the ```Platform``` section. Just provide the bucket name, the full path is not required
4. For a sample BWA run  we will be using the platinum fastq genomes which are publicly hosted on the ```genomics-public-data/platinum-genomes``` cloud storage bucket. The two fastq files we will be using are ERR194159_1.fastq.gz and ERR194159_2.fastq.gz. In the ```input``` field under ```Downsample``` add ```gs://genomics-public-data/platinum-genomes/fastq/ERR194159_1.fastq.gz``` to INPUT_R1 and ```gs://genomics-public-data/platinum-genomes/fastq/ERR194159_2.fastq.gz``` to INPUT_R2. For your own input files provide the full Google cloud bucket path including gs://. The inputs need to be specified in key-value format. The key will be used to interpret the value later on. For example in your command you can refer to the first fastq file as ${INPUT_R1}
5. ```fractions``` represents the extent to which the whole input will be downsampled. You can keep the values as is, or tinker around with it to get different results.
6. After this we will be adding the output and logging bucket names to the configuration file. The output and logging buckets will be created under the bucket created in step 3. You will only need to provide the bucket path relative to the bucket created in step 3. For example if you created a bucket called bwa-example in step 3 and then you created bwa under bwa-example, and then under bwa you created bwa-logging and bwa-output, then you need to provide bwa/bwa-output to the output field and bwa/bwa-logging to the logging field.
7. The ```fullrun``` field indicates whether the input will be downsampled or not. Keep it as false to run the downsampled input on the pipeline. Setting it to true would mean that the entire pipeline will be run on the whole input
8. In the ```image``` field under ```Profiling``` provide the contianer image that contains the pipeline you wish to execute Hummingbird on
9. In the ```logging``` field provide a bucket where Hummingbird will write the log files that are generated during the profiling step. This should be different from the logging bucket you provided under the Downsample field. It should be relative to the bucket created in step 3
10. In the ```result``` field provide a bucket that will store the profiling results. It should be relative to the bucket created in step 3
11. threads?
12. ```input-recursive``` is where you need to provide any additional files that will be needed during execution. For example if you have your reference files under the ```references/GRCh37lite``` bucket(relative to the bucket created in step 3) then you can mention it in the ```input-recursive``` field with a key such as ```REF```
13. In the ```command``` field provide the command that is to be executed in the container. Use the keys that were mentioned in the input field and the input-recursive field(if any).
14. The output file name and path can be mentioned in the ```output``` field. It should be relative to the bucket created in step 3
15. Once the configuration file is created you can execute Hummingbird by executing ```hummingbird <path to conf file>```

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
    - `fractions` (optional) A list of decimals representing the downsample size. The default list is `[0.0001, 0.01, 0.1]` which means the sample will be downsized to 0.1%, 1%, 10% of its original size.
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
