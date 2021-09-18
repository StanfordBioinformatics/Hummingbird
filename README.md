<img src="https://github.com/StanfordBioinformatics/Hummingbird/blob/main/docs/hummingbird_2_2x.png" width="200" align="right">

# Hummingbird: Efficient Performance Prediction for Executing Genomic Applications in the Cloud #

## Overview

Hummingbird is a Python framework that gives a variety of optimum instance configurations to run your favorite genomics pipeline on cloud platforms.

The input for this framework is the necessary information required to run a cloud job and it generates different instance configurations that the user can use to run the pipeline on the cloud. The user can choose from a variety of instance configurations, such as the fastest, the cheapest, and the most efficient. The detailed explanation on these configurations can be found in the latter section of this README.

The unique feature about Hummingbird is that it takes the input files, downsamples them, runs the whole computational pipeline on these downwsampled files and subsequently provides the user with different optimum instance configurations. Therefore, the users obtain the resulting configurations in a short amount of time compared to a run on the entire pipeline with the whole input file(s) for different instance configurations.

Currently, Hummingbird supports Google Cloud (GCP), Amazon Web Service (AWS) and Microsoft Azure, and we hope to add other cloud providers in the future.

## Installation Instructions

Hummingbird can be installed using
```
pip install CloudHummingbird
```

It is recommended to use the ```--install-option="--prefix=$PREFIX_PATH"``` along with pip while installing Hummingbird. This would give users easy access to the sample configuration files located in conf/examples which the users might need to refer to while writing their own configuration file(s) for their own computational pipeline. Alternatively, the configuration files can be found here: ```<virtualenv_name>/lib/<python_ver>/site-packages/Hummingbird/conf/examples```

Hummingbird requires pip and python 3 as prerequesites for installation.

It is highly recommended to use a [virtual environment](https://packaging.python.org/guides/installing-using-pip-and-virtual-environments/) to isolate the execution environment. Please follow the instructions from the above link to create a virtual environment, and then activate it:
```
source <virtual-environment-name>/bin/activate
```

## Section 1: [Getting Started on Google Cloud, AWS Batch, Azure Batch](./docs/GettingStarted.md)
This section explains how to get started on Google Cloud, AWS and Azure.

## Section 2: [Sample Run on Google Cloud](./docs/SampleRun.md)
This section provides instructions to execute a sample run of BWA on Google Cloud using Hummingbird

## Section 3: [Editing the Configuration File](./docs/EditConf.md)
This section provides information about the configuration file and how to edit it

## Section 4: [Executing Hummingbird](./docs/ExecHummingbird.md)
This section provides information about how to execute Hummingbird

## Section 5: [Hummingbird Result](./docs/HummingbirdResult.md)
This section provides a guide to interpret the results provided by Hummingbird

## Section 6: [Using Different Input File Formats and Tools for Format Conversions](./docs/FormatConv.md)
This section provides a guide for users who want to leverage the downsampling step in Hummingbird but have input files in formats different than BAM or fastq/fastq.gz

## Section 7: [Alternative Downsampling Methods](./docs/AltDownsampling.md)
This section provides users a guide to alternative downsampling techniques other than the ones supported by Hummingbird

## Section 8: [Workflow Parser](./docs/WorkflowParser.md)
This section explains how Hummingbird parses workflows provided by the user

## Section 9: [Container Technology](./docs/ContainerTech.md)
This section explains how Hummingbird takes advantage of the container technology for execution

## Section 10: [I/O Profiling](./docs/IOProfiling.md)
This section explains how future versions of Hummingbird will profile I/O throughput as well

## Section 11: [Fault Tolerance](./docs/FaultTolerance.md)
This section describes the fault tolerant capabilities of Hummingbird

## Section 12: [Requirements for Running Hummingbird on Cloud Platform](./docs/CloudProviderRequirements.md)
This section lists all required components for running Hummingbird on a Cloud Platform provider.

* Logo Credit: Camille Berry
