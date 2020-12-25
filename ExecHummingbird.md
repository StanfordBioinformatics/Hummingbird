## Section 4: Executing Hummingbird

Once the pre-requisites have been installed and the chosen configuration file has been modified properly, you can execute Hummingbird.

To execute Hummingbird run the following command:
```
hummingbird [options] <path to your configuration file>
```
For example:
```
hummingbird conf/bwa.conf.json
```
Hummingbird has two options:

1. `--fa_downsample` (optional) specifies the tool used to downsample the input files. Choose between seqtk and zless, default is seqtk.

1. `-p` or `--profiler` (optional) specifies the profiling tool used to monitor memory and runtime information. Default is `time` which uses /usr/bin/time on local backend.

During execution, Hummingbird will first downsample the input file(s) and place the downsampled input file(s) in the bucket. At this point, Hummingbird will ask for your input in order to continue. Enter 'N' to stop Hummingbird so you can configure the input json files for your pipeline using the newly downsampled files. You will need to write a separate input file for each thread (change cpu count to match each thread) and downsample size. Then upload these files to the Google cloud path specified in the `json_input` section of your configuration file.

Now start Hummingbird again using the same command. The previously downsampled files will be saved, so it should be quick. Type 'y' when the same prompt appears, and Hummingbird will begin profiling. No further user input is required. The expected runtime will depend on your pipeline and downsample sizes chosen.

Please keep in mind that Hummingbird will download your input files, downsample them, and then upload the downsampled files to the Google bucket folder mentioned in the `output` field of `Downsample` in the configuration file, so make sure that there is enough space locally in the machine started by dsub to download the input files. In addition, ensure that the boot disk is large enough to support your docker image. The default boot-disk size is 50GB and disk size is 1000GB.
