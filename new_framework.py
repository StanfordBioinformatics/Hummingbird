import os
import time

#vcpu=[4,8,16,32,64,4,8,16,32,64,4,8,16,32,64]
#memory=[3.6,7.2,14.4,28.8,57.6,15,30,60,120,240,26,52,104,208,416]
vcpu=[8,16,32,8,16,32,8,16,32]
memory=[7.2,14.4,28.8,30,60,120,52,104,208]
#on the basis of us-east
#cost=[0.00003954,0.000078931,0.000157709,0.000315263,0.000630373,0.00005293,0.000105708,0.000211263,0.000422374,0.000844597,0.000065928,0.000131709,0.000263265,0.000526374,0.001052595]
cost=[0.000078931,0.000157709,0.000315263,0.000105708,0.000211263,0.000422374,0.000131709,0.000263265,0.000526374]
timing_high_cpu=[]
timing_standard=[]
timing_high_mem=[]
cost_high_cpu=[]
cost_standard=[]
cost_high_mem=[]
speedup_high_cpu=[]
speedup_standard=[]
speedup_high_mem=[]
normalized_speedup_high_cpu=[]
normalized_speedup_standard=[]
normalized_speedup_high_mem=[]
input_file=[]
errors=[]
timing_info=[]
total_cost=[]
job_name=[]
logging_info=[]
commands=[]
multi_threading=[]
dsub=[]
dsub_temp=""
num_files=0
flag=0
flag1=0

#parse the input path
f=open("./User_Provided_Input.txt","r")
for line in f:
	if "INPUT_PATH" in line:
		split1 = line.split('=')
                second_split1 = split1[1]
                second_split1 = second_split1[:-3]
                input_path = second_split1
                break
f.close()

#parse the project name
f=open("./User_Provided_Input.txt","r")
for line in f:
	if "project " in line:
		splitpr = line.split("project ")
                proj = splitpr[1]
                proje = proj.split(' ')
                project = proje[0]
                break
f.close()

#parse the logging bucket
f=open("./User_Provided_Input.txt","r")
for line in f:
        if "--logging" in line:
                splitlg = line.split("--logging ")
                lg = splitlg[1]
                lgg = lg.split(' ')
                log_bucket = lgg[0]
                break
f.close()

#parse the name of the input files
f=open("./User_Provided_Input.txt","r")
for line in f:
        if "Input_File" in line:
                num_files+=1
                equal2 = line.split('=')
                equals2 = equal2[1]
                equals2 = equals2[1:]
                equals2 = equals2[:-2]
                input_file.append(equals2)
		if ".fastq.gz" in input_file[0]:
			file_extension=".fastq.gz"
		elif ".fastq.gz.bz2" in input_file[0]:
			file_extension=".fastq.gz.bz2"
		else:
			file_extension=".bam"
f.close()

#parse the multi-threading option
f=open("./User_Provided_Input.txt","r")
for line in f:
	if "Multi-threaded" in line:
		multi_t=line.split('=')
                multi_thread=multi_t[1]
                multi_thread=multi_thread.strip('\n')
		if ',' in multi_thread:
			multi_threading=multi_thread.split(',')
		else:
			multi_threading.append(multi_thread)
f.close()

#parse the command line(s)
f=open("./User_Provided_Input.txt","r")
for line in f:
	if "command" in line:
		split_command = line.split('command')
                splits1=split_command[1]
                splits1 = splits1[:2] + "time "  + splits1[2:]
                line = split_command[0]+"command"+splits1[:-1]
		commands.append(line)
f.close()

#parse the dsub command excluding the command line(s)
f=open("./User_Provided_Input.txt","r")
for line in f:
        if "command" in line:
		break
	else:
		dsub_temp=dsub_temp+line[:-2]
f.close()

#add the command line(s) to the dsub command
for i in range(len(commands)):
	dsub.append(dsub_temp+commands[i])

'''
#download the input files, downsample them according to the file extension and upload them to the google cloud bucket
for i in range(num_files):
        os.system('gsutil cp '+input_path+'/'+input_file[i]+' .')
	if file_extension == ".fastq.gz":
		os.system('zless '+input_file[i]+' | head -n 1000000 > '+str(i+1)+file_extension)
	elif file_extension == ".fastq.gz.bz2":
		print "DOWNSAMPLE FOR BZ2"
	else:
		print "DOWNSAMPLE FOR BAM"
        os.system('gsutil cp '+str(i+1)+file_extension+' '+input_path)

	#move the input files one directory above
	a=input_path.split('/')
	y=""
	for j in range(len(a)):
        	if j == len(a)-1:
                	break
        	elif j == len(a)-2:
                	y=y+a[j]
        	else:
                	y=y+a[j]+"/"
	os.system('gsutil mv '+input_path+'/'+input_file[i]+' '+y)
	#print ('gsutil mv '+input_path+'/'+input_file[i]+' '+y)
'''
#replace the input files in the first dsub command with the downsampled files
for i in range(num_files):
	if input_file[i] in dsub[0]:
		dsub[0]=dsub[0].replace(input_file[i], str(i+1)+file_extension)

#execute dsub jobs for each command line given in the input text file
for i in range(len(commands)):

	job_name=[]
	logging_info=[]
	timing_info=[]
	total_cost=[]
	errors=[]
	flag=0
	timing_high_cpu=[]
	timing_standard=[]
	timing_high_mem=[]
	cost_high_cpu=[]
	cost_standard=[]
	cost_high_mem=[]
	speedup_high_cpu=[]
	speedup_standard=[]
	speedup_high_mem=[]
	normalized_speedup_high_cpu=[]
	normalized_speedup_standard=[]
	normalized_speedup_high_mem=[]

	#execute the dsub job	
	for j in range(len(memory)):
		find_ram=dsub[i].split('min-ram ')
	        ram_split=find_ram[1]
       		ram_space=ram_split.split(' ')
	        ram=ram_space[0]
	        dsub[i]=dsub[i].replace("min-ram "+str(ram), "min-ram "+str(memory[j]))

	        find_cores=dsub[i].split('min-cores ')
        	cores_split=find_cores[1]
	        cores_space=cores_split.split(' ')
	        cores=cores_space[0]
        	dsub[i]=dsub[i].replace("min-cores "+str(cores), "min-cores "+str(vcpu[j]))

		#if command has multi-threading, replace the multi-threading argument by the number of VCPUs in the instance
	        if (multi_threading[i])!='NO':
			multit=dsub[i].split(' '+multi_threading[i]+' ')
			multith=multit[1]
			multithr=multith.split(' ')
			multithre=multithr[0]
			dsub[i]=dsub[i].replace(multi_threading[i]+" "+str(multithre), multi_threading[i]+" "+str(vcpu[j]))
		
		if i==(len(commands)-1):
			os.system(dsub[i]+"'"+" > dsub_information.txt")
		else:
			os.system(dsub[i]+" > dsub_information.txt")

		#get the name of the job and the log file
	        f=open("./dsub_information.txt","r")
	        for line in f:
	                job_name.append(line.strip('\n'))
	                log=line.split("--")
	                log_info=log[2]
       	 	        log_info=log_info.strip('\n')
               		logging_info.append(log_info)
		f.close()
	
	time.sleep(240)

	#check whether the job has finished successfully or not
	while True:
	        os.system("dstat --provider google --project "+project+" --jobs '"+job_name[flag]+"' --status '*' > dstat_info.txt")
       		f=open("./dstat_info.txt","r")
        	for line in f:
                	if "Success" in line:
                        	flag=flag+1
				#print i
			if "fail" in line:
				flag=flag+1
				print job_name[flag-1]+" was unsuccessful"
        	f.close()

        	if flag==(len(vcpu)):
                	break
	
		time.sleep(30)

	#once the job is done, save the time taken to a file and check the log file for errors
	for j in range(len(vcpu)):
		os.system("gsutil cp "+log_bucket+"/"+job_name[j]+"-stderr.log .")
	        os.system("tail -n 3 "+job_name[j]+"-stderr.log | head -n 1 > timing.txt")
		os.system("tail -n 10 "+job_name[j]+"-stderr.log > check_error.txt")	
		
		#f=open("./"+job_name[j]+"-stderr.log")
		f=open("./check_error.txt")
		for line in f:
                	if (("fail" in line.lower()) and ("memory" in line.lower())) or ("error" in line.lower()):
				#commonly occuring error present in all log files, but doesn't affect the job
                        	if "GPG error" in line:
                                	continue
                        	print line
                        	print "Log file "+job_name[j]+"-stderr.log seems to have an error"
                        	print "Skipping this configuration while providing the best configuration"
                        	if j in errors:
                                	continue
                        	else:
                                	errors.append(j)
	        f.close()
	
		#calculate the time taken by the machine to execute the job
		f=open("./timing.txt")
        	for line in f:
                	timin=line.split('\t')
                	timin[1]=timin[1].strip('\t')
                	timin[1]=timin[1].strip('\n')
                	timing=timin[1]
                	timing=timing[:-1]
                	timing_i=timing.split('m')
                	m_to_s=(float(timing_i[0]))*60
                	s=float(timing_i[1])
                	total_time=m_to_s+s
                	x=total_time*float(cost[j])
			
			if "-Xms" in dsub[i]:
				temp_mem=dsub[i].split("-Xms")
				java_mem=temp_mem[1][0]
				if (int(java_mem)) > (memory[j]-2):
					total_cost.append(0)
                                	timing_info.append(0)
					break
			#if there are errors in the log file, do not take that instance into account during cost and timing calculation
                	if j in errors:
                        	total_cost.append(0)
                        	timing_info.append(0)
                	else:
                        	total_cost.append(x)
                        	timing_info.append(total_time)
        	f.close()

	#calculate the timing and cost for different machine types
	timing_high_cpu=timing_info[:int(len(timing_info))/3]
	while True:
        	if 0 in timing_high_cpu:
                	timing_high_cpu.remove(0)
        	else:
                	break
	#print "Timing High-CPU: ",timing_high_cpu

	cost_high_cpu=total_cost[:int(len(total_cost))/3]
	while True:
        	if 0 in cost_high_cpu:
                	cost_high_cpu.remove(0)
	        else:
        	        break
	#print "Cost High-CPU: ",cost_high_cpu

	timing_standard=timing_info[int(len(timing_info))/3:int(len(timing_info))/3+int(len(timing_info))/3]
	while True:
        	if 0 in timing_standard:
                	timing_standard.remove(0)
	        else:
        	        break
	#print "Timing Standard: ",timing_standard

	cost_standard=total_cost[int(len(total_cost))/3:int(len(total_cost))/3+int(len(total_cost))/3]
	while True:
        	if 0 in cost_standard:
                	cost_standard.remove(0)
	        else:
        	        break
	#print "Cost Standard: ",cost_standard

	timing_high_mem=timing_info[int(len(timing_info))/3+int(len(timing_info))/3:]
	while True:
        	if 0 in timing_high_mem:
                	timing_high_mem.remove(0)
	        else:
        	        break
	#print "Timing High-Mem: ",timing_high_mem

	cost_high_mem=total_cost[int(len(total_cost))/3+int(len(total_cost))/3:]
	while True:
        	if 0 in cost_high_mem:
                	cost_high_mem.remove(0)
	        else:
	                break
	#print "Cost High-Mem: ",cost_high_mem
	
	#calculate normalized speedup for different instance types
	flag1=0
	for j in range(len(timing_high_cpu)):
	        x=timing_high_cpu[0]/timing_high_cpu[j]
	        speedup_high_cpu.append(x)
	        if j is 0:
	                normalized_speedup_high_cpu.append(1)
	        else:
	                x=speedup_high_cpu[j]/(2**flag1)
	                normalized_speedup_high_cpu.append(x)
	                flag1=flag1+1
	
	flag1=0
	for j in range(len(timing_standard)):
	        x=timing_standard[0]/timing_standard[j]
	        speedup_standard.append(x)
	        if j is 0:
	                normalized_speedup_standard.append(1)
	        else:
	                x=speedup_standard[j]/(2**flag1)
	                normalized_speedup_standard.append(x)
	                flag1=flag1+1
	
	flag1=0
	for j in range(len(timing_high_mem)):
	        x=timing_high_mem[0]/timing_high_mem[j]
	        speedup_high_mem.append(x)
        	if j is 0:
                	normalized_speedup_high_mem.append(1)
	        else:
        	        x=speedup_high_mem[j]/(2**flag1)
                	normalized_speedup_high_mem.append(x)
	                flag1=flag1+1

	#calculate maximum normalized speedup amongst different instances of the same machine type
	max_normalized_speedup_high_cpu=max(normalized_speedup_high_cpu)
	max_normalized_speedup_standard=max(normalized_speedup_standard)
	max_normalized_speedup_high_mem=max(normalized_speedup_high_mem)

	'''
	print "normalized_speedup_high_cpu: ",normalized_speedup_high_cpu
	print "normalized_speedup_standard: ",normalized_speedup_standard
	print "normalized_speedup_high_mem: ",normalized_speedup_high_mem
	'''

	max_normalized_speedup=max(max_normalized_speedup_high_cpu,max_normalized_speedup_standard,max_normalized_speedup_high_mem)

	#retrieve the actual configuration from the value of max normalized speedup
	if max_normalized_speedup == max_normalized_speedup_high_cpu:
        	idx1=normalized_speedup_high_cpu.index(max_normalized_speedup)
        	idx2=timing_info.index(timing_high_cpu[idx1])

	elif max_normalized_speedup == max_normalized_speedup_standard:
        	idx1=normalized_speedup_standard.index(max_normalized_speedup)
        	idx2=timing_info.index(timing_standard[idx1])
	else:
        	idx1=normalized_speedup_high_mem.index(max_normalized_speedup)
        	idx2=timing_info.index(timing_high_mem[idx1])

	idx1_high_cpu=normalized_speedup_high_cpu.index(max_normalized_speedup_high_cpu)
	idx2_high_cpu=timing_info.index(timing_high_cpu[idx1_high_cpu])
	
	idx1_standard=normalized_speedup_standard.index(max_normalized_speedup_standard)
	idx2_standard=timing_info.index(timing_standard[idx1_standard])
	
	idx1_high_mem=normalized_speedup_high_mem.index(max_normalized_speedup_high_mem)
	idx2_high_mem=timing_info.index(timing_high_mem[idx1_high_mem])

	#print str(total_cost[idx2_high_cpu])+" "+str(total_cost[idx2_standard])+" "+str(total_cost[idx2_high_mem])
	cheapest_among_fastest=min(total_cost[idx2_high_cpu],total_cost[idx2_standard],total_cost[idx2_high_mem])
	idx_cheapest_among_fastest=total_cost.index(cheapest_among_fastest)
	
	fastest_high_cpu=min(timing_high_cpu)
	fastest_standard=min(timing_standard)
	fastest_high_mem=min(timing_high_mem)
	
	fastest=min(fastest_high_cpu,fastest_standard,fastest_high_mem)
	fastest_index=timing_info.index(fastest)
	
	cheapest_high_cpu=min(cost_high_cpu)
	cheapest_standard=min(cost_standard)
	cheapest_high_mem=min(cost_high_mem)
	
	cheapest=min(cheapest_high_cpu,cheapest_standard,cheapest_high_mem)
	
	cheapest_index=total_cost.index(cheapest)
	'''
	print "Speedup High-CPU: ",speedup_high_cpu
	print "Speedup Standard: ",speedup_standard
	print "Speedup High-Mem: ",speedup_high_mem
	'''
	#print "Max Normalized Speedup: ",max_normalized_speedup
	print "STAGE: "+str(i+1)
	print "Cheapest Configuration: "+str(vcpu[cheapest_index])+" cores & "+str(memory[cheapest_index])+" GB RAM. Time: "+str(timing_info[cheapest_index])
	print "Fastest Configuration in terms of normalized speedup: "+str(vcpu[idx2])+" cores & "+str(memory[idx2])+" GB RAM"
	print "Fastest Configuration in terms of execution time: "+str(vcpu[fastest_index])+" cores & "+str(memory[fastest_index])+" GB RAM. Time: "+str(timing_info[fastest_index])
	print "Fastest and Cheapest Configuration: "+str(vcpu[idx_cheapest_among_fastest])+" cores & "+str(memory[idx_cheapest_among_fastest])+" GB RAM. Time: "+str(timing_info[idx_cheapest_among_fastest])	
