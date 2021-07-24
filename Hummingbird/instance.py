#!/usr/bin/env python3

import sqlite3
import subprocess

from retry import retry


class Instance:
    @staticmethod
    def get_machine_types(conf, min_mem):
        service = conf['Platform']['service']
        service = service.lower()
        cpu_list = conf['Profiling'].get('thread', [2])
        if service in ['google', 'gcp']:
            region = conf['Platform']['regions']
            return GCPInstance.get_machine_types(region, cpu_list, min_mem)
        elif service == 'aws':
            return AWSInstance.get_machine_types(AWSInstance.get_instance_families(), cpu_list, min_mem)
        elif service in ['azure', 'az']:
            return AzureInstance.get_machine_types(conf, cpu_list, min_mem)

    def __init__(self, name, cpu, mem):
        self.name = name
        self.cpu = int(cpu)
        if mem:
            self.mem = float(mem)
        else:
            self.mem = None

    def __eq__(self, other):
        if self is other:
            return True
        if type(self) != type(other):
            return False
        if self.name == other.name and self.cpu == other.cpu and self.mem == other.mem:
            return True
        else:
            return False

    def __hash__(self):
        return hash((self.name, self.cpu, self.mem))

    def __repr__(self):
        return self.name

    def get_core(self):
        return str(self.cpu)

    @staticmethod
    def check_threads_supported(threads, service):
        if service in ['google', 'gcp']:
            return GCPInstance.check_threads_supported(threads)
        elif service == 'aws':
            return AWSInstance.check_threads_supported(threads)
        elif service in ['azure', 'az']:
            return AzureInstance.check_threads_supported(threads)

    @staticmethod
    def get_supported_threads(service):
        if service in ['google', 'gcp']:
            return GCPInstance.get_supported_threads()
        elif service == 'aws':
            return AWSInstance.get_supported_threads()
        elif service in ['azure', 'az']:
            return AzureInstance.get_supported_threads()


class GCPInstance(Instance):
    pricing = {
        'n1-standard-1': 0.0475,
        'n1-standard-2': 0.0950,
        'n1-standard-4': 0.1900,
        'n1-standard-8': 0.3800,
        'n1-standard-16': 0.7600,
        'n1-standard-32': 1.5200,
        'n1-standard-64': 3.0400,
        'n1-highmem-2': 0.1184,
        'n1-highmem-4': 0.2368,
        'n1-highmem-8': 0.4736,
        'n1-highmem-16': 0.9472,
        'n1-highmem-32': 1.8944,
        'n1-highmem-64': 3.7888,
        'n1-highcpu-2': 0.0709,
        'n1-highcpu-4': 0.1418,
        'n1-highcpu-8': 0.2836,
        'n1-highcpu-16': 0.5672,
        'n1-highcpu-32': 1.1344,
        'n1-highcpu-64': 2.2688
    }

    @staticmethod
    def get_machine_types(region, cpu_list, min_mem):
        valid, invalid = [], []
        output = subprocess.check_output('gcloud compute machine-types list', shell=True)
        conn = sqlite3.connect(':memory:')
        cur = conn.cursor()
        cur.execute('''CREATE TABLE instance
                    (NAME text, ZONE text, CPUS INTEGER, MEMORY_GB real, DEPRECATED text)''')
        machine_types = output.decode("utf-8").splitlines()  # decode to convert bytes array to string in python3
        machine_types = [line.split() for line in machine_types]
        # Deprecated filed is empty for returned output
        machine_types = [line + [None] for line in machine_types if len(line) == 4]
        cur.executemany("INSERT INTO instance VALUES (?,?,?,?,?)", machine_types)
        conn.commit()
        sql_statement = '''SELECT DISTINCT NAME, CPUS, MEMORY_GB FROM instance
WHERE NAME LIKE ?
AND ZONE LIKE ?
AND MEMORY_GB >= ?
AND CPUS = ? '''
        for cpu, mem in zip(cpu_list, min_mem):
            cur.execute(sql_statement, ['n1%', region + '%', mem, cpu])
            entries = cur.fetchall()
            for entry in entries:
                valid.append(GCPInstance(entry[0], entry[1], entry[2]))
        sql_statement = sql_statement.replace('>=', '<')
        for cpu, mem in zip(cpu_list, min_mem):
            cur.execute(sql_statement, ['n1%', region + '%', mem, cpu])
            entries = cur.fetchall()
            for entry in entries:
                invalid.append(GCPInstance(entry[0], entry[1], entry[2]))
        conn.close()
        return valid, invalid

    def set_price(self, price=None):
        if price:
            self.price = price
        elif self.name in GCPInstance.pricing:
            self.price = GCPInstance.pricing[self.name]
        else:
            raise Exception('Fail to set price.')

    def __init__(self, name=None, cpu=None, mem=None):
        if name is None and (cpu is None or mem is None):
            Instance.__init__(self, 'n1-standard-1', 1, 3.75)
        elif name:
            _, family, vcpu = name.split('-')
            if family == 'standard':
                multiplier = 3.75
            elif family == 'highmem':
                multiplier = 7.5
            elif family == 'highcpu':
                multiplier = 0.9
            else:
                raise Exception('Unsupported instance family')
            Instance.__init__(self, name, vcpu, int(vcpu) * multiplier)
        else:
            Instance.__init__(self, 'custom', cpu, mem)

    @staticmethod
    def check_threads_supported(threads):
        return set(threads).issubset(GCPInstance.get_supported_threads())

    @staticmethod
    def get_supported_threads():
        return set([int(key.rpartition('-')[-1]) for key in GCPInstance.pricing.keys()])


class AWSInstance(Instance):
    thread_suffix = {1: '.micro', 2: '.large', 4: '.xlarge', 8: '.2xlarge', 16: '.4xlarge', 32: '.8xlarge',
                     36: '.9xlarge', 64: '.16xlarge'}
    # TODO dynamically fetch instance pricing. e.g.:
    #   aws pricing get-products --service-code AmazonEC2 \
    #       --filters "Type=TERM_MATCH,Field=instanceType,Value=m5.xlarge" \
    #       "Type=TERM_MATCH,Field=location,Value=US East (N. Virginia)"
    pricing = {
        't2.micro': 0.0116,
        'r5.large': 0.126,
        'r5.xlarge': 0.252,
        'r5.2xlarge': 0.504,
        'r5.4xlarge': 1.008,
        'r5.8xlarge': 2.016,
        'r5.16xlarge': 4.032,
        'm5.large': 0.096,
        'm5.xlarge': 0.192,
        'm5.2xlarge': 0.384,
        'm5.4xlarge': 0.768,
        'm5.8xlarge': 1.536,
        'm5.16xlarge': 3.072,
        'c5.large': 0.085,
        'c5.xlarge': 0.17,
        'c5.2xlarge': 0.34,
        'c5.4xlarge': 0.68,
        'c5.9xlarge': 1.53,
        'c5a.large': 0.077,
        'c5a.xlarge': 0.154,
        'c5a.2xlarge': 0.308,
        'c5a.4xlarge': 0.616,
        'c5a.8xlarge': 1.232,
        'c5a.16xlarge': 2.464,
        'c5d.large': 0.096,
        'c5d.xlarge': 0.192,
        'c5d.2xlarge': 0.384,
        'c5d.4xlarge': 0.768,
        'c5d.9xlarge': 1.728,
        'c5ad.large': 0.086,
        'c5ad.xlarge': 0.172,
        'c5ad.2xlarge': 0.344,
        'c5ad.4xlarge': 0.688,
        'c5ad.8xlarge': 1.376,
        'c5ad.16xlarge': 2.752,
        'i3.large': 0.156,
        'i3.xlarge': 0.312,
        'i3.2xlarge': 0.624,
        'i3.4xlarge': 1.248,
        'i3.8xlarge': 2.496,
        'i3.16xlarge': 4.992,
    }

    def __init__(self, name=None, cpu=None, mem=None):
        if name is None and (cpu is None or mem is None):
            Instance.__init__(self, 't2.micro', 1, 1)
        elif name:
            vcpu, mem = AWSInstance.desc_instance(name)
            Instance.__init__(self, name, vcpu, mem)
        else:
            Instance.__init__(self, 'custom', cpu, mem)

    def set_price(self, price=None):
        if price:
            self.price = price
        elif self.name in AWSInstance.pricing:
            self.price = AWSInstance.pricing[self.name]
        else:
            raise Exception('Fail to set price.')

    @staticmethod
    def desc_instance(name):
        import boto3
        ec2_client = boto3.client('ec2')
        response = ec2_client.describe_instance_types(InstanceTypes=[name])
        desc = response['InstanceTypes'][0]
        return desc['VCpuInfo']['DefaultVCpus'], desc['MemoryInfo']['SizeInMiB'] / 1024

    @staticmethod
    def get_machine_types(families, cpu_list, min_mem):
        valid, invalid = [], []
        for family in families:
            for cpu, mem in zip(cpu_list, min_mem):
                instance_type = AWSInstance.match_instance_type(family, cpu)
                if not AWSInstance.is_instance_valid(instance_type):
                    continue

                ins = AWSInstance(instance_type)
                if ins.mem >= mem:
                    valid.append(ins)
                else:
                    invalid.append(ins)
        return valid, invalid

    @staticmethod
    def match_instance_type(family, cpu):
        instance_type = family + AWSInstance.thread_suffix.get(cpu, '')
        # if the current instance is not supported for the thread specified, get the next one
        # this can happen for instances such as c5d.9xlarge where there is no c5d.8xlarge
        if not AWSInstance.is_instance_valid(instance_type) and cpu < 36:
            instance_type = AWSInstance.match_instance_type(family, cpu + 1)

        return instance_type

    @staticmethod
    def check_threads_supported(threads):
        return set(threads).issubset(AWSInstance.get_supported_threads())

    @staticmethod
    def get_supported_threads():
        return set(AWSInstance.thread_suffix.keys())

    @staticmethod
    def get_instance_families():
        return set([ins.split('.')[0] for ins in AWSInstance.pricing.keys() if ins != 't2.micro'])

    @staticmethod
    def is_instance_valid(instance_type):
        return instance_type and instance_type in AWSInstance.pricing


class AzureInstance(Instance):
    pricing = {
        'Standard_B1s': 0.0113,
        'Standard_E2_v3': 0.126,
        'Standard_E4_v3': 0.252,
        'Standard_E8_v3': 0.504,
        'Standard_E16_v3': 1.008,
        'Standard_E32_v3': 2.016,
        'Standard_D2_v3': 0.096,
        'Standard_D4_v3': 0.192,
        'Standard_D8_v3': 0.384,
        'Standard_D16_v3': 0.768,
        'Standard_D32_v3': 1.536,
    }

    machine_thread_mapping = {
        1: ['Standard_B1s'],
        2: ['Standard_E2_v3', 'Standard_D2_v3'],
        4: ['Standard_E4_v3', 'Standard_D4_v3'],
        8: ['Standard_E8_v3', 'Standard_D8_v3'],
        16: ['Standard_E16_v3', 'Standard_D16_v3'],
        32: ['Standard_E32_v3', 'Standard_D32_v3'],
    }

    def __init__(self, conf, machine=None, name=None, cpu=None, mem=None):
        self.conf = conf
        if machine:
            vcpu, mem = self.get_machine_specs(machine)
            super(AzureInstance, self).__init__(name, vcpu, mem)
        elif name:
            vcpu, mem = self.desc_instance(name, conf['Platform']['location'])
            super(AzureInstance, self).__init__(name, vcpu, mem)
        else:
            super(AzureInstance, self).__init__('Standard_E4_v3', cpu or 4, mem or 32)

    def set_price(self, price=None):
        if price:
            self.price = price
        elif self.name in self.pricing:
            self.price = self.pricing[self.name]
        else:
            raise Exception('Fail to set price.')

    def desc_instance(self, name, location):
        desc = self.filter_machine(self.conf, location, name)
        return self.get_machine_specs(desc)

    @staticmethod
    def get_machine_specs(machine):
        return int(machine['numberOfCores']), int(machine['memoryInMB']) / 1024

    @staticmethod
    def get_machine_types(conf, cpu_list, min_mem):
        valid, invalid = [], []
        for cpu, mem in zip(cpu_list, min_mem):
            cpu = int(cpu)
            if cpu not in AzureInstance.machine_thread_mapping:
                continue

            for name in AzureInstance.machine_thread_mapping[cpu]:
                ins = AzureInstance(conf, name=name)
                if ins.mem >= mem:
                    valid.append(ins)
                else:
                    invalid.append(ins)
        return valid, invalid

    @staticmethod
    @retry(tries=3, delay=1)
    def filter_machines(conf, location, machine_names):
        from azure.identity import AzureCliCredential
        from azure.mgmt.compute import ComputeManagementClient

        machine_names = set(machine_names)
        credential = AzureCliCredential()
        subscription_id = conf['Platform']['subscription']
        compute_client = ComputeManagementClient(credential, subscription_id)

        machines = []
        for vm in compute_client.virtual_machine_sizes.list(location):
            if vm.name in machine_names:
                machines.append(vm.serialize())
        return machines

    @staticmethod
    def filter_machine(conf, location, machine_name):
        return AzureInstance.filter_machines(conf, location, [machine_name])[0]

    @staticmethod
    def check_threads_supported(threads):
        return set(threads).issubset(AzureInstance.get_supported_threads())

    @staticmethod
    def get_supported_threads():
        return set(AzureInstance.machine_thread_mapping.keys())
