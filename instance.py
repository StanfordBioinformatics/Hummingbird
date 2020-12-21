from __future__ import print_function
from __future__ import unicode_literals
import json
import subprocess
import sqlite3

from retry import retry


class Instance:
    @staticmethod
    def get_machine_types(conf, min_mem):
        service = conf['Platform']['service']
        service = service.lower()
        cpu_list = conf['Profiling'].get('thread', [4])
        if service in ['google', 'gcp']:
            region = conf['Platform']['regions']
            return GCPInstance.get_machine_types(region, cpu_list, min_mem)
        elif service == 'aws':
            return AWSInstance.get_machine_types(['r4', 'r5'], cpu_list, min_mem)
        elif service in ['azure', 'az']:
            location = conf['Platform']['location']
            return AzureInstance.get_machine_types(conf, location, cpu_list, min_mem)

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
        machine_types = output.decode("utf-8").splitlines() # decode to convert bytes array to string in python3
        machine_types = [line.split() for line in machine_types]
        # Deprecated filed is empty for returned output
        machine_types = [line + [None] for line in machine_types if len(line) == 4]
        cur.executemany("INSERT INTO instance VALUES (?,?,?,?,?)", machine_types)
        conn.commit()
        SQL = '''SELECT DISTINCT NAME, CPUS, MEMORY_GB FROM instance
WHERE NAME LIKE ?
AND ZONE LIKE ?
AND MEMORY_GB >= ?
AND CPUS = ? '''
        for cpu, mem in zip(cpu_list, min_mem):
            cur.execute(SQL, ['n1%', region + '%', mem, cpu])
            entries = cur.fetchall()
            for entry in entries:
                valid.append(GCPInstance(entry[0], entry[1], entry[2]))
        SQL = SQL.replace('>=', '<')
        for cpu, mem in zip(cpu_list, min_mem):
            cur.execute(SQL, ['n1%', region + '%', mem, cpu])
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


class AWSInstance(Instance):
    thread_suffix = {2: '.large', 4: '.xlarge', 8: '.2xlarge', 16: '.4xlarge', 32: '.8xlarge'}
    pricing = {'r4.large': 0.133,
               'r4.xlarge': 0.266,
               'r4.2xlarge': 0.532,
               'r4.4xlarge': 1.064,
               'r4.8xlarge': 2.128,
               'r5.large': 0.126,
               'r5.xlarge': 0.252,
               'r5.2xlarge': 0.504,
               'r5.4xlarge': 1.008,
               'r5.8xlarge': 2.016,
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
        output = subprocess.check_output(['aws', 'ec2', 'describe-instance-types', '--instance-types', name])
        desc = json.loads(output)['InstanceTypes'][0]
        return desc['VCpuInfo']['DefaultVCpus'], desc['MemoryInfo']['SizeInMiB'] / 1024

    @staticmethod
    def get_machine_types(families, cpu_list, min_mem):
        valid, invalid = [], []
        for family in families:
            for cpu, mem in zip(cpu_list, min_mem):
                ins = AWSInstance(family + AWSInstance.thread_suffix[cpu])
                if ins.mem >= mem:
                    valid.append(ins)
                else:
                    invalid.append(ins)
        return valid, invalid


AZURE_INSTANCE_PRICES = {
    'Standard_E2_v3': 0.148,
    'Standard_E4_v3': 0.296,
    'Standard_E8_v3': 0.56,
    'Standard_E16_v3': 1.12,
    'Standard_E32_v3': 2.24,
}


class AzureInstance(Instance):
    pricing = AZURE_INSTANCE_PRICES

    def __init__(self, conf, machine=None, name=None, cpu=None, mem=None):
        self.conf = conf
        if name is None:
            if not mem:
                mem = 16
            if not cpu:
                cpu = 2
            super(AzureInstance, self).__init__('Standard_E2_v3', cpu, mem)
        elif machine:
            vcpu, mem = self.get_machine_specs(machine)
            super(AzureInstance, self).__init__(name, vcpu, mem)
        elif name:
            vcpu, mem = self.desc_instance(name)
            super(AzureInstance, self).__init__(name, vcpu, mem)

    def set_price(self, price=None):
        if price:
            self.price = price
        elif self.name in self.pricing:
            self.price = self.pricing[self.name]
        else:
            raise Exception('Fail to set price.')

    def desc_instance(self, name, region):
        desc = self.filter_machines(self.conf, region, [name])[0]
        return self.get_machine_specs(desc)

    @staticmethod
    def get_machine_specs(machine):
        return int(machine['numberOfCores']), int(machine['memoryInMB']) / 1024

    @staticmethod
    def get_machine_types(conf, location, cpu_list, min_mem):
        valid, invalid = [], []
        all_machines = AzureInstance.filter_machines(conf, location, AZURE_INSTANCE_PRICES.keys())
        for machine in all_machines:
            for cpu, mem in zip(cpu_list, min_mem):
                ins = AzureInstance(conf, machine=machine, name=machine['name'], cpu=cpu, mem=mem)
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
