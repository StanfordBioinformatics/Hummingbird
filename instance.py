from __future__ import print_function
from __future__ import unicode_literals
import json
import subprocess
import sqlite3

class Instance:
    @staticmethod
    def get_machine_types(conf, min_mem):
        service = conf['Platform']['service']
        service = service.lower()
        region = conf['Platform']['regions']
        cpu_list = conf['Profiling'].get('thread', [4])
        if service == 'google' or service == 'gcp':
            return GCP_Instance.get_machine_types(region, cpu_list, min_mem)
        elif service == 'aws':
            return AWS_Instance.get_machine_types(['r4', 'r5'], cpu_list, min_mem)

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

class GCP_Instance(Instance):
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
                valid.append(GCP_Instance(entry[0], entry[1], entry[2]))
        SQL = SQL.replace('>=', '<')
        for cpu, mem in zip(cpu_list, min_mem):
            cur.execute(SQL, ['n1%', region + '%', mem, cpu])
            entries = cur.fetchall()
            for entry in entries:
                invalid.append(GCP_Instance(entry[0], entry[1], entry[2]))
        conn.close()
        return valid, invalid

    def set_price(self, price=None):
        if price:
            self.price = price
        elif self.name in GCP_Instance.pricing:
            self.price = GCP_Instance.pricing[self.name]
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

class AWS_Instance(Instance):
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
            vcpu, mem = AWS_Instance.desc_instance(name)
            Instance.__init__(self, name, vcpu, mem)
        else:
            Instance.__init__(self, 'custom', cpu, mem)

    def set_price(self, price=None):
        if price:
            self.price = price
        elif self.name in AWS_Instance.pricing:
            self.price = AWS_Instance.pricing[self.name]
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
                ins = AWS_Instance(family+AWS_Instance.thread_suffix[cpu])
                if ins.mem >= mem:
                    valid.append(ins)
                else:
                    invalid.append(ins)
        return valid, invalid

if __name__ == "__main__":
    ins = AWS_Instance('t2.micro')
    print(ins.mem)
