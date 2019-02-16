from __future__ import print_function
import subprocess
import sqlite3

class Instance:
    @staticmethod
    def get_machine_types(conf, min_mem):
        service = conf['Platform']['service']
        service = service.lower()
        region = conf['Platform']['regions']
        cpu_list = conf['Profiling']['thread'].split(',')
        if service == 'google' or service == 'gcp':
            return GCP_Instance.get_machine_types(region, cpu_list, min_mem)

    def __init__(self, cpu, mem):
        self.cpu = cpu
        self.mem = mem

class GCP_Instance(Instance):
    @staticmethod
    def get_machine_types(region, cpu_list, min_mem):
        output = subprocess.check_output('gcloud compute machine-types list', shell=True)
        conn = sqlite3.connect(':memory:')
        cur = conn.cursor()
        cur.execute('''CREATE TABLE instance
                    (NAME text, ZONE text, CPUS INTEGER, MEMORY_GB real, DEPRECATED text)''')
        machine_types = output.splitlines()
        machine_types = [line.split() for line in machine_types]
        # Deprecated filed is empty for returned output
        machine_types = [line + [None] for line in machine_types if len(line) == 4]
        cur.executemany("INSERT INTO instance VALUES (?,?,?,?,?)", machine_types)
        SQL = '''SELECT DISTINCT NAME, MEMORY_GB FROM instance
WHERE ZONE LIKE ?
AND MEMORY_GB >= ?
AND CPUS = ? '''
        for cpu, mem in zip(cpu_list, min_mem):
            cur.execute(SQL, [region + '%', mem, cpu])
        valid = cur.fetchall()
        SQL.replace('>=', '<')
        for cpu, mem in zip(cpu_list, min_mem):
            cur.execute(SQL, [region + '%', mem, cpu])
        invalid = cur.fetchall()
        return valid, invalid

    def __init__(self, cpu, mem):
        Instance.__init__(self, cpu, mem)
