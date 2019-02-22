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

    def __init__(self, name, cpu, mem):
        self.name = name
        self.cpu = int(cpu)
        if mem:
            self.mem = float(mem)
        else:
            self.mem = None

    def get_core(self):
        return str(self.cpu)

class GCP_Instance(Instance):
    @staticmethod
    def get_machine_types(region, cpu_list, min_mem):
        valid, invalid = [], []
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
        conn.commit()
        SQL = '''SELECT DISTINCT NAME, CPUS, MEMORY_GB FROM instance
WHERE ZONE LIKE ?
AND MEMORY_GB >= ?
AND CPUS = ? '''
        for cpu, mem in zip(cpu_list, min_mem):
            cur.execute(SQL, [region + '%', mem, cpu])
            entries = cur.fetchall()
            for entry in entries:
                valid.append(GCP_Instance(entry[0], entry[1], entry[2]))
        SQL = SQL.replace('>=', '<')
        for cpu, mem in zip(cpu_list, min_mem):
            cur.execute(SQL, [region + '%', mem, cpu])
            entries = cur.fetchall()
            for entry in entries:
                invalid.append(GCP_Instance(entry[0], entry[1], entry[2]))
        conn.close()
        return valid, invalid

    def __init__(self, name, cpu, mem):
        Instance.__init__(self, name, cpu, mem)
