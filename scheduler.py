import subprocess

def humanize(num):
    """A utility function to help generate human readable number string"""
    if not isinstance(num, int) and num.isdigit():
        num = int(num)
    for unit in ['', 'K', 'M']:
        if num % 1000:
            return '%d%s' % (num, unit)
        else:
            num /= 1000
    return "%d%s" % (num, 'G')

class Scheduler(object):
    """Dsub scheduler construction and execution."""
    NO_JOB = 'NO_JOB'
    job_pool = []
    def __init__(self, tool, conf):
        self.tool = tool
        if self.tool == 'dsub':
            self.cmd = 'dsub'
            self.add_argument('--provider', 'google-v2')
            self.add_argument('--project', conf['Platform']['project'])
            if conf.get('Platform', 'zones', fallback=None):
                self.add_argument('--zones', conf['Platform']['zones'])
            else:
                self.add_argument('--regions', conf['Platform']['regions'])

    def add_argument(self, argname, value=None):
        """Add one argument to cmd string."""
        if value:
            self.cmd += ' ' + argname + ' ' + value
        else:
            self.cmd += ' ' + argname

    def run(self):
        """Run cmd as subprocess and add it in job_pool."""
        print self.cmd
        job_id = subprocess.check_output(self.cmd, shell=True)
        job_id = job_id.strip()
        if job_id != Scheduler.NO_JOB:
            Scheduler.job_pool.append(job_id)

    def run_after(self):
        """Run cmd as subprocess after processes in job_pool have finished.

        This subprocess will not be blocked, and allows multiple subprocesses
        run in parallel.

        Returns:
            A process desc itself.
        """
        if Scheduler.job_pool:
            jobs = ' '.join(Scheduler.job_pool)
            self.add_argument('--after', jobs)
        return subprocess.Popen(self.cmd, shell=True)
