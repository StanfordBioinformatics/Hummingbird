import subprocess
import logging

class Scheduler(object):
    """Dsub scheduler construction and execution."""
    def __init__(self, tool, conf):
        self.tool = tool
        if self.tool == 'dsub':
            self.cmd = 'dsub'
            self.add_argument('--provider', 'google-v2')
            self.add_argument('--project', conf['Platform']['project'])
            if 'zones' in conf['Platform']:
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
        """Run cmd as subprocess and return the Popen object."""
        logging.debug(self.cmd)
        return subprocess.Popen(self.cmd, shell=True)
