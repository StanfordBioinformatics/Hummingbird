from Hummingbird.hummingbird_utils import bcolors


class SchedulerException(Exception):
    def __init__(self, message):
        self.message = message
        super().__init__(self.message)

    def __str__(self):
        return bcolors.FAIL + self.message + bcolors.ENDC
