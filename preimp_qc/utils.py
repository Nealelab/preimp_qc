import time
import secrets
import humanize


def gen_uid(n=5):
    return secrets.token_urlsafe(n)


class TimeLogger:
    def __init__(self, label, verbose=True):
        self.label = label
        self.verbose = verbose
        self.start = None

    def __aenter__(self):
        self.start = time.time()
        if self.verbose:
            print(f"Started '{self.label}'")
        return self

    def __aexit__(self, exc_type, exc_val, exc_tb):
        assert self.start is not None
        duration = humanize.time.naturaltime(time.time() - self.start)
        if self.verbose:
            print(f"Finished '{self.label}' in {duration} seconds")
