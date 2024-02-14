#! python
"""Created on Sun Mar 30 13:25:29 2014.

@author: Mao
"""
import contextlib

if __name__ == "__main__":
    import os
    import sys

    logfile = sys.argv[1] if len(sys.argv) > 1 else "install-vaspy.txt"
    file = open(logfile)
    for line in file:
        line = line.rstrip("\n")
        with contextlib.suppress(FileNotFoundError):
            os.remove(line)
        with contextlib.suppress(OSError):
            os.removedirs(os.path.dirname(line))
