#! python
"""
Created on Sun Mar 30 13:25:29 2014

@author: Mao
"""

if __name__ == "__main__":
    import os
    import sys

    logfile = sys.argv[1] if len(sys.argv) > 1 else "install-vaspy.txt"
    file = open(logfile)
    for line in file:
        line = line.rstrip("\n")
        try:
            os.remove(line)
        except FileNotFoundError:
            pass
        try:
            os.removedirs(os.path.dirname(line))
        except OSError:
            pass
