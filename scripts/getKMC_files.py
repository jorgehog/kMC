from os import system
from sys import argv, exit


if len(argv) < 5:
    print "Usage: ssh_id fromPath toPath r0:r1 ..."
    exit(1)

ssh_id = argv[1]
fromPath = argv[2]
toPath = argv[3]
ranges = " ".join(argv[4:])

subranges = ranges.split()

for subrange in subranges:

    bottom, top = [int(x) for x in subrange.split(":")]

    for n in range(bottom, top + 1):
        call = "scp %s:%s/kMC%d.lmp %s" % (ssh_id, fromPath, n, toPath)
        
        system(call)
        
