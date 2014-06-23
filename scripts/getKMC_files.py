from os import system
from sys import argv, exit


if len(argv) < 5:
    print "Usage: ssh_id fromPath toPath r0:r1 ..."
    exit(1)

ssh_id = argv[1]
fromPath = argv[2]
toPath = argv[3]

for subrange in argv[4:]:
    
    splits = subrange.split(":")

    if len(splits) == 2:
        bottom, top = [int(x) for x in splits]
        step = 1
    else:
        bottom, top, step = [int(x) for x in splits]

    for n in range(bottom, top + 1, step):
        call = "scp %s:%s/kMC%d.lmp %s" % (ssh_id, fromPath, n, toPath)
        
        system(call)
        
