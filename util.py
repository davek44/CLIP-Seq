#!/usr/bin/env python
import operator, os, subprocess, time

############################################################
# util
#
# Helpful methods that are difficult to categorize.
############################################################

############################################################
# condorify
############################################################
def condorify(cmds):
    return ['runCmd -c "%s"' % c for c in cmds]


############################################################
# exec_par
#
# Execute the commands in the list 'cmds' in parallel, but
# only running 'max_proc' at a time.
############################################################
def exec_par(cmds, max_proc, print_cmd=False):
    total = len(cmds)
    finished = 0
    running = 0
    p = []

    if max_proc == 1:
        while finished < total:
            op = subprocess.Popen(cmds[finished], shell=True)
            os.waitpid(op.pid, 0)
            finished += 1

    else:
        while finished + running < total:
            # launch jobs up to max
            while running < max_proc and finished+running < total:
                if print_cmd:
                    print cmds[finished+running]
                p.append(subprocess.Popen(cmds[finished+running], shell=True))
                #print 'Running %d' % p[running].pid
                running += 1

            # are any jobs finished
            new_p = []
            for i in range(len(p)):
                if p[i].poll() != None:
                    running -= 1
                    finished += 1
                else:
                    new_p.append(p[i])

            # if none finished, sleep
            if len(new_p) == len(p):
                time.sleep(1)
            p = new_p

        # wait for all to finish
        for i in range(len(p)):
            p[i].wait()


############################################################
# sort_dict
#
# Sort a dict by the values, returning a list of tuples
############################################################
def sort_dict(hash, reverse=False):
    return sorted(hash.items(), key=operator.itemgetter(1), reverse=reverse)

