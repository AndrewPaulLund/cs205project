#!/usr/bin/env python
"""
This code perform MPI parallization on the Harvard
Medical School compute cluster. 

This code was adapted and modified based on the sample
example code demonstrating parallization from:
https://github.com/jbornschein/mpi4py-examples

This code starts off a master node which will then compute
the bins needed for parallelization. Then the master node
will send instruction to the child nodes to tell it
what are the jobs it will run. Following processing,
the child nodes will send a signal back to the master
node to get it to send another job to it. Once
all tasks are completed, the MPI run will stop.

Author: Kar-Tong Tan
Email: ktan@g.harvard.edu
"""
from __future__ import print_function
from mpi4py import MPI
import sys
import os
import datetime

def enum(*sequential, **named):
    """
    Handy way to fake an enumerated type in Python
    http://stackoverflow.com/questions/36932/how-can-i-represent-an-enum-in-python
    """
    enums = dict(zip(sequential, range(len(sequential))), **named)
    return type('Enum', (), enums)

# Define MPI message tags
tags = enum('READY', 'DONE', 'EXIT', 'START')



def calculateBinsForJobs(chrom, chromSize, binSize, outputPrefix):
    '''
    Calculate the bins needed for splitting the jobs
    into the different workers. 
    '''
    binsList = []
    numberBins = int(chromSize) / int(binSize) + 1;
    
    for i in range(1, numberBins + 1):
        start       = 1 + (i - 1) * int(binSize)
        end         = i * int(binSize);
        region      = str(chrom) + ":" + str(start) + "-" + str(end);
        outputFile  = str(outputPrefix) + "_" + str(region) + ".bcf";
        bins = [start, end, region, outputFile]
        binsList.append(bins)
    
    return binsList
    


'''
Initializations and preliminaries. This 
'''
comm = MPI.COMM_WORLD   # get MPI communicator object
size = comm.size        # total number of processes
rank = comm.rank        # rank of this process
status = MPI.Status()   # get MPI status object



def master_function():
    '''
    Master process executes code below
    '''

    # Get the chromsize and the binsize from the command line arguments
    binSize     = sys.argv[1]
    bamFile     = sys.argv[2]
    refGenome   = sys.argv[3]
    numCores    = int(sys.argv[4])
    outputPrefix = sys.argv[5];
    chrom       = sys.argv[6];
    chromSize   = sys.argv[7]
    
	
    # Calculate the bins sizes and the jobs that need to be done
    # given the paramenters that were defined here.
    taskList = calculateBinsForJobs(chrom, chromSize, binSize, outputPrefix)
	
	
	#tasks = range(2*size)
	#tasks       = numCores
	#tasks       = len(taskList)
    task_index  = 0
    num_workers = numCores - 1
    closed_workers = 0
    allFiles = []
	
    print("Master starting with %d workers" % num_workers)
    startTime = datetime.datetime.now()
    

    '''
    Check if all workers have closed. This while loop 
    would constantly check if all workers have finished.
    If there are spare workers and spare task, start them.
    '''
    while closed_workers < num_workers:
        data = comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)
        source = status.Get_source()
        tag = status.Get_tag()
        if tag == tags.READY:
            # Worker is ready, so send it a task
            if task_index < len(taskList):
                taskToSend = taskList[task_index] + [refGenome, bamFile]
                print(taskToSend)

                currFile = taskList[task_index][3]
                allFiles.append(currFile)
				

                comm.send(taskToSend, dest=source, tag=tags.START)
                print("Sending task %d to worker %d" % (task_index, source))
                task_index += 1
            else:
                # We send a signal to the worker to exit if there is
                # no more tasks to work on.
                comm.send(None, dest=source, tag=tags.EXIT)
        elif tag == tags.DONE:
            results = data
            print(results)
            print("Got data from worker %d" % source)
        elif tag == tags.EXIT:
            print("Worker %d exited." % source)
            closed_workers += 1
	
	# Merge results from all the workers back into a single file
	outputFilesFil = []
	for file in allFiles:
        # Check if the file is empty. Add to arr if not.
		if os.path.exists("%s" %file):
			if(os.stat("%s" % file).st_size != 0):
				outputFilesFil.append(file)

    filesStrFil = " ".join(outputFilesFil)

    outputBcf = outputPrefix + ".combined.bcf";
    #os.system("bcftools concat %s > %s" %(filesStrFil, outputBcf))

    print("Master finishing")
    endTime = datetime.datetime.now()
    timeDiff = endTime - startTime
    finalTimeStr = timeDiff.seconds + float(timeDiff.microseconds) / 1000000
    print("Finished %s cores %s time" %(num_workers, finalTimeStr))



def worker_function():
    '''
    Worker processes execute code below
    '''
    name = MPI.Get_processor_name()
    print("I am a worker with rank %d on %s." % (rank, name))
    while True:
        comm.send(None, dest=0, tag=tags.READY)
        task = comm.recv(source=0, tag=MPI.ANY_TAG, status=status)
        tag = status.Get_tag()

        if tag == tags.START:
            # Do the work here
            #result = task**2
            start       = task[0]
            end         = task[1]
            region      = task[2]
            outputFile  = task[3]
            refGenome   = task[4]
            bamFile     = task[5]
            print("Starting analysis of: %s" %region)

            # Unpack the data sent from the master node and perform the task
            os.system("samtools mpileup -d 100000000000 -Buf %s -r %s %s | bcftools call -O b -c -v > %s" 
                %(refGenome, region, bamFile, outputFile))
            
            result = 'Finished run of %s' %region

            # Send a signal to the master node telling it
            # that it has finished the jobs.
            comm.send(result, dest=0, tag=tags.DONE)

        # We stop the MPI run for this worker if we receive the EXIT tag
        elif tag == tags.EXIT:
            break

    # Send no data back to the master node and tell
    # it that it is going to exit. As a fall back.
    # (Should not run)        
    comm.send(None, dest=0, tag=tags.EXIT)


if rank == 0:
	'''
	Master node work
	'''
	master_function()
else:
	'''
	Worker node work
	'''
	worker_function()

