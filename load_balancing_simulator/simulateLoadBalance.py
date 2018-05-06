#!/usr/bin/python

'''
Takes a list of runtime and simulates a queue
for a parallel jobs. Identify what is the best,
worst and randomize runtimes given a heterogeneous
input data.

Author : Kar-Tong Tan
Email  : ktan@g.harvard.edu
Date   : 6th May 2018
'''

import random
import sys


def readFile(file):
	'''
	Read a file and collect the runtimes
	from the file
	'''
	f = open(file, "r")
	timeList = []
	for line in f:
		lineArr = line.strip().split()
		time = lineArr[2]
		timeList.append(time)
	# Remove last two elements which gives the overall times
	timeList.pop()
	timeList.pop()
	
	return(timeList)


class worker:
	'''
	The worker class that simulates the single core that is
	doing the work itself.
	'''
	def __init__(self, jobDuration):
		self.jobDuration = jobDuration
		self.currentTime = 0

	def moveTime(self, timeToShift):
		'''
		Move time and check if job has finished yet
		'''
		self.currentTime += timeToShift
		
		

	def checkJobFinish(self):
		'''
		Check if this particular worker has finished
		it's job. Worker is deemed to have finished its
		job if the worker currentTime is longer than
		the total job duration defined for this job.
		'''

		if float(self.currentTime) > float(self.jobDuration):
			# Job has finished:
			#print "orange"
			return 1
		else:
			return 0


class manager:
	'''
	Overall manager managing the global time as well
	as the workers.
	'''
	def __init__(self, jobRuntimes, numProc, timeShift):
		self.currTime = 0
		self.numProc  = numProc
		self.workerList  = []
		self.jobRuntimes = jobRuntimes
		self.idleProcTime = 0
		self.totalProcTime = 0
		self.totalTime = 0
		self.timeShift = timeShift
		
		# Create collection of workers based on number
		# of processors
		for i in range(self.numProc):
			if len(self.jobRuntimes) > 0:
				runTimeNeededForJob = self.jobRuntimes.pop(0)
				newWorker = worker(runTimeNeededForJob)
				self.workerList.append(newWorker)
	
	def tickAndProceed(self):
		'''
		Central Manager that does the global processing of
		all the files.
		'''
		timeShift = 0.01
		self.moveTime(timeShift)

		# Check if workers have finished their job and
		# remove it if it has finished.
		self.checkWorkersFinishedJob()


		# Check if there is spare capacity to run new jobs and
		# if there is still jobs to do. If so, start new job
		while(self.checkSpareJobsTodo() and self.checkNumSpareWorkersAvailable() > 0):
			self.startNewjob()

		# Count proc time and unused time
		self.totalTime += timeShift
		self.totalProcTime += timeShift * int(self.numProc)
		self.idleProcTime += timeShift * (int(self.numProc) - len(self.workerList))

	
	def moveTime(self, timeshift):
		'''
		Move the global time. For both the manager
		as well as for the worker
		'''
		self.currTime += timeshift
		self.moveWorkerTimes(timeshift)
		
		
	def moveWorkerTimes(self, timeshift):
		'''
		Move the time of each of the workers that is
		being managed by the manager.
		'''
		for worker in self.workerList:
			worker.moveTime(timeshift)
	
	
	def checkWorkersFinishedJob(self):
		'''
		Check if worker has finished job. If worker
        has finished job. Remove it.
		'''
		workerListNew = []

		# Check which jobs have completed. Only add
		# uncompleted jobs to new worker list
		for i in range(len(self.workerList)):
			worker = self.workerList[i]
			finishStatus = worker.checkJobFinish()
			#print finishStatus
			if finishStatus == 1:
				pass
			else:
				workerListNew.append(worker)

		self.workerList = workerListNew
		
		
	def checkNumSpareWorkersAvailable(self):
		'''
		Check number of spare workers available to do work.
		'''
		spareWorkers = self.numProc - len(self.workerList)
		return spareWorkers
	
	
	def checkSpareJobsTodo(self):
		'''
		Check if there are still spare jobs that the manager
		need to complete.
		'''
		if len(self.jobRuntimes) > 0:
			return 1
		else:
			return 0
			

	def startNewjob(self):
		'''
		Start a new job and append it to the list of workers
		that is currently working.
		'''
		runTimeNeededForJob = self.jobRuntimes.pop()
		newWorker = worker(runTimeNeededForJob)
		self.workerList.append(newWorker)


	def checkAllJobsFinished(self):
		'''
		Check that there are no more jobs to do
		and that all workers have finished their jobs
		as well.
		'''
		#print self.checkSpareJobsTodo()
		if(self.checkSpareJobsTodo() == 0 and len(self.workerList) == 0):
			return 1
		else:
			return 0



'''
Read the input file with the timing for all the runs
'''
#time = readFile("runRNA2.cores1.binSize1000000.runtime.txt")
time = readFile(sys.argv[1])




def startAndRunSimulation(cores, time):
	'''
	Start running the actual simulation. Start
	a new manager and perform the simulation give the 
	different runtime list and different number of cores
	available.
	'''
	timeshift = 0.01
	count = 0

	theManager = manager(time, cores, timeshift)

	while(theManager.checkAllJobsFinished() != 1):
		theManager.tickAndProceed()
		count += 1

	pctIdleTime = float(theManager.idleProcTime / theManager.totalProcTime)

	return (theManager.idleProcTime, theManager.totalProcTime, theManager.totalTime, pctIdleTime)



def simulateOriginal(cores, time):
	'''
	Simulate run with the original order of 
	run times.
	'''
	result = startAndRunSimulation(cores, time)
	return result

def simulateBestCase(cores, time):
	'''
	Simulate run where the biggest tasks are fed in and
	run first in the analysis. We regard this method of
	load balancingas a type of 'Best case'
	'''
	timeSortBest = sorted(time)
	result = startAndRunSimulation(cores, timeSortBest)
	return result

def simulateWorstCase(cores, time):
	'''
	Simulate run where the smallest tasks are fed in and
	run first in the analysis. We regard this method
	of load balancing as a type of 'Worst case' scenario
	'''
	timeSortWorst = sorted(time, reverse=True)
	result = startAndRunSimulation(cores, timeSortWorst)
	return result

def simulateRandomCase(cores, time):
	'''
	We randomize the order in which the tasks are being run
	in the simulation. 
	'''
	timeSortRandom = time
	random.shuffle(timeSortRandom)
	result = startAndRunSimulation(cores, timeSortRandom)
	return result

 	

# Simulate for different number of cores
cores = [1,2,4,8,16,20,32,64,128,256]

for core in cores:
	'''
	Pass in the copied version of time and see what
	happens in each of the cases.
	'''
	original = simulateOriginal(core, time[:])
	bestCase = simulateBestCase(core, time[:])
	worstCase = simulateWorstCase(core, time[:])
	random1 = simulateRandomCase(core, time[:])
	random2 = simulateRandomCase(core, time[:])
	random3 = simulateRandomCase(core, time[:])
	ave1 = float(random1[0] + random2[0] +random3[0]) / 3
	ave2 = float(random1[1] + random2[1] +random3[1]) / 3
	ave3 = float(random1[2] + random2[2] +random3[2]) / 3
	ave4 = float(random1[3] + random2[3] +random3[3]) / 3
	randomData = [ave1, ave2, ave3, ave4]

	# Print out the results
	print "\t".join(map(str, ("original", core) + original))
	print "\t".join(map(str, ("bestCase", core) + bestCase))
	print "\t".join(map(str, ("worstCase", core) + worstCase))
	print "\t".join(map(str, ["random", core] + randomData))
