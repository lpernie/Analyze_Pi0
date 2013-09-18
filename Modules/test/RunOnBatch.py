#!/usr/bin/env python

import subprocess, time, sys, os
from Method import *

pwd         = os.getcwd()
eosPath     = '/store/caf/user/lpernie' 
#queue    = 'cmscaf1nd'
queue    = '8nh'
OnlyContCorr = False
nInter = -1
useES = 'True'
Gamma_MVA = 'False'

#cut
cuts4s9 = 0.8
isGun    = True

if isGun:
   isaGun       = 'True'
   Gamma_MVA    = False #!
   OnlyContCorr = False  #!
   dirname      = 'Prova_01'#!
   if OnlyContCorr : 
      inputlist_n = pwd + '/InputPi0Gun_group.txt'
      ijobmax     = 15
   else:
      # Se vuoi il TTree dell'MVA
      if Gamma_MVA:
         inputlist_n = pwd + '/InputPi0Gun_groupGAMMA.txt'
      else:
         # Se vuoi girare su tutti i Pi0
         inputlist_n = pwd + '/InputPi0Gun_group.txt'#!
         #inputlist_n = pwd + '/InputPi0Gun_group_short.txt'#!
      ijobmax     = 25#! 7 per MVA solo
else:
   inputlist_n = pwd + '/InputPi0Alca_Short.txt'
   isaGun      = 'False'
   dirname     = 'NormaleAlca_02_noES'
   useES       = 'True'
   ijobmax     = 1
   nInter      = 1000000

workdir  = pwd + '/' + dirname
srcPath  = workdir + '/src/'
cfgPath  = workdir + '/cfg/'

#-------- check if you have right access to queues --------#
checkAccessToQueues = subprocess.Popen(['bjobs'], stderr=subprocess.PIPE, shell=True);
output = checkAccessToQueues.communicate()[1]
if(output.find('command not found')==-1):
    print "[calib] Correct setup for batch submission"
else:
    print "[calib] Missing access to queues"
    sys.exit(1)

#-------- create folders --------#

print "[calib] Creating local folders (" + dirname + ")"

folderCreation = subprocess.Popen(['mkdir -p ' + workdir ], stdout=subprocess.PIPE, shell=True);
folderCreation.communicate()
folderCreation = subprocess.Popen(['mkdir -p ' + srcPath ], stdout=subprocess.PIPE, shell=True);
folderCreation.communicate()
folderCreation = subprocess.Popen(['mkdir -p ' + cfgPath ], stdout=subprocess.PIPE, shell=True);
folderCreation.communicate()
print "[calib] Creating folders on EOS"
folderCreation = subprocess.Popen(['cmsMkdir ' + eosPath + '/' + dirname ], stdout=subprocess.PIPE, shell=True);
folderCreation.communicate()

#-------- fill cfg files --------#

# open list of input files
inputlist_f = open( inputlist_n )
# read the list containing all the input files
inputlistbase_v = inputlist_f.readlines()

inputlist_v = inputlistbase_v[:]
ijob=0
# Creating different list for hadd
NrelJob = float(len(inputlist_v)) / float(ijobmax)
if( float(int(NrelJob) - NrelJob) < 0. ):
    NrelJob = int(NrelJob) + 1

print "[calib] Total number of files to be processed: " , len(inputlistbase_v)
print "[calib] Creating cfg Files"

hhaddSrc_Gun = workdir + "/finalHadd_Gun.list"
hhaddSrc_n = workdir + "/finalHadd.list"
FinalHaddList = open( hhaddSrc_n, 'w')
FinalHaddList_Gun = open( hhaddSrc_Gun, 'w')
while (len(inputlist_v) > 0):
    # List for Hadd of Conteinment Correction
    if isGun:
       FinalHaddList.write("root://eoscms//eos/cms" + eosPath + "/" + dirname + "/ContCorr_"  + str(ijob) + ".root\n")
       FinalHaddList_Gun.write("root://eoscms//eos/cms" + eosPath + "/" + dirname + "/LocalPi0Gun_"  + str(ijob) + ".root\n")
    else:
       FinalHaddList.write("root://eoscms//eos/cms" + eosPath + "/" + dirname + "/LocalPi0Alca_"  + str(ijob) + ".root\n")
    # create CFG file
    fill_cfg_n = cfgPath + "config_" + str(ijob) + ".py"
    fill_cfg_f = open( fill_cfg_n, 'w' )
    printFillCfg( Gamma_MVA, fill_cfg_f, str(ijob), workdir, isaGun, OnlyContCorr, str(nInter), useES, str(cuts4s9) )
    # loop over the names of the input files to be put in a single cfg
    lastline = min(ijobmax,len(inputlist_v)) - 1
    for line in range(min(ijobmax,len(inputlist_v))):
        ntpfile = inputlist_v.pop(0)
        ntpfile = ntpfile.rstrip()
        if ntpfile != '':
            if(line != lastline):
                fill_cfg_f.write("        '" + ntpfile + "',\n")
            else:
                fill_cfg_f.write("        '" + ntpfile + "'\n")
    fill_cfg_f.write("    )\n")
    fill_cfg_f.write(")\n")
    fill_cfg_f.close()

    # print SRC file
    fillSrc_n = srcPath + "config_" + str(ijob) + ".sh"
    fillSrc_f = open( fillSrc_n, 'w')
    if (isGun):
       source_s1 = "/tmp/LocalPi0Gun_" + str(ijob) + ".root"
    else:
       source_s1 = "/tmp/LocalPi0Alca_" + str(ijob) + ".root"
    source_s2 = "/tmp/ContCorr_" + str(ijob) + ".root"
    destination_s = eosPath + '/' + dirname + '/'
    printSubmitSrc(fillSrc_f, fill_cfg_n, source_s1, source_s2, destination_s, OnlyContCorr, pwd, isGun)
    fillSrc_f.close()

    # make the source file executable
    changePermission = subprocess.Popen(['chmod 777 ' + fillSrc_n], stdout=subprocess.PIPE, shell=True);
    debugout = changePermission.communicate()

    env_script_n = 'source ' + srcPath + "config_" + str(ijob) + ".sh"
    submit_s = "bsub -q " + queue + " -o " + workdir + "/prova_" + str(ijob) + ".log " + env_script_n 
    submitJobs = subprocess.Popen([submit_s], stdout=subprocess.PIPE, shell=True);
    output = (submitJobs.communicate()[0]).splitlines()
    print "[calib]  '-- " + output[0]

    ijob = ijob+1

FinalHaddList.close()
FinalHaddList_Gun.close()

# checking number of running/pending jobs
checkJobs = subprocess.Popen(['bjobs -q ' + queue], stdout=subprocess.PIPE, shell=True);
datalines = (checkJobs.communicate()[0]).splitlines()

print 'Waiting for jobs to be finished...'
nData = 999
#Daemon cheking running jobs
while len(datalines)>2 :#>= stessa coda
    for entry in datalines:
        entry = entry.rstrip()
        entry = entry.split()[0]
        #print entry
        if(entry.find('JOBID')!=-1): continue
        i = int(entry)

    time.sleep(1)
    
    checkJobs = subprocess.Popen(['bjobs -q ' + queue], stdout=subprocess.PIPE, shell=True);
    datalines = (checkJobs.communicate()[0]).splitlines()
    if( datalines<nData ):
       nData = datalines
       print "N datalines: " + str(nData)

print "Done with all jobs! Now merge'em all!!!"
## Now The final Hadd  
if( isGun ):  
   if not ( OnlyContCorr ):
       hadd_s1 = 'hadd -f /tmp/LocalPi0Gun_TOT.root @' + hhaddSrc_Gun
       print '[hadd] :: ' + hadd_s1
       addFiles1 = subprocess.Popen([hadd_s1],stdout=subprocess.PIPE, shell=True)
       nFilesAdded1 = 0
       filesAdded1 = (addFiles1.communicate()[0]).splitlines()
       print 'Now staging LocalPi0Gun_TOT.root on EOS'
       stage_s1 = 'cmsStage -f /tmp/LocalPi0Gun_TOT.root ' + eosPath + '/' + dirname
       print stage_s1
       stageEpsilonFile1 = subprocess.Popen([stage_s1], stdout=subprocess.PIPE, shell=True);
       print stageEpsilonFile1.communicate()

   else:
      hadd_s = 'hadd -f /tmp/ContCorr_TOT.root @' + hhaddSrc_n
      print '[hadd] :: ' + hadd_s
      addFiles = subprocess.Popen([hadd_s],stdout=subprocess.PIPE, shell=True)
      nFilesAdded = 0
      filesAdded = (addFiles.communicate()[0]).splitlines()
      print 'Now staging ContCorr_TOT.root on EOS'
      stage_s = 'cmsStage -f /tmp/ContCorr_TOT.root ' + eosPath + '/' + dirname
      print stage_s
      stageEpsilonFile = subprocess.Popen([stage_s], stdout=subprocess.PIPE, shell=True);
      print stageEpsilonFile.communicate()
else:
   hadd_s = 'hadd -f /tmp/LocalPi0Alca_TOT.root @' + hhaddSrc_n
   print '[hadd] :: ' + hadd_s
   addFiles = subprocess.Popen([hadd_s],stdout=subprocess.PIPE, shell=True)
   nFilesAdded = 0
   filesAdded = (addFiles.communicate()[0]).splitlines()
   print 'Now staging ContCorr_TOT.root on EOS'
   stage_s = 'cmsStage -f /tmp/LocalPi0Alca_TOT.root ' + eosPath + '/' + dirname
   print stage_s
   stageEpsilonFile = subprocess.Popen([stage_s], stdout=subprocess.PIPE, shell=True);
   print stageEpsilonFile.communicate()
print 'Congrats... This is The End!'
