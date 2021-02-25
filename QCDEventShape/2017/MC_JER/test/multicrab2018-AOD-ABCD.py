#!/usr/bin/env python
"""
This is a small script that does the equivalent of multicrab.
"""
import os
from optparse import OptionParser

from CRABAPI.RawCommand import crabCommand
from CRABClient.ClientExceptions import ClientException
from httplib import HTTPException

def getOptions():
	"""
	Parse and return the arguments provided by the user.
	"""
	usage = ("Usage: %prog --crabCmd CMD [--workArea WAD --crabCmdOpts OPTS]"
		 "\nThe multicrab command executes 'crab CMD OPTS' for each project directory contained in WAD"
		 "\nUse multicrab -h for help")
	parser = OptionParser(usage=usage)

	parser.add_option('-c', '--crabCmd',
		dest = 'crabCmd',
         	default = '',
        	help = "crab command",
         	metavar = 'CMD')

	parser.add_option('-w', '--workArea',
		dest = 'workArea',
		default = '',
		help = "work area directory (only if CMD != 'submit')",
		metavar = 'WAD')

	parser.add_option('-m', '--crabCmdOpts',
		dest = 'crabCmdOpts',
                default = '',
                help = "options for crab command CMD",
                metavar = 'OPTS')



	(options, arguments) = parser.parse_args()



	if arguments:
		parser.error("Found positional argument(s): %s." % (arguments))
	if not options.crabCmd:
		parser.error("(-c CMD, --crabCmd=CMD) option not provided.")
	if options.crabCmd != 'submit':
		if not options.workArea:
			parser.error("(-w WAR, --workArea=WAR) option not provided.")
		if not os.path.isdir(options.workArea):
			parser.error("'%s' is not a valid directory." % (options.workArea))
	return options

def main():
	options = getOptions()
# The submit command needs special treatment.
	if options.crabCmd == 'submit':
	#--------------------------------------------------------
	# This is the base config:
	#--------------------------------------------------------
		from CRABClient.UserUtilities import config
		config = config()
	
		config.General.requestName = None
		config.General.workArea = 'CRAB_L1_HCAL_2018ABCD'
		config.General.transferOutputs = True
		config.General.transferLogs = True
		config.JobType.pluginName = 'Analysis'
		config.JobType.psetName = 'jetStudy.py'
#		config.Data.ignoreLocality = True
		config.Data.inputDBS = 'global'
		config.JobType.maxMemoryMB = 2500
#		config.JobType.priority = 9999
		config.Data.splitting = 'EventAwareLumiBased'
#		config.Data.splitting = 'LumiBased'
		config.Data.unitsPerJob = 80000
		config.Data.useParent = True
		config.Data.inputDataset = None
#		config.Data.splitting = 'LumiBased'
#		config.Data.splitting = 'Automatic'
#		config.Data.unitsPerJob = 20
#		config.Data.totalUnits = 30
		config.Data.outputDatasetTag = None
		config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/ReReco/Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt'
#		config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
		config.Data.publication = False
                config.JobType.allowUndistributedCMSSW = True
		config.Site.storageSite = 'T2_TW_NCHC' 



	# Will submit one task for each of these input datasets.
		inputDatasets = [
			'/SingleMuon/Run2018A-17Sep2018-v2/AOD',
			'/SingleMuon/Run2018B-17Sep2018-v1/AOD',
			'/SingleMuon/Run2018C-17Sep2018-v1/AOD',
			'/SingleMuon/Run2018D-22Jan2019-v2/AOD',
				]
		for inDS in inputDatasets:
			config.General.requestName = inDS.split('/')[2]
			config.Data.inputDataset = inDS
			config.Data.outputDatasetTag = '%s_%s' % (config.General.workArea, config.General.requestName)
			# Submit.
			try:
				print "Submitting for input dataset %s" % (inDS)
				crabCommand(options.crabCmd, config = config, *options.crabCmdOpts.split())
			except HTTPException as hte:
				print "Submission for input dataset %s failed: %s" % (inDS, hte.headers)
			except ClientException as cle:
				print "Submission for input dataset %s failed: %s" % (inDS, cle)
	elif options.workArea:
		for dir in os.listdir(options.workArea):
			projDir = os.path.join(options.workArea, dir)
			if not os.path.isdir(projDir):
				continue
			# Execute the crab command.
			msg = "Executing (the equivalent of): crab %s --dir %s %s" % (options.crabCmd, projDir, options.crabCmdOpts)
			print "-"*len(msg)
			print msg
			print "-"*len(msg)
			try:
				crabCommand(options.crabCmd, dir = projDir, *options.crabCmdOpts.split())
			except HTTPException as hte:
				print "Failed executing command %s for task %s: %s" % (options.crabCmd, projDir, hte.headers)
			except ClientException as cle:
				print "Failed executing command %s for task %s: %s" % (options.crabCmd, projDir, cle)
if __name__ == '__main__':
	main()

