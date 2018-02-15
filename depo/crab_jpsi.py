import sys
import os

jsonFile="Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON_MuonPhys.txt"
run = "F"

from WMCore.Configuration import Configuration
config = Configuration()

#print("Test = " + str(skipevt))

jobdir = 'X4140_' + run

if not os.path.exists(jobdir):
    os.makedirs(jobdir)

datasetnames = {
"F" : '/MuOnia/Run2017F-PromptReco-v1/MINIAOD',
}

runNumber = [
#'274094-274240',
'305388-307000', #0
'305388-306000', #1
'305388-305700', #2
'305700-306000', #3
'306000-306100', #4
'306100-306200', #5
'306200-306300', #6
'306300-306500', #7
'' #8
#'273158',
]


datasetName = datasetnames[run]
runNum = runNumber[8]
#lumi = jsonfile[jNum]
lumi = jsonFile
#HLT = HLTPath[0]

import datetime
timestamp = datetime.datetime.now().strftime("_%Y%m%d_%H%M%S")

dataset = filter(None, datasetName.split('/'))

config.section_('General')
config.General.transferOutputs  = True
config.General.workArea         = jobdir
#config.General.requestName     = 'JetHT_Run2015D_PromptReco_v4_RECO'+timestamp
#config.General.requestName             = dataset[0]+'_'+dataset[1]+'_'+dataset[2]+'_'+runNum+'_'+HLT+timestamp
config.General.requestName      = dataset[0]+'_'+dataset[1]+'_'+dataset[2]+'_'+runNum+'_'+timestamp #+'_split_'+ jsonFile.split('_')[-1].split('.')[0]
config.General.transferLogs     = False

config.section_('JobType')
config.JobType.psetName         = '/lustre/home/adrianodif/CMSSW_9_2_2/src/mmkk/mmkk/test/run-jpsi-cand.py'
config.JobType.pluginName       = 'Analysis'
config.JobType.maxMemoryMB      = 2500
config.JobType.maxJobRuntimeMin = 2750
config.JobType.allowUndistributedCMSSW = True
#config.JobType.inputFiles      = ['Run_time.txt']
#config.JobType.outputFiles     = ['EventList.txt','EventListClean.txt']

config.section_('Data')
config.Data.inputDataset        = datasetName
config.Data.inputDBS            = 'global'
config.Data.totalUnits          = -1
config.Data.unitsPerJob         = 1
config.Data.splitting           = 'LumiBased'
config.Data.runRange            = runNum
config.Data.lumiMask            = lumi
config.Data.outLFNDirBase       = '/store/user/adiflori/'
config.Data.publication         = False
config.Data.ignoreLocality      = True


config.section_('Site')
#config.Site.storageSite        = 'T2_CH_CERN'
config.Site.storageSite         = 'T2_IT_Bari'
#config.Site.blacklist          = ['T2_IN_TIFR','T2_US_Vanderbilt']
config.Site.blacklist           = ['T1*']
#config.Site.whitelist           = ['T2_IT_Legnaro']
#config.Site.whitelist          = ['T2_IT_Bari','T2_IT_Pisa','T2_IT_Rome','T2_IT_Legnaro']
