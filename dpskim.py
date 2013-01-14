import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32( 300 ) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(

      'file:/local/cms/phedex/store/data/Run2012C/SinglePhoton/RECO/EXODisplacedPhoton-PromptSkim-v3/000/200/190/00000/18BB0794-8CDF-E111-B9B0-0025B31E3D3C.root'
 
    ),
    # explicitly drop photons resident in AOD/RECO, to make sure only those locally re-made (uncleaned photons) are used
    inputCommands = cms.untracked.vstring('keep *'
                                          #,'drop  *_photonCore_*_RECO' # drop hfRecoEcalCandidate as remade in this process
                                          #, 'drop *_photons_*_RECO' # drop photons as remade in this process
                                          )

)

process.options   = cms.untracked.PSet(
                    wantSummary = cms.untracked.bool(True),  
                    SkipEvent = cms.untracked.vstring('ProductNotFound')
)   


process.skim = cms.EDFilter('DPSkim',

    triggerName  = cms.vstring('HLT_Photon50_CaloIdVL_IsoL', 'HLT_DisplacedPhoton65_CaloIdVL_IsoL_PFMET25'),
    trigSource   = cms.InputTag("TriggerResults","","HLT"),
    jetSource    = cms.InputTag("ak5PFJets"),
    metSource    = cms.InputTag("pfMet"),
    photonSource = cms.InputTag("photons"),

    # photon cuts                pt   eta  Num  leadingPt 
    photonCuts    = cms.vdouble( 15,  2.4,   1 ,     30   ),
    # photon isolation           trk,  ecalSumEt, ecalR, hcalSumEt, hcalR 
    photonIso     = cms.vdouble(  0.2,       99.0,   0.2,      99.0,   0.2 ),
    #photonIso     = cms.vdouble(  0.2,       4.5,   0.1,       4.0,   0.1 ),
    # jet cuts                   pt    eta  nJets
    jetCuts       = cms.vdouble( 35. , 2.4,     0 ),
    metCuts       = cms.vdouble(  30. )

)

process.p = cms.Path(
                     process.skim
                    )

from Configuration.EventContent.EventContent_cff import AODEventContent

process.out = cms.OutputModule("PoolOutputModule"
		, fileName = cms.untracked.string( 'skim_data.root' )
	 	, SelectEvents = cms.untracked.PSet(
		        SelectEvents = cms.vstring('p')
			    )
		, outputCommands = cms.untracked.vstring(
			'drop *'
			)
		)

process.out.outputCommands.extend (AODEventContent.outputCommands )
process.out.outputCommands.extend( ['keep *_cscSegments_*_*'] )
process.outpath = cms.EndPath(process.out)

