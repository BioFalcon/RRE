#!/usr/bin/env nextflow

def helpMessage() {
  log.info"""
  ==============================================================
    ▄▄▄▄▄▄▄   ▄▄▄▄▄▄▄    ▄▄▄▄▄▄▄ 
    ███▀▀███▄ ███▀▀███▄ ███▀▀▀▀▀ 
    ███▄▄███▀ ███▄▄███▀ ███▄▄    
    ███▀▀██▄  ███▀▀██▄  ███      
    ███  ▀███ ███  ▀███ ▀███████ 
    Recursive Repeat    Extension
  ==============================================================

  Usage: nextflow run RRE.nf [PARAMETERS]

  Parameters:
  --Genome               Genome to be analyzed in fasta format.
  --consensus            Consensus sequences in fasta format.
  --outDir               Output directory (default:./Results).
  --extension            Number of bases to extend the consensus (default:250)
  --workflow             Workflow to use for the extension [HEEA|RRE]
  --maxExtensionSize     Maximum size for the extension to be considered  (default:25000)
                         NOTE:In recursive extension, this limit only applies to each side separately
  --hyperT               Activate Hyperthreading (off by default)
  --consensusAln         Alignment file of the consensus sequences in stockholm format.
  --maxRounds           Maximum number of rounds for the recursive extension (default:30)
  


  Help:
  --h,--help            Show this help message and exit.
  """.stripIndent()
  return ""
}


/************************************************
****                                         ****
****         Variable Initialization         ****
****                                         ****
************************************************/

//Initialize variables
params.Genome                      = ""
params.consensus                   = ""
params.outDir                      = "./Results"
params.extension                   = 250
params.workflow                    = ""
params.help                        = ""
params.h                           = ""
params.maxExtensionSize            = 25000
params.hyperT                      = ""
params.consensusAln                = ""
params.maxRounds                   = 30  
params.maxFamilies                 = 10
params.SampleSeqs                  = 25

//Filtering parameters
params.hmmEValue                   = 1e-5

//Initialize Curation parameters
params.noiseThreshold              = 0.2

//Variables for debugging and advanced usage
params.maxCentralRounds            = 1 

//Other parameters
params.prevRoundCoverage           = 0.5

//Parameter for HMMER Results (in case of available)
params.hmmResults                  = ""

//Parameter for AncientMode
params.AncientMode                   = ""

//Parameter for Horizontal cleaning
params.percentVertical             = 0.3

//Parameters for Polishing
params.polishminHits                = 10
params.polishLengthLimitMultiplier  = 1.5
params.polishWindow                 = 10000
params.polishCoverage               = 50

/************************************************
****                                         ****
****               Variable Check            ****
****                                         ****
************************************************/

//Check if help needed
if ( params.help || params.h ){
    exit 0, helpMessage()
}

//Check if required files/arguments to run are present
if ( !params.Genome){
    exit 1, helpMessage() + "--Genome not specified !!!"
} else if (!params.consensus){
    exit 1, helpMessage() + "--consensus not specified !!!"
} else if (!params.consensusAln){
    exit 1, helpMessage() + "--consensusAln not specified !!!"
} else if ( !params.workflow ){
    exit 1, helpMessage() + "--workflow not specified !!!"
} else if ( params.workflow != "RRE" && params.workflow != "HEEA" ){
    exit 1, helpMessage() + "--workflow must be either 'RRE' or 'HEEA' !!!"
}

/************************************************
****                                         ****
****              Variable Setup             ****
****                                         ****
************************************************/

//Setup hyperT parameter
if (params.hyperT){
    hyperT = "True"
} else {
    hyperT = "False"
}


if (params.hmmResults){
    hmmResults = file(params.hmmResults)
}

//Obtain the absolute paths for all the required files
outDir             = file(params.outDir)
Genome             = file(params.Genome)
consensus          = file(params.consensus)
consensusAln       = file(params.consensusAln)

CurateScript       = file("./utils/MSACuration.py")
CoverageScript     = file("./utils/CoverageSelection.py")
RREScript          = file("./utils/RecursiveExtension.sh")
RREOLDScript       = file("./utils/RecursiveOldExt.sh")
CoordExtractScript = file("./utils/CoordinateExtraction.py")
OutlierScript      = file("./utils/FamilyDetection.py")
VerticalScript     = file("./utils/VerticalCleaning.py")

/************************************************
****                                         ****
****                Print log                ****
****                                         ****
************************************************/
//Print log with used paths
log.info"""
  ==============================================================
    ▄▄▄▄▄▄▄   ▄▄▄▄▄▄▄    ▄▄▄▄▄▄▄ 
    ███▀▀███▄ ███▀▀███▄ ███▀▀▀▀▀ 
    ███▄▄███▀ ███▄▄███▀ ███▄▄    
    ███▀▀██▄  ███▀▀██▄  ███      
    ███  ▀███ ███  ▀███ ▀███████ 
    Recursive Repeat    Extension
  ==============================================================   
""".stripIndent()
log.info "Running the pipeline with the following parameters: "
log.info "============================================================="
log.info "Genome path                : ${params.Genome}"
log.info "Consensus path             : ${params.consensus}"
log.info "Consensus alignment file   : ${params.consensusAln}"
log.info "Output path                : ${params.outDir}"
log.info "Extension size             : ${params.extension}"
log.info "Output path                : ${params.outDir}"
log.info "Maximum extension size     : ${params.maxExtensionSize}"
log.info "============================================================="
log.info "Running the pipeline with the following hardcoded parameters:"
log.info "Curation script path       : ${CurateScript}"
log.info "Coverage script path       : ${CoverageScript}"
log.info "Outlier detection script   : ${OutlierScript}"
log.info "Coordinate extraction scr. : ${CoordExtractScript}"
log.info "Vertical cleaning script   : ${VerticalScript}"
log.info "============================================================="

/************************************************
****                                         ****
****            Include modules              ****
****                                         ****
************************************************/

//Preprocessing
include {WorkDirPrep    } from './modules/01_Preprocessing/01_WorkDirPrep'
include {HMMDBBuild     } from './modules/01_Preprocessing/02_HMMDBBuild'
include {CreateNSplitHMM} from './modules/01_Preprocessing/03_CreateHMM'
include {HMMER_Run      } from './modules/01_Preprocessing/04_HMMRun'
include {HMMSelection   } from './modules/01_Preprocessing/05_HMMSelection'
include {FAIndex        } from './modules/01_Preprocessing/06_FAIndex'

//RRE
////Extension Central
if ( params.AncientMode ){
    include { RRE_CentralExtension } from './modules/04_OldExt/01_CentralExtension/main.nf'
} else {
    include { RRE_CentralExtension } from './modules/02_RRE/01_CentralExtension/main.nf'
}

////Recursive Extension
include {RRE_RecursiveExtension } from './modules/02_RRE/02_RecursiveExtension/main.nf'

//HEEA
////Exntension HEEA
include { HEEA_Extension } from './modules/03_HEEA/01_HEEAExtension/main.nf'

//Polishing
include { Polishing      } from './modules/05_Polishing/01_Polishing/main.nf'

//CDHit
include { CDHit          } from './modules/06_CDHit/01_CDHit/main.nf'
include { CDHitHEEA      } from './modules/06_CDHit/02_CDHitHEEA/main.nf'

/************************************************
*************************************************
****                Workflow                 ****
*************************************************
************************************************/

workflow HEEA{
    take:
        Consensus
        HMMRDB
        ConsensusFasta
        GenomeIndex
        HMMROut
        outDir
        Genome
        CurateScript
        OutlierScript
        CoordExtractScript

    main:
    HEEA_Extension(Consensus,
                        HMMRDB,
                        ConsensusFasta,
                        GenomeIndex,
                        HMMROut,
                        outDir,
                        Genome,
                        CurateScript,
                        OutlierScript,
                        CoordExtractScript)

    CDHitHEEA( HEEA_Extension.out.flatten() ,
               outDir,
               ConsensusFasta,
               consensusAln )
}

workflow RRE{
    take:
        ConsensusSelection
        HMMRDB
        ConsensusFasta
        FAIndex
        HMMROut
        outDir
        Genome
        CurateScript
        OutlierScript
        CoordExtractScript
        RREScript
        RREOLDScript
        VerticalScript
        CoverageScript

    main:
    //////////////////////////////////////////////
    //Central Extension
    //////////////////////////////////////////////
    RRE_CentralExtension( ConsensusSelection,
                          HMMRDB,
                          ConsensusFasta,
                          FAIndex,
                          HMMROut,
                          outDir,
                          Genome,
                          CurateScript,
                          OutlierScript,
                          CoordExtractScript
                        )

    //////////////////////////////////////////////
    //Recursive Extension
    //////////////////////////////////////////////
    //In case of AncientMode script
    if ( params.AncientMode ){
        RREScript = RREOLDScript
    }

    RRE_RecursiveExtension( RRE_CentralExtension.out.AllSides,
                        HMMRDB,
                        FAIndex,
                        outDir,
                        CurateScript,
                        OutlierScript,
                        Genome,
                        RREScript,
                        CoordExtractScript,
                        VerticalScript )
    
    //////////////////////////////////////////////
    //Polishing
    //////////////////////////////////////////////
    Polishing( RRE_RecursiveExtension.out ,
                outDir,
                CoverageScript,
                CurateScript,
                Genome,
                HMMRDB )

    //////////////////////////////////////////////
    //Merging
    //////////////////////////////////////////////
    CDHit( Polishing.out.flatMap{ n -> [n[1]]}.collect(),
           outDir,
           ConsensusFasta,
           consensusAln )
    
}

//Main workflow
workflow{
    //Create Workdir and result dir
    if (!outDir.exists()){
        outDir.mkdirs()
    }
    outDir=Channel.fromPath(outDir)
    WorkDirPrep(outDir)

    //Make Channels 
    Genome             = Channel.fromPath(Genome)
    consensus          = Channel.fromPath(consensus)
    consensusAln       = Channel.fromPath(consensusAln)
    CurateScript       = Channel.fromPath(CurateScript)
    CoverageScript     = Channel.fromPath(CoverageScript)
    OutlierScript      = Channel.fromPath(OutlierScript)
    RREScript          = Channel.fromPath(RREScript)
    RREOLDScript       = Channel.fromPath(RREOLDScript)
    CoordExtractScript = Channel.fromPath(CoordExtractScript)
    VerticalScript     = Channel.fromPath(VerticalScript)
    hyperTCh           = Channel.from(hyperT)

    /************************
    ***** Preprocessing ***** 
    *************************/
    //////////////////////////////////////////////
    //Preprocessing steps
    //////////////////////////////////////////////
    
    ////Create HMM database
    HMMDBBuild(Genome)

    ////Create HMM models and split them
    CreateNSplitHMM(consensusAln)

    if ( params.hmmResults ){
        HMMER_Results = Channel.fromPath(hmmResults)
    } else {
        ////Run HMMERsearch
        HMMER_Run( HMMDBBuild.out.collect(),
                   CreateNSplitHMM.out.SplitHMM.flatMap() ) 

        ////Select the consensi to extend
        HMMER_Results = HMMER_Run.out.collectFile(name:"HMMER_Results.txt")
    }

    HMMSelection(HMMER_Results,
                 consensus.collect() )
    
    ////Create fai index for posterior runs
    FAIndex(Genome)

    if ( params.workflow == "RRE" ){
        RRE( HMMSelection.out.splitText(by:1, file:true),
                           HMMDBBuild.out.collect(),
                           consensus.collect(),
                           FAIndex.out.collect(),
                           HMMER_Results.collect(),
                           outDir.collect(),
                           Genome.collect(),
                           CurateScript.collect(),
                           OutlierScript.collect(),
                           CoordExtractScript.collect(),
                           RREScript.collect(),
                           RREOLDScript.collect(),
                           VerticalScript.collect(),
                           CoverageScript.collect() )

    } else if ( params.workflow == "HEEA" ){
        HEEA( HMMSelection.out.splitText(by:1, file:true),
                           HMMDBBuild.out.collect(),
                           consensus.collect(),
                           FAIndex.out.collect(),
                           HMMER_Results.collect(),
                           outDir.collect(),
                           Genome.collect(),
                           CurateScript.collect(),
                           OutlierScript.collect(),
                           CoordExtractScript.collect() )
    } 
}