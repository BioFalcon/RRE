#!/bin/bash

#Defaults for the AncientMode-specific options
maxRetries=3
minConservation=0.45
supportWindow=10
minBits=30
rfMinOcc=0.5
minCopies=10
maxIdentity=auto
refineIterations=8
curationIdentityExp=1.4
curationGapExp=0.3
curationMaxThresh=0.1
anchorLength=500
anchorSpan=300
anchorMinRatio=0.2
anchorMinCopies=5
anchorMaxJump=3

while getopts "i:l:r:c:t:e:g:h:x:o:C:O:s:m:p:S:R:F:T:P:V:M:K:W:B:Q:N:D:G:A:Y:X:U:I:J:Z:L:" opt
do
   case "$opt" in
      i ) RepeatID="$OPTARG" ;;
      t ) CPU="$OPTARG" ;;
      e ) Extension="$OPTARG" ;;
      g ) Genome="$OPTARG" ;; 
      h ) HMMRDB="$OPTARG" ;;
      x ) GenomeIndex="$OPTARG" ;;
      o ) outDir="$OPTARG" ;;
      C ) CurateScript="$OPTARG" ;;
      O ) OutlierScript="$OPTARG" ;;
      s ) CoordExtractScript="$OPTARG" ;;
      m ) maxExtensionSize="$OPTARG" ;;
      p ) prevRoundCoverage="$OPTARG" ;;
      S ) SampleSeqs="$OPTARG" ;;
      R ) maxRounds="$OPTARG" ;;
      F ) maxFamilies="$OPTARG" ;;
      T ) noiseThreshold="$OPTARG" ;;
      P ) percentVertical="$OPTARG" ;;
      V ) VerticalScript="$OPTARG" ;;
      M ) maxRetries="$OPTARG" ;;
      K ) minConservation="$OPTARG" ;;
      W ) supportWindow="$OPTARG" ;;
      B ) minBits="$OPTARG" ;;
      Q ) rfMinOcc="$OPTARG" ;;
      N ) minCopies="$OPTARG" ;;
      D ) maxIdentity="$OPTARG" ;;
      G ) refineIterations="$OPTARG" ;;
      A ) curationIdentityExp="$OPTARG" ;;
      Y ) curationGapExp="$OPTARG" ;;
      X ) curationMaxThresh="$OPTARG" ;;
      U ) anchorLength="$OPTARG" ;;
      I ) anchorSpan="$OPTARG" ;;
      J ) anchorMinRatio="$OPTARG" ;;
      Z ) anchorMinCopies="$OPTARG" ;;
      L ) anchorMaxJump="$OPTARG" ;;
   esac
done

#The AncientMode helpers live next to the other utils (resolved here, before changing directory)
UtilsDir=$(dirname $(realpath ${VerticalScript}))
HandRFScript=${UtilsDir}/HandRF.py
StitchScript=${UtilsDir}/StitchMerge.py
if [ ! -f ${HandRFScript} ] || [ ! -f ${StitchScript} ];then
    echo "HandRF.py or StitchMerge.py not found next to ${VerticalScript}, exiting"
    exit 1
fi
    
#Save nf directory into a variable
oldWorkDir=$(pwd)

#Move to new WorkDir 
cd ${outDir}/WorkDir/${RepeatID}

#Backup new path 
NewPath=$( pwd )

#Initialize some some variables
ExtensionSize=${Extension}
#AnchorSize=$(echo ${ExtensionSize} | awk '{printf "%.0f\n", $1*0.5}') 
AnchorSize=$(  echo ${ExtensionSize} | awk '{print int($1*1.5)}' )
Increase=$(  echo ${ExtensionSize} | awk '{$0=$0*0.125;printf("%d\n",$0+=$0<0?0:0.9)}' )
LongConsensus=False

######################
#####  Function  #####
######################

#Build an HMM whose match states are the RF columns defined by HandRF.py,
#and its consensus (one consensus base per match state, so HMM and consensus coordinates agree)
#Usage: BuildHandHMM <aln.fa> <OutPrefix> <HMMName> <TrimSide: Left|Right|None>
#Writes <OutPrefix>.stk, .hmm, .Consensus.fa and .support.tsv; replaces <aln.fa> by the trimmed alignment
function BuildHandHMM() {
    local AlnIn=$1
    local Prefix=$2
    local HMMName=$3
    local TrimSide=$4

    python3 ${HandRFScript} \
        --input ${AlnIn} \
        --output ${Prefix}.HandRF \
        --minOcc ${rfMinOcc} \
        --trimSide ${TrimSide} \
        --minCons ${minConservation} \
        --window ${supportWindow} \
        --minSeqs ${minCopies} \
        --name ${HMMName} > ${Prefix}.HandRF.log || return 1

    #Not enough match columns left to build a model
    if [ $(awk '$1=="MATCH_COLUMNS"{print $2}' ${Prefix}.HandRF.log) -lt ${supportWindow} ];then
        return 1
    fi

    cp ${Prefix}.HandRF.aln.fa ${AlnIn}
    mv ${Prefix}.HandRF.stk ${Prefix}.stk
    mv ${Prefix}.HandRF.support.tsv ${Prefix}.support.tsv
    mv ${Prefix}.HandRF.match.aln.fa ${Prefix}.match.aln.fa

    hmmbuild \
        --hand \
        --dna \
        --cpu ${CPU} \
        --informat stockholm \
        --seed 1992 \
        --fragthresh 1.0 \
        --maxinsertlen 10 \
        --eentexp \
        -n ${HMMName} \
        ${Prefix}.hmm \
        ${Prefix}.stk > /dev/null || return 1

    hmmemit -c ${Prefix}.hmm > ${Prefix}.Consensus.fa || return 1
}

#Write the coordinates of the outermost <Size> bases of a consensus, on the side being extended
#Usage: SetOuterCoord <Consensus.fa> <Side> <Size> <OutCoord>
function SetOuterCoord() {
    local LenCons=$(seqkit fx2tab --length --name $1 | cut -f2)
    if [[ $2 == "Left" ]];then
        echo -e "1\t$(( $3 < LenCons ? $3 : LenCons ))" > $4
    else
        echo -e "$(( LenCons - $3 > 1 ? LenCons - $3 : 1 ))\t${LenCons}" > $4
    fi
}

#Search and sampling settings for each retry level. A round that fails to extend is retried
#from the last good model, searching from the very edge and, at higher levels, with a wider flank,
#more sampled copies and a lower bit-score threshold. After maxRetries failed retries the side stops.
function RetrySettings() {
    local Level=$1
    case ${Level} in
        0|1) RoundExt=${ExtensionSize}
             RoundSample=${SampleSeqs}
             RoundMinBits=${minBits} ;;
        2)   RoundExt=$(( ExtensionSize * 2 ))
             RoundSample=$(( SampleSeqs * 3 / 2 ))
             RoundMinBits=$(( minBits - 5 )) ;;
        *)   RoundExt=$(( ExtensionSize * 2 ))
             RoundSample=$(( SampleSeqs * 2 ))
             RoundMinBits=$(( minBits - 10 )) ;;
    esac
}

function HMMExtensionSides() {
#Define starting variables
local CurrRound=01
local BreakNow=False
local CentralOK=True
local SkipRounds=No
local Family=00
local NumFamilies=1

#Define directory
local CentralDir=$(realpath ../Central/)
local RightDir=$(realpath ../Right/)
local LeftDir=$(realpath ../Left/)

#Side Variables
if [[ $1 == "Left" ]];then
        local LowerCase=Left
        local UpperCase=LEFT
        local LetterCase=L
else
        local LowerCase=Right
        local UpperCase=RIGHT
        local LetterCase=R
fi

#Check if Staging dir exist
if [ ! -d Staging ] || [ ! -f ./Staging/GoodCHECK ] ;then
        mkdir -p Staging
        ln -s -f ${oldWorkDir}/* ./Staging/

        #Make check for Staging
        > ./Staging/GoodCHECK
fi

#Update number of families
if [[ $1 == "Left" ]];then
    NumFamilies=$( ls -d ${CentralDir}/Family_*/ | wc -l )
    for FamilyDir in $( ls -d ${CentralDir}/Family_*/ );do
        mkdir -p Family_$(echo ${FamilyDir} | sed 's/.*Family_//;s/\///')
    done
else
    NumFamilies=$( ls -d ${LeftDir}/Family_*/ | wc -l )
    for FamilyDir in $( ls -d ${LeftDir}/Family_*/ );do
        mkdir -p Family_$(echo ${FamilyDir} | sed 's/.*Family_//;s/\///')
    done
fi

while [ ${Family} -lt ${NumFamilies} ];do
    cd Family_${Family}

    #Reset variables
    BreakNow=False
    SkipRounds=No
    CurrRound=01

    #Check if should continue
    if [ -f ${CentralDir}/Family_${Family}/Final/${LowerCase}PackCentral/${UpperCase}.BAD.CHECK ];then
        CentralOK=False

        > ./BAD.SIDE.CHECK

    else
        #Check if Round_00 exists
        if [ ! -d Round_00 ] || [ ! -f ./Round_00/GoodCHECK ]; then
            ##Copy Files for this round
            mkdir -p Round_00
            if [[ $1 == "Left" ]];then
                cp ${CentralDir}/Family_${Family}/Final/FinalConsensi.aln.fa \
                    ./Round_00/07_CurrentConsensi.${LowerCase}.Round00.Extended.Curated.aln.fa
                cp ${CentralDir}/Family_${Family}/Final/FinalConsensi.consensus.fa \
                    ./Round_00/07_CurrentConsensi.${LowerCase}.Round00.Extended.Curated.Consensus.fa
                cp ${CentralDir}/Family_${Family}/Final/Final.Consensus.${LowerCase}Side.Coord \
                    ./Round_00/08_CurrentConsensi.${LowerCase}.Round00.ExtendedSide.Coord
                cp ${CentralDir}/Family_${Family}/Final/FinalConsensi.Extension.Coord.bed \
                    ./Round_00/02_CurrentConsensi.${LowerCase}.Round00.Extended.bed
                cp ${CentralDir}/Family_${Family}/Final/FinalConsensi.aln.fa \
                    ./Round_00/05_CurrentConsensi.${LowerCase}.Round00.Extended.aln.fa
                cp ./Round_00/07_CurrentConsensi.${LowerCase}.Round00.Extended.Curated.Consensus.fa \
                    ./Round_00/05_CurrentConsensi.${LowerCase}.Round00.Extended.Consensus.fa
                #Add Round00 IDs
                awk '{ 
                    if($1 ~ "^>"){
                        gsub("$","___Central",$1);
                    };
                    print
                    }' ./Round_00/07_CurrentConsensi.${LowerCase}.Round00.Extended.Curated.aln.fa \
                > ./Round_00/07_TEMP.aln.fa

                mv ./Round_00/07_TEMP.aln.fa ./Round_00/07_CurrentConsensi.${LowerCase}.Round00.Extended.Curated.aln.fa

                cp ./Round_00/07_CurrentConsensi.${LowerCase}.Round00.Extended.Curated.aln.fa \
                    ./Round_00/05_CurrentConsensi.${LowerCase}.Round00.Extended.aln.fa
            else
                #Check if left finished good
                if [ -f ${LeftDir}/Family_${Family}/GOOD.FAM.CHECK ];then
                        #Determine highest good left round
                        HighestRound=$( ls ${LeftDir}/Family_${Family}/Round_*/GoodCHECK | sed "s/\/GoodCHECK//;s/Left\///;s/.*\///"| tail -1 )
                        tac -b -r -s "^>"  ${LeftDir}/Family_${Family}/${HighestRound}/07_CurrentConsensi.Left*Extended.Curated.aln.fa \
							> ./Round_00/07_CurrentConsensi.${LowerCase}.Round00.Extended.Curated.aln.fa
                        cp ${LeftDir}/Family_${Family}/${HighestRound}/07_CurrentConsensi.Left*Extended.Curated.Consensus.fa \
                            ./Round_00/07_CurrentConsensi.${LowerCase}.Round00.Extended.Curated.Consensus.fa
                        #Column calls of the merged model (StitchMerge.py), absent if the left side never extended
                        for File in ${LeftDir}/Family_${Family}/${HighestRound}/07_CurrentConsensi.Left*Extended.Curated.Columns.txt;do
                            [ -f ${File} ] && cp ${File} ./Round_00/07_CurrentConsensi.${LowerCase}.Round00.Extended.Curated.Columns.txt
                        done
                        cp ${LeftDir}/Family_${Family}/${HighestRound}/02_CurrentConsensi.Left*Extended.bed \
                            ./Round_00/02_CurrentConsensi.${LowerCase}.Round00.Extended.bed                            
                        cp ${CentralDir}/Family_${Family}/Final/Final.Consensus.RightSide.Coord \
                            ./Round_00/08_CurrentConsensi.${LowerCase}.Round00.ExtendedSide.Coord
                        cp ${CentralDir}/Family_${Family}/Final/FinalConsensi.consensus.fa \
                            ./Round_00/05_CurrentConsensi.${LowerCase}.Round00.Extended.Consensus.fa
                        cp ${CentralDir}/Family_${Family}/Final/FinalConsensi.aln.fa \
                            ./Round_00/05_CurrentConsensi.${LowerCase}.Round00.Extended.aln.fa


                        #Generate Coord for Right
                        #echo -e "$(( ${ModelLen} - ${ExtSize} ))\t${ModelLen}" > ./Round_00/08_CurrentConsensi.${LowerCase}.Round00.ExtendedSide.Coord
                else
                    cp ${CentralDir}/Family_${Family}/Final/FinalConsensi.aln.fa \
                        ./Round_00/07_CurrentConsensi.${LowerCase}.Round00.Extended.Curated.aln.fa
                    cp ${CentralDir}/Family_${Family}/Final/FinalConsensi.consensus.fa \
                        ./Round_00/07_CurrentConsensi.${LowerCase}.Round00.Extended.Curated.Consensus.fa
                    cp ${CentralDir}/Family_${Family}/Final/Final.Consensus.${LowerCase}Side.Coord \
                        ./Round_00/08_CurrentConsensi.${LowerCase}.Round00.ExtendedSide.Coord
                    cp ${CentralDir}/Family_${Family}/Final/FinalConsensi.Extension.Coord.bed \
                        ./Round_00/02_CurrentConsensi.${LowerCase}.Round00.Extended.bed
                    cp ${CentralDir}/Family_${Family}/Final/FinalConsensi.aln.fa \
                        ./Round_00/05_CurrentConsensi.${LowerCase}.Round00.Extended.aln.fa
                    cp ${CentralDir}/Family_${Family}/Final/FinalConsensi.consensus.fa \
                        ./Round_00/05_CurrentConsensi.${LowerCase}.Round00.Extended.Consensus.fa
                fi
                
                #Add Round00 IDs
                awk '{ 
                    if($1 ~ "^>"){
                        gsub("$","___Central",$1);
                    };
                    print
                    }' ./Round_00/05_CurrentConsensi.${LowerCase}.Round00.Extended.aln.fa \
                > ./Round_00/05_TEMP.aln.fa

                mv ./Round_00/05_TEMP.aln.fa ./Round_00/05_CurrentConsensi.${LowerCase}.Round00.Extended.aln.fa
            fi
        
            if ( [[ $1 == "Left" ]] || ( [[ $1 == "Right" ]] && [ ! -f ${LeftDir}/Family_${Family}/GOOD.FAM.CHECK ] ) ) && \
                [ -f ${CentralDir}/Family_${Family}/Final/FinalConsensi.SeedModel ];then
                ##The seed is the starting model: use its HMM and match columns as they are
                cp ${CentralDir}/Family_${Family}/Final/FinalConsensi.hmm \
                    ./Round_00/05_CurrentConsensi.${LowerCase}.Round00.Extended.hmm
                cp ./Round_00/07_CurrentConsensi.${LowerCase}.Round00.Extended.Curated.Consensus.fa \
                    ./Round_00/05_CurrentConsensi.${LowerCase}.Round00.Extended.Consensus.fa
                #Every seed position is a consensus position of the merged model (StitchMerge.py column calls)
                grep -v ">" ./Round_00/07_CurrentConsensi.${LowerCase}.Round00.Extended.Curated.Consensus.fa | tr -d '\n' | \
                    tr 'acgtn' 'ACGTN' > ./Round_00/07_CurrentConsensi.${LowerCase}.Round00.Extended.Curated.Columns.txt
                echo >> ./Round_00/07_CurrentConsensi.${LowerCase}.Round00.Extended.Curated.Columns.txt
            elif [[ $1 == "Left" ]] || ( [[ $1 == "Right" ]] && [ ! -f ${LeftDir}/Family_${Family}/GOOD.FAM.CHECK ] );then
                ##Make HMM model with match columns defined by RF (hmmbuild --hand)
                cp ./Round_00/07_CurrentConsensi.${LowerCase}.Round00.Extended.Curated.aln.fa \
                    ./Round_00/08_CurrentConsensi.${LowerCase}.Round00.Model.aln.fa
                BuildHandHMM \
                    ./Round_00/08_CurrentConsensi.${LowerCase}.Round00.Model.aln.fa \
                    ./Round_00/08_CurrentConsensi.${LowerCase}.Round00.Extended.Curated \
                    ${RepeatID} \
                    None || { echo -e "\nCould not build the Round00 model, exiting"; exit 1; }

                cp ./Round_00/08_CurrentConsensi.${LowerCase}.Round00.Extended.Curated.hmm \
                    ./Round_00/05_CurrentConsensi.${LowerCase}.Round00.Extended.hmm
                cp ./Round_00/08_CurrentConsensi.${LowerCase}.Round00.Extended.Curated.Consensus.fa \
                    ./Round_00/05_CurrentConsensi.${LowerCase}.Round00.Extended.Consensus.fa

                ##Merged model: one column per consensus position, extended by stitching each round
                cp ./Round_00/08_CurrentConsensi.${LowerCase}.Round00.Extended.Curated.match.aln.fa \
                    ./Round_00/07_CurrentConsensi.${LowerCase}.Round00.Extended.Curated.aln.fa
                cp ./Round_00/08_CurrentConsensi.${LowerCase}.Round00.Extended.Curated.Consensus.fa \
                    ./Round_00/07_CurrentConsensi.${LowerCase}.Round00.Extended.Curated.Consensus.fa
            else
                ##Copy HMM model from left
                cp ${LeftDir}/Family_${Family}/Round_00/05_CurrentConsensi.Left*Extended.hmm \
                    ./Round_00/05_CurrentConsensi.${LowerCase}.Round00.Extended.hmm
                cp ${LeftDir}/Family_${Family}/Round_00/05_CurrentConsensi.Left*Extended.Consensus.fa \
                    ./Round_00/05_CurrentConsensi.${LowerCase}.Round00.Extended.Consensus.fa
            fi

            ##Search window: the outer edge of the model, in the model's own coordinates
            SetOuterCoord \
                ./Round_00/05_CurrentConsensi.${LowerCase}.Round00.Extended.Consensus.fa \
                ${LowerCase} \
                ${ExtensionSize} \
                ./Round_00/08_CurrentConsensi.${LowerCase}.Round00.ExtendedSide.Coord

            echo 0 > ./Round_00/NextRetryLevel

            ##Make check
            > ./Round_00/GoodCHECK
        fi

        #Check if other Rounds exist and backup if needed
        if [ $( ls -d ./Round_*| wc -l ) -gt 1 ];then
            #Check if Final dir is there
            if [ -f ./GOOD.FAM.CHECK ] ;then
                #Set variable to not do rounds
                SkipRounds=Yes
            fi
            #Get which is the last round, and back it up
            CurrRound=$( ls -d ./Round_* | tail -n 1 | sed 's/.\/Round_//' )

            #Check if it is complete and increase round
            if [ -f ./Round_${CurrRound}/GoodCHECK ];then
                CurrRound=$(echo ${CurrRound} | awk '{printf "%02d\n", $1+1}')
            else
                #Make backup of failed round 
                BackupNum=01
                Back=True
                while [ $Back == "True" ];do
                    if [ -d ./BackUp_Round${CurrRound}_${BackupNum} ];then
                        BackupNum=$(echo ${BackupNum} | awk '{printf "%02d\n", $1+1}')
                    else
                        mv ./Round_${CurrRound} ./BackUp_Round${CurrRound}_${BackupNum}
                        Back=False
                    fi
                done
            fi

            #Check if round has surpassed max rounds
            if [ ${CurrRound} -ge ${maxRounds} ];then
                SkipRounds=Yes
            fi 
        fi

        #Rounds in a row whose new bases are not all supported by copies of the model edge,
        #and the number of copies of the model edge in each round
        AnchoredLowStreak=0
        AnchoredHistory=""

        #Main loop
        while [ ${BreakNow} == "False" ] && [ ${SkipRounds} == "No" ];do
            #Make dir for this round
            mkdir -p Round_${CurrRound}

            #Set previous round
            PrevRound=$(echo ${CurrRound} | awk '{printf "%02d\n", $1-1}')

            #Retry level set by the previous round (0 after a successful extension)
            RetryLevel=$(cat ./Round_${PrevRound}/NextRetryLevel 2>/dev/null || echo 0)
            RetrySettings ${RetryLevel}
            RoundFailed=False
            echo -e "\n${LowerCase} round ${CurrRound}: retry level ${RetryLevel} (flank ${RoundExt}, ${RoundSample} copies, bits > ${RoundMinBits})"

            #HMMER Search
            nhmmer \
            --cpu ${CPU} \
            --tblout ./Round_${CurrRound}/00_CurrentConsensi.${LowerCase}.Round${CurrRound}.hmmout \
            --notextw \
            --noali \
            --seed 1992 \
            ./Round_${PrevRound}/05_CurrentConsensi.${LowerCase}.Round${PrevRound}.Extended.hmm \
            ../Staging/HMMDB_Genome  > /dev/null

            ErrorCode=$?
            if [ ${ErrorCode} -ne 0 ];then
                echo -e "\nHMMER search failed, exiting"
                exit ${ErrorCode}
            fi
            
            #Format HMMER output
            sed '/#/d; s/\s\{1,\}/\t/g' ./Round_${CurrRound}/00_CurrentConsensi.${LowerCase}.Round${CurrRound}.hmmout | \
            awk -v MinBits=${RoundMinBits} '$14 > MinBits {print}' \
            > ./Round_${CurrRound}/00_CurrentConsensi.${LowerCase}.Round${CurrRound}.hmmout.tmp
            
            cp ./Round_${CurrRound}/00_CurrentConsensi.${LowerCase}.Round${CurrRound}.hmmout.tmp \
                ./Round_${CurrRound}/00_CurrentConsensi.${LowerCase}.Round${CurrRound}.hmmout

            #Filter hits with high bias
            #awk '{OFS="\t"} ($15/$14) <= 0.20{print}' ./Round_${CurrRound}/00_CurrentConsensi.${LowerCase}.Round${CurrRound}.hmmout > ./Round_${CurrRound}/00_Temp.hmmout
            #mv ./Round_${CurrRound}/00_Temp.hmmout ./Round_${CurrRound}/00_CurrentConsensi.${LowerCase}.Round${CurrRound}.hmmout

            #Filter those that cover extended part
            LeftCoord=$( cat ./Round_${PrevRound}/08_CurrentConsensi.${LowerCase}.Round${PrevRound}.ExtendedSide.Coord  | cut -f1 )
            RightCoord=$( cat ./Round_${PrevRound}/08_CurrentConsensi.${LowerCase}.Round${PrevRound}.ExtendedSide.Coord | cut -f2 )

            bedtools intersect \
            -a <( cat ./Round_${PrevRound}/08_CurrentConsensi.${LowerCase}.Round${PrevRound}.ExtendedSide.Coord | sed 's/^/Target\t/' ) \
            -b <( awk '{OFS="\t"}{print "Target",$5,$6,$1,$9,$10,$12,$14}' ./Round_${CurrRound}/00_CurrentConsensi.${LowerCase}.Round${CurrRound}.hmmout ) \
            -wb \
            -f ${prevRoundCoverage} | \
            awk -v Side=${LowerCase} '{OFS="\t"}{ 
            FragLen=int( 1.10 * ($3 - $2) );
            if (Side == "Left"){
                    if( $10 == "+"){print $7,$8,$8 + FragLen,".",$11,$10}else{print $7,$8 - FragLen,$8,".",$11,$10}
            } else {
                    if( $10 == "+"){print $7,$9 - FragLen,$9,".",$11,$10}else{print $7,$9,$9 + FragLen,".",$11,$10}
            }
            }'  > ./Round_${CurrRound}/00_CurrentConsensi.${LowerCase}.Round${CurrRound}.Raw.bed

            #Extend the sequences to be aligned
            if [[ $1 == "Left" ]];then
                bedtools slop \
                    -s \
                    -l ${RoundExt} \
                    -r ${AnchorSize} \
                    -i ./Round_${CurrRound}/00_CurrentConsensi.${LowerCase}.Round${CurrRound}.Raw.bed \
                    -g <(cut -f1,2  ../Staging/${GenomeIndex}) > ./Round_${CurrRound}/01_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.bed
            else
                bedtools slop \
                    -s \
                    -r ${RoundExt} \
                    -l ${AnchorSize} \
                    -i ./Round_${CurrRound}/00_CurrentConsensi.${LowerCase}.Round${CurrRound}.Raw.bed \
                    -g <(cut -f1,2  ../Staging/${GenomeIndex}) > ./Round_${CurrRound}/01_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.bed
            fi

            #Select Sequences to be used
            CounterSeqs=0
            FileLine=1
            TotalLines=$(wc -l < ./Round_${CurrRound}/01_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.bed )
            > ./Round_${CurrRound}/02_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.bed

            while [ $CounterSeqs -lt ${RoundSample} ];do
                #Extract current line
                head -n ${FileLine} ./Round_${CurrRound}/01_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.bed | tail -n 1  > ./Round_${CurrRound}/Temp.bed

                #Check if there is no overlap
                if [[ $(bedtools intersect -a ./Round_${CurrRound}/Temp.bed -b  ./Round_${CurrRound}/02_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.bed | wc -l ) -eq 0 ]];then
                        cat ./Round_${CurrRound}/Temp.bed >>  ./Round_${CurrRound}/02_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.bed
                        CounterSeqs=$((${CounterSeqs}+1))
                fi

                #Increase FileLine
                FileLine=$((${FileLine}+1))

                if [[ ${FileLine} -gt ${TotalLines} ]];then
                        break
                fi
            done

            #Perform Checks, save log just in case
            ##Check if it has enough sequences to continue 
            #if [ $( wc -l < ./Round_${CurrRound}/02_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.bed ) -lt 10 ];then
            #    echo -e "\nNot enough sequences to continue"
            #    break
            #fi

            #Cross-check with previous round 
            if [[ ${LowerCase} == "Left" ]];then
                bedtools intersect \
                    -a ./Round_${CurrRound}/02_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.bed \
                    -b ./Round_${PrevRound}/02_CurrentConsensi.${LowerCase}.Round${PrevRound}.Extended.bed \
                    -wb | \
                awk -v Increase=${RoundExt} -v OFMT='%f' '{OFS="\t"}{
                    if($6=="+"){
                        if($2>$8 || sqrt(($8-$2)*($8-$2))< Increase ){
                            $2=$8-Increase;
                            print            
                        }else{
                            print
                        }
                    }else{
                        if($3<$9 || sqrt(($9-$3)*($9-$3))<Increase ){
                            $3=$9+Increase;
                            print
                        }else{
                            print
                        }
                    }}' | \
                cut -f1-6 > ./Round_${CurrRound}/02_Temp.bed
            else
                bedtools intersect \
                    -a ./Round_${CurrRound}/02_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.bed \
                    -b ./Round_${PrevRound}/02_CurrentConsensi.${LowerCase}.Round${PrevRound}.Extended.bed \
                    -wb | \
                awk -v Increase=${RoundExt} -v OFMT='%f' '{OFS="\t"}{
                    if($6=="+"){
                        if($3<$9 || sqrt(($9-$3)*($9-$3))< Increase ){
                            $3=$9+Increase;
                            print            
                        }else{
                            print
                        }
                    }else{
                        if($2>$8 || sqrt(($8-$2)*($8-$2))<Increase ){
                            $2=$8-Increase;
                            print
                        }else{
                            print
                        }
                    }}' | \
                cut -f1-6 > ./Round_${CurrRound}/02_Temp.bed
            fi

            #Print the rest
            bedtools intersect \
                -a ./Round_${CurrRound}/02_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.bed \
                -b ./Round_${PrevRound}/02_CurrentConsensi.${LowerCase}.Round${PrevRound}.Extended.bed \
                -v >> ./Round_${CurrRound}/02_Temp.bed

            cp ./Round_${CurrRound}/02_Temp.bed ./Round_${CurrRound}/02_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.bed
            
            #Get fasta sequences
            bedtools getfasta \
                    -fi ../Staging/${Genome} \
                    -fo  ./Round_${CurrRound}/03_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.fa \
                    -bed ./Round_${CurrRound}/02_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.bed \
                    -s
            #Transform Ns into gaps
            awk '{if($0 ~ "^>"){print $0}else{gsub("N","-",$0); print}}' \
                ./Round_${CurrRound}/03_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.fa \
            > ./Round_${CurrRound}/03_TEMP.fa

            cp ./Round_${CurrRound}/03_TEMP.fa ./Round_${CurrRound}/03_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.fa

            #Add ID about Round
            awk -v RoundL=${CurrRound} \
            -v Side=${LowerCase} \
            '{ 
                    if($1 ~ "^>"){
                    gsub("$","___"Side"R"RoundL,$1);
                    print
                    } else{print}
            }' \
            ./Round_${CurrRound}/03_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.fa > ./Round_${CurrRound}/Temp.fa
        
            cp ./Round_${CurrRound}/Temp.fa ./Round_${CurrRound}/03_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.fa

            #One representative per group of near-identical copies (segmental duplications)
            bash ${UtilsDir}/CollapseNearIdentical.sh \
                ./Round_${CurrRound}/03_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.fa \
                ${CPU} \
                ${maxIdentity}

            #Check if there is enough sequences to continue
            if [ $( grep -c ">" ./Round_${CurrRound}/03_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.fa ) -lt 5 ];then
                echo -e "\nNot enough sequences to continue"
                RoundFailed=True
            fi

            if [[ ${RoundFailed} == "False" ]];then
            #Make alignment of only selected sequences
            mafft \
                --genafpair \
                --thread ${CPU} \
                --maxiterate 1000 \
                ./Round_${CurrRound}/03_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.fa \
            > ./Round_${CurrRound}/04_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.fa.aln \
            2> ./Round_${CurrRound}/04_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.MAFFTlog

            #Curate 
            ##Calculate MinCoverage threshold
            MinCoverageThresh=$( grep -c ">" ./Round_${CurrRound}/04_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.fa.aln | awk '{if($1>10){print int($1*0.30)}else{print 2}}')

            python3 ../Staging/${CurateScript} \
                    --input ./Round_${CurrRound}/04_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.fa.aln \
                    --output ./Round_${CurrRound}/05_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFamily} \
                    --SeqID ${RepeatID}#___Round${LetterCase}${CurrRound} \
                    --MinCoverage ${MinCoverageThresh} \
                    --StepWindow 25 \
                    --identityExp ${curationIdentityExp} \
                    --gapExp ${curationGapExp} \
                    --MaxThresh ${curationMaxThresh} \
                    --firstPeak \
                    --Plots

            #Rename output 
            mv ./Round_${CurrRound}/05_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFamily}__Peak0.aln.fa \
                ./Round_${CurrRound}/05_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.aln.fa

            #Make all sequences uppercase
            seqkit seq \
                    -u \
                    ./Round_${CurrRound}/05_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.aln.fa \
                    > ./Round_${CurrRound}/05_TEMP.fa.aln

            cp ./Round_${CurrRound}/05_TEMP.fa.aln \
                ./Round_${CurrRound}/05_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.aln.fa

            #Trim the edge of the extension that is not conserved across copies (random flank),
            #then make the hmmer model with match columns defined by RF, and its consensus
            BuildHandHMM \
                ./Round_${CurrRound}/05_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.aln.fa \
                ./Round_${CurrRound}/05_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended \
                ${RepeatID}__${CurrRound} \
                ${LowerCase} || RoundFailed=True
            fi

            if [[ ${RoundFailed} == "False" ]];then
            ##Make mafft alignment of the previous consensus
            mafft \
                --localpair \
                --maxiterate 1000 \
                --thread ${CPU} \
                <( cat ./Round_${PrevRound}/05_CurrentConsensi.${LowerCase}.Round${PrevRound}.Extended.Consensus.fa \
                    ./Round_${CurrRound}/05_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.Consensus.fa ) \
                > ./Round_${CurrRound}/08_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.aln

            ../Staging/${CoordExtractScript} \
                --input ./Round_${CurrRound}/08_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.aln \
                --output ./Round_${CurrRound}/08_CurrentConsensi.${LowerCase}.Round${CurrRound}.ExtendedSide.Coord \
                --Side ${LowerCase} \
                --BasesOverl 20

            #Get extended size: new consensus bases beyond the previous model, after trimming the unsupported edge
            #(bases of the previous model before the anchor were already in the model: not counted, not stitched twice)
            ExtendedSize=$( python3 ${UtilsDir}/ExtensionSize.py \
                --input ./Round_${CurrRound}/08_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.aln \
                --side ${LowerCase} \
                --minRun 20 )
            if ! [[ ${ExtendedSize} =~ ^[0-9]+$ ]];then
                echo -e "\nCould not measure the extension of ${LowerCase} round ${CurrRound}, exiting"
                exit 1
            fi

            #Genome-anchored support: genomic copies of the model's outer anchorLength bases must continue colinearly
            #into the new bases. A real end shows as a drop in the fraction of anchored copies that continue,
            #drift into unique flank or a neighbouring repeat as new bases reached by few anchored copies.
            #A round is low when not all its new bases are supported: keep the supported bases (or fail the round if
            #too few), and stop this side at the second low round in a row (single low rounds happen inside repeats).
            AnchoredEnd=False
            AnchoredStop=False
            if [[ ${ExtendedSize} -ge ${Increase} ]] && awk -v r=${anchorMinRatio} 'BEGIN{exit !(r>0)}';then
                AnchoredPrefix=./Round_${CurrRound}/09_CurrentConsensi.${LowerCase}.Round${CurrRound}.Anchored
                python3 ${UtilsDir}/AnchoredSupport.py build \
                    --round ./Round_${CurrRound}/05_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.Consensus.fa \
                    --merged ./Round_${PrevRound}/07_CurrentConsensi.${LowerCase}.Round${PrevRound}.Extended.Curated.Consensus.fa \
                    --side ${LowerCase} \
                    --extension ${ExtendedSize} \
                    --anchor ${anchorLength} \
                    --output ${AnchoredPrefix}.fa
                nhmmer \
                    --cpu ${CPU} \
                    --dna \
                    --notextw \
                    --noali \
                    -E 1e-5 \
                    --tblout ${AnchoredPrefix}.tblout \
                    ${AnchoredPrefix}.fa \
                    ../Staging/HMMDB_Genome > /dev/null
                python3 ${UtilsDir}/AnchoredSupport.py score \
                    --query ${AnchoredPrefix}.fa \
                    --tbl ${AnchoredPrefix}.tblout \
                    --side ${LowerCase} \
                    --extension ${ExtendedSize} \
                    --span ${anchorSpan} \
                    --minCopies ${anchorMinCopies} \
                    --minRatio ${anchorMinRatio} > ${AnchoredPrefix}.txt
                AnchoredSupported=$( awk '$1=="SUPPORTED"{print $2}' ${AnchoredPrefix}.txt )
                AnchoredCopies=$( awk '$1=="ANCHORED"{print $2}' ${AnchoredPrefix}.txt )
                echo -e "\nAnchored support: $( tr '\n' ' ' < ${AnchoredPrefix}.txt)of ${ExtendedSize} new bases"
                #A jump in the copies of the model edge (more than anchorMaxJump times the lower count of the
                #two previous rounds) means the edge now lies in a more abundant neighbouring repeat
                AnchoredLowest=$( echo ${AnchoredHistory} | tr ' ' '\n' | tail -2 | sort -n | head -1 )
                AnchoredHistory="${AnchoredHistory} ${AnchoredCopies}"
                if [[ $( echo ${AnchoredHistory} | wc -w ) -ge 3 ]] && [[ ${AnchoredCopies} -ge 20 ]] && \
                    awk -v c=${AnchoredCopies} -v l=${AnchoredLowest} -v j=${anchorMaxJump} 'BEGIN{exit !(c > j*l)}';then
                    echo -e "\nCopies of the model edge jump from ${AnchoredLowest} to ${AnchoredCopies}: a more abundant repeat next to this one"
                    AnchoredStop=True
                    RoundFailed=True
                elif [[ ${AnchoredSupported} -ge $(( ExtendedSize - 10 )) ]];then
                    AnchoredLowStreak=0
                else
                    AnchoredLowStreak=$(( AnchoredLowStreak + 1 ))
                    if [[ ${AnchoredLowStreak} -ge 2 ]];then
                        echo -e "\nSecond round in a row not supported by copies of the model edge: end of the repeat"
                        AnchoredStop=True
                        RoundFailed=True
                    elif [[ ${AnchoredSupported} -lt ${Increase} ]];then
                        echo -e "\nThe new bases are not supported by copies of the model edge"
                        RoundFailed=True
                    else
                        AnchoredEnd=True
                    fi
                fi
            fi

            if [[ ${ExtendedSize} -lt ${Increase} ]];then
                echo -e "\nOnly ${ExtendedSize} supported bases added (minimum ${Increase})"
                RoundFailed=True
            elif [[ ${RoundFailed} == "False" ]];then
                #Add this round to the merged model along the alignment of its consensus to the merged consensus
                python3 ${StitchScript} \
                    --merged ./Round_${PrevRound}/07_CurrentConsensi.${LowerCase}.Round${PrevRound}.Extended.Curated.aln.fa \
                    --round ./Round_${CurrRound}/05_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.match.aln.fa \
                    --extension ${ExtendedSize} \
                    --side ${LowerCase} \
                    --output ./Round_${CurrRound}/07_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.Curated \
                    --SeqID ${RepeatID}#___Round${LetterCase}${CurrRound} | tee ./Round_${CurrRound}/07_Stitch.log
                [ ${PIPESTATUS[0]} -eq 0 ] || RoundFailed=True

                #Low anchored support: drop the new bases past it
                if [[ ${AnchoredEnd} == "True" ]] && [[ ${RoundFailed} == "False" ]];then
                    NewPositions=$( awk '$1=="STITCH_NEW_POSITIONS"{print $2}' ./Round_${CurrRound}/07_Stitch.log )
                    echo -e "\nLow anchored support: keeping ${AnchoredSupported} of ${NewPositions:-?} new positions"
                    python3 ${UtilsDir}/TrimEdge.py \
                        --prefix ./Round_${CurrRound}/07_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.Curated \
                        --side ${LowerCase} \
                        --positions $(( ${NewPositions:-${ExtendedSize}} - AnchoredSupported ))
                fi

                #Next search window: the new bases, and at least ExtensionSize in from the edge,
                #so the next round is not left with only a handful of copies after a small extension
                SetOuterCoord \
                    ./Round_${CurrRound}/05_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.Consensus.fa \
                    ${LowerCase} \
                    $(( ExtendedSize > ExtensionSize ? ExtendedSize : ExtensionSize )) \
                    ./Round_${CurrRound}/08_CurrentConsensi.${LowerCase}.Round${CurrRound}.ExtendedSide.Coord
            fi
            fi

            if [[ ${RoundFailed} == "True" ]];then
                NextLevel=$(( RetryLevel + 1 ))

                #Keep the failed attempt for inspection
                FailedDir=./Failed_Round${CurrRound}_Retry${RetryLevel}
                while [ -d ${FailedDir} ];do
                    FailedDir=${FailedDir}_x
                done
                mv ./Round_${CurrRound} ${FailedDir}

                if [[ ${AnchoredStop} == "True" ]];then
                    echo -e "\nStopping the ${LowerCase} side at the end of the repeat"
                    BreakNow=True
                    break
                fi

                if [[ ${NextLevel} -gt ${maxRetries} ]];then
                    echo -e "\nNo supported extension after ${maxRetries} retries, stopping the ${LowerCase} side"
                    BreakNow=True
                    break
                fi

                #Restart this round from the last good model, searching from its outer edge
                mkdir -p ./Round_${CurrRound}
                for File in 02_CurrentConsensi.${LowerCase}.RoundXX.Extended.bed \
                            05_CurrentConsensi.${LowerCase}.RoundXX.Extended.aln.fa \
                            05_CurrentConsensi.${LowerCase}.RoundXX.Extended.hmm \
                            05_CurrentConsensi.${LowerCase}.RoundXX.Extended.Consensus.fa \
                            07_CurrentConsensi.${LowerCase}.RoundXX.Extended.Curated.aln.fa \
                            07_CurrentConsensi.${LowerCase}.RoundXX.Extended.Curated.Consensus.fa \
                            07_CurrentConsensi.${LowerCase}.RoundXX.Extended.Curated.Columns.txt;do
                    [ -f ./Round_${PrevRound}/${File/XX/${PrevRound}} ] && \
                        cp ./Round_${PrevRound}/${File/XX/${PrevRound}} ./Round_${CurrRound}/${File/XX/${CurrRound}}
                done

                RetrySettings ${NextLevel}
                SetOuterCoord \
                    ./Round_${CurrRound}/05_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.Consensus.fa \
                    ${LowerCase} \
                    ${RoundExt} \
                    ./Round_${CurrRound}/08_CurrentConsensi.${LowerCase}.Round${CurrRound}.ExtendedSide.Coord
                echo ${NextLevel} > ./Round_${CurrRound}/NextRetryLevel
            else
                echo 0 > ./Round_${CurrRound}/NextRetryLevel
            fi

            #Check if it has reached max rounds
            if [ ${CurrRound} -ge ${maxRounds} ];then
                BreakNow=True
                break
            fi

            #Make check for this round
            > ./Round_${CurrRound}/GoodCHECK

            #Increase round number
            CurrRound=$(echo ${CurrRound} | awk '{printf "%02d\n", $1+1}')

        done

        #Wrap-up family
        if [ $(ls Round_*/GoodCHECK | wc -l ) -ge 2 ] && [ ! -f ./GOOD.FAM.CHECK ];then
            #Make Final dir
            ##Find which round was last
            HighestRound=$(ls Round_*/GoodCHECK | sed 's/\/GoodCHECK//;s/Round_//;s/.*\///' | sort -n | tail -1)

            ##Make symbolic link
            ln -s ./Round_${HighestRound}/ ./Final

            ##Make check
            > ./GOOD.FAM.CHECK
        elif [ ! -f ./GOOD.FAM.CHECK ];then
            > ./BAD.FAM.CHECK
        fi 
    fi

    #Increase Family Number
    Family=$( echo ${Family} | awk '{printf "%02d\n", $1+1}' )

    NumFamilies=$(ls -d ../Family_* | wc -l)

    #Return to the families dir
    cd ..
done
}

#############
# LEFT SIDE  #

#Create directory for repeat
mkdir -p ./Left
cd ./Left

#Run the extension
if [ ! -f ./LEFT.DONE.CHECK ];then
    HMMExtensionSides Left

    > ./LEFT.DONE.CHECK
fi

# RIGHT SIDE  #

#Create directory for repeat
cd ${NewPath}
mkdir -p ./Right
cd ./Right

#Run the extension
if [ ! -f ./RIGHT.DONE.CHECK ];then
    
    HMMExtensionSides Right

    > ./RIGHT.DONE.CHECK
fi 

######################
##### MERGE SIDE  ####
######################
#Create directory for results
cd ${NewPath}
mkdir -p ./Merge
cd ./Merge

if [ ! -f ./MERGE.GOOD.CHECK ];then
    if [ -f ../Left/Family_00/GOOD.FAM.CHECK ];then
        if [ -f ../Right/Family_00/GOOD.FAM.CHECK ];then
            HighestRound=$( ls ../Right/Family_00/Round_*/GoodCHECK | sed "s/\/GoodCHECK//;s/Right\///;s/.*\///"| tail -1 )
            cp ../Right/Family_00/${HighestRound}/07_CurrentConsensi.Right*Extended.Curated.aln.fa \
                ./Merge.aln.fa
            for File in ../Right/Family_00/${HighestRound}/07_CurrentConsensi.Right*Extended.Curated.Columns.txt;do [ -f ${File} ] && cp ${File} ./Merge.Columns.txt; done
            cp ../Right/Family_00/${HighestRound}/07_CurrentConsensi.Right*Extended.Curated.Consensus.fa \
                ./Merge.Consensus.fa
        else
            HighestRound=$( ls ../Left/Family_00/Round_*/GoodCHECK | sed "s/\/GoodCHECK//;s/Left\///;s/.*\///"| tail -1 )
            cp ../Left/Family_00/${HighestRound}/07_CurrentConsensi.Left*Extended.Curated.aln.fa \
                ./Merge.aln.fa
            for File in ../Left/Family_00/${HighestRound}/07_CurrentConsensi.Left*Extended.Curated.Columns.txt;do [ -f ${File} ] && cp ${File} ./Merge.Columns.txt; done
            cp ../Left/Family_00/${HighestRound}/07_CurrentConsensi.Left*Extended.Curated.Consensus.fa \
                ./Merge.Consensus.fa
        fi
    elif [ -f ../Right/Family_00/GOOD.FAM.CHECK ];then
        HighestRound=$( ls ../Right/Family_00/Round_*/GoodCHECK | sed "s/\/GoodCHECK//;s/Right\///;s/.*\///"| tail -1 )
        cp ../Right/Family_00/${HighestRound}/07_CurrentConsensi.Right*Extended.Curated.aln.fa \
            ./Merge.aln.fa
        for File in ../Right/Family_00/${HighestRound}/07_CurrentConsensi.Right*Extended.Curated.Columns.txt;do [ -f ${File} ] && cp ${File} ./Merge.Columns.txt; done
        cp ../Right/Family_00/${HighestRound}/07_CurrentConsensi.Right*Extended.Curated.Consensus.fa \
            ./Merge.Consensus.fa
    fi

    #Genome-wide refinement of the final model: hits aligned with hmmalign, match columns by occupancy (--hand)
    if [ ${refineIterations} -gt 0 ] && [ -f ./Merge.aln.fa ];then
        echo -e "\nRefining the final model genome-wide (${refineIterations} iterations at most)"
        bash ${UtilsDir}/RefineModel.sh \
            -a ./Merge.aln.fa \
            -d ${oldWorkDir}/${HMMRDB} \
            -g ${oldWorkDir}/${Genome} \
            -o ./${RepeatID}.Refined \
            -n ${RepeatID} \
            -t ${CPU} \
            -I ${refineIterations} \
            -Q ${rfMinOcc} \
            -D ${maxIdentity} || echo "Refinement failed, the unrefined model is in Merge.aln.fa"
    fi

    #Make check
    > ./MERGE.GOOD.CHECK
fi

#Make check
> ./MERGE.GOOD.CHECK


#Make outputs
cd ${NewPath}

#Copy to original dir
ln -s ${PWD}/Left/LEFT.DONE.CHECK ${oldWorkDir}/
ln -s ${PWD}/Right/RIGHT.DONE.CHECK ${oldWorkDir}/
ln -s ${PWD}/Merge/MERGE.GOOD.CHECK ${oldWorkDir}/

#Go back to original dir
cd ${oldWorkDir} 