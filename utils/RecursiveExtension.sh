#!/bin/bash

while getopts "i:l:r:c:t:e:g:h:x:o:C:O:s:m:p:S:R:F:T:P:V:" opt
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
   esac
done
    
#Save nf directory into a variable
oldWorkDir=$(pwd)

#Move to new WorkDir 
cd ${outDir}/WorkDir/${RepeatID}

#Backup new path 
NewPath=$( pwd )

#Initialize some some variables
ExtensionSize=${Extension}
AnchorSize=$(echo ${ExtensionSize} | awk '{printf "%.0f\n", $1*0.5}') 
Increase=$(  echo ${ExtensionSize} | awk '{$0=$0*0.25;printf("%d\n",$0+=$0<0?0:0.9)}' )
LongConsensus=False

######################
#####  Function  #####
######################

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
            else
                #Check if left finished good
                if [ -f ${LeftDir}/Family_${Family}/GOOD.FAM.CHECK ];then
                        #Determine highest good left round
                        HighestRound=$( ls ${LeftDir}/Family_${Family}/Round_*/GoodCHECK | sed "s/\/GoodCHECK//;s/Left\///;s/.*\///"| tail -1 )
                        cp ${LeftDir}/Family_${Family}/${HighestRound}/07_CurrentConsensi.Left*Extended.Curated.aln.fa \
                            ./Round_00/07_CurrentConsensi.${LowerCase}.Round00.Extended.Curated.aln.fa
                        cp ${LeftDir}/Family_${Family}/${HighestRound}/07_CurrentConsensi.Left*Extended.Curated.Consensus.fa \
                            ./Round_00/07_CurrentConsensi.${LowerCase}.Round00.Extended.Curated.Consensus.fa

                        #Calculate extension size
                        ExtSize=$( cat ${LeftDir}/Family_${Family}/Round_00/08_CurrentConsensi.Left.Round00.ExtendedSide.Coord | \
                        awk '{print $2 - $1}' )    

                        #Check size of LEFT extension
                        ModelLen=$( seqkit fx2tab --length --name ${LeftDir}/Family_${Family}/${HighestRound}/07_CurrentConsensi.Left*Extended.Curated.Consensus.fa | cut -f2)

                        #Generate Coord for Right
                        echo -e "$(( ${ModelLen} - ${ExtSize} ))\t${ModelLen}" > ./Round_00/08_CurrentConsensi.${LowerCase}.Round00.ExtendedSide.Coord
                else
                    cp ${CentralDir}/Family_${Family}/Final/FinalConsensi.aln.fa \
                        ./Round_00/07_CurrentConsensi.${LowerCase}.Round00.Extended.Curated.aln.fa
                    cp ${CentralDir}/Family_${Family}/Final/FinalConsensi.consensus.fa \
                        ./Round_00/07_CurrentConsensi.${LowerCase}.Round00.Extended.Curated.Consensus.fa
                    cp ${CentralDir}/Family_${Family}/Final/Final.Consensus.${LowerCase}Side.Coord \
                        ./Round_00/08_CurrentConsensi.${LowerCase}.Round00.ExtendedSide.Coord
                fi

            fi
        
            ##Make HMM model 
            hmmbuild \
            --symfrac 0 \
            --dna \
            --cpu ${CPU} \
            --informat afa \
            --seed 1992 \
            --fragthresh 1.0 \
            -n ${RepeatID} \
            --wnone \
            ./Round_00/08_CurrentConsensi.${LowerCase}.Round00.Extended.Curated.hmm \
            ./Round_00/07_CurrentConsensi.${LowerCase}.Round00.Extended.Curated.aln.fa

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

        #Main loop
        while [ ${BreakNow} == "False" ] && [ ${SkipRounds} == "No" ];do
            #Make dir for this round
            mkdir -p Round_${CurrRound}

            #Set previous round
            PrevRound=$(echo ${CurrRound} | awk '{printf "%02d\n", $1-1}')

            #HMMER Search
            nhmmer \
            --cpu ${CPU} \
            --tblout ./Round_${CurrRound}/00_CurrentConsensi.${LowerCase}.Round${CurrRound}.hmmout \
            --notextw \
            --noali \
            --seed 1992 \
            ./Round_${PrevRound}/08_CurrentConsensi.${LowerCase}.Round${PrevRound}.Extended.Curated.hmm \
            ../Staging/HMMDB_Genome  > /dev/null

            ErrorCode=$?
            if [ ${ErrorCode} -ne 0 ];then
                echo -e "\nHMMER search failed, exiting"
                exit ${ErrorCode}
            fi
            
            #Format HMMER output
            sed '/#/d; s/\s\{1,\}/\t/g' ./Round_${CurrRound}/00_CurrentConsensi.${LowerCase}.Round${CurrRound}.hmmout > ./Round_${CurrRound}/00_CurrentConsensi.${LowerCase}.Round${CurrRound}.hmmout.tmp
            mv ./Round_${CurrRound}/00_CurrentConsensi.${LowerCase}.Round${CurrRound}.hmmout.tmp \
                ./Round_${CurrRound}/00_CurrentConsensi.${LowerCase}.Round${CurrRound}.hmmout

            #Filter hits with high bias
            awk '{OFS="\t"} ($15/$14) <= 0.20{print}' ./Round_${CurrRound}/00_CurrentConsensi.${LowerCase}.Round${CurrRound}.hmmout > ./Round_${CurrRound}/00_Temp.hmmout
            mv ./Round_${CurrRound}/00_Temp.hmmout ./Round_${CurrRound}/00_CurrentConsensi.${LowerCase}.Round${CurrRound}.hmmout

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
                    -l ${ExtensionSize} \
                    -r ${AnchorSize} \
                    -i ./Round_${CurrRound}/00_CurrentConsensi.${LowerCase}.Round${CurrRound}.Raw.bed \
                    -g <(cut -f1,2  ../Staging/${GenomeIndex}) > ./Round_${CurrRound}/01_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.bed
            else
                bedtools slop \
                    -s \
                    -r ${ExtensionSize} \
                    -l ${AnchorSize} \
                    -i ./Round_${CurrRound}/00_CurrentConsensi.${LowerCase}.Round${CurrRound}.Raw.bed \
                    -g <(cut -f1,2  ../Staging/${GenomeIndex}) > ./Round_${CurrRound}/01_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.bed
            fi

            #Select Sequences to be used
            CounterSeqs=0
            FileLine=1
            TotalLines=$(wc -l < ./Round_${CurrRound}/01_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.bed )
            > ./Round_${CurrRound}/02_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.bed

            while [ $CounterSeqs -lt ${SampleSeqs} ];do
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
            if [ $( wc -l < ./Round_${CurrRound}/02_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.bed ) -lt 10 ];then
                echo -e "\nNot enough sequences to continue"
                break
            fi

            #Get fasta sequences
            bedtools getfasta \
                    -fi ../Staging/${Genome} \
                    -fo  ./Round_${CurrRound}/03_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.fa \
                    -bed ./Round_${CurrRound}/02_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.bed \
                    -s
            
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
        
            mv ./Round_${CurrRound}/Temp.fa ./Round_${CurrRound}/03_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.fa

            #Make alignment of only selected sequences
            mafft \
            --localpair \
            --maxiterate 1000 \
            --thread ${CPU} \
            ./Round_${CurrRound}/03_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.fa \
            > ./Round_${CurrRound}/04_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.aln 

            #Check for families
            ../Staging/${OutlierScript} \
            --input ./Round_${CurrRound}/04_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.aln \
            --output ./Round_${CurrRound}/04_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.Fams

            #Check if there are multiple families
            if [ $(ls ./Round_${CurrRound}/04_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.Fams.Cluster*fa.aln | wc -l ) -gt 1 ];then
                #Remove those that only contain less than N sequences
                CurrFamily=00
                for FamilyAln in $(ls ./Round_${CurrRound}/04_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.Fams.Cluster*fa.aln );do
                        if [[ $(grep -c ">" ${FamilyAln}) -ge 7 ]];then
                            cp ${FamilyAln} ./Round_${CurrRound}/04_00_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFamily}.aln.fa
                        fi
                        CurrFamily=$( echo ${CurrFamily} | awk '{printf "%02d\n", $1+1}' )
                done
                    
                #Make each family uppercase
                for FamilyAln in $(ls ./Round_${CurrRound}/04_00_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family*aln.fa );do
                        #Make all sequences uppercase
                        seqkit seq \
                        -u \
                        ${FamilyAln} \
                        > ./Round_${CurrRound}/04_Temp.aln
                        
                        mv ./Round_${CurrRound}/04_Temp.aln \
                        ${FamilyAln}
                done

                TotalFamilies=$(ls ../Family_* | wc -l)
                if [[ ${TotalFamilies} -ge ${maxFamilies} ]];then
                    #Move the families that wont be used
                    for NotUsedFam in $(seqkit stats -T ./Round${CurrRound}/04_00_*aln.fa| sed '1d' | sort -k4| sed '1d'| cut -f1);do
                        mv ${NotUsedFam} $( echo ${NotUsedFam} | sed 's/04_00/04_99/' )
                    done
                fi

                #Clean and expand each family
                CurrFamily_Num=1
                for FamilyAln in $( ls ./Round_${CurrRound}/04_00_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family*aln.fa );do
                    #Get family number
                    CurrFamily=$( echo ${FamilyAln} | \
                            sed 's/.*Family//;s/.aln.fa//' )

                    #Clean alignment
                    ../Staging/${CurateScript} \
                    --input ${FamilyAln} \
                    --output ./Round_${CurrRound}/04_01_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFamily} \
                    --SeqID ${RepeatID}#___Round${LetterCase}${CurrRound} \
                    --MinCoverage 1 \
                    --Plots

                    PeakNumLoop=0
                    for PeakAln in $( ls ./Round_${CurrRound}/04_01_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFamily}__Peak*aln.fa );do
                        if [[ $(seqkit stats -T ${PeakAln} | cut -f8 | tail -1) -gt 0 ]];then 
                            #Create hmmer model
                            hmmbuild \
                            --symfrac 0 \
                            --dna \
                            --informat afa \
                            --seed 1992 \
                            --fragthresh 1.0 \
                            --wnone \
                            -n Peak${PeakNumLoop} \
                            ./Round_${CurrRound}/04_02_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFamily}_Peak${PeakNumLoop}.hmm \
                            ${PeakAln}

                            #Search for the model in the genome
                            nhmmer \
                            --cpu ${CPU} \
                            --tblout ./Round_${CurrRound}/04_03_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFamily}_Peak${PeakNumLoop}.tblout \
                            --notextw \
                            --noali \
                            --seed 1992 \
                            -E 1e-5 \
                            ./Round_${CurrRound}/04_02_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFamily}_Peak${PeakNumLoop}.hmm \
                            ../Staging/HMMDB_Genome > /dev/null

                            #Cleanup and remove high bias
                            sed '/#/d; s/\s\{1,\}/\t/g' ./Round_${CurrRound}/04_03_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFamily}_Peak${PeakNumLoop}.tblout | \
                            awk '{OFS="\t"} ($15/$14) <= 0.20{print}' > ./Round_${CurrRound}/04_03_Temp.tblout 
                            mv ./Round_${CurrRound}/04_03_Temp.tblout ./Round_${CurrRound}/04_03_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFamily}_Peak${PeakNumLoop}.tblout
                        else
                            > ./Round_${CurrRound}/04_03_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFamily}_Peak${PeakNumLoop}.tblout
                            > ./Round_${CurrRound}/04_02_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFamily}_Peak${PeakNumLoop}.hmm
                        fi

                        #Increase peak number
                        PeakNumLoop=$((${PeakNumLoop}+1))
                    done

                    #Select Peak with best hits
                    PeakNumLoop=0
                    HighestMedian=0
                    SelectedPeak=0
                    for PeakHMMOut in $(ls ./Round_${CurrRound}/04_03_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFamily}_Peak*.tblout );do
                        #Calculate median of bitscores
                        MedianScores=$( sed 's/\s\{1,\}/\t/g' ${PeakHMMOut} | \
                                        cut -f14 | \
                                        sort -n | \
                                        awk '{a[i++]=$1} END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]}' )

                            #Check if median is higher than the highest median
                        if (( $(echo "${MedianScores} > ${HighestMedian}" | bc -l) )) ;then
                            HighestMedian=${MedianScores}
                            SelectedPeak=${PeakNumLoop}
                        fi
                        PeakNumLoop=$((${PeakNumLoop}+1))
                    done

                    #Select the best peak and hits
                    cp ./Round_${CurrRound}/04_01_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFamily}__Peak${SelectedPeak}.aln.fa \
                    ./Round_${CurrRound}/04_04_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFamily}.aln.fa

                    cp ./Round_${CurrRound}/04_03_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFamily}_Peak${SelectedPeak}.tblout \
                    ./Round_${CurrRound}/04_04_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFamily}.tblout

                    cp ./Round_${CurrRound}/04_01_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFamily}.Consensus__Peak${SelectedPeak}.fa \
                    ./Round_${CurrRound}/04_04_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFamily}.Consensus.fa

                    #Select sequences that cover 80% of the model
                    ModelLen=$(seqkit stats -T ./Round_${CurrRound}/04_04_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFamily}.aln.fa | \
                            cut -f7 | \
                            tail -1 | \
                            sed 's/\..*//' )

                    bedtools intersect \
                    -a <(echo -e "Target\t0\t${ModelLen}") \
                    -b <( sed '/#/d' ./Round_${CurrRound}/04_04_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFamily}.tblout |  awk '{OFS="\t"}{print "Target",$5,$6,$1,$9,$10,$12,$14}' ) \
                    -wb \
                    -f 0.8 | \
                    awk '{OFS="\t"}{ 
                    if( $10 == "+"){print $7,$8,$9 ,".",$11,$10}else{print $7,$9,$8,".",$11,$10}
                    }' > ./Round_${CurrRound}/04_05_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFamily}.bed

                    #Select sequences to use
                    CounterSeqs=0
                    FileLine=1
                    TotalLines=$(wc -l < ./Round_${CurrRound}/04_05_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFamily}.bed )
                    > ./Round_${CurrRound}/04_06_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFamily}.bed

                    while [[ $CounterSeqs -lt ${SampleSeqs} ]];do
                        #Extract current line
                        head -n ${FileLine} ./Round_${CurrRound}/04_05_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFamily}.bed | tail -n 1  > ./Round_${CurrRound}/04_Temp.bed

                            #Check if there is no overlap
                        if [[ $(bedtools intersect -a ./Round_${CurrRound}/04_Temp.bed -b ./Round_${CurrRound}/04_06_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFamily}.bed | wc -l ) -eq 0 ]];then
                            cat ./Round_${CurrRound}/04_Temp.bed >>  ./Round_${CurrRound}/04_06_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFamily}.bed
                            CounterSeqs=$((${CounterSeqs}+1))
                        fi

                        #Increase FileLine
                        FileLine=$((${FileLine}+1))

                        if [[ ${FileLine} -gt ${TotalLines} ]];then
                            break
                        fi
                    done

                    #Extract sequences
                    bedtools getfasta \
                    -fi ../Staging/${Genome} \
                    -fo  ./Round_${CurrRound}/04_07_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFamily}.fa \
                    -bed ./Round_${CurrRound}/04_06_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFamily}.bed \
                    -s
                    
                    #Add ID about Round
                    awk -v RoundL=${CurrRound} \
                    -v Side=${LowerCase} \
                    '{ 
                            if($1 ~ "^>"){
                            gsub("$","___"Side"R"RoundL,$1);
                            print
                            } else{print}
                    }' \
                    ./Round_${CurrRound}/04_07_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFamily}.fa > ./Round_${CurrRound}/04_Temp.fa
                    
                    mv ./Round_${CurrRound}/04_Temp.fa ./Round_${CurrRound}/04_07_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFamily}.fa

                    #Make alignment for this round
                    mafft \
                    --localpair \
                    --maxiterate 1000 \
                    --thread ${CPU} \
                    ./Round_${CurrRound}/04_07_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFamily}.fa \
                    > ./Round_${CurrRound}/04_08_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFamily}.aln

                    #Make all sequences uppercase
                    seqkit seq \
                    -u \
                    ./Round_${CurrRound}/04_08_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFamily}.aln \
                    > ./Round_${CurrRound}/04_08_Temp.aln

                    mv ./Round_${CurrRound}/04_08_Temp.aln ./Round_${CurrRound}/04_08_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFamily}.aln

                    #Curate alignment and make consensus
                    ../Staging/${CurateScript} \
                    --input ./Round_${CurrRound}/04_08_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFamily}.aln \
                    --output ./Round_${CurrRound}/04_09_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFamily} \
                    --SeqID ${RepeatID}#___Round${LetterCase}${CurrRound} \
                    --NoiseThresh ${noiseThreshold} \
                    --Plots

                    #Search each peak in genome
                    PeakNumLoop=0

                    for PeakAln in $(ls ./Round_${CurrRound}/04_09_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFamily}__Peak*aln.fa );do
                        if [[ $(seqkit stats -T ${PeakAln} | cut -f8 | tail -1) -gt 0 ]];then 
                            #Create hmmer model
                            hmmbuild \
                            --symfrac 0 \
                            --dna \
                            --informat afa \
                            --seed 1992 \
                            --fragthresh 1.0 \
                            --wnone \
                            -n Peak${PeakNumLoop} \
                            ./Round_${CurrRound}/04_10_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFamily}_Peak${PeakNumLoop}.hmm \
                            ${PeakAln}

                            #Search for the model in the genome
                            nhmmer \
                            --cpu ${CPU} \
                            --tblout ./Round_${CurrRound}/04_11_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFamily}_Peak${PeakNumLoop}.tblout \
                            --notextw \
                            --noali \
                            --seed 1992 \
                            -E 1e-5 \
                            ./Round_${CurrRound}/04_10_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFamily}_Peak${PeakNumLoop}.hmm \
                            ../Staging/HMMDB_Genome > /dev/null

                            #Cleanup and remove high bias
                            sed '/#/d; s/\s\{1,\}/\t/g' ./Round_${CurrRound}/04_11_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFamily}_Peak${PeakNumLoop}.tblout | \
                            awk '{OFS="\t"} ($15/$14) <= 0.20{print}' > ./Round_${CurrRound}/04_11_Temp.tblout 
                            mv ./Round_${CurrRound}/04_11_Temp.tblout ./Round_${CurrRound}/04_11_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFamily}_Peak${PeakNumLoop}.tblout
                        else
                            > ./Round_${CurrRound}/04_11_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFamily}_Peak${PeakNumLoop}.tblout
                            > ./Round_${CurrRound}/04_10_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFamily}_Peak${PeakNumLoop}.hmm
                        fi

                        #Select Peak with best hits
                        PeakNumLoop=0
                        HighestMedian=0
                        SelectedPeak=0

                        for PeakHMMOut in $(ls ./Round_${CurrRound}/04_11_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFamily}_Peak*.tblout );do
                            #Calculate median of bitscores
                            MedianScores=$( sed 's/\s\{1,\}/\t/g' ${PeakHMMOut} | \
                                            cut -f14 | \
                                            sort -n | \
                                            awk '{a[i++]=$1} END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]}' )

                            #Check if median is higher than the highest median
                            if (( $(echo "${MedianScores} > ${HighestMedian}" | bc -l) )) ;then
                                HighestMedian=${MedianScores}
                                SelectedPeak=${PeakNumLoop}
                            fi
                            PeakNumLoop=$((${PeakNumLoop}+1))
                        done
                        
                        #Select the best peak and hits
                        cp ./Round_${CurrRound}/04_09_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFamily}__Peak${SelectedPeak}.aln.fa \
                        ./Round_${CurrRound}/04_12_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFamily}.aln.fa

                        cp ./Round_${CurrRound}/04_11_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFamily}_Peak${SelectedPeak}.tblout \
                        ./Round_${CurrRound}/04_12_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFamily}.tblout

                        cp ./Round_${CurrRound}/04_09_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFamily}.Consensus__Peak${SelectedPeak}.fa \
                        ./Round_${CurrRound}/04_12_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFamily}.Consensus.fa

                    done

                    #Select those families that have homology with the previous round
                    #Make a db from previous consensus
                    makeblastdb \
                            -in ./Round_${PrevRound}/07_CurrentConsensi.${LowerCase}.Round${PrevRound}.Extended.Curated.Consensus.fa \
                            -dbtype nucl \
                            -out ./Round_${CurrRound}/04_13_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFamily}.db
                    
                    #Only do blast if there is sequence left
                    if [[ $( seqkit stat -T ./Round_${CurrRound}/04_12_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFamily}.Consensus.fa| cut -f7| tail -1| sed 's/\..*//' ) -gt 0 ]];then
                        #Perform blast
                        blastn \
                        -query ./Round_${CurrRound}/04_12_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFamily}.Consensus.fa \
                        -db ./Round_${CurrRound}/04_13_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFamily}.db \
                        -outfmt 6 \
                        -out ./Round_${CurrRound}/04_13_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFamily}.blast \
                        -num_threads ${CPU}

                        #if there are hits, then keep the family
                        if [[ $(wc -l < ./Round_${CurrRound}/04_13_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFamily}.blast) -gt 0 ]];then
                            cp ${FamilyAln} ./Round_${CurrRound}/04_14_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFamily}.aln.fa
                        fi
                    fi

                    #Increase family number
                    CurrFamily_Num=$( echo ${CurrFamily} | awk '{printf "%02d\n", $1+1}' )
                done

                #Lets sort out the new families
                if [ $(ls ./Round_${CurrRound}/04_14_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family*aln.fa | wc -l ) -gt 0 ];then
                    #Loop through families
                    for CurrFam in $(ls ./Round_${CurrRound}/04_14_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family*aln.fa | sed 's/.*Family//;s/.aln.fa//');do
                        #Do all the steps for this one before sorting
                        ##Cleanup alignment
                        ../Staging/${CurateScript} \
                            --input ./Round_${CurrRound}/04_14_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFam}.aln.fa \
                            --output ./Round_${CurrRound}/04_15_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFam} \
                            --SeqID ${RepeatID}#___Round${LetterCase}${CurrRound} \
                            --NoiseThresh ${noiseThreshold} \
                            --Plots
                        
                        ##Select the best peak
                        PeakNumLoop=0
                        if [[ $(ls ./Round_${CurrRound}/04_15_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFam}__Peak*aln.fa | wc -l ) -gt 1 ]];then
                            for PeakAln in $(ls ./Round_${CurrRound}/04_15_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFam}__Peak*aln.fa );do
                                if [[ $(seqkit stats -T ${PeakAln} | cut -f8 | tail -1) -gt 0 ]];then 
                                    #Create hmmer model
                                    hmmbuild \
                                    --symfrac 0 \
                                    --dna \
                                    --informat afa \
                                    --seed 1992 \
                                    --fragthresh 1.0 \
                                    --wnone \
                                    -n Peak${PeakNumLoop} \
                                    ./Round_${CurrRound}/04_15_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFam}_Peak${PeakNumLoop}.hmm \
                                    ${PeakAln}

                                    #Search for the model in the genome
                                    nhmmer \
                                    --cpu ${CPU} \
                                    --tblout ./Round_${CurrRound}/04_15_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFam}_Peak${PeakNumLoop}.tblout \
                                    --notextw \
                                    --noali \
                                    --seed 1992 \
                                    -E 1e-5 \
                                    ./Round_${CurrRound}/04_15_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFam}_Peak${PeakNumLoop}.hmm \
                                    ../Staging/HMMDB_Genome > /dev/null

                                    #Cleanup and remove high bias
                                    sed '/#/d; s/\s\{1,\}/\t/g' ./Round_${CurrRound}/04_15_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFam}_Peak${PeakNumLoop}.tblout | \
                                    awk '{OFS="\t"} ($15/$14) <= 0.20{print}' > ./Round_${CurrRound}/04_15_Temp.tblout 
                                    mv ./Round_${CurrRound}/04_15_Temp.tblout ./Round_${CurrRound}/04_15_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFam}_Peak${PeakNumLoop}.tblout
                                else
                                    > ./Round_${CurrRound}/04_15_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFam}_Peak${PeakNumLoop}.tblout
                                    > ./Round_${CurrRound}/04_15_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFam}_Peak${PeakNumLoop}.hmm
                                fi
                                #Increase peak number
                                PeakNumLoop=$((${PeakNumLoop}+1))
                            done

                            #Select Peak with best hits
                            PeakNumLoop=0
                            HighestMedian=0
                            SelectedPeak=0

                            for PeakHMMOut in $(ls ./Round_${CurrRound}/04_15_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFam}_Peak*.tblout );do
                                #Calculate median of bitscores
                                MedianScores=$( sed 's/\s\{1,\}/\t/g' ${PeakHMMOut} | \
                                                cut -f14 | \
                                                sort -n | \
                                                awk '{a[i++]=$1} END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]}' )

                                #Check if median is higher than the highest median
                                if (( $(echo "${MedianScores} > ${HighestMedian}" | bc -l) )) ;then
                                    HighestMedian=${MedianScores}
                                    SelectedPeak=${PeakNumLoop}
                                fi
                                PeakNumLoop=$((${PeakNumLoop}+1))
                            done
                            
                            #Select the best peak and hits
                            cp ./Round_${CurrRound}/04_15_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFam}__Peak${SelectedPeak}.aln.fa \
                            ./Round_${CurrRound}/04_16_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFam}.aln.fa

                            cp ./Round_${CurrRound}/04_15_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFam}_Peak${SelectedPeak}.tblout \
                            ./Round_${CurrRound}/04_16_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFam}.tblout

                            cp ./Round_${CurrRound}/04_15_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFam}.Consensus__Peak${SelectedPeak}.fa \
                            ./Round_${CurrRound}/04_16_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFam}.Consensus.fa
                        else
                            cp ./Round_${CurrRound}/04_15_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFam}__Peak0.aln.fa \
                            ./Round_${CurrRound}/04_16_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFam}.aln.fa

                            cp ./Round_${CurrRound}/04_15_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFam}.Consensus__Peak0.fa \
                            ./Round_${CurrRound}/04_16_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFam}.Consensus.fa
                        fi

                        ##Add alignment to previous alignment
                        mafft \
                            --localpair \
                            --maxiterate 1000 \
                            --thread ${CPU} \
                            --lexp -1.5 \
                            --lop 0.5 \
                            --add ./Round_${CurrRound}/04_16_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFam}.aln.fa \
                            ./Round_${PrevRound}/07_CurrentConsensi.${LowerCase}.Round${PrevRound}.Extended.Curated.aln.fa \
                            > ./Round_${CurrRound}/04_17_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFam}.aln

                        #Clean alignment
                        #Make zero consensus and cleaning
                        ../Staging/${CurateScript} \
                        --input ./Round_${CurrRound}/04_17_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFam}.aln \
                        --output ./Round_${CurrRound}/04_18_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFam}.Curated \
                        --SeqID ${RepeatID}#___Round${LetterCase}${CurrRound} \
                        --Plots \
                        --zeroOnly
                        
                        #Remove peak suffix
                        mv ./Round_${CurrRound}/04_18_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFam}.Curated.Consensus__Peak0.fa ./Round_${CurrRound}/04_18_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFam}.Curated.Consensus.fa
                        mv ./Round_${CurrRound}/04_18_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFam}.Curated__Peak0.aln.fa ./Round_${CurrRound}/04_18_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFam}.Curated.aln.fa
                        
                        #Check if consensus has grown too much
                        CurrConsensLength=$( seqkit fx2tab --length --name ./Round_${CurrRound}/04_18_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFam}.Curated.Consensus.fa | cut -f2 )

                        if [ ${CurrConsensLength} -gt ${maxExtensionSize} ];then
                            BreakNow=True
                            LongConsensus=True
                            continue
                        elif [ ${CurrConsensLength} -lt 1 ];then
                            BreakNow=True
                            continue
                        fi

                        ##Previous consensus length
                        PrevConsensLength=$( seqkit fx2tab --length --name ./Round_${PrevRound}/07_CurrentConsensi.${LowerCase}.Round${PrevRound}.Extended.Curated.Consensus.fa | cut -f2 )

                        ##Make mafft alignment of the previous consensus
                        mafft \
                            --localpair \
                            --maxiterate 1000 \
                            --thread ${CPU} \
                            <( cat ./Round_${PrevRound}/07_CurrentConsensi.${LowerCase}.Round${PrevRound}.Extended.Curated.Consensus.fa \
                                ./Round_${CurrRound}/04_18_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFam}.Curated.Consensus.fa ) \
                            > ./Round_${CurrRound}/04_19_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFam}.Curated.aln
                                
                        ../Staging/${CoordExtractScript} \
                            --input ./Round_${CurrRound}/04_19_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFam}.Curated.aln \
                            --output ./Round_${CurrRound}/04_19_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFam}.Curated.Coord \
                            --Side ${LowerCase}

                        #Get extended size
                        ExtendedSize=$( cat ./Round_${CurrRound}/04_19_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFam}.Curated.Coord | \
                                        awk '{print $2-$1}' )
                        
                        if [ ${ExtendedSize} -lt ${Increase} ];then
                            if [ ${ExtendedSize} -ne 0 ];then
                                > ./Round_${CurrRound}/04_20.GOODCHECK_Family${CurrFam}
                                > ./Round_${CurrRound}/04_20_BreakNow_Family${CurrFam}
                            fi
                        else
                            > ./Round_${CurrRound}/04_20.GOODCHECK_Family${CurrFam}
                            #Make model for next round
                            hmmbuild \
                            --symfrac 0 \
                            --dna \
                            --cpu ${CPU} \
                            --informat afa \
                            --seed 1992 \
                            --fragthresh 1.0 \
                            --wnone \
                            -n ${RepeatID} \
                            ./Round_${CurrRound}/04_21_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFam}.Curated.hmm \
                            ./Round_${CurrRound}/04_18_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFam}.Curated.aln.fa
                        fi

                        #Lets sort families into new dirs if they exist
                        FamNum=0
                        if [ $( ls ./Round_${CurrRound}/04_20.GOODCHECK_Family* | wc -l ) -gt 0 ];then
                                for CurrFamily in $( ls ./Round_${CurrRound}/04_20.GOODCHECK_Family* | sed 's/.*Family//' );do
                                    #First family is reference
                                    if [ ${FamNum} -eq 0 ];then
                                        #Fragment alignment 05
                                        cp ./Round_${CurrRound}/04_12_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFamily}.aln.fa \
                                        ./Round_${CurrRound}/05_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.Curated.aln.fa
                                        cp ./Round_${CurrRound}/04_12_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFamily}.Consensus.fa \
                                        ./Round_${CurrRound}/05_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.Curated.Consensus.fa
                                        
                                        #Add Alignment 06
                                        cp ./Round_${CurrRound}/04_16_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFamily}.aln.fa \
                                        ./Round_${CurrRound}/06_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.Curated.aln.fa

                                        #Clean Alignment 07
                                        cp ./Round_${CurrRound}/04_18_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFamily}.Curated.aln.fa \
                                        ./Round_${CurrRound}/07_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.Curated.aln.fa
                                        cp ./Round_${CurrRound}/04_18_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFamily}.Curated.Consensus.fa \
                                        ./Round_${CurrRound}/07_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.Curated.Consensus.fa

                                        #Coords and model 08
                                        cp ./Round_${CurrRound}/04_19_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFamily}.Curated.Coord \
                                        ./Round_${CurrRound}/08_CurrentConsensi.${LowerCase}.Round${CurrRound}.ExtendedSide.Coord
                                        cp ./Round_${CurrRound}/04_21_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFamily}.Curated.hmm \
                                        ./Round_${CurrRound}/08_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.Curated.hmm

                                        #CHECK 
                                        cp ./Round_${CurrRound}/04_20.GOODCHECK_Family${CurrFam} \
                                        ./Round_${CurrRound}/GoodCHECK
                                    else
                                        #Make new family directory
                                        ##Last family number
                                        HighestCurFam=$(ls ../Family_* | sed 's/.*Family_//' | sort -n | tail -1)

                                        ##Increase by one the family number
                                        HighestCurFam=$( echo ${HighestCurFam} | awk '{printf "%02d\n", $1+1}' )

                                        #Make dir
                                        mkdir ../Family_${HighestCurFam}

                                        #Copy all previous rounds to new family
                                        cp -r ./Round_* ../Family_${HighestCurFam}/

                                        #Copy that family into round
                                        ##Fragment alignment 05
                                        cp ./Round_${CurrRound}/04_12_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFamily}.aln.fa \
                                        ../Family_${HighestCurFam}/Round_${CurrRound}/05_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.Curated.aln.fa
                                        cp ./Round_${CurrRound}/04_12_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFamily}.Consensus.fa \
                                        ../Family_${HighestCurFam}/Round_${CurrRound}/05_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.Curated.Consensus.fa
                                        
                                        ##Add Alignment 06
                                        cp ./Round_${CurrRound}/04_16_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFamily}.aln.fa \
                                        ../Family_${HighestCurFam}/Round_${CurrRound}/06_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.Curated.aln.fa

                                        ##Clean Alignment 07
                                        cp ./Round_${CurrRound}/04_18_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFamily}.Curated.aln.fa \
                                        ../Family_${HighestCurFam}/Round_${CurrRound}/07_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.Curated.aln.fa
                                        cp ./Round_${CurrRound}/04_18_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFamily}.Curated.Consensus.fa \
                                        ../Family_${HighestCurFam}/Round_${CurrRound}/07_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.Curated.Consensus.fa

                                        ##Coords and model 08
                                        cp ./Round_${CurrRound}/04_19_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFamily}.Curated.Coord \
                                        ../Family_${HighestCurFam}/Round_${CurrRound}/08_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.Curated.Coords
                                        cp ./Round_${CurrRound}/04_21_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFamily}.Curated.hmm \
                                        ../Family_${HighestCurFam}/Round_${CurrRound}/08_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.Curated.hmm

                                        #CHECK 
                                        cp ./Round_${CurrRound}/04_20.GOODCHECK_Family${CurrFamily} \
                                        ../Family_${HighestCurFam}/Round_${CurrRound}/GoodCHECK
                                        cp ./Round_${CurrRound}/04_20_BreakNow_Family${CurrFamily} \
                                        ../Family_${HighestCurFam}/Round_${CurrRound}/BreakNow
                                    fi
                                    FamNum=$((${FamNum}+1))
                                done
                        else 
                                #Round over
                                BreakNow=True
                        fi
                    done
                else
                    break
                fi  
            else
                #Business as usual
                #Make all sequences uppercase
                seqkit seq \
                        -u \
                        ./Round_${CurrRound}/04_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.Fams.Cluster_0.fa.aln \
                        > ./Round_${CurrRound}/04_TEMP.fa.aln

                mv ./Round_${CurrRound}/04_TEMP.fa.aln \
                        ./Round_${CurrRound}/04_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.Fams.Cluster_0.fa.aln

                #Clean the alignment
                ../Staging/${CurateScript} \
                        --input ./Round_${CurrRound}/04_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.Fams.Cluster_0.fa.aln \
                        --output ./Round_${CurrRound}/04_00_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.Curated \
                        --SeqID ${RepeatID}#___Round${LetterCase}${CurrRound} \
                        --NoiseThresh ${noiseThreshold} \
                        --Plots

                #Select Peak with best hits
                ##Dont look for the best peak if there is only one
                if [ $(ls ./Round_${CurrRound}/04_00_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.Curated__Peak*aln.fa | wc -l ) -eq 1 ];then
                    cp ./Round_${CurrRound}/04_00_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.Curated__Peak0.aln.fa \
                    ./Round_${CurrRound}/05_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.Curated.aln.fa
                    cp ./Round_${CurrRound}/04_00_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.Curated.Consensus__Peak0.fa \
                    ./Round_${CurrRound}/05_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.Curated.Consensus.fa
                else
                    PeakNumLoop=0
                    for PeakAln in $(ls ./Round_${CurrRound}/04_00_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.Curated__Peak*aln.fa);do
                        if [[ $(seqkit stats -T ${PeakAln} | cut -f8 | tail -1) -gt 0 ]];then 
                                #Create hmmer model
                                hmmbuild \
                                --symfrac 0 \
                                --dna \
                                --informat afa \
                                --seed 1992 \
                                --fragthresh 1.0 \
                                --wnone \
                                -n Peak${PeakNumLoop} \
                                ./Round_${CurrRound}/04_01_CurrentConsensi__Peak${PeakNumLoop}.hmm \
                                ${PeakAln}

                                #Search for the model in the genome
                                nhmmer \
                                --cpu ${CPU} \
                                --tblout ./Round_${CurrRound}/04_02_CurrentConsensi__Peak${PeakNumLoop}.tblout \
                                --notextw \
                                --noali \
                                --seed 1992 \
                                -E 1e-5 \
                                ./Round_${CurrRound}/04_01_CurrentConsensi__Peak${PeakNumLoop}.hmm \
                                ../Staging/HMMDB_Genome > /dev/null

                                #Cleanup and remove high bias
                                sed '/#/d; s/\s\{1,\}/\t/g' ./Round_${CurrRound}/04_02_CurrentConsensi__Peak${PeakNumLoop}.tblout | \
                                awk '{OFS="\t"} ($15/$14) <= 0.20{print}' > ./Round_${CurrRound}/04_02_Temp.tblout
                        else
                                > ./Round_${CurrRound}/04_02_CurrentConsensi__Peak${PeakNumLoop}.tblout
                                > ./Round_${CurrRound}/04_01_CurrentConsensi__Peak${PeakNumLoop}.hmm
                        fi
                        #Increase Peak Number
                        PeakNumLoop=$((${PeakNumLoop}+1))
                    done
                fi

                #Select Peak with best hits
                PeakNumLoop=0
                HighestMedian=0
                SelectedPeak=0
                for PeakHMMOut in $(ls ./Round_${CurrRound}/04_02_CurrentConsensi__Peak*.tblout );do
                    #Calculate median of bitscores
                    MedianScores=$( sed 's/\s\{1,\}/\t/g' ${PeakHMMOut} | \
                                    cut -f14 | \
                                    sort -n | \
                                    awk '{a[i++]=$1} END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]}' )

                    #Check if median is higher than the highest median
                    if (( $(echo "${MedianScores} > ${HighestMedian}" | bc -l) )) ;then
                        HighestMedian=${MedianScores}
                        SelectedPeak=${PeakNumLoop}
                    fi
                    PeakNumLoop=$((${PeakNumLoop}+1))
                done

                #Select the best peak and hits
                cp ./Round_${CurrRound}/04_00_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.Curated__Peak${SelectedPeak}.aln.fa \
                ./Round_${CurrRound}/05_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.Curated.aln.fa
                cp ./Round_${CurrRound}/04_00_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.Curated.Consensus__Peak${SelectedPeak}.fa \
                ./Round_${CurrRound}/05_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.Curated.Consensus.fa

                #Merge alignments
                mafft \
                --localpair \
                --maxiterate 1000 \
                --lexp -1.5 \
                --lop 0.5 \
                --thread ${CPU} \
                --add ./Round_${CurrRound}/05_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.Curated.aln.fa \
                ./Round_${PrevRound}/07_CurrentConsensi.${LowerCase}.Round${PrevRound}.Extended.Curated.aln.fa \
                > ./Round_${CurrRound}/06_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.aln \
                2> ./Round_${CurrRound}/06_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.MAFFTlog

                #Make zero consensus and cleaning
                ../Staging/${CurateScript} \
                --input ./Round_${CurrRound}/06_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.aln \
                --output ./Round_${CurrRound}/07_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.Curated \
                --SeqID ${RepeatID}#___Round${LetterCase}${CurrRound} \
                --Plots \
                --zeroOnly
                
                #Remove peak suffix
                mv ./Round_${CurrRound}/07_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.Curated.Consensus__Peak0.fa ./Round_${CurrRound}/07_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.Curated.Consensus.fa
                mv ./Round_${CurrRound}/07_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.Curated__Peak0.aln.fa ./Round_${CurrRound}/07_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.Curated.aln.fa
                
                #Check if consensus has grown too much
                CurrConsensLength=$( seqkit fx2tab --length --name ./Round_${CurrRound}/07_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.Curated.Consensus.fa | cut -f2 )

                if [ ${CurrConsensLength} -gt ${maxExtensionSize} ];then
                    BreakNow=True
                    LongConsensus=True
                    break
                elif [ ${CurrConsensLength} -lt 1 ];then
                    BreakNow=True
                    break
                fi

                ##Previous consensus length
                PrevConsensLength=$( seqkit fx2tab --length --name ./Round_${PrevRound}/07_CurrentConsensi.${LowerCase}.Round${PrevRound}.Extended.Curated.Consensus.fa | cut -f2 )

                ##Make mafft alignment of the previous consensus
                mafft \
                    --localpair \
                    --maxiterate 1000 \
                    --thread ${CPU} \
                    <( cat ./Round_${PrevRound}/07_CurrentConsensi.${LowerCase}.Round${PrevRound}.Extended.Curated.Consensus.fa \
                        ./Round_${CurrRound}/07_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.Curated.Consensus.fa ) \
                    > ./Round_${CurrRound}/08_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.aln

                ../Staging/${CoordExtractScript} \
                    --input ./Round_${CurrRound}/08_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.aln \
                    --output ./Round_${CurrRound}/08_CurrentConsensi.${LowerCase}.Round${CurrRound}.ExtendedSide.Coord \
                    --Side ${LowerCase}
                        
                #Get extended size
                ExtendedSize=$( cat ./Round_${CurrRound}/08_CurrentConsensi.${LowerCase}.Round${CurrRound}.ExtendedSide.Coord | \
                                awk '{print $2-$1}' )

                if [ ${ExtendedSize} -lt ${Increase} ];then
                    if [ ${ExtendedSize} -ne 0 ];then
                        > ./Round_${CurrRound}/GoodCHECK
                    fi
                    BreakNow=True
                    break
                else
                    > ./Round_${CurrRound}/GoodCHECK
                    #Make model for next round
                    hmmbuild \
                    --symfrac 0 \
                    --dna \
                    --cpu ${CPU} \
                    --informat afa \
                    --seed 1992 \
                    --fragthresh 1.0 \
                    --wnone \
                    -n ${RepeatID} \
                    ./Round_${CurrRound}/08_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.Curated.hmm \
                    ./Round_${CurrRound}/07_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.Curated.aln.fa
                fi

                #Check if it has reached max rounds
                if [ ${CurrRound} -ge ${maxRounds} ];then
                    BreakNow=True
                    break
                fi
            fi

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
    #Determine the number of families on each side
    NumFamiliesLeft=$(ls -d ../Left/Family_* | wc -l)
    NumFamiliesRight=$(ls -d ../Right/Family_* | wc -l)
    
    if [ ${NumFamiliesLeft} -gt ${NumFamiliesRight} ];then
        NumFamilies=${NumFamiliesLeft}
    elif [ ${NumFamiliesRight} -gt 0 ];then
        NumFamilies=${NumFamiliesRight}
    else
        NumFamilies=1
    fi

    for CurrFam in $( seq -w 00 1 $(( ${NumFamilies} - 1 )) );do
        #Make Family directory
        mkdir -p ./Family_${CurrFam}

        ExtendSucc=True
        #Check which side to use for final output
        #Bothsides or RIGHT OK
        if [ -f ../Right/Family_${CurrFam}/GOOD.FAM.CHECK ];then
            SideMerge=Right
        elif [ -f ../Left/Family_${CurrFam}/GOOD.FAM.CHECK ];then
            SideMerge=Left
        elif [ -f ../Right/Family_${CurrFam}/BAD.FAM.CHECK ] && [ -f ../Left/Family_${CurrFam}/BAD.FAM.CHECK ] && [ -f ../Central/Family_${CurrFam}/Central.GOOD.CHECK ];then
            SideMerge=Central
        elif ( [ -f ../Right/Family_${CurrFam}/BAD.SIDE.CHECK ] || [ -f ../Right/Family_${CurrFam}/BAD.FAM.CHECK ] ) && ( [ -f ../Left/Family_${CurrFam}/BAD.SIDE.CHECK ] || [ -f ../Left/Family_${CurrFam}/BAD.FAM.CHECK ] ) && [ -f ../Central/Family_${CurrFam}/Central.GOOD.CHECK ];then
            SideMerge=Central
        else
            ExtendSucc=False
        fi

        if [ ${ExtendSucc} == "True" ] && [ ${SideMerge} != "Central" ] ;then
            HighestRound=$( ls ../${SideMerge}/Family_${CurrFam}/Round_*/GoodCHECK | sed "s/\/GoodCHECK//;s/${SideMerge}\///;s/.*\///"| tail -1 )
            cp ../${SideMerge}/Family_${CurrFam}/${HighestRound}/07_CurrentConsensi*Extended.Curated.aln.fa ./Family_${CurrFam}/FinalConsensi.aln.fa
            cp ../${SideMerge}/Family_${CurrFam}/${HighestRound}/07_CurrentConsensi*Extended.Curated.Consensus.fa ./Family_${CurrFam}/FinalConsensi.consensus.fa
            esl-reformat stockholm ./Family_${CurrFam}/FinalConsensi.aln.fa  >  ./Family_${CurrFam}/FinalConsensi.aln.stk
        elif [ ${ExtendSucc} == "True" ] && [ ${SideMerge} == "Central" ];then
            HighestRound=$( ls ../${SideMerge}/Family_${CurrFam}/Round_*/GoodCHECK | sed "s/\/GoodCHECK//;s/${SideMerge}\///;s/.*\///"| tail -1 )
            cp ../${SideMerge}/Family_${CurrFam}/${HighestRound}/07_CurrentConsensi*Extended.Curated.aln.fa \
                ./Family_${CurrFam}/FinalConsensi.aln.fa
            cp ../${SideMerge}/Family_${CurrFam}/${HighestRound}/07_CurrentConsensi*Extended.Curated.Consensus.fa \
                ./Family_${CurrFam}/FinalConsensi.consensus.fa
            esl-reformat stockholm ./Family_${CurrFam}/FinalConsensi.aln.fa  >  ./Family_${CurrFam}/FinalConsensi.aln.stk
        else
            echo ">${RepeatID}#BAD" > ./Family_${CurrFam}/FinalConsensi.aln.fa
            echo "A" >> ./Family_${CurrFam}/FinalConsensi.aln.fa
            cp ./Family_${CurrFam}/FinalConsensi.aln.fa ./Family_${CurrFam}/FinalConsensi.consensus.fa
            > ./Family_${CurrFam}/FinalConsensi.aln.stk
        fi

        #Make Final Check
        > ./Family_${CurrFam}/GoodCHECK
    done
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