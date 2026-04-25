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
#AnchorSize=$(echo ${ExtensionSize} | awk '{printf "%.0f\n", $1*0.5}') 
AnchorSize=$(  echo ${ExtensionSize} | awk '{print int($1*1.5)}' )
Increase=$(  echo ${ExtensionSize} | awk '{$0=$0*0.125;printf("%d\n",$0+=$0<0?0:0.9)}' )
LongConsensus=False
SecondGo=False

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
                        cp ${LeftDir}/Family_${Family}/${HighestRound}/07_CurrentConsensi.Left*Extended.Curated.aln.fa \
                            ./Round_00/07_CurrentConsensi.${LowerCase}.Round00.Extended.Curated.aln.fa
                        cp ${LeftDir}/Family_${Family}/${HighestRound}/07_CurrentConsensi.Left*Extended.Curated.Consensus.fa \
                            ./Round_00/07_CurrentConsensi.${LowerCase}.Round00.Extended.Curated.Consensus.fa
                        cp ${LeftDir}/Family_${Family}/${HighestRound}/02_CurrentConsensi.Left*Extended.bed \
                            ./Round_00/02_CurrentConsensi.${LowerCase}.Round00.Extended.bed
                        cp ${LeftDir}/Family_${Family}/${HighestRound}/07_CurrentConsensi.Left*Extended.Curated.Consensus.fa \
                            ./Round_00/05_CurrentConsensi.${LowerCase}.Round00.Extended.Consensus.fa

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
                    cp ${CentralDir}/Family_${Family}/Final/FinalConsensi.Extension.Coord.bed \
                        ./Round_00/02_CurrentConsensi.${LowerCase}.Round00.Extended.bed
                fi

            fi
        
            #Cutoff for symfrac based on number of sequences
            #Symfrac_Thresh=$(grep -c ">" ./Round_00/07_CurrentConsensi.${LowerCase}.Round00.Extended.Curated.aln.fa | awk '{print 50/$1}')
            Symfrac_Thresh=0.3 

            ##Make HMM model 
            hmmbuild \
            --symfrac ${Symfrac_Thresh} \
            --dna \
            --cpu ${CPU} \
            --informat afa \
            --seed 1992 \
            --fragthresh 1.0 \
            -n ${RepeatID} \
            --wnone \
            ./Round_00/08_CurrentConsensi.${LowerCase}.Round00.Extended.Curated.hmm \
            ./Round_00/07_CurrentConsensi.${LowerCase}.Round00.Extended.Curated.aln.fa

            cp ./Round_00/08_CurrentConsensi.${LowerCase}.Round00.Extended.Curated.hmm \
                ./Round_00/05_CurrentConsensi.${LowerCase}.Round00.Extended.hmm
            
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
            ./Round_${PrevRound}/05_CurrentConsensi.${LowerCase}.Round${PrevRound}.Extended.hmm \
            ../Staging/HMMDB_Genome  > /dev/null

            ErrorCode=$?
            if [ ${ErrorCode} -ne 0 ];then
                echo -e "\nHMMER search failed, exiting"
                exit ${ErrorCode}
            fi
            
            #Format HMMER output
            sed '/#/d; s/\s\{1,\}/\t/g' ./Round_${CurrRound}/00_CurrentConsensi.${LowerCase}.Round${CurrRound}.hmmout | \
            awk '$14 > 30 {print}' \
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
                awk -v Increase=${ExtensionSize} -v OFMT='%f' '{OFS="\t"}{
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
                awk -v Increase=${ExtensionSize} -v OFMT='%f' '{OFS="\t"}{
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
                    --firstPeak \
                    --Plots

            #Rename output 
            mv ./Round_${CurrRound}/05_CurrentConsensi.${LowerCase}.Round${CurrRound}.Family${CurrFamily}__Peak0.aln.fa \
                ./Round_${CurrRound}/05_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.aln.fa

            #Make hmmer model out of this alignment
            hmmbuild \
                --fragthresh 1.0 \
                --dna \
                --cpu ${CPU} \
                --seed 1992 \
                -n ${RepeatID}__${CurrRound} \
                --wnone \
                --symfrac ${percentVertical} \
            ./Round_${CurrRound}/05_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.hmm \
            ./Round_${CurrRound}/05_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.aln.fa

            #Get sequence
            hmmemit \
                -c \
                ./Round_${CurrRound}/05_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.hmm \
            > ./Round_${CurrRound}/05_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.Consensus.fa

            #Make all sequences uppercase
            seqkit seq \
                    -u \
                    ./Round_${CurrRound}/05_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.aln.fa \
                    > ./Round_${CurrRound}/05_TEMP.fa.aln
            
            cp ./Round_${CurrRound}/05_TEMP.fa.aln \
                ./Round_${CurrRound}/05_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.aln.fa
            
            #Reconcile previous alignments
            ##Make alignment with previous round
            mafft \
                --localpair \
                --maxiterate 1000 \
                --thread ${CPU} \
                --lexp -1.5 \
                --lop 0.5 \
                --add ./Round_${PrevRound}/05_CurrentConsensi.${LowerCase}.Round${PrevRound}.Extended.aln.fa \
                ./Round_${CurrRound}/05_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.aln.fa | \
            sed "s/${LowerCase}R${PrevRound}/&DUP/" | \
            sed "s/Central/&DUP/" \
            > ./Round_${CurrRound}/06_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.tmp1.aln.fa

            #Clean alignment with vertical script
            python3 ../Staging/${VerticalScript} \
                --input ./Round_${CurrRound}/06_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.tmp1.aln.fa \
                --output ./Round_${CurrRound}/06_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.tmp2 \
                --perc  ${percentVertical}

            #Merge alignment with previous rounds
            mafft \
                --seed ./Round_${CurrRound}/06_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.tmp2.aln.fa \
                --seed ./Round_${PrevRound}/07_CurrentConsensi.${LowerCase}.Round${PrevRound}.Extended.Curated.aln.fa \
				/dev/null | \
            seqkit grep -v -r -p "DUP" \
            > ./Round_${CurrRound}/06_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.aln.fa

            #Remove gappy columns
            python3 ../Staging/${CurateScript} \
                --input ./Round_${CurrRound}/06_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.aln.fa \
                --output ./Round_${CurrRound}/07_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.Curated \
                --SeqID ${RepeatID}#___Round${LetterCase}${CurrRound} \
                --zeroOnly \
                --Plots

            #Remove peak suffix
            mv ./Round_${CurrRound}/07_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.Curated.Consensus__Peak0.fa ./Round_${CurrRound}/07_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.Curated.Consensus.fa
            mv ./Round_${CurrRound}/07_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.Curated__Peak0.aln.fa ./Round_${CurrRound}/07_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.Curated.aln.fa                

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

            #Get extended size
            ExtendedSize=$( cat ./Round_${CurrRound}/08_CurrentConsensi.${LowerCase}.Round${CurrRound}.ExtendedSide.Coord | \
                            awk '{print $2-$1}' )

            if [[ ${ExtendedSize} -lt ${Increase} ]];then
                if [[ ${SecondGo} == "False" ]];then
                    #Remake Coord file to git a second go
                    if [[ $1 == "Left" ]];then
                        echo -e "0\t${ExtensionSize}" > ./Round_${CurrRound}/08_CurrentConsensi.${LowerCase}.Round${CurrRound}.ExtendedSide.Coord
                    else
                        LengthConsensus=$(seqkit fx2tab --length --name ./Round_${CurrRound}/07_CurrentConsensi.${LowerCase}.Round${CurrRound}.Extended.Curated.Consensus.fa | cut -f2)
                        echo | awk -v LenCons=${LengthConsensus} -v Ext=${ExtensionSize} '{print LenCons-Ext"\t"LenCons}' > ./Round_${CurrRound}/08_CurrentConsensi.${LowerCase}.Round${CurrRound}.ExtendedSide.Coord
                    fi
                    SecondGo=True
                else 
                    if [ ${ExtendedSize} -ne 0 ];then
                        > ./Round_${CurrRound}/GoodCHECK
                    fi
                    BreakNow=True
                    break
                fi
            else
                SecondGo=False
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
        mkdir -p Left
        cp ../Left/Family_00/Round_00/05_CurrentConsensi.Left.Round00.Extended.aln.fa \
            ./Left/Round_00.aln.fa
        for CurrRound in $(ls -d ../Left/Family_00/Round_*| sed 1d| sed 's/.*_//');do
            PrevRound=$(echo ${CurrRound} | awk '{printf "%02d\n", $1-1}')
            mafft \
                --seed ../Left/Family_00/Round_${CurrRound}/06_CurrentConsensi.Left.Round${CurrRound}.Extended.tmp2.aln.fa \
                --seed ./Left/Round_${PrevRound}.aln.fa \
                /dev/null | \
            seqkit grep -v -r -p "DUP" > ./Left/TEMP_Round${CurrRound}.aln.fa

            python ../Left/Staging/${VerticalScript} \
                --input ./Left/TEMP_Round${CurrRound}.aln.fa \
                --output ./Left/Round_${CurrRound} \
                --perc 0.3
        done

        if [ -f ../Right/Family_00/GOOD.FAM.CHECK ];then
            mkdir -p Right

            cp ./Left/MergedAln_Round${CurrRound}.aln.fa \
                ./Right/Round00.aln.fa

            for CurrRound in $(ls -d ../Right/Family_00/Round_*| sed 1d| sed 's/.*_//');do
                PrevRound=$(echo ${CurrRound} | awk '{printf "%02d\n", $1-1}')
                mafft \
                    --seed ../Right/Family_00/Round_${CurrRound}/06_CurrentConsensi.Right.Round${CurrRound}.Extended.tmp2.aln.fa \
                    --seed ./Right/Round_${PrevRound}.aln.fa \
                    /dev/null | \
                seqkit grep -v -r -p "DUP" > ./Right/TEMP_Round${RoundNum}.aln.fa

                python ../Right/Staging/${VerticalScript} \
                    --input ./Right/TEMP_Round${RoundNum}.aln.fa \
                    --output Round_${CurrRound} \
                    --perc 0.3
            done
        fi
    elif [ ../Right/Family_00/GOOD.FAM.CHECK ];then
        mkdir -p Right

        cp ../Right/Family_00/Round_00/05_CurrentConsensi.Right.Round00.Extended.aln.fa \
            ./Right/Round_00.aln.fa

        for CurrRound in $(ls -d ../Right/Family_00/Round_*| sed 1d| sed 's/.*_//');do
            PrevRound=$(echo ${CurrRound} | awk '{printf "%02d\n", $1-1}')
            mafft \
                --seed ../Right/Family_00/Round_${CurrRound}/06_CurrentConsensi.Right.Round${CurrRound}.Extended.tmp2.aln.fa \
                --seed ./Right/Round_${PrevRound}.aln.fa \
                /dev/null | \
            seqkit grep -v -r -p "DUP" > ./Right/TEMP_Round${CurrRound}.aln.fa

            python ../Right/Staging/${VerticalScript} \
                --input ./Right/TEMP_Round${CurrRound}.aln.fa \
                --output ./Right/Round_${CurrRound} \
                --perc 0.3
        done
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