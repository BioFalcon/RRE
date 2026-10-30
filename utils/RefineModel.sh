#!/bin/bash
#Refine the final model of a repeat genome-wide, Dfam style: search the genome with the model, align the hits
#back to the model with hmmalign, define match columns (#=GC RF) from the occupancy of the hits, rebuild the
#model with hmmbuild --hand, and iterate until the model is stable: few new hits and an almost unchanged consensus.
#The exact consensus never settles for an old repeat (near-tied columns flip as the sample of copies shifts).
#Iterating can also slowly recruit a related family, so refinement stops if the consensus drifts too far from
#that of the first iteration, keeping the last model within the limit.
#Models are built as Dfam builds them: --maxinsertlen 10 (without it, hits bridge long insertions, often younger
#TEs nested in the copies, and those count as copies of the repeat) and exponential entropy weighting (--eentexp,
#to the default target of 0.62 bits/position), which down-weights densely covered columns more than sparse ones,
#so regions with few copies (e.g. eroded ends) are not flattened to the prior.
#Model ends supported by only a few copies (e.g. flank picked up during the extension) are trimmed by coverage.
#Blocks of a younger, more abundant family the extension ran into (many hits, much more similar to each other
#than the copies of the rest of the model) are cut off at the depth valley between them and the rest (ForeignBlocks.py).
#
#Usage: RefineModel.sh -a <aln.fa> -d <HMMER db of the genome> -g <genome.fa> -o <OutPrefix> [Options]
#   -a  Starting alignment in fasta format (the merged model). If <aln>.Columns.txt (StitchMerge.py) is
#       next to it, its column calls are the match columns of the first model, otherwise HandRF.py decides
#   -d  Genome as searched by nhmmer (fasta or makehmmerdb database)
#   -g  Genome fasta, to extract the hits (indexed, .fai)
#   -o  Output prefix
#   -n  Model name (Default: basename of the output prefix)
#   -t  CPUs (Default: 1)
#   -I  Max. number of iterations (Default: 8)
#   -E  Max. E-value of a hit (Default: 1e-5)
#   -X  Max. number of hits used to build the model, best scoring first (Default: 2000, as Dfam seeds)
#   -L  Min. length of a hit (Default: 50)
#   -Q  Min. occupancy of a column, among the hits spanning it, for it to be a match column (Default: 0.5)
#   -T  Trim model ends with fewer hits per column than this fraction of the median (Default: 0.1)
#   -D  Collapse hits more similar than this into one representative, or auto (Default: auto)
#   -G  Converged when hits grow by less than this fraction, together with -C (Default: 0.02)
#   -C  Converged when the consensus changes at less than this fraction of positions, together with -G (Default: 0.03;
#       near-tied columns of an old repeat keep flipping between iterations, ~1-2% of positions)
#   -R  Stop when the consensus differs from that of iteration 1 at more than this fraction of positions (Default: 0.10)
#
#Writes:
#   <OutPrefix>.hmm                 final model
#   <OutPrefix>.Consensus.fa        consensus of the final model (one base per match state)
#   <OutPrefix>.stk                 final alignment in Stockholm, with #=GC RF, match and insert columns
#   <OutPrefix>.aln.fa              the same alignment in fasta, with all the indels
#   <OutPrefix>.RFclean.aln.fa      the same alignment with the match (consensus) columns only
#   <OutPrefix>.support.tsv         residues and conservation per column of the final alignment
#   <OutPrefix>.iterations.tsv      hits, copies, match columns and consensus change per iteration

CPU=1
Iterations=8
MaxE=1e-5
MaxSeqs=2000
MinLen=50
MinOcc=0.5
TrimCov=0.1
MaxIdentity=auto
MinGain=0.02
MaxChange=0.03
MaxDrift=0.10
Name=""

while getopts "a:d:g:o:n:t:I:E:X:L:Q:T:D:G:C:R:" opt;do
    case "$opt" in
        a ) Aln="$OPTARG" ;;
        d ) Database="$OPTARG" ;;
        g ) Genome="$OPTARG" ;;
        o ) Out="$OPTARG" ;;
        n ) Name="$OPTARG" ;;
        t ) CPU="$OPTARG" ;;
        I ) Iterations="$OPTARG" ;;
        E ) MaxE="$OPTARG" ;;
        X ) MaxSeqs="$OPTARG" ;;
        L ) MinLen="$OPTARG" ;;
        Q ) MinOcc="$OPTARG" ;;
        T ) TrimCov="$OPTARG" ;;
        D ) MaxIdentity="$OPTARG" ;;
        G ) MinGain="$OPTARG" ;;
        C ) MaxChange="$OPTARG" ;;
        R ) MaxDrift="$OPTARG" ;;
    esac
done

if [ -z "${Aln}" ] || [ -z "${Database}" ] || [ -z "${Genome}" ] || [ -z "${Out}" ];then
    sed -n '2,/^$/p' $0 | sed 's/^#//'
    exit 1
fi
UtilsDir=$(dirname $(realpath $0))
[ -z "${Name}" ] && Name=$(basename ${Out})
Work=${Out}.Refine
mkdir -p ${Work}

#Stockholm of a fasta alignment with the given column calls ('-' = not a match column) as RF
function FastaToStockholm() {
    python3 - "$1" "$2" "$3" "$4" <<'EOF'
import sys
Aln, Calls, Out, Name = sys.argv[1:5]
Names, Seqs = [], []
for line in open(Aln):
    line = line.strip()
    if line.startswith(">"):
        Names.append(line[1:].split()[0]); Seqs.append([])
    elif line:
        Seqs[-1].append(line.upper().replace(".", "-"))
Seqs = ["".join(s) for s in Seqs]
RF = "".join("." if c == "-" else "x" for c in open(Calls).read().strip())
W = max(len(n) for n in Names + ["#=GC RF"]) + 1
with open(Out, "w") as fh:
    fh.write("# STOCKHOLM 1.0\n#=GF ID " + Name + "\n\n")
    for n, s in zip(Names, Seqs):
        if s.strip("-"):
            fh.write(n.ljust(W) + s + "\n")
    fh.write("#=GC RF".ljust(W) + RF + "\n//\n")
EOF
}

#Fraction of positions at which two consensi differ (mismatches and indels of their global alignment)
function ConsensusChange() {
    (sed "1s/.*/>a/" $1; sed "1s/.*/>b/" $2) > ${Work}/Pair.fa
    mafft --quiet --globalpair --maxiterate 100 ${Work}/Pair.fa 2> /dev/null | \
        awk '/^>/{n++; next}{s[n]=s[n] toupper($0)}
             END{for(i=1;i<=length(s[1]);i++){a=substr(s[1],i,1); b=substr(s[2],i,1); if(a!=b) d++; if(a!="-") la++}
                 printf "%.4f\n", d/la}'
}

#Model of iteration 0, from the starting alignment
ColumnsFile=${Aln%.aln.fa}.Columns.txt
if [ -f ${ColumnsFile} ] && [ $(tr -d '\n' < ${ColumnsFile} | wc -c) -eq $(awk 'NR==2{print length($0)}' ${Aln}) ];then
    FastaToStockholm ${Aln} ${ColumnsFile} ${Work}/Iter00.stk ${Name}
    echo "Iteration 00: match columns from ${ColumnsFile}"
else
    python3 ${UtilsDir}/HandRF.py --input ${Aln} --output ${Work}/Iter00 --minOcc ${MinOcc} --name ${Name} > ${Work}/Iter00.HandRF.log || exit 1
    echo "Iteration 00: match columns from HandRF.py"
fi
hmmbuild --hand --dna --cpu ${CPU} --seed 1992 --fragthresh 1.0 --maxinsertlen 10 --eentexp -n ${Name} \
    ${Work}/Iter00.hmm ${Work}/Iter00.stk > /dev/null || { echo "hmmbuild failed on the starting alignment"; exit 1; }
hmmemit -c ${Work}/Iter00.hmm > ${Work}/Iter00.Consensus.fa
echo -e "iteration\thits\tused\tcopies\tmatch_columns\tcoverage_trimmed\tconsensus_length\thit_gain\tconsensus_change\tdrift_from_iter1" > ${Out}.iterations.tsv
echo -e "0\tNA\tNA\tNA\t$(grep "^LENG" ${Work}/Iter00.hmm | awk '{print $2}')\tNA\t$(grep -v ">" ${Work}/Iter00.Consensus.fa | tr -d '\n' | wc -c)\tNA\tNA\tNA" >> ${Out}.iterations.tsv
PrevHits=0

Prev=00
for (( k=1; k<=Iterations; k++ ));do
    It=$(printf "%02d" ${k})

    #Genome-wide search with the current model
    nhmmer --cpu ${CPU} --seed 1992 --noali --notextw -E ${MaxE} \
        --tblout ${Work}/Iter${It}.tbl ${Work}/Iter${Prev}.hmm ${Database} > /dev/null || { echo "nhmmer failed"; exit 1; }

    #Blocks of a younger abundant family: only the model between them is kept
    Len=$(grep "^LENG" ${Work}/Iter${Prev}.hmm | awk '{print $2}')
    python3 ${UtilsDir}/ForeignBlocks.py -t ${Work}/Iter${It}.tbl -L ${Len} -g ${Genome} -E ${MaxE} -c ${CPU} \
        > ${Work}/Iter${It}.foreign.log || echo "KEEP 1 ${Len}" > ${Work}/Iter${It}.foreign.log
    read KeepS KeepE < <(awk '$1=="KEEP"{print $2, $3}' ${Work}/Iter${It}.foreign.log)
    if [ -z "${KeepS}" ];then KeepS=1; KeepE=${Len}; fi
    if [ ${KeepS} -gt 1 ] || [ ${KeepE} -lt ${Len} ];then
        echo "Iteration ${It}: blocks of a younger abundant family found, keeping model positions ${KeepS}-${KeepE} of ${Len}"
    fi

    #Hits as bed (envelope coordinates), best first, at least MinLen long, overlapping the kept model positions
    grep -v "^#" ${Work}/Iter${It}.tbl | \
        awk -v MinLen=${MinLen} -v MaxE=${MaxE} -v S=${KeepS} -v E=${KeepE} '{OFS="\t"}
            $13 <= MaxE && $5 <= E && $6 >= S {
                if($12 == "+"){s=$9-1; e=$10}else{s=$10-1; e=$9}
                if(e - s >= MinLen) print $1, s, e, $1":"s"-"e, $14, $12
            }' | sort -k5,5gr > ${Work}/Iter${It}.hits.bed
    Hits=$(wc -l < ${Work}/Iter${It}.hits.bed)
    if [ ${Hits} -lt 5 ];then
        echo "Iteration ${It}: only ${Hits} hits, keeping the model of iteration ${Prev}"
        break
    fi
    head -n ${MaxSeqs} ${Work}/Iter${It}.hits.bed > ${Work}/Iter${It}.used.bed
    bedtools getfasta -s -nameOnly -fi ${Genome} -bed ${Work}/Iter${It}.used.bed -fo ${Work}/Iter${It}.copies.fa
    sed -i 's/[()]//g' ${Work}/Iter${It}.copies.fa

    #Segmental duplications are not independent copies
    bash ${UtilsDir}/CollapseNearIdentical.sh ${Work}/Iter${It}.copies.fa ${CPU} ${MaxIdentity} > ${Work}/Iter${It}.collapse.log
    Copies=$(grep -c ">" ${Work}/Iter${It}.copies.fa)

    #Align the copies to the model: match columns plus insert columns, as in Dfam seed alignments
    hmmalign --dna --trim --outformat afa ${Work}/Iter${Prev}.hmm ${Work}/Iter${It}.copies.fa \
        > ${Work}/Iter${It}.hmmalign.afa || { echo "hmmalign failed"; exit 1; }
    #Columns of the kept match states (and the inserts between them), as an uppercase fasta alignment
    python3 - ${Work}/Iter${It}.hmmalign.afa ${KeepS} ${KeepE} > ${Work}/Iter${It}.aln.fa <<'EOF'
import sys
Afa, S, E = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
Names, Seqs = [], []
for line in open(Afa):
    line = line.strip()
    if line.startswith(">"):
        Names.append(line); Seqs.append([])
    elif line:
        Seqs[-1].append(line)
Seqs = ["".join(x) for x in Seqs]
#Insert columns hold lowercase residues or '.'; match columns uppercase residues or '-'
Insert = [any(q[j] == "." or q[j].islower() for q in Seqs) for j in range(len(Seqs[0]))]
Match = [j for j, i in enumerate(Insert) if not i]
First, Last = Match[S - 1], Match[min(E, len(Match)) - 1]
for n, q in zip(Names, Seqs):
    Row = q[First:Last + 1].replace(".", "-").upper()
    if Row.strip("-"):
        print(n + "\n" + Row)
EOF

    #New match columns from the occupancy of the hits; ends with few hits are trimmed
    python3 ${UtilsDir}/HandRF.py --input ${Work}/Iter${It}.aln.fa --output ${Work}/Iter${It} \
        --minOcc ${MinOcc} --trimCoverage ${TrimCov} --name ${Name} > ${Work}/Iter${It}.HandRF.log || exit 1
    Match=$(awk '$1=="MATCH_COLUMNS"{print $2}' ${Work}/Iter${It}.HandRF.log)
    CovTrim=$(awk '$1=="COVERAGE_TRIMMED_MATCH_COLUMNS"{print $2}' ${Work}/Iter${It}.HandRF.log)

    hmmbuild --hand --dna --cpu ${CPU} --seed 1992 --fragthresh 1.0 --maxinsertlen 10 --eentexp -n ${Name} \
        ${Work}/Iter${It}.hmm ${Work}/Iter${It}.stk > /dev/null || { echo "hmmbuild failed"; exit 1; }
    hmmemit -c ${Work}/Iter${It}.hmm > ${Work}/Iter${It}.Consensus.fa

    #Change from the previous model, gain in hits, and drift from the first refined model
    Change=$(ConsensusChange ${Work}/Iter${Prev}.Consensus.fa ${Work}/Iter${It}.Consensus.fa)
    Gain=$(awk -v h=${Hits} -v p=${PrevHits} 'BEGIN{if(p>0) printf "%.4f", (h-p)/p; else print "NA"}')
    Drift=$( [ ${k} -gt 1 ] && ConsensusChange ${Work}/Iter01.Consensus.fa ${Work}/Iter${It}.Consensus.fa || echo 0 )
    echo -e "${k}\t${Hits}\t$(wc -l < ${Work}/Iter${It}.used.bed)\t${Copies}\t${Match}\t${CovTrim}\t$(grep -v ">" ${Work}/Iter${It}.Consensus.fa | tr -d '\n' | wc -c)\t${Gain}\t${Change}\t${Drift}" >> ${Out}.iterations.tsv
    echo "Iteration ${It}: ${Hits} hits (gain ${Gain}), ${Copies} copies after collapsing, ${Match} match columns (${CovTrim} trimmed at low coverage), consensus change ${Change}, drift ${Drift}"

    if awk -v d=${Drift} -v m=${MaxDrift} 'BEGIN{exit !(d > m)}';then
        echo "Iteration ${It}: the consensus drifted more than ${MaxDrift} from iteration 1, keeping the model of iteration ${Prev}"
        break
    fi
    Prev=${It}
    PrevHits=${Hits}
    if [ ${Gain} != "NA" ] && awk -v g=${Gain} -v c=${Change} -v mg=${MinGain} -v mc=${MaxChange} 'BEGIN{exit !(g < mg && c < mc)}';then
        echo "Iteration ${It}: converged"
        break
    fi
done

if [ ${Prev} == "00" ];then
    echo "No genome-wide refinement was possible"
    exit 1
fi

#Final outputs
cp ${Work}/Iter${Prev}.hmm ${Out}.hmm
sed "1s/.*/>${Name}/" ${Work}/Iter${Prev}.Consensus.fa > ${Out}.Consensus.fa
cp ${Work}/Iter${Prev}.stk ${Out}.stk
cp ${Work}/Iter${Prev}.aln.fa ${Out}.aln.fa
cp ${Work}/Iter${Prev}.match.aln.fa ${Out}.RFclean.aln.fa
cp ${Work}/Iter${Prev}.support.tsv ${Out}.support.tsv
echo "Refined model: ${Out}.hmm ($(grep "^LENG" ${Out}.hmm | awk '{print $2}') match states, $(grep -c ">" ${Out}.aln.fa) sequences)"
