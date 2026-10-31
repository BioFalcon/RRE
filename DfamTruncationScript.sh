#!/bin/bash
# Usage: bash Script_TruncateSTK.sh <input.stk> <start> <end> [min_freq]
#   min_freq: majority-rule cutoff (0-1, default 0.5); columns below it become N

if [ $# -lt 3 ] || [ $# -gt 4 ]; then
  echo "Usage: $0 <input.stk> <start> <end> [min_freq]" >&2
  exit 1
fi

File="$1"
DesStart="$2"
DesEnd="$3"
MinFreq="${4:-0.5}"

[ -f "${File}" ] || { echo "ERROR: ${File} not found" >&2; exit 1; }
[[ "${DesStart}" =~ ^[0-9]+$ && "${DesEnd}" =~ ^[0-9]+$ ]] || { echo "ERROR: start and end must be positive integers" >&2; exit 1; }
[ "${DesStart}" -lt "${DesEnd}" ] || { echo "ERROR: start must be smaller than end" >&2; exit 1; }
[[ "${MinFreq}" =~ ^(0(\.[0-9]+)?|1(\.0+)?)$ ]] || { echo "ERROR: min_freq must be between 0 and 1" >&2; exit 1; }

# Repeat name from "#=GF ID" line; fall back to the file name
RepName=$(awk '$1=="#=GF" && $2=="ID"{print $3; exit}' "${File}")
if [ -n "${RepName}" ]; then HasID=1; else HasID=0; RepName=$(basename "${File}" .stk); fi

OutFile="${RepName}.${DesStart}-${DesEnd}.stk"
OutCons="${RepName}.${DesStart}-${DesEnd}.consensus.fa"
StkName="${RepName}_${DesStart}-${DesEnd}"   # STK ID (no class)
Header="${StkName}#Unknown"                  # FASTA header (with class)

awk -v DesStart=${DesStart} -v DesEnd=${DesEnd} -v ConsFile="${OutCons}" \
    -v StkName="${StkName}" -v Name="${Header}" -v HasID=${HasID} -v MinFreq=${MinFreq} '
{
  if ($0 ~ "#=GC RF") {
    split($3,A,""); Start=0; End=0; Counter=0
    for (i=1; i<=length(A); ++i) {
      if (A[i] ~ "[XATCG]") { Counter+=1; if (Counter==DesStart) Start=i; if (Counter==DesEnd) End=i }
    }
    RF=substr($3,Start,End-Start)
    print $1, $2, RF
    # remember which truncated columns are consensus (RF letter) columns
    n=split(RF,R,""); NC=0
    for (j=1; j<=n; j++) if (R[j] ~ /[A-Za-z]/) { NC++; CCol[NC]=j }
  }
  else if ($1=="#=GF" && $2=="ID") { print "#=GF ID " StkName }
  else if ($0 ~ "^# STOCKHOLM") { print; if (!HasID) print "#=GF ID " StkName }
  else if ($0 !~ "^#") {
    S=substr($2,Start,End-Start)
    print $1"\t"S
    # count bases at each consensus column
    for (k=1; k<=NC; k++) {
      b=toupper(substr(S,CCol[k],1))
      if (b ~ /^[ACGT]$/) Cnt[k,b]++
    }
  }
  else if ($0 ~ "^#" && $2 != "RF") { print }
}
END {
  # majority rule: most frequent base among non-gap bases; N if below MinFreq or no bases
  Cons=""
  split("A C G T", Bases, " ")
  for (k=1; k<=NC; k++) {
    best="N"; bc=0; tot=0
    for (m=1; m<=4; m++) { c=Cnt[k,Bases[m]]+0; tot+=c; if (c>bc) { bc=c; best=Bases[m] } }
    if (tot==0 || bc/tot < MinFreq) best="N"
    Cons=Cons best
  }
  print ">"Name > ConsFile; print Cons > ConsFile
}' "${File}" | sed '/\t\.\{1,\}$/d' > "${OutFile}"

echo "Wrote ${OutFile} (ID: ${StkName}) and ${OutCons} (>${Header}, min_freq=${MinFreq})"
