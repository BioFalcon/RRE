DesStart=2000
DesEnd=2400
File="OldhAT1.stk"
OutFile="OldhAT1.3800-5800.stk"

awk -v DesStart=${DesStart} -v DesEnd=${DesEnd} '{if($0~"#=GC RF"){split($3,A,"");Start=0;End=0;Counter=0;for(i=1;i<=length(A);++i){if(A[i]~"[XATCG]"){Counter+=1;if(Counter==DesStart){Start=i};if(Counter==DesEnd){End=i}}};print $1,$2,substr($3,Start,End-Start)}else if($0 !~ "^#"){print $1"\t"substr($2,Start,End-Start)}else if($0 ~ "^#" && $2 != "RF"){print}}' ${File}| sed '/\t\.\{1,\}$/d' > ${OutFile}
