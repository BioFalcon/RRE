process HMMSelection{
    label "HMMSelection"

    input:
    file(HMMOutput)
    file(consensusFasta)

    output:
    file("02_Selection.txt")

    script:
    """
    #Select those copies that have at least 10 hits
    cut -f3 ${HMMOutput} | \\
    uniq -c | \\
    sed 's/^\\s\\{1,\\}//;s/\\s\\{1,\\}/\\t/' | \\
    awk '\$1 >= 10 { print \$2}' > 01_Selection.txt

    #Remove LTRs 
    grep -f <(grep ">" ${consensusFasta} | grep -v "Type=LTR"| sed 's/#.*/\\\$/;s/>//') \\
    01_Selection.txt > 02_Selection.txt
    """
}