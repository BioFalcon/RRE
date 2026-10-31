#!/usr/bin/env python3

import getopt
import sys

def usage():
        print("Script to prepare the seed alignment of AncientMode: sets the match columns (#=GC RF) of each family")
        print("for hmmbuild --hand and writes the families' names and classes\n")
        print("Usage: SeedPrep.py -i <Alignment.stk> -o <OutputPrefix> [Options]")
        print("\nArguments:")
        print("MANDATORY:")
        print("--input           | -i\t Seed alignment(s) in Stockholm format")
        print("--output          | -o\t Output prefix (writes .stk and .classes.tsv)")
        print("OPTIONAL (without them, every family must have a #=GC RF line):")
        print("--rfAll           | -a\t Every column with a residue is a match column")
        print("--rfFromOccupancy | -c\t Columns where at least --rfMinOcc of the sequences spanning them have a base are match columns")
        print("--rfMinOcc        | -m\t Min. fraction for --rfFromOccupancy (Default: 0.5)")
        print("--help            | -h\t Show this help message and exit")
        exit()

#Characters that mark a gap in a sequence, or a non-match column in #=GC RF
GapChars = ".-~_"

def ReadStockholm(input):
    #Families in file order: ID, class (#=GF TP), sequences (order kept, interleaved blocks joined) and RF
    Families = []
    Cur = None
    with open(input) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("# STOCKHOLM"):
                Cur = {"ID": None, "TP": None, "Order": [], "Seqs": {}, "RF": ""}
            elif line.startswith("//"):
                if Cur is not None:
                    Families.append(Cur)
                Cur = None
            elif Cur is None or not line.strip():
                continue
            elif line.startswith("#=GF"):
                Fields = line.split(None, 2)
                if len(Fields) == 3 and Fields[1] == "ID":
                    Cur["ID"] = Fields[2].strip()
                elif len(Fields) == 3 and Fields[1] == "TP":
                    Cur["TP"] = Fields[2].strip()
            elif line.startswith("#=GC"):
                Fields = line.split()
                if len(Fields) == 3 and Fields[1] == "RF":
                    Cur["RF"] += Fields[2]
            elif not line.startswith("#"):
                Name, Seq = line.split()[:2]
                if Name not in Cur["Seqs"]:
                    Cur["Order"].append(Name)
                    Cur["Seqs"][Name] = []
                Cur["Seqs"][Name].append(Seq)
    return Families

def OccupancyRF(Seqs, minOcc, minRes):
    #Same rule as HandRF.py: bases / sequences spanning the column (first to last residue of each sequence)
    Width = len(Seqs[0])
    Spanning = [0] * Width
    NBase = [0] * Width
    for s in Seqs:
        Res = [j for j, c in enumerate(s) if c not in GapChars]
        if not Res:
            continue
        for j in range(Res[0], Res[-1] + 1):
            Spanning[j] += 1
        for j in Res:
            if s[j].upper() in "ACGT":
                NBase[j] += 1
    return "".join("x" if Spanning[j] and NBase[j] / Spanning[j] >= minOcc and NBase[j] >= minRes else "." for j in range(Width))

def main():
    try:
        options, remainder = getopt.getopt(sys.argv[1:], 'i:o:acm:h', ['input=','output=','rfAll','rfFromOccupancy','rfMinOcc=','help'])
    except getopt.GetoptError as err:
        print(str(err))
        usage()

    #Default values
    RFAll = False
    RFOcc = False
    minOcc = 0.5
    Input = None
    Output = None

    for opt, arg in options:
        if opt in ('-i', '--input'):
            Input = arg
        elif opt in ('-o', '--output'):
            Output = arg
        elif opt in ('-a', '--rfAll'):
            RFAll = True
        elif opt in ('-c', '--rfFromOccupancy'):
            RFOcc = True
        elif opt in ('-m', '--rfMinOcc'):
            minOcc = float(arg)
        elif opt in ('-h', '--help'):
            usage()
    if Input is None or Output is None:
        usage()
    if RFAll and RFOcc:
        sys.exit("ERROR: --rfAll and --rfFromOccupancy are alternatives, use only one")

    Families = ReadStockholm(Input)
    if not Families:
        sys.exit("ERROR: no Stockholm alignment found in " + Input)

    Errors = []
    IDs = [F["ID"] for F in Families if F["ID"]]
    for i in sorted(set(x for x in IDs if IDs.count(x) > 1)):
        Errors.append(i + ": the same ID is used by %d alignments" % IDs.count(i))
    Out = open(Output + ".stk", "w")
    Classes = open(Output + ".classes.tsv", "w")
    for k, F in enumerate(Families):
        Label = F["ID"] if F["ID"] else "alignment %d" % (k + 1)
        if not F["ID"]:
            Errors.append(Label + ": no #=GF ID line")
            continue
        if any(c.isspace() or c == "#" for c in F["ID"]):
            Errors.append(Label + ": the ID cannot contain spaces or '#'")
            continue
        Seqs = ["".join(F["Seqs"][n]) for n in F["Order"]]
        if not Seqs:
            Errors.append(Label + ": no sequences")
            continue
        Width = len(Seqs[0])
        if any(len(s) != Width for s in Seqs):
            Errors.append(Label + ": sequences have different lengths")
            continue

        #Match columns: the alignment's own RF line, unless --rfAll or --rfFromOccupancy define them
        if RFAll:
            RF = "".join("x" if any(s[j] not in GapChars for s in Seqs) else "." for j in range(Width))
            How = "all columns with a residue"
        elif RFOcc:
            RF = OccupancyRF(Seqs, minOcc, min(3, len(Seqs)))
            How = "occupancy >= %g" % minOcc
        elif F["RF"]:
            RF = F["RF"]
            How = "#=GC RF line"
            if len(RF) != Width:
                Errors.append(Label + ": the #=GC RF line has %d columns, the alignment %d" % (len(RF), Width))
                continue
        else:
            Errors.append(Label + ": no #=GC RF line (add one, or run with --rfAll or --rfFromOccupancy)")
            continue
        NMatch = sum(c not in GapChars for c in RF)
        if NMatch == 0:
            Errors.append(Label + ": no match columns (" + How + ")")
            continue

        Class = F["TP"].replace(" ", "_") if F["TP"] else "Unknown"
        Pad = max(len(n) for n in F["Order"] + ["#=GC RF"]) + 1
        Out.write("# STOCKHOLM 1.0\n#=GF ID %s\n" % F["ID"])
        if F["TP"]:
            Out.write("#=GF TP %s\n" % F["TP"])
        for n, s in zip(F["Order"], Seqs):
            Out.write(n.ljust(Pad) + s + "\n")
        Out.write("#=GC RF".ljust(Pad) + RF + "\n//\n")
        Classes.write("%s\t%s\n" % (F["ID"], Class))
        print("%s: %d sequences, %d columns, %d match columns (%s), class %s" % (F["ID"], len(Seqs), Width, NMatch, How, Class))
    Out.close()
    Classes.close()

    if Errors:
        sys.exit("ERROR in the seed alignment " + Input + ":\n  " + "\n  ".join(Errors))

if __name__ == "__main__":
    main()
