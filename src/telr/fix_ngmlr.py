import sys

for line in sys.stdin:
    line = line.replace("\n","")
    entry = line.split("\t")
    if len(entry) > 5:
        if int(entry[4]) < 0:
            entry[4] = "0"
            print("\t".join(entry))
        else:
            print(line)
    else:
        print(line)