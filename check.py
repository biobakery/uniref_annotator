import sys
from collections import Counter

total = 0
counts = Counter( )

with open( sys.argv[1] ) as fh:
    for line in fh:
        line = line.strip( )
        if line[0] == ">":
            items = line.split( "|" )
            total += 1
            a90, a50, b90, b50 = items[-5], items[-4], items[-2], items[-1]

            if "unknown" not in a90 + b90:
                if a90 == b90:
                    counts["uniref90 same"] += 1
                else:
                    counts["uniref90 diff"] += 1
            elif "unknown" in a90 and "unknown" in b90:
                counts["uniref90 ?->?"] += 1
            elif "unknown" in a90:
                counts["uniref90 ?->*"] += 1
            elif "unknown" in b90:
                counts["uniref90 *->?"] += 1
                #print line
            if "unknown" not in a50 + b50:
                if a50 == b50:
                    counts["uniref50 same"] += 1
                else:
                    counts["uniref50 diff"] += 1
            elif "unknown" in a50 and "unknown" in b50:
                counts["uniref50 ?->?"] += 1
            elif "unknown" in a50:
                counts["uniref50 ?->*"] += 1
            elif "unknown" in b50:
                counts["uniref50 *->?"] += 1

for k in sorted( counts ):
    print k, "{:.3f}".format( 1.0 * counts[k] / total ), counts[k]
