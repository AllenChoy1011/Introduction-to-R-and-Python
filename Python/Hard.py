mutation = {134: 'A', 443: 'G', 769: 'C', 955: 'C', 990: 'C',
            1051: 'A', 1078: 'T', 1941: 'A', 2138: 'C', 2638: 'T', 3003: 'T'}

print ("            T134A   A443G   G769C   G955C   A990C  G1051A  G1078T  T1941A  T2138C   G2638T  A3003T")

with open('/Users/allenchoy/Desktop/stu_10.dat') as file:
    lines = file.readlines()

    patient = 1
    for lineIdx, lineVal in enumerate(lines):
        if 'BLASTN 2.2.31+' in lineVal:     # Check if it's a new patient
            print('pat'+str(patient).zfill(2), end='')
            patient += 1

        if 'Query_' in lineVal:     # Check if it's the query line
            queryLine = lineVal.split()
            hapA = lines[lineIdx+1].split()     # Store hapA below query line
            hapB = lines[lineIdx+2].split()     # Store hapB below query line

            index = 0
            position = int(queryLine[1])
            while position <= int(queryLine[3]):        # Iterate through hap and find mutations
                count = 0
                if position in mutation:        # Check if it's on the wanted positions
                    if(hapA[2][index] != '.' and hapA[2][index] == mutation[position]):     # Check if there is a mutation on hapA
                        count += 1
                    if(hapB[2][index] != '.' and hapB[2][index] == mutation[position]):     # Check if there is a mutation on hapB
                        count += 1
                    print("       " + str(count), end='')

                    if(position == 3003):   
                        print()     # Print a new line for the next patient

                position += 1
                index += 1
