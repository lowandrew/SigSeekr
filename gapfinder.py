__author__ = 'mikeknowles'

gapfile = "/home/blais/Downloads/P719KYB7114-Alignment.txt"
gapdict = {}
count = 0
for line in open("/home/blais/Downloads/P719KYB7114-Alignment.txt"):
    array = line.split("\t")
    start = int(array[6])
    stop = int(array[7])
    if array[1] not in gapdict:
        gapdict[array[1]] = [[start, stop]]
    for coor in gapdict[array[1]]:
        if start <= coor[0]:
            coor[0] = start
            if stop > coor[1]:
                coor[1] = stop
        elif start > coor[1]:
            gapdict[array[1]].append([start, stop])
    # print gapdict[array[1]]
    # count += 1
    if count > 12:
        break
    # for coor in gapdict:
# print gapdict
space = {}
for alignment in gapdict:
    space[alignment] = []
    previous = [1, 1]
    for coor in gapdict[alignment]:

        # print coor, previous
        if coor[0] > previous[1]:
            space[alignment].append((coor[0], previous[1], coor[0] - previous[1]))
        previous = coor
print space



    # print array[1], array[6]
