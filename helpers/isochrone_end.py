import csv
from numpy import median
from cluster_isochrone import find_data_in_files_beta

# B-V V
# R-I I


csvFile = 'R-I-I.csv'
writefile = open(csvFile, 'a', newline='', encoding='utf-8')
writer = csv.writer(writefile, delimiter=',', quotechar='\'')


age = 6.60
while age < 10.25:
    metallicity = -2.2
    while metallicity < 0.75:
        data = find_data_in_files_beta(age, metallicity, ['R', 'I', 'I'])

        dataArray = [[tup[0], tup[1]] for tup in data]

        deltas = [None]
        deltaXs = [None]
        deltaYs = [None]
        maxDeltaIndex = 0
        count = [0]
        for i in range(len(dataArray)):
            x_i = dataArray[i][0]
            y_i = dataArray[i][1]
            if i > 0:
                deltaX = abs(x_i - dataArray[i - 1][0])
                deltaY = abs(y_i - dataArray[i - 1][1])
                deltaXs.append(deltaX)
                deltaYs.append(deltaY)

        xMedianValue = median(deltaXs[1:])
        yMedianValue = median(deltaYs[1:])
        # From the beginning of delta_i, find the nth = 1st i such that delta_i < sqrt(2).
        # Call it iBeg. KEEP all points before iBeg.
        N = len(deltaXs)
        c1 = 0
        c2 = 0
        c3 = 0
        for i in range(N):
            delta = 0
            if i > 0:
                delta = ((deltaXs[i] / xMedianValue) ** 2 + (deltaYs[i] / yMedianValue) ** 2) ** 0.5
                if delta < (2 ** 0.5):
                    count.append(count[i - 1] - 1)
                else:
                    count.append(count[i - 1] + 1)

            sqaure = (i ** 2) - N * i
            c1 += count[i] * sqaure
            c2 += i * sqaure
            c3 += sqaure ** 2
            deltas.append(delta)

        countN = count[N - 1]
        c2 *= countN / N
        c = (c1 - c2) / c3
        b = countN / N - c * N

        # From the end of delta_i, find the nth = 1st i such that delta_i < sqrt(2).
        # Call it iEnd. REMOVE all points after iEnd.
        iBeg = 0  # init iBeg to be 0
        iEnd = N  # init iEnd as the last count
        minTemp = 0
        maxTemp = 0
        deltas.pop(0)

        for i in range(1, N):
            temp = count[i] - (b * i + c * (i ** 2))
            if temp < minTemp:
                minTemp = temp
                iEnd = i
            elif temp >= maxTemp:
                maxTemp = temp
                iBeg = i
        if iBeg > iEnd:
            iBeg = 0
            maxTemp = 0
            for i in range(1, iEnd):
                temp = count[i] - i / (N) * count[N - 1]
                if temp >= maxTemp:
                    maxTemp = temp
                iBeg = i

        # maxDeltaIndex = deltas.indexOf(Math.max.apply(null, deltas.slice(iBeg-1, iEnd+1)));
        # if type(deltas) == 'float': print(deltas)
        iSkip = round(-33.185*age+456.77)
        # iSkip = deltas.index(max(deltas[iBeg-1:iEnd+1]))
        iBeg_skip = iSkip - 20 if iSkip >= 20 else 0
        iEnd_skip = iSkip + 40 if iSkip < round(-25.84*age+451.77) else round(-25.84*age+451.77)
        iEnd_skip = iEnd_skip if iEnd_skip < N else N
        iSkip = deltas.index(max(deltas[iBeg_skip:iEnd_skip]))
        writer.writerow([format(age, ".2f"), format(metallicity, ".2f"), iSkip])
        metallicity += 0.05
    age += 0.05


writefile.close()
