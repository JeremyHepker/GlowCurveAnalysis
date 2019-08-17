import csv
import matplotlib.pylab as plt
import numpy as np
import sys
import csv
import os

def smooth(x,window_len=11,window='hanning'):
    if window_len<3:
        return x
    s = np.r_[x[window_len-1:0:-1],x,x[-2:-window_len-1:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')
    
    y=np.convolve(w/w.sum(),s,mode='valid')
    return y

def main(directory = sys.argv[1]):
    # open the file in universal line ending mode
    os.makedirs(directory+"/clean_data", exist_ok=True);
    for filename in os.listdir(directory):
        header = []
        value = []
        newValue = []
        count = 0
        if filename == "clean_data" or filename == ".DS_Store":
            continue
        with open(directory+filename) as f:
            for line in f:
                if count > 8:
                    print(line)
                    try:
                        sampleNo = line.split(',')[0]
                        barcode = line.split(',')[1]
                        time = line.split(',')[2]
                        header.append(float(line.split(',')[3]))
                        value.append(float(line.split(',')[4]))
                    except:
                        pass
                count += 1
        newValue.append(value[0])
        for v in value:
            if v > (newValue[-1] + 30) or v < (newValue[-1] - 30):
                newValue.append(newValue[-1])
            else:
                newValue.append(v)
        array = np.array(newValue)
        y = smooth(array)
        y = smooth(y)
        with open(directory+"clean_data/"+filename, "w+") as outfile:
            outfile.write("barcode,temp,Count\n")
            for i in range(len(header)):
                outfile.write("barcode"+","+str(header[i])+","+str(y[i])+"\n")

if __name__ == "__main__":
    # pylint: disable=no-value-for-parameter
    main()

import numpy

