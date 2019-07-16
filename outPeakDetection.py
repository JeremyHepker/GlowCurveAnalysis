import csv
import matplotlib.pyplot as plt
import numpy as np
import sys
import csv
import os

def main(directory = sys.argv[1]):
    # open the file in universal line ending mode
    fig = plt.figure()
    ax = fig.add_subplot(111);
    for filename in os.listdir(directory):
        if filename == "maxima.csv":
            x = []
            y = []
            with open(directory+"/"+filename, 'r') as infile:
                reader = csv.reader(infile, delimiter=',')
                for row in reader:
                    x.append(float(row[0]))
                    y.append(float(row[1]))
                infile.close()
            ax.scatter(x, y, marker = 'x',label="maxima", alpha=1.0)
        elif filename == "inflection.csv":
            x = []
            y = []
            with open(directory+"/"+filename, 'r') as infile:
                reader = csv.reader(infile, delimiter=',')
                for row in reader:
                    x.append(float(row[0]))
                    y.append(float(row[1]))
                infile.close()
            ax.scatter(x, y, marker = '8',label="inflections", alpha=1.0)
        elif filename == "minimum.csv":
            x = []
            y = []
            with open(directory+"/"+filename, 'r') as infile:
                reader = csv.reader(infile, delimiter=',')
                for row in reader:
                    x.append(float(row[0]))
                    y.append(float(row[1]))
                infile.close()
            ax.scatter(x, y, marker = 's',label="minimum", alpha=1.0)
        elif filename == "curve.csv":
            x = []
            y = []
            with open(directory+"/"+filename, 'r') as infile:
                reader = csv.reader(infile, delimiter=',')
                for row in reader:
                    x.append(float(row[0]))
                    y.append(float(row[1]))
                infile.close()
            ax.plot(x, y)
        elif filename == "peakTest_output.csv":
            with open(directory+"/"+filename, 'r') as infile:
        # read the file as a dictionary for each row ({header : value})
                reader = csv.DictReader(infile)
                data = {}
                for row in reader:
                    for header, value in row.items():
                        try:
                            try:
                                value = float(value)
                            except:
                                value = 0.0
                            data[header].append(value)
                        except KeyError:
                            data[header] = [value]
                infile.close()
            for d in data:
                if d != 'temp':
                    ax.plot(data['temp'],data[d])
        elif filename == "subtraction_output.csv":
            with open(directory+"/"+filename, 'r') as infile:
        # read the file as a dictionary for each row ({header : value})
                reader = csv.DictReader(infile)
                data = {}
                for row in reader:
                    for header, value in row.items():
                        try:
                            try:
                                value = float(value)
                            except:
                                value = 0.0
                            data[header].append(value)
                        except KeyError:
                            data[header] = [value]
                infile.close()
            for d in data:
                if d != 'temp':
                    ax.plot(data['temp'],data[d])
    plt.xlabel('temp')
    plt.ylabel('counts')
    plt.legend()
    fig.savefig(directory+"/plot.png",format='png', dpi=1000);
    #plt.show()

if __name__ == "__main__":
    # pylint: disable=no-value-for-parameter
    main()

        