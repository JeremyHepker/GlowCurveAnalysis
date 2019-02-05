import csv
import matplotlib.pylab as plt

import csv

# open the file in universal line ending mode
with open('output.csv', 'rU') as infile:
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

for d in data:
    if d != 'temp':
        plt.plot(data['temp'],data[d])
plt.show()
