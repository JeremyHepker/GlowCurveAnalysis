import csv
import matplotlib.pylab as plt
import numpy as np
import sys
import csv
import os

def main(directory = sys.argv[1]):
    # open the file in universal line ending mode
    os.makedirs(directory+"/output_plots", exist_ok=True);
    for filename in os.listdir(directory):
        try:
            if(filename.startswith('.')):
                continue
            if(filename == "temp_v_FOM.csv" or filename == "barcode_v_FOM.csv"):
                header = []
                value = []
                with open(directory+"/"+filename) as file:
                    for line in file:
                        d = line.split(',')
                        header.append(float(line.split(',')[0]))
                        value.append(float(line.split(',')[1]))

                fig, ax = plt.subplots()
                ax.scatter(header, value)
                ax.axis([0, 12, 0, 0.10])
                ax.set(xlabel='Time Temperature Profile (C/s)', ylabel='Figure of Merit (%)')
                #ax.grid()
                ax.tick_params(labelsize=12)
                z = np.polyfit(header, value, 1)
                p = np.poly1d(z)
                ax.plot(header, p(header), 'm-')
                fig.savefig(directory+"/output_plots/"+filename+"_plot.png",format='png', dpi=1200)
                plt.close()       
            else:
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

                fig, ax = plt.subplots()
                for d in data:
                    if d != 'temp':
                        ax.plot(data['temp'],data[d])
                ax.set(xlabel='Temperature (K)', ylabel='Light counts (arbitrary units)')
                ax.xaxis.label.set_size(11)
                ax.yaxis.label.set_size(11)
                ax.tick_params(labelsize=12)
                ax.axis('off')
                fig.savefig(directory+"/output_plots/"+filename+"_plot.png",format='png', dpi=1000);
                plt.close()
        except IsADirectoryError:
            pass


if __name__ == "__main__":
    # pylint: disable=no-value-for-parameter
    main()
