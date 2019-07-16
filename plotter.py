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
            if(filename.startswith('.') or filename == "statistics.csv"):
                continue

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
#                ax.axis('off')
            fig.savefig(directory+"/output_plots/"+filename+"_plot.png",format='png', dpi=1000);
            plt.close()
        except IsADirectoryError:
            pass


if __name__ == "__main__":
    # pylint: disable=no-value-for-parameter
    main()
