import csv
import matplotlib.pylab as plt
import sys
import csv
def main(filename = sys.argv[1]):
    # open the file in universal line ending mode
    with open(filename, 'r') as infile:
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
    plt.title(filename)
    plt.subplots_adjust(left=0.0, right=1.0, bottom=0.0, top=1.0)
    plt.show(block=False)
    input("Hit Enter To Close")
    plt.close()


if __name__ == "__main__":
    # pylint: disable=no-value-for-parameter
    main()
