
import csv
import matplotlib.pylab as plt
import numpy as np
import sys
import csv
import os
one = (1133,1095,1105,1068,1042,1041,1044,1090,1108,1057)
two = (1126,1080,1103,1051,1064,1047,1088,1053,1139,1136)
three = (1034,1056,1073,1027,1098,1065,1059,1092,1046,1093)
four = (1138,1067,1077,1124,1125,1132,1061,1048,1029,1084)
five = (1128,1107,1113,1116,1123,1085,1079,1104,1109,1084)
six = (1094,1091,1086,1137,1039,1120,1057,1127,1131,1083)
seven = (1032,1038,1030,1060,1070,1089,1111,1135,1082,1087)
eight = (1036,1096,1121,1136,1115,1043,1075,1045,1114,1037)
nine = (1100,1063,1055,1050,1122,1054,1081,1074,1066,1099)
ten = (1133,1130,1100,1126,1034,1138,1128,1094,1032,1036)
def main(file = sys.argv[1]):
    # open the file in universal line ending mode
    header = []
    value = []
    valueNew = []
    with open(file+ "/dose_vs_fom.csv") as f:
        for line in f:
            d = line.split(',')
            header.append(float(line.split(',')[0]))
            value.append(float(line.split(',')[1]))
    fig, ax = plt.subplots()
#    for i in range( len(header)):
#        if header[i] in one:
#            new = 2.43
#            valueNew.append(new)
#        elif header[i] in two:
#            new = 4.853
#            valueNew.append(new)
#        elif header[i] in three:
#            new = 7.05
#            valueNew.append(new)
#
#        elif header[i] in four:
#            new = 9.705
#            valueNew.append(new)
#
#        elif header[i] in five:
#            new = 14.558
#            valueNew.append(new)
#
#        elif header[i] in six:
#            new = 15
#            valueNew.append(new)
#
#        elif header[i] in seven:
#            new = 19.41
#            valueNew.append(new)
#
#        elif header[i] in eight:
#            new = 24.26
#            valueNew.append(new)
#
#        elif header[i] in nine:
#            new = 29.117
#            valueNew.append(new)
#
#        elif header[i] in ten:
#            new = 30
#            valueNew.append(new)
#        else:
#            valueNew.append(7.05)
    ax.scatter(header, value)
    ax.axis([0, 35, 0, 0.2])
    ax.set(xlabel='Dose (mGy)', ylabel='Figure of Merit (%)')
    #ax.grid()
    ax.tick_params(labelsize=12)
    z = np.polyfit(header, value, 1)
    p = np.poly1d(z)
    ax.plot(header, p(header), 'm-')
    fig.savefig(file + "/output_plots/dose_vs_fom_plot.png",format='png', dpi=1200)
#    with open("dose_vs_fom.csv", "w+") as f:
#        for i in range(len(value)):
#            f.write(str(valueNew[i])+","+str(value[i])+"\n")
    plt.close()

if __name__ == "__main__":
    # pylint: disable=no-value-for-parameter
    main()

