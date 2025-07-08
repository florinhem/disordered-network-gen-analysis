import matplotlib.pyplot as plt
import numpy as np

data_measured=[
    ["ctn", 28, 0.25, 180.0, 0.1, 0.01, 0.13, 0.14],
    ["ctn", 224, 0.25, 180.0, 0.1, 0.01, 0.40, 0.47],
    ["ctn", 756, 0.25, 180.0, 0.1, 0.01, 1.0, 1.15],

    ["dia", 8, 0.25, 180.0, 0.1, 0.01, 0.17, 0.18],
    ["dia", 64, 0.25, 180.0, 0.1, 0.01, 0.17, 0.18],
    ["dia", 216, 0.25, 180.0, 0.1, 0.01, 0.38, 0.44],

    ["pto", 14, 0.25, 180.0, 0.1, 0.01, 0.07, 0.08],
    ["pto", 112, 0.25, 180.0, 0.1, 0.01, 0.30, 0.35],
    ["pto", 378, 0.25, 180.0, 0.1, 0.01, 1.10, 1.17],

    ["srd", 80, 0.25, 180.0, 0.1, 0.01, 0.21, 0.22],
    ["srd", 270, 0.25, 180.0, 0.1, 0.01, 0.25, 0.29],

    ["srs", 64, 0.25, 180.0, 0.1, 0.01, 0.10, 0.11],
    ["srs", 216, 0.25, 180.0, 0.1, 0.01, 0.20, 0.21],
]

data = [
    ["ctn", 28, 0.25, 180.0, 0.001, 0.08225647682891701 ],
    ["ctn", 224, 0.25, 180.0, 0.001, 0.5004858817153227 ],
    ["ctn", 756, 0.25, 180.0, 0.001, 1.709489436721156 ],
    ["dia", 8, 0.25, 180.0, 0.001, 0.07433115173976412 ],
    ["dia", 64, 0.25, 180.0, 0.001, 0.12629940504665044 ],
    ["dia", 216, 0.25, 180.0, 0.001, 0.35359802369526516 ],
    ["pto", 14, 0.25, 180.0, 0.001, 0.02641560785080463 ],
    ["pto", 112, 0.25, 180.0, 0.001, 0.34237940371427256 ],
    ["pto", 378, 0.25, 180.0, 0.001, 1.4211384360754205 ],
    ["srd", 10, 0.25, 180.0, 0.001, 0.042627801941970835 ],
    ["srd", 80, 0.25, 180.0, 0.001, 0.21951131653347236 ],
    ["srd", 270, 0.25, 180.0, 0.001, 0.9274487670570966 ],
    ["srs", 8, 0.25, 180.0, 0.001, 0.06701158786879931 ],
    ["srs", 64, 0.25, 180.0, 0.001, 0.0948350052111373 ],
    ["srs", 216, 0.25, 180.0, 0.001, 0.4498312032774132 ],
]

networks = ["dia", "srs", "srd", "pto", "ctn"]
colors = np.array([
    (46,37,133),
    (51,117,56),
    (93,168,153),
    (148,203,236),
    (220,205,125),
    (194,106,119),
    (159,74,150),
    (126,41,84),
    (221,221,221)
])/255

plt.axhline(0, color="black", linewidth=1.5, linestyle='-')

for i, net in enumerate(networks):
    x = [row[1] for row in data_measured if row[0] == net]
    y1 = [row[6] for row in data_measured if row[0] == net]
    y2= [row[7] for row in data_measured if row[0] == net]
    print(x)
    print(y1)
    print(y2)
    plt.fill_between(x, y1, y2, alpha=0.25, color=colors[i], label=net)


for i, net in enumerate(networks):
    x = [row[1] for row in data if row[0] == net]
    y = [row[5] for row in data if row[0] == net]
    plt.plot(x, y, "-o", color=colors[i], label=net)


fontsize=11
plt.xlabel(r'$\beta / \alpha$', fontsize=fontsize)
plt.ylabel(r'$T_\mathrm{melt}$',fontsize=fontsize)
plt.legend(fontsize=fontsize)
all_x = sorted(set(row[1] for row in data_measured + data))
plt.xticks(all_x, fontsize=fontsize)
plt.yticks(fontsize=fontsize)
plt.xscale("log")

plt.show()