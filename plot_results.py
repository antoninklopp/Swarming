import sys
import matplotlib.pyplot as plt

key_word = "steps"

filename = sys.argv[1]
filename2 = sys.argv[2]

variable = [[], []]
variable2 = [[], []]

with open(filename) as f:
    for line in f:
        if key_word in line:
            variable[0].append([int(line.split(" ")[3]))
            variable[1].append(float(line.split(" ")[5])])

with open(filename2) as f:
    for line in f:
        if key_word in line:
            variable2[0].append([int(line.split(" ")[3]))
            variable2[1].append(float(line.split(" ")[5])])

plt.plot(variable[0], variable[1])
plt.plot(variable2[0], variable2[1])
plt.show()
