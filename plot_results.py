filename = ""

steps = []
agents = []

with open(filename) as f:
    for line in f:
        if "steps" in line:
            steps.append([int(line.split(" ")[3]), float(line.split(" ")[5])])
        if "agent" in line:
            agents.append([int(line.split(" ")[3]), float(line.split(" ")[5])])
