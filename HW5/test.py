x = []
y = []
file = open("input.data")
for line in file:
    tmp = line.split(' ')
    x.append(float(tmp[0]))
    y.append(float(tmp[1].split('\n')[0]))
print(type(x[0]))
print(x)
print(y)
