# reads the ancestral states from the output of badger

file_arrangements = open("../data_cooked/silene.res.prm","r").readlines()
infile = open("../data_cooked/silene.infile","r").readline().split(",")[2:-1]
correction = []
for i in infile:
    if int(i) < 0:
        correction.append(abs(int(i)))

half_length = int(len(file_arrangements)/2)
file_arrangements = file_arrangements[half_length:]


arrangements = {}
for line in file_arrangements:
    value = line.split("(")[0]
    permutation = line.split("(")[1].split(")")[0]
    if not permutation in arrangements:
        arrangements[permutation] = 0
    arrangements[permutation] = arrangements[permutation] + 1
    
keys = list(arrangements.keys())

keys = sorted(arrangements.keys(), key=lambda x: -arrangements[x])

a = keys[0]
chaine = ""
permutation = a.split(",")[:-1]
for p in permutation:
        valeur = int(p)
        if abs(valeur) in correction:
            valeur = -valeur
        chaine = chaine + str(valeur)+" "
print(chaine)
        
output = open("../data_cooked/silene.res.ancestral_arrangements","w")
for a in keys:
    output.write(str(arrangements[a])+" : ")
    permutation = a.split(",")[:-1]
    for p in permutation:
        valeur = int(p)
        if abs(valeur) in correction:
            valeur = -valeur
        output.write(str(valeur)+" ")
    output.write("\n")
