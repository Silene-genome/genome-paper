# input : latifolia annotations, and conica annotations, for the coordinates of the genes, and the orthologs and gametologs

# output : set of orthologs and gametologs with conserved syntenies

#gene_positions = open("../data_raw/S.latifolia_v3.0.chr.gff3","r").readlines()
gene_positions_latifolia = open("../data_raw/S.latifolia_v4.0.gff3_polished","r").readlines()
gene_positions_conica = open("../data_raw/S_conica_v.3.gff","r").readlines()
gene_positions_vulgaris = open("../data_raw/vulgaris.for.eric","r").readlines()
file_orthologs = open("../data_raw/single.copy.ortho.txt","r").readlines()
file_strata = open("../data_raw/strata_RCRV.csv","r").readlines()
#file_strata = open("../data_raw/strata_genes_coordinates.txt","r").readlines()

# read strata
s = {}
for line in file_strata:
    words = line.strip().split()
    name_strata = words[0]
    min_strata = int(words[1])
    max_strata = int(words[2])
    s[name_strata] = [min_strata,max_strata]

# definition of strata, 
def strata(g):

    if genes[g][0] == "chrX":
        for name_strata in s:

            if int(genes[g][1]) >= s[name_strata][0] and int(genes[g][1]) <= s[name_strata][1]:

                return name_strata

        print("Warning, no strata found")
    else:
        return 0

# read gene coordinates : latifolia (X and Y)
genes = {"NA":["0",0,0,"0"]} # assignation name -> [chromosome,start,stop,strand]
for line in gene_positions_latifolia:
    words = line.strip().split("\t")
    if len(words) > 8 and words[2] == "mRNA":
        chromosome = words[0]
        start = int(words[3])
        stop = int(words[4])
        strand = words[6]
        name = words[8].split(";")[0].split("=")[1]
        genes[name] = [chromosome,start,stop,strand]


# read gene coordinates : conica
for line in gene_positions_conica:
    words = line.strip().split("\t")
    if len(words) > 8 and words[2] == "mRNA":
        chromosome = words[0]
        start = int(words[3])
        stop = int(words[4])
        strand = words[6]
        name = words[8].split(";")[0].split("=")[1]
        genes[name] = [chromosome,start,stop,strand]

for line in gene_positions_vulgaris:
    words = line.strip().split("\t")
    chromosome = words[0]
    start = int(words[1])
    stop = int(words[2])
    strand = words[3]
    name = words[4].split(";")[0].split("=")[1]
    genes[name] = [chromosome,start,stop,strand]




# read gametology and orthology
orthologs = []
for line in file_orthologs[1:]:
    words = line.strip().split()
    gene_X_name = words[2]
    gene_Y_name = words[3]
    gene_conica_name = words[1]
    gene_vulgaris_name = words[4]
    
    #print(gene_conica_name,gene_conica_name in genes)
    if (gene_X_name in genes and
        gene_Y_name in genes and
        gene_conica_name in genes and 
        gene_vulgaris_name in genes):
        #genes[gene_X_name][0] == "chrX" and 
        #genes[gene_Y_name][0] == "chrY" and 
        #genes[gene_conica_name][0] == "Chr05" and 
        #genes[gene_vulgaris_name][0] == "scaffold_1"):
        orthologs.append([gene_X_name]+genes[gene_X_name]+
                         [gene_Y_name]+genes[gene_Y_name]+
                         [gene_conica_name]+genes[gene_conica_name]+
                         [gene_vulgaris_name]+genes[gene_vulgaris_name])


# orthologs: [0-4] X
#            [5-9] Y
#            [10-14] conica
#            [15-19] vulgaris
#            [20] number
#            [21] size synteny XY
#            [22] size synteny Xconica
#            [23] size synteny Xvulgaris
#            [24] strata

#print(len(orthologs))
  
# compute numbers and signs
orthologs.sort(key=lambda x:(x[1], x[2])) # sort according to X latifolia
print(len(orthologs),"orthologs")

for i in range(len(orthologs)):
    # signed version
    #if orthologs[i][4] == "+" and orthologs[i][14] == "+":
        #sign = 1
    #elif orthologs[i][4] == "-" and orthologs[i][14] == "-":
        #sign = 1
    #elif orthologs[i][4] == "-" and orthologs[i][14] == "+":
        #sign = -1
    #elif orthologs[i][4] == "+" and orthologs[i][14] == "-":
        #sign = -1
    #orthologs[i].append((i+1)*sign)
    orthologs[i].append(i+1)



# compute synteny blocs Y
orthologs.sort(key=lambda x:(x[6], x[7])) # sort according to Y latifolia


#synteny_output = open("../data_cooked/synteny_strata","w")
size = 0
i = 0
depart = 0
while i < len(orthologs)+1:
    if ((i > 0) and 
        ((i == len(orthologs)) or
         (orthologs[i][1] != orthologs[i-1][1]) or
         (orthologs[i][6] != orthologs[i-1][6]) or
         # signed version
         #(orthologs[i][20] != orthologs[i-1][20] + 1))):
         # unsigned version
         (abs(abs(orthologs[i][20]) - abs(orthologs[i-1][20])) != 1))):
            for j in range(depart,i):
                orthologs[j].append(size)
            size = 1
            depart = i
    else:
        size = size + 1
    i = i + 1

# compute synteny blocs con
orthologs.sort(key=lambda x:(x[11], x[12])) # sort according to conica

size = 0
i = 0
depart = 0
while i < len(orthologs)+1:
    if ((i > 0) and 
        ((i == len(orthologs)) or
         (orthologs[i][1] != orthologs[i-1][1]) or
         (orthologs[i][11] != orthologs[i-1][11]) or
         # signed version
         #(orthologs[i][20] != orthologs[i-1][20] + 1))):
         # unsigned version
         (abs(abs(orthologs[i][20]) - abs(orthologs[i-1][20])) != 1))):
            for j in range(depart,i):
                orthologs[j].append(size)
            size = 1
            depart = i
    else:
        size = size + 1
    i = i + 1
    

# compute synteny blocs vulg
orthologs.sort(key=lambda x:(x[16], x[17])) # sort according to vulgaris

size = 0
i = 0
depart = 0
while i < len(orthologs)+1:
    if ((i > 0) and 
        ((i == len(orthologs)) or
         (orthologs[i][1] != orthologs[i-1][1]) or
         (orthologs[i][16] != orthologs[i-1][16]) or
         # signed version
         #(orthologs[i][20] != orthologs[i-1][20] + 1))):
         # unsigned version
         (abs(abs(orthologs[i][20]) - abs(orthologs[i-1][20])) != 1))):
            for j in range(depart,i):
                orthologs[j].append(size)
            size = 1
            depart = i
    else:
        size = size + 1
    i = i + 1
    

# write summarized data
orthologs.sort(key=lambda x:(x[1], x[2])) # sort according to X latifolia
output = open("../data_cooked/orthogametologs_init.csv","w")
output.write("name_gene_X chr_X start_X stop_X strand_X name_gene_Y chr_Y start_Y stop_Y strand_Y name_gene_conica chr_conica start_conica stop_conica strand_conica name_gene_vulgaris chr_vulgaris start_vulgaris stop_vulgaris strand_vulgaris gene_number synteny_Y synteny_con synteny_vulg strata\n")
for i in range(len(orthologs)):
    for j in orthologs[i]:
        output.write(str(j)+" ")
    #output.write("strata0\n")
    output.write(str(strata(orthologs[i][0]))+"\n")


# orthologs: [0-4] X
#            [5-9] Y
#            [10-14] conica
#            [15-19] vulgaris
#            [20] number
#            [21] size synteny XY
#            [22] size synteny Xconica
#            [23] size synteny Xvulgaris
#            [24] strata

# cured file: filter all isolated genes, remove
i = 0
while i < len(orthologs):
    #print(orthologs[i][11],orthologs[i][6],orthologs[i][16],orthologs[i][17])
    if (
        
        # FILTERS BY SYNTENY, CAN BE REMOVED FOR THE S2 analysis
        (orthologs[i][21] >= 2) and # with Y
        (orthologs[i][22] >= 2) and # with con
        (orthologs[i][23] >= 2) and # with vulg
        
        ## FILTERS: S3 is always removed, and there is an "only S2" mode
        #(strata(orthologs[i][0]) == "S2") and ###### only for rearrangements of S2
        (strata(orthologs[i][0]) != "S3") and
        
        (orthologs[i][1] == "chrX") and 
        (orthologs[i][6] == "chrY") and 
        (orthologs[i][11] == "Chr05") and 
        (orthologs[i][16] == "scaffold_1") and
        (True)):
        del orthologs[i][20:]
        i = i + 1
    else:
        del orthologs[i]

print(len(orthologs),"orthologs after filtering by synteny")

output_badger = open("../data_cooked/silene.infile","w")

orthologs.sort(key=lambda x:(x[1], x[2]))

for i in range(len(orthologs)):
    orthologs[i].append((i+1))
    
output = open("../data_cooked/chr_X","w")
orthologs.sort(key=lambda x:(x[1], x[2]))
output_badger.write("latifolia_X,start,")
#output.write("0 ")
for i in range(len(orthologs)):
    number = orthologs[i][20]
    if orthologs[i][4] == "-":
        number = -number
    output.write(str(number) + " ")

    output_badger.write(str(number) + ",")
output_badger.write("PAR\n")
output.write(str(len(orthologs)+1)+"\n")
 
# write Y : achtung : opposite order because PAR is LEFT while in X PAR is right
output = open("../data_cooked/chr_Y","w")
orthologs.sort(key=lambda x:(x[6], x[7]))
output_badger.write("latifolia_Y,-PAR,")
chaine = ""
for i in range(len(orthologs)):
    number = orthologs[i][20]
    if orthologs[i][9] != "-":
        number = -number
    chaine = str(number) + " " + chaine
    output_badger.write(str(-number) + ",")
output_badger.write("-start\n")
output.write(chaine + str(len(orthologs)+1)+"\n")

output = open("../data_cooked/chr_conica","w")
orthologs.sort(key=lambda x:(x[11], x[12]))
output_badger.write("conica_5,start,")
for i in range(len(orthologs)):
    number = orthologs[i][20]
    if orthologs[i][14] == "-":
        number = -number
    output.write(str(number) + " ")
    output_badger.write(str(number) + ",")
output_badger.write("PAR\n")
output.write(str(len(orthologs)+1)+"\n")

output = open("../data_cooked/chr_vulgaris","w")
orthologs.sort(key=lambda x:(x[16], x[17]))
output_badger.write("vulgaris_1,start,")
for i in range(len(orthologs)):
    number = orthologs[i][20]
    if orthologs[i][19] == "-":
        number = -number
    output.write(str(number) + " ")
    output_badger.write(str(number) + ",")
output_badger.write("PAR\n")
output.write(str(len(orthologs)+1)+"\n")

orthologs.sort(key=lambda x:(x[1], x[2]))
output = open("../data_cooked/orthogametologs_cured.csv","w")
output.write("name_gene_X chr_X start_X stop_X strand_X name_gene_Y chr_Y start_Y stop_Y strand_Y name_gene_conica chr_conica start_conica stop_conica strand_conica name_gene_vulgaris chr_vulgaris start_vulgaris stop_vulgaris strand_vulgaris gene_number strata\n")
for i in range(len(orthologs)):
    #print(orthologs[i])
    for j in orthologs[i]:
        output.write(str(j)+" ")
    output.write(str(strata(orthologs[i][0]))+"\n")


