import sys
import re

gffs = sys.argv[1:]

classif_d = {"45S_rDNA_Slat":"rDNA/45S","CL99_5SrDNA+spacer_364bp":"rDNA/5S",
             "CL119_LINE_2075bp":"LINE",
             "CL12_TR1_74bp":"Satellites/TR1","CL12_X43.1_312bp":"Satellites/X43.1","CL16_STAR_43bp":"Satellites/STAR-C","CL55_LINE_3458bp":"LINE","CL68_TRAYC_165bp":"Satellites/TRAYC",
             "CL84_LINE_3865bp":"LINE","STAR-Y":"Satellites/STAR-Y"}

regex=r'Motif:(.*)\"'
for gff in gffs:
    out_name = gff.replace(".gff","_rmOutMod.gff")
    with open(out_name,"w") as out:
        with open(gff) as gf:
            for l in gf:
                if "Motif" in l and re.search(regex,l):
                    match = (re.search(regex,l).group(1))
                    ll =l.rstrip().split("\t")
                    pref="\t".join(ll[:-1])
                    old_att = "|".join(ll[-1].split())
                    attr = f"Classification={classif_d[match]};RM_attribute={old_att}\n"
                    out.write(pref + "\t" + attr)
