"""
1. program goes through DANTE annotation TABLE and look for 
corresponding TE id in EDTA GFF file and adds annotation
to attributes
2. while parsing table and gff file the script collects
information on number of full-length TEs and length of their
fragmented siblings and generates two tables:
(i) Full-length LTR RTs:
chr, TEid, len, annot, LTRsim, LTRlen
(ii) Rate of full-length to fragmented for each TE (LTR/DNA):
chr, TEid, fl_cnt, fl_len, fr_len
"""

import sys
from collections import namedtuple

gff, dante_tab = sys.argv[1], sys.argv[2]

# define namedtuple for each track line
GffLine = namedtuple('GffLine', ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes'])

def attributes2dict(attributes):
    """transform attributes to dictionary"""
    return {rec.split("=")[0]:rec.split("=")[1] for rec in attributes.split(";")}

def annot_line(gfl_obj):
    # fragmentation status
    fr_status = "full_length"
    if "TE_homo_" in gfl_obj.attributes['ID']:
        fr_status = "fragmented"

    # DANTE annotation
    teid = "_".join(gfl_obj.attributes['Name'].split("_")[:2])
    annot = "NA"
    with open(dante_tab) as tb:
        for rec in tb:
            if rec.startswith(teid+"|") or rec.startswith(teid+"_"):
                annot = rec.rstrip().split("\t")[2]
    return annot, fr_status


out_name = "_".join(gff.split(".")[:-2])+"_reAnnot.gff3"

with open(out_name,"w") as out:
    with open(gff) as gff_in:
        for line in gff_in:
            if not line.startswith("#"):
                #print(line)
                line_list = line.rstrip().split("\t")
                line_list[-1] = attributes2dict(line_list[-1])
                new_attributes = ""
                #print(*line_list)
                gl = GffLine(*line_list)
                if 'LTR' in gl.attributes['Classification']:
                    annot, fr_status = annot_line(gl)
                    if 'LTR' in annot:                
                        if gl.attributes['Classification'] == 'LTR/unknown' and annot != 'NA':
                            #print(annot)
                            gl.attributes['Classification'] = 'LTR/' + annot.split("|")[2].split("/")[1].capitalize()
                            new_attributes = ";".join([k+"="+v for k,v in gl.attributes.items()]) + f";dante_annotation={annot};frag_status={fr_status}"
                        else:
                            new_attributes = ";".join([k+"="+v for k,v in gl.attributes.items()]) + f";dante_annotation={annot};frag_status={fr_status}"
                        #print(new_attributes)
                if new_attributes:
                    line_list[-1] = new_attributes
                    line = "\t".join(line_list)+"\n"
            out.write(line)
