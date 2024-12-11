import sys
import os
from os.path import exists
from collections import namedtuple
import subprocess
import re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

gff, fa = sys.argv[1], sys.argv[2]

# define namedtuple for each track line
GffLine = namedtuple('GffLine', ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes'])

def attributes2dict(attributes):
    """transform attributes to dictionary"""
    return {rec.split("=")[0]:rec.split("=")[1] for rec in attributes.split(";")}
## call: GffLine(*ll)

def rep_region_list():
    rrl = []
    with open(gff) as gf:
        for l in gf:
            if "\trepeat_region\t" in l:
                ll = l.rstrip().split("\t")
                ll[-1] = attributes2dict(ll[-1])
                rrl.append(GffLine(*ll))
    return rrl

def find_ltrs(ch, rep_reg_id):
    ltr_tracks = []
    with open(gff) as gf:
        for l in gf:
            if l.startswith(ch+"\t"):
                ll = l.rstrip().split("\t")
                ll[-1] = attributes2dict(ll[-1])
                track = GffLine(*ll)
                if track.type == "long_terminal_repeat" and track.attributes['Parent'] == rep_reg_id:
                    ltr_tracks.append(track)
    if len(ltr_tracks) == 2:
        return ltr_tracks
    else:
        print(f"Not able to find LTRs for repeat_region ID: {rep_reg_id}")
                
def get_ltr_fa_pair(ltrs_tracks):
    ltr_d = {'L':'ltr5','R':'ltr3'}
    ltr_len = []
    for ltr in ltrs_tracks:
        start, end = int(ltr.start)-1, int(ltr.end)
        ltr_len.append(end-start)
        head = ltr_d[ltr.attributes['ID'][0].upper()]
        out_name = head+".fa"
        bed_file = head+".bed"
        with open(bed_file,"w") as bed:
            bed.write(f"{ltr.seqid}\t{start}\t{end}\n")
        cmd = f"bedtools getfasta -fi {fa} -bed {bed_file} > {out_name}"
        process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
        process.wait()
        print(process.returncode)
        os.remove(bed_file)            
    if exists('ltr5.fa') and exists('ltr3.fa'):
        print("FASTA files are written.")
        return ltr_len

def div_from_distmat(dm_file):
    with open(dm_file) as df:
        for l in df:
            if l.rstrip()[-2:] == " 1":
                return float(l.split()[1])
            
def ltr_div(ltr_len):
    files2remove = ['ltr5.fa','ltr3.fa']
    # run stretcher
    cmd = "stretcher " + "ltr5.fa ltr3.fa -outfile ltrIdent.txt"
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    process.wait()
    print(process.returncode)
    # parse ltr identity
    identReg = r' \(\s?(\d+\.\d+)%'
    ident = 0
    avgLtrLen = (ltr_len[0] + ltr_len[1]) / 2
    with open("ltrIdent.txt") as identFile:
        for l in identFile:
            if "# Identity:" in l:
                print(l)
                ident = float(re.search(identReg, l.rstrip()).group(1))
    files2remove.append("ltrIdent.txt")

    # get Kimura distance:
    ## run clustalw to get .phy file
    k80 = ""
    k2p_g = "NA"
    k2p_l = "NA"
    if not os.path.exists("outfile"):
        with open("outfile", "w") as of:
            of.write("")
    if ltr_len[0] > 0 and ltr_len[1] > 0:
        # merge two fa files:
        seq_List = []
        for r in SeqIO.parse("ltr5.fa","fasta"):
            seq_List.append(str(r.seq))
        for s in SeqIO.parse("ltr3.fa","fasta"):
            seq_List.append(str(s.seq))
        with open("infile.fa", "w") as infOut:
            infOut.write(">ltr5\n")
            infOut.write(seq_List[0] + "\n")
            infOut.write(">ltr3\n")
            infOut.write(seq_List[1] + "\n")                
        cmd1 = "clustalw infile.fa -output=PHYLIP"
        process = subprocess.Popen(cmd1, shell=True, stdout=subprocess.PIPE)
        process.wait()
        print(process.returncode)

        os.rename('infile.phy', 'infile')
        with open("opt_dnadist", "w") as out:
            out.write("infile\nR\nD\nY\n")

        cmd1 = "cat opt_dnadist | phylip dnadist"
        process = subprocess.Popen(cmd1, shell=True, stdout=subprocess.PIPE)
        process.wait()
        print(process.returncode)

        with open("outfile") as outfile:
            next(outfile)
            for l in outfile:
                lList = l.split()
                k80 = lList[2]  # return K80
                break
                files2remove.append("infile.opt_dnadist")
        files2remove.append("infile")
        files2remove.append("infile.fa")
        files2remove.append("infile.dnd")
        files2remove.append("outfile")
        files2remove.append("infile")
        # get Kimura 2-parameters distance:
        ## get MAFFT alignment    
        cmd_mg1 = "mafft --globalpair --maxiterate 1000 infile.fa > mafft_global.txt"
        process = subprocess.Popen(cmd_mg1, shell=True, stdout=subprocess.PIPE)
        process.wait()
        print(process.returncode)
        
        with open("opt_distmat_g", "w") as out:
            out.write("mafft_global.txt\n2\n\n")
        cmd_mg2 = "cat opt_distmat_g | distmat"
        process = subprocess.Popen(cmd_mg2, shell=True, stdout=subprocess.PIPE)
        process.wait()
        print(process.returncode)
        if div_from_distmat("mafft_global.distmat"):
            k2p_g = round(div_from_distmat("mafft_global.distmat") / 100,4)

        cmd_ml1 = "mafft --localpair --maxiterate 1000 infile.fa > mafft_local.txt"
        process = subprocess.Popen(cmd_ml1, shell=True, stdout=subprocess.PIPE)
        process.wait()
        print(process.returncode)
        
        with open("opt_distmat_l", "w") as out:
            out.write("mafft_local.txt\n2\n\n")
        cmd_ml2 = "cat opt_distmat_l | distmat"
        process = subprocess.Popen(cmd_ml2, shell=True, stdout=subprocess.PIPE)
        process.wait()
        print(process.returncode)

        if div_from_distmat("mafft_local.distmat"):
            k2p_l = round(div_from_distmat("mafft_local.distmat") / 100, 4)

        files2remove += ["mafft_global.txt","mafft_local.txt","opt_distmat_g","opt_distmat_l","mafft_global.distmat","mafft_local.distmat"]
                
    if not k80:
        k80 = "NA"
        
    # remove redundant files
    files2remove+=['opt_dnadist','outfile']
# remove redundant files
    for file in files2remove:
        if os.path.exists(file):
            os.remove(file)                
    return [str(ident), str(avgLtrLen), str(ltr_len[0]), str(ltr_len[1]), str(k80.rstrip()), str(k2p_l), str(k2p_g)]
    #return [str(k2p_l), str(k2p_g)]

def main():
    # collect all repeat_region tracks
    rep_regs = rep_region_list()

    # 1. get both LTRs tracks for each repeat region
    # 2. gen fasta sequences -> get their identity and K80 divergence
    out_name = "Slat_v4_FL_annot_ltrDiv_ident_k2p_k80.tab"
    with open(out_name,'w') as out:
        for rr in rep_regs:
            rr_id = rr.attributes['ID']
            ltrs_tracks = find_ltrs(rr.seqid, rr_id)
            ltr_len = get_ltr_fa_pair(ltrs_tracks)
            fam = rr.attributes['dante_annotation'].split("|")[2].split("/")[1] + "_" + rr.attributes['dante_annotation'].split("|")[-1]
            ltr_div_ident = "\t".join(ltr_div(ltr_len))
            out.write(f"{rr.seqid}\t{rr.start}\t{rr.end}\t{fam}\t{rr.attributes['ltr_identity']}\t{ltr_div_ident}\n")
            
if __name__ == "__main__":
    main()
