import subprocess
from pathlib import Path
import shutil
import difflib

OLD_JSCB_CMD = "python3 ~/HGT/jscb/JSCB/run_jscb.py --genbank {infile} -o ~/HGT/jscb/JSCB/"
NEW_JSCB_CMD = "python3 src/jscb_new.py {infile} -o outputs"

genbanks = list(Path("genbanks").glob('*.gbk'))

OLD_JSCB_DIR = Path("~/HGT/jscb/JSCB/").expanduser()

gi_diff = open("res/gi_diff.tsv", 'w')
hgt_diff = open("res/hgt_diff.tsv", 'w')

for gbk in genbanks:
    subprocess.run(OLD_JSCB_CMD.format(gbk), shell=True, check=True)

    
    subprocess.run(NEW_JSCB_CMD.format(gbk), shell=True, check=True)
    
    with open(OLD_JSCB_DIR / "genomic_islands_summary.tsv") as old, open("outputs/genomic_islands_summary.tsv") as new:
        old_lines = old.readlines()
        new_lines = new.readlines()
        
    diff = difflib.unified_diff(
        old_lines,
        new_lines,
        fromfile= "OLD",
        tofile= "NEW"
    )
    gi_diff.writelines(diff)
    
    with open(OLD_JSCB_DIR / "hgt_genes_summary.tsv") as old, open("outputs/hgt_genes_summary.tsv") as new:
        old_lines = old.readlines()
        new_lines = new.readlines()
        
    diff = difflib.unified_diff(
        old_lines,
        new_lines,
        fromfile= "OLD",
        tofile= "NEW"
    )
    hgt_diff.writelines(diff)
    
    shutil.rmtree("outputs")
    
    (OLD_JSCB_DIR / "combined.gbk").unlink()
    (OLD_JSCB_DIR / "JSCB_output.gi").unlink()
    (OLD_JSCB_DIR / "JSCB_output.tsv").unlink()
    (OLD_JSCB_DIR / "hgt_genes_summary.tsv").unlink()
    (OLD_JSCB_DIR / "genomic_islands_summary.tsv").unlink()

gi_diff.close()
hgt_diff.close()

    
    
