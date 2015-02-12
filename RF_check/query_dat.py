import Bio
import sqlite3 as sql
from Bio import SeqIO
from Bio.Blast import NCBIWWW, NCBIXML



dv_dat = SeqIO.to_dict(SeqIO.parse(open('dv4_36.fasta', 'rU'), 'fasta'))

blast_desc = list()
blast_lens = list()
blast_e = list()

resdb = sql.connect('out_data.db')
cursor = resdb.cursor()
cursor.execute('CREATE TABLE IF NOT EXISTS denv4(name TEXT UNIQUE, desc TEXT, aln_len INTEGER, eval REAL)')


for i in dv_dat.keys():
    print 'querying NCBI for %s' %i
    temp_res = None
    max_iter = 10
    while (temp_res is None) and max_iter > 1: 
        try:
            temp_res = NCBIWWW.qblast('blastn', 'nt', dv_dat[i].seq, alignments = 1, descriptions = 1, megablast = True)

            print 'downloaded data for %s' %i
            temp_proc = NCBIXML.read(temp_res)
#            blast_desc.append(temp_proc.alignments[0].title)
#            blast_lens.append(temp_proc.alignments[0].length)
#            blast_e.append(temp_proc.expect)
            cursor.execute('INSERT INTO denv4 VALUES(?, ?, ?, ?)', (i, temp_proc.alignments[0].title, temp_proc.alignments[0].length, temp_proc.expect))
            resdb.commit()
        except:
            print 'Error. %s not found' %i
            max_iter -= 1
            continue

resdb.close()
