#!/usr/bin/env python
#coding:utf-8
"""
Translate RNA in 3-frames to ORFs
"""
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

if __name__ == '__main__':
    if len(sys.argv) != 3:
        sys.stderr.write('Usage: program input.fa output.fa\n')
        sys.exit(1)

    infile, outfile = sys.argv[1:3]
    output_handle = open(outfile, 'w')

    numer,count=0,0
    for r in SeqIO.parse(open(infile), "fasta"):
        numer+=1
        cdna = Seq(str(r.seq), IUPAC.unambiguous_dna)
        rna = cdna.upper().transcribe()
        print numer
        for frame in xrange(3):
            newrna = rna[frame:]
            newrna = newrna[:int(len(newrna)/3)*3]
            end=frame
            try:
                cds=newrna.translate()
                block_list = cds.split('*')
                for block in block_list[:-1]:
                    start=end
                    end+=len(block)*3+3
                    Pm = block.find('M')
                    start += Pm*3
                    block = block[Pm:]
                    if len(block) >= 20: 
                        count+=1
                        proid = r.id + ":" + str(start) + '-' + str(end)
                        proseq = block
                        SeqIO.write([SeqRecord(proseq, id=proid, description='')],output_handle, "fasta")
            except: 
                   break   
    print 'total orfs:%s'%count  
    output_handle.close()

