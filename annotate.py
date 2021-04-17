#!/usr/bin/env python

import subprocess
import os

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Blast import NCBIXML

def get_predicted_protein():
	#Use getAnnoFastaFromJoingenes.py to extract prot sequences from main_assembly.fasta
	subprocess.run(["getAnnoFastaFromJoingenes.py",
					"-g","main_assembly.fasta",
					"-f","./braker/braker.gtf",
					"-o","main_assembly_prot.fa"
					]
	)

def makeblastdb(dbseqs):
	'''Create blast database with makeblastdb command from blast+.

	Arguments:
	dbseqs: name of uniprot file to build database.

	Returns: makeblastdb output files in the current directory
	'''
	command = [
		'makeblastdb','-hash_index','-in',dbseqs,'-dbtype','prot',
		'-out','blast.db',"-parse_seqids"
	]

	subprocess.run(command)

def run_blastp(thread):
	'''Run blastp command. The query is predicted proteins or other proteins and the database
	is the refseq database, unless a different database is supplied.

	Argument:
	query: protein set to search against blast database
	dbname: blast database name

	Returns:
	XML file of blast results.
	'''
	command = [
		"blastp","-db","blast.db","-query","main_assembly_prot.fa.aa",
		'-evalue','1e-16','-max_target_seqs','2',
		'-out',"./blastp_result.xml","-outfmt","5",
		'-num_threads',thread
	]

	print("Running blastp..")
	subprocess.run(command)
	print("Done")

def seq_lookup_table(fasta_file):  #Extract name and sequence from fasta, ie predicted protein from braker
	'''
	'''
	lookup_table = {}
	for record in SeqIO.parse(fasta_file,"fasta"):
		lookup_table[record.id] = record.seq
	return lookup_table

def go_through(blast_record):  #Take the first blast match which is supposed to best and extract the protein description from title, separated by OS and a space " "

	prot_functions = []
	for alignment in blast_record.alignments:
		title = alignment.title 
		index = title.find("sp") if "sp" in title else 0 #Start reading title from sp, begin extracting at first space " " and end extraction at "OS" e.g: sp|P39707|SEN34_YEAST <start>tRNA-splicing endonuclease subunit SEN34 <end>OS=Saccharomyces cerevisiae
		prot_function = title[title.find(" ",index):title.find('OS=')] #OS= in case hitid has OS like OSH1
		prot_functions.append(prot_function) #list of lists, each piece is a list of all proteins functions matched with one query
	return prot_functions

def hits_from_blast_results(result_file): #

	with open(result_file) as blast_file:
		blast_records = NCBIXML.parse(blast_file)
		hits = {} #hits is a dictionary, each element is mapped to another value
		for blast_record in blast_records:
			query = blast_record.query				#get each query in blast_records
			#query = query[:query.find("#")].strip(" ")
			protein_functions = go_through(blast_record)
			if protein_functions:
				hits[query] = protein_functions[0] #attach first [0] (best match) protein function to the query name and put them in hits dictionary
	return hits

def label_proteins(predicted_proteins_file,blast_result_file,outfile):
	lookup_table = seq_lookup_table(predicted_proteins_file)
	hits = hits_from_blast_results(blast_result_file)
	annotations = []
	for i,item in enumerate(hits.items()): 	#enumerate: make loop counter through all item in list
											#items(): convert dictionary to list of lists, each element is one object and its matched value
		fasta_id,predicted_function = item
		seq = lookup_table[fasta_id]
		new_record = SeqRecord(	id="fungy{}".format(i),
								description="fungy{} {}".format(i,predicted_function),
								seq=seq
		)
		annotations.append(new_record)
	SeqIO.write(annotations,outfile,"fasta")

def annotate_proteins(dbseqs,thread):
	get_predicted_protein()
	makeblastdb(dbseqs)
	run_blastp(thread)
	label_proteins(	"main_assembly_prot.fa.aa",
					"blastp_result.xml",
					"./annotated_prot.faa"
	)

def cleanup(seqfile):
	if not os.path.isdir("blast_out"):
		os.mkdir("blast_out")
	subprocess.run(["mv","./blastp_result.xml","./blast_out/"])
	subprocess.run(["mv","./blast_out","./"+seqfile+"_outdir/"])
	subprocess.run(["mv","./main_assembly_prot.fa.aa","./braker/"])
	subprocess.run(["mv","./main_assembly_prot.fa.codingseq","./braker/"])
	subprocess.run(["mv","./braker","./"+seqfile+"_outdir/"])
	subprocess.run(["mv","./main_assembly.fasta","./"+seqfile+"_outdir/"])
	subprocess.run(["mv",seqfile+".correctedReads.fasta.gz","./"+seqfile+"_outdir/"])
	subprocess.run(["mv","./prothint_out/prothint_augustus.gff","./"+seqfile+"_outdir/"])
	subprocess.run(["rm","-r","./prothint_out"])
	subprocess.run(["mv","./annotated_prot.faa","./"+seqfile+"_outdir/"])
	subprocess.run(["mv","./repeatmasker_out","./"+seqfile+"_outdir/"])
	subprocess.run(["mv","./Assembly_chosen.txt","./"+seqfile+"_outdir/"])
	subprocess.run(["rm","blast.db.*","gmes.log","run.cfg"])
