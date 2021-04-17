#!/usr/bin/env python

import os
import sys
import argparse
import subprocess
import assembly
import annotate
import multiprocessing

thread=str(multiprocessing.cpu_count())

def main():
	parser = argparse.ArgumentParser(description="Assembles and annotates genome given a set of proteins.")
	parser.add_argument("-s",
						type=str,
						required=True,
						metavar="<Nanopore file>",
						help="Single read of genome sequencing"
						)
	parser.add_argument("-gs",
						type=str,
						required=True,
						metavar="<estimated genome size e.g 5m>",
						help="Approximate genome size, e.g 5m, 4.6g"
						)
	parser.add_argument("-l",
						type=str,
						required=True,
						choices=['bacteria','archaea','fungi'],
						metavar="<Lineage e.g: bacteria, archaea>",
						help="Kingdom, acceptable input: archaea, bacteria, fungi"
						)
	parser.add_argument("-db",
						type=str,
						required=True,
						metavar="<Swissprot database in fasta>",
						help="Protein fasta to be used in blast and annotation"
						)
	args = parser.parse_args()
	
	assembly.assemble_genome(
		args.s,
		args.gs,
		thread
	)
	
	assembly.run_busco(
		args.s,
		args.l
	)
	
	compare=assembly.busco_result(
		args.s,
		args.l
	)

	with open("Assembly_chosen.txt","w") as file:
		file.write("Choosing assembly from "+compare)
	
	assembly.repeatmask(thread)
	
	assembly.run_braker(args.l,thread)
	
	annotate.annotate_proteins(args.db,thread)
	annotate.cleanup(args.s)

if __name__ == "__main__": #script is being executed, not imported
	main()