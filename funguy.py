#!/usr/bin/env python

import os
import sys
import argparse
import subprocess
import assembly


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
	args = parser.parse_args()
	
	assembly.assemble_genome(
		args.s,
		args.gs
	)
	
	assembly.run_busco(
		args.s,
		args.l
	)
	
	compare=assembly.busco_result(
		args.s,
		args.l
	)

	print(f"Choosing result from {compare} assembler")


if __name__ == "__main__": #script is being executed, not imported
	main()