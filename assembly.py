import glob, os, shutil, sys
import subprocess
import multiprocessing

def run_flye(seqfile,gsize): #run flye assembler and move the result folder to <seqfile>_outdir
	if os.path.exists('flye_out'):
		shutil.rmtree('flye_out')
	path=os.getcwd()
	thread=str(multiprocessing.cpu_count())
	seqcorr=seqfile+".correctedReads.fasta.gz"
	print(f"Running flye on {seqcorr}, output to assembly.fasta in dir canu_out")
	command=["flye","--nano-corr",seqcorr,"-o","flye_out","-g",gsize,"-t",thread]
	subprocess.run(command)
	shutil.move(path+'/flye_out',path+'/'+seqfile+'_outdir')
	
def run_canu(seqfile,gsize): #run canu assembler and move the result folder to <seqfile>_outdir
	subprocess.run(["canu",  #Generate corrected reads
		"-correct",
		"-p",seqfile,
		"-d","corrected_canu",
		"genomeSize="+gsize,
		"-nanopore",
		seqfile]
	)
	subprocess.run(["mv","./corrected_canu/"+seqfile+".correctedReads.fasta.gz","./"])
	if os.path.exists('corrected_canu'):
		shutil.rmtree('corrected_canu')
	if os.path.exists('canu_out'):
		shutil.rmtree('canu_out')
	path=os.getcwd()
	seqcorr=seqfile+".correctedReads.fasta.gz"
	print(f"Running canu on {seqcorr}, output to {seqfile}.contigs.fasta in dir canu_out")
	command=["canu","-p",seqfile,"-d","canu_out","genomeSize="+gsize,"-nanopore","-corrected",seqcorr]
	subprocess.run(["canu","-p",seqfile,"-d","canu_out","genomeSize="+gsize,"-nanopore","-corrected",seqcorr])
	shutil.move(path+'/canu_out',path+'/'+seqfile+'_outdir')
	
def run_busco_canu(seqfile,lineage): #run busco and move result to <seqfile>_outdir
	path=os.getcwd()
	if os.path.exists('busco_out_canu'):
		shutil.rmtree('busco_out_canu')
	subprocess.run(["busco",  #busco on canu
		"-i","./"+seqfile+"_outdir/"+"canu_out/"+seqfile+".contigs.fasta",
		"-o","busco_out_canu","-m","genome",
		"-l",lineage+"_odb10"]
	)
	shutil.move(path+'/busco_out_canu',path+'/'+seqfile+'_outdir')
	
def run_busco_flye(seqfile,lineage): #run busco and move result to <seqfile>_outdir
	path=os.getcwd()
	if os.path.exists('busco_out_flye'):
		shutil.rmtree('busco_out_flye')
	subprocess.run(["busco",  #busco on flye
		"-i","./"+seqfile+"_outdir/"+"flye_out/assembly.fasta",
		"-o","busco_out_flye","-m","genome",
		"-l",lineage+"_odb10"]
	)
	shutil.move(path+'/busco_out_flye',path+'/'+seqfile+'_outdir')
	
def run_busco(seqfile,lineage):  #parallelize busco
	p1= multiprocessing.Process(target=run_busco_canu, args=(seqfile,lineage, ))
	p2= multiprocessing.Process(target=run_busco_flye, args=(seqfile,lineage, ))
	
	p1.start()
	p2.start()
	p1.join()
	p2.join()

def assemble_genome(seqfile,gsize):
	outdir=(f"{seqfile}_outdir")
	path=os.getcwd()
	try:
		os.mkdir(outdir)
	except FileExistsError as error:
		print("Folder {} already exists.".format('outdir'))
		val=input("Clear out content(y/n)? ")
		if val == 'y':
			shutil.rmtree(outdir)
			os.mkdir(outdir)
			run_canu(seqfile,gsize)
			run_flye(seqfile,gsize)
		elif val == 'n':
			sys.exit('Exiting')
		else:
			print("Unrecognize input, exiting function")
			sys.exit()
	else:
		run_canu(seqfile,gsize)
		run_flye(seqfile,gsize)
		return
	
def busco_result(seqfile,lineage):
#Parse the result file from busco using input of kingdom type then return total complete busco
	canupath=os.path.relpath(seqfile+"_outdir/busco_out_canu/short_summary.specific."+lineage+"_odb10.busco_out_canu.txt")
	flyepath=os.path.relpath(seqfile+"_outdir/busco_out_flye/short_summary.specific."+lineage+"_odb10.busco_out_flye.txt")
	with open(canupath,"r") as file:
		lines=file.readlines()
		canu_score=lines[8].split('%')[0].split(':')[1] #extract the total complete score 
	with open(flyepath,"r") as file:
		lines=file.readlines()
		flye_score=lines[8].split('%')[0].split(':')[1]
	if canu_score>flye_score:
		return 'canu'
	elif canu_score<flye_score:
		return 'flye'
	else:
		return 'flye'