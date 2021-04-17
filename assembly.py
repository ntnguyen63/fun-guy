import glob, os, shutil, sys
import subprocess
import multiprocessing

##TODO check what happen if seqfile is not in the same location of cmd eg seqfile = ./test/final/mushroom/oxford.fasta
##TODO maybe convert ALL relpath to os.getcwd for more compatibility?
##TODO


def run_flye(seqfile,gsize,thread): #run flye assembler and move the result folder to <seqfile>_outdir on the corrected reads from canu
	if os.path.exists('flye_out'):
		shutil.rmtree('flye_out')
	path=os.getcwd()
	seqcorr=seqfile+".correctedReads.fasta.gz"
	print(f"Running flye on {seqcorr}, output to assembly.fasta in dir flye_out")
	command=["flye",
			"--nano-corr",seqcorr,
			"-o","flye_out",
			"-g",gsize,
			"-t",thread]
	subprocess.run(command)
	shutil.move(path+'/flye_out',path+'/'+seqfile+'_outdir')
	
def run_canu(seqfile,gsize): #run canu assembler and move the result folder to <seqfile>_outdir
	subprocess.run(["canu",
		"-correct",
		"-p",seqfile,
		"-d","corrected_canu", ##output to corrected_canu/ dir because canu can output other files not just the corrected reads, this is cleaner
		"genomeSize="+gsize,
		"-nanopore",
		seqfile]				#Correct the raw data, output corrected reads to <seqfile>.correctedReads.fasta.gz
	)
	
	subprocess.run(["mv",
					"./corrected_canu/"+seqfile+".correctedReads.fasta.gz",
					"./"]) #move it to working directory
	
	if os.path.exists('corrected_canu'):
		shutil.rmtree('corrected_canu')
	if os.path.exists('canu_out'):
		shutil.rmtree('canu_out')
	path=os.getcwd()
	seqcorr=seqfile+".correctedReads.fasta.gz"
	print(f"Running canu on {seqcorr}, output to {seqfile}.contigs.fasta in dir canu_out")
	#command=["canu","-p",seqfile,"-d","canu_out","genomeSize="+gsize,"-nanopore","-corrected",seqcorr]
	subprocess.run(["canu",
					"-p",seqfile,
					"-d","canu_out",
					"genomeSize="+gsize,
					"-nanopore","-corrected",seqcorr]
	)
	shutil.move(path+'/canu_out',path+'/'+seqfile+'_outdir')
	
def run_busco_canu(seqfile,lineage): #run busco and move result to <seqfile>_outdir
	path=os.getcwd()
	if os.path.exists('busco_out_canu'):
		shutil.rmtree('busco_out_canu')
	subprocess.run(["busco",
					"-i","./"+seqfile+"_outdir/"+"canu_out/"+seqfile+".contigs.fasta",
					"-o","busco_out_canu","-m","genome",
					"-l","./busco_downloads/"+lineage+"_odb10",
					"-c","2"]
	)
	shutil.move(path+'/busco_out_canu',path+'/'+seqfile+'_outdir')
	
def run_busco_flye(seqfile,lineage): #run busco and move result to <seqfile>_outdir
	path=os.getcwd()
	if os.path.exists('busco_out_flye'):
		shutil.rmtree('busco_out_flye')
	subprocess.run(["busco",
		"-i","./"+seqfile+"_outdir/flye_out/assembly.fasta",
		"-o","busco_out_flye","-m","genome",
		"-l","./busco_downloads/"+lineage+"_odb10",
		"-c","2"]
	)
	shutil.move(path+'/busco_out_flye',path+'/'+seqfile+'_outdir')
	
def download_busco(lineage): #Make a folder {lineage}_odb10 with data for busco
	if os.path.exists('busco_downloads')==False:
		os.mkdir("busco_downloads")
	if lineage=="fungi":
		subprocess.run(["wget","-nc","-O","./busco_downloads/"+lineage+"_odb10.tar.gz",
						"https://busco-data.ezlab.org/v5/data/lineages/fungi_odb10.2020-09-10.tar.gz"])
	elif lineage=="eukaryota":
		subprocess.run(["wget","-nc","-O","./busco_downloads/"+lineage+"_odb10.tar.gz",
						"https://busco-data.ezlab.org/v5/data/lineages/eukaryota_odb10.2020-09-10.tar.gz"])
	elif lineage=="bacteria":
		subprocess.run(["wget","-nc","-O","./busco_downloads/"+lineage+"_odb10.tar.gz",
						"https://busco-data.ezlab.org/v5/data/lineages/bacteria_odb10.2020-03-06.tar.gz"])
	else:
		subprocess.run(["wget","-nc","-O","./busco_downloads/"+lineage+"_odb10.tar.gz",
						"https://busco-data.ezlab.org/v5/data/lineages/archaea_odb10.2021-02-23.tar.gz"])
	subprocess.run(["tar","-xzf","./busco_downloads/"+lineage+"_odb10.tar.gz","-C","./busco_downloads/"])
	#Leave busco_downloads dir outside for subsequent runs
	
def download_db(lineage): #produce {lineage}_proteins.fasta for prothint. Can be predownloaded and leave in working dir to save time
	if lineage=="archaea" and os.path.exists('archeae_proteins.fasta')==False: 	###Database link mispelled archaea as archea
		subprocess.run(["wget","-nc","https://www.orthodb.org/v9.1/download/odb9v1_archea_fasta.tar.gz"])
		subprocess.run(["tar","-xzf","odb9v1_archea_fasta.tar.gz"])
		subprocess.run(["cat","archea/Rawdata/*>"+lineage+"_proteins.fasta"])
		subprocess.run(["rm","./odb9v1_archea_fasta.tar.gz"])
		shutil.rmtree('archea')
	elif os.path.exists(lineage+"_proteins.fasta")==False:
		subprocess.run(["wget","-nc","https://v100.orthodb.org/download/odb10_"+lineage+"_fasta.tar.gz"])
		subprocess.run(["tar","-xzf","odb10_"+lineage+"fasta.tar.gz"])
		subprocess.run(["cat",lineage+"/Rawdata/*>"+lineage+"_proteins.fasta"])
		subprocess.run(["rm","-r","odb10_"+lineage+"fasta.tar.gz"])
		shutil.rmtree(lineage)

	
def run_busco(seqfile,lineage):  #parallelize 
	p0= multiprocessing.Process(target=download_db, args=(lineage, ))
	p1= multiprocessing.Process(target=download_busco, args=(lineage, ))
	p2= multiprocessing.Process(target=run_busco_canu, args=(seqfile,lineage, ))
	p3= multiprocessing.Process(target=run_busco_flye, args=(seqfile,lineage, ))
	p0.start()
	p1.start()
	p1.join() #Wait for busco download to finish
	p2.start()
	p3.start()
	p0.join()
	p3.join()
	p2.join()

def assemble_genome(seqfile,gsize,thread): #check existing output folder and ask if we can clear it. ##TODO maybe implement custom output folder?
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
			run_canu(seqfile,gsize) #run canu first for the corrected reads
			run_flye(seqfile,gsize,thread)
		elif val == 'n':
			sys.exit('Exiting')
		else:
			print("Unrecognize input, exiting function")
			sys.exit()
	else:
		run_canu(seqfile,gsize)
		run_flye(seqfile,gsize,thread)
		return
	
def busco_result(seqfile,lineage):
#Parse the result file from busco then return total complete busco. Then move the corresponding assembly to working dir, rename it to main_assembly.fasta
	canupath=os.path.relpath(seqfile+"_outdir/busco_out_canu/short_summary.specific."+lineage+"_odb10.busco_out_canu.txt")
	flyepath=os.path.relpath(seqfile+"_outdir/busco_out_flye/short_summary.specific."+lineage+"_odb10.busco_out_flye.txt")
	with open(canupath,"r") as file:
		lines=file.readlines()
		canu_score=float(lines[8].split('%')[0].split(':')[1]) #extract the total complete score as float
	with open(flyepath,"r") as file:
		lines=file.readlines()
		flye_score=float(lines[8].split('%')[0].split(':')[1])
	if canu_score>flye_score:
		compare="canu"
		subprocess.run(["cp","./"+seqfile+"_outdir/"+"canu_out/"+seqfile+".contigs.fasta","."])
		subprocess.run(["mv","./"+seqfile+".contigs.fasta","main_assembly.fasta"])
	elif canu_score<flye_score:
		compare='flye'
		subprocess.run(["cp","./"+seqfile+"_outdir/flye_out/assembly.fasta","."])
		subprocess.run(["mv","assembly.fasta","main_assembly.fasta"])
	else:
		compare='flye'
		subprocess.run(["cp","./"+seqfile+"_outdir/flye_out/assembly.fasta","."])
		subprocess.run(["mv","assembly.fasta","main_assembly.fasta"])
	return compare
		
def repeatmask(thread):
	subprocess.run(["BuildDatabase","-name","temporary_db","./main_assembly.fasta"]) #produce temporary_db.* files as database
	subprocess.run(["RepeatModeler","-database","temporary_db","-pa",thread]) #produce temporary-families.fa and .stk and RM_* folders
	
	subprocess.run(["RepeatMasker",
					"-lib","temporary_db-families.fa", #use *families.fa for repeatmasker
					"main_assembly.fasta",
					"-dir","repeatmasker_out",
					"-xsmall"]
					)
	subprocess.run(["rm","-rf","RM_?????.*"])
	subprocess.run(["rm","-f","temporary_db.*","temporary_db-families.stk"])
	
##need to move repeatmasker_out to {seqfile}_outdir later

def run_braker(lineage,thread): #produce prothint_out and braker directories
	subprocess.run(["prothint.py",
					"main_assembly.fasta",
					lineage+"_proteins.fasta",
					"--workdir","prothint_out"]
					)
	command=[	"braker.pl",
				"--genome=./repeatmasker_out/main_assembly.fasta.masked",
				"--hints=./prothint_out/prothint_augustus.gff", #hint file is generated from running prothint with orthodb fungi v10
				"--softmasking","--cores="+thread
			]
	if lineage=="fungi":
		command.append("--fungus")
	subprocess.run(command)
	#need to move the 2 dir to {seqfile}_outdir