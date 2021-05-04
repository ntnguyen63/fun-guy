# fun-guy

Program to assemble fungal genome from nanopore read using canu and flye. Then using busco to compare the assemblies and proceed to predict and annotate the proteins frm the better assembly.

Usage: funguy.py -l <lineage> -s <nanopore read> -gs <approx genome size> -db <fasta to build blast database for annotation>

May also work with archaeal and bacterial genome but functionality for those were not tested.

Dependencies:
Canu, Flye, Busco v5.1.2, braker2, RepeatMasker, RepeatModeler, blast, GenemarkES, Prothint, Augustus

Most dependencies can be setup with miniconda

  conda install -c bioconda -c conda-forge -c thiesgehrmann busco=5.1.2 canu flye blast repeatmodeler repeatmasker braker

repeatmodeler and repeatmasker need to be setup for the type of organism whose genome is being assembled

conda will not fully setup prothint, augustus, and braker for you. These will need to be installed manually. You need to make sure the path for PERL5LIB is pointed to your perl installation where you have the dependencies, and GENEMARK_PATH points to a working directory of genemark
