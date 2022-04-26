**Streaming workflow for evolutionary and comparative genomics**

**Contributors**

Qingxiang Guo

**About**

The following contents contains source code for the analyses and plots in BMC Biology paper, 2022, “A myxozoan genome reveals mosaic evolution in a parasitic cnidarian”

**Abstract**

Parasite evolution has been conceptualized as a process of genetic loss and simplification. Contrary to this model, there is evidence of expansion and conservation of gene families related to essential functions of parasitism in some parasite genomes, reminiscent of widespread mosaic evolution where subregions of a genome have different rates of evolutionary change. We found evidence of mosaic genome evolution in the cnidarian *Myxobolus honghuensis*, a myxozoan parasite of fish, with extremely simple morphology.

**Notes**

Written by Qingxiang Guo, qingxiang.guo@outlook.com, distributed without any guarantees or restrictions.

**Codes**

**1. Use Blobtools to remove the contamination in myxozoan genome**

**1.1 Filter the genome to 200 bp**

filter\_fasta\_by\_length.pl genome.fasta 200 200000 genome\_200.fasta

**1.2 Get the coverage info**

\# Build a index for the genome

nohup bowtie2-build M.honghuensis\_200\_scaf.fasta index --threads 8 &

\# Mapping

nohup bowtie2 -p 24 -x index -1 M.honghuensis\_1.fq -2 M.honghuensis\_2.fq -k 1 --very-fast-local -S out.sam &

\# Convert the style

samtools view -bS out.sam > out.bam

**1.3 Get the megablast results**

nohup blastn -task megablast -query M.honghuensis\_200\_scaf.fasta -db nt -culling\_limit 5 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle staxids sscinames sskingdoms' -num\_threads 48 -evalue 1e-25 -out assembly\_megablast\_25.out &

\# Then add species information to the last line using the following method

ln -s /home/train/public\_database/accession2taxid/acc2tax\_nucl\_all.txt ./

nohup perl -lne '

`  `BEGIN{open UT, "<acc2tax\_nucl\_all.txt" or die $!; while (<UT>) { $ut{$2}=$3 if /^(\S+)\t(\S+)\t(\S+)/ } }

`  `{print "$\_\t$ut{$1}" if /^\S+\t(\S+)/ and exists $ut{$1}}' \

< assembly\_megablast\_25.out \

\> assembly\_megablast\_25.out\_taxid &

\# Change the output style

awk -v OFS="\t" -F"\t" '{print $1,$17,$12}' assembly\_megablast\_25.out\_taxid  > mega.out

**1.4 Get the diamond-blastx results**

\# Make a index

diamond makedb --in uniref90.aa -d uniref90  

nohup diamond blastx -q ../../2\_assembly/M.honghuensis\_200\_scaf.fasta --sensitive -k 20 -c 1 --evalue 1e-10 --threads 48 --db /home/train/public\_database/uniref90/uniref90.dmnd --out assembly\_diamond\_10.out &

\# Add species information to the last line

ln -s /home/train/public\_database/uniref90/uniref90.taxlist ./

perl -lne '

`  `BEGIN{open UT, "<uniref90.taxlist" or die $!; while (<UT>) { $ut{$1}=$2 if /^(\S+)\t(\S+)$/ } }

`  `{print "$\_\t$ut{$1}" if /^\S+\t(\S+)/ and exists $ut{$1}}' \

< assembly\_diamond\_10.out \

\> assembly\_diamond\_10.out\_taxid

\# Change the ouput style

awk -v OFS="\t" -F"\t" '{print $1,$13,$12}' assembly\_diamond\_10.out\_taxid  > diamond.out

**1.5 A database file was then generated based on the genome sequence, blast results, maps results**

\# File preparation

cd /opt/biosoft/blobtools/M.honghuensis.plot\_1

ln -s /home/gqx/genome/M.honghuensis/3\_contam\_remove/diamond/diamond.out ./

ln -s /home/gqx/genome/M.honghuensis/3\_contam\_remove/megablast/mega.out ./

ln -s ~/genome/M.honghuensis/2\_1\_assembly/M.honghuensis\_200\_scaf.fasta ./

ln -s /home/gqx/genome/M.honghuensis/3\_contam\_remove/mapping\_results/out.bam ./

\# Merging and ordering blast results

cat mega.out diamond.out > blast.out

sort\_blast\_by\_query\_name.pl blast.out

mv sorted\_output blast.out

python2.7 ../blobtools create -i M.honghuensis\_200\_scaf.fasta -b out.bam -t blast.out -o M.honghuensis\_1\_blob --names ../names.dmp --nodes ../nodes.dmp

**1.6 Statistics performed on the database**

python2.7 ../blobtools view -i M.honghuensis\_1\_blob.blobDB.json -o ./

grep '^##' M.honghuensis\_1\_blob.blobDB.table.txt ; 

grep -v '^##' M.honghuensis\_1\_blob.blobDB.table.txt | \

column -t -s $'\t'

**1.7 Visualization**

python2.7 ../blobtools blobplot -i M.wulii\_1\_blob.blobDB.json -o ./ --format pdf --colours colours.txt

\# colour.txt contains the file with the assigned color in the following format:

no-hit,#c6c6c6

other,#ffffff

Nematoda,#dd0000

Arthropoda,#3c78d8

Cnidaria,# 6aa84f

Chordata,#fdfd99

Firmicutes,#fca55d

Proteobacteria,#47a0b3

**1.8 Contamination was removed according to blobtools visualization**

\# Extract the abnormal sequences (score > 200) from the M.honghuensis\_1\_blob.blobDB.table.txt file. Blast it with nt database and remove those aligned to bacteria and host database

\# Adjust the output format

\# format.sh

#!/bin/bash

grep '^##' M.honghuensis\_1\_blob.blobDB.table.txt;\

grep -v '^##' M.honghuensis\_1\_blob.blobDB.table.txt | \

column -t -s $'\t'

./format.sh > result

**1.9 Processing the output file and get the blast results**

blob\_result\_seq\_extract.pl M.honghuensis\_1\_blob.blobDB.table.txt

\# Extract the potential contamination reads

extract\_seq\_from\_fasta.pl M.honghuensis\_200\_scaf.fasta seq\_for\_blast

mv extracted.fasta contam\_candidate.fa

**1.10** **NT alignment of suspected contaminant sequences**

mkdir -p /home/gqx/genome/M.wulii/3\_contam\_remove/contam\_blast

cd /home/gqx/genome/M.wulii/3\_contam\_remove/contam\_blast

cp /home/gqx/genome/M.wulii/3\_contam\_remove/M.wulii.plot\_1/contam\_candidate.fa ./

nohup blastn -query contam\_candidate.fa -db /home/train/public\_database/nt\_160828/nt -evalue 1e-5 -max\_target\_seqs 20 -num\_threads 24 -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle frames sscinames sskingdoms" -out nt\_result &

**1.11 Process alignments and add species information**

\# Extract the best hits

export LANG=C; export LC\_ALL=C; sort -k1,1 -k12,12gr -k11,11g -k3,3gr nt\_result | sort -u -k1,1 --merge > bestHits

cat bestHits | cut -f2 > acc

perl -p -i -e 's/\.(\d)//g' acc

acc2tax -i acc -o result -d /home/train/public\_database/accession2taxid/

\# Extract the hits with key words (teleostomi, bacteria) and get the accession number

./Teleostomi\_Bacteria\_extract.pl result

\# Get the seq ID from the accession number

cat bestHits | cut -f1,2 > header

perl -p -i -e 's/\.(\d)//g' header

./get\_contam\_from\_accesion.pl contam\_accession header

\# Also remove those reads with coverage lower than 20

./another\_blob\_result\_seq\_extract.pl M.honghuensis\_1\_blob.blobDB.table.txt

\# Get the seq\_remove\_by\_bam0

cat ../M.honghuensis.plot\_1/seq\_remove\_by\_bam0 final\_contam\_header > true\_bad\_list

remove\_duplicate.pl true\_bad\_list

mv duplicate\_remove true\_bad\_list

\# true\_bad\_list contains the final bad list

**1.12 Remove the contamination by ID**

remove\_contaminant\_by\_ID.pl 

~/genome/M.honghuensis/2\_1\_assembly/M.honghuensis\_200\_scaf.fasta  ../3\_contam\_remove/contam\_blast/ true\_bad\_list

mv survive.fasta M.wulii\_200\_clean.fasta

\# The contam removal process is finished

**2. Fourfold synonymous third-codon transversion (4dTV) analyses**

\# Please check the 4DTv\_calculation.pl in attached Scripts

**3. Plot the genome KEGG annotation results**

\# Extract the KEGG path way information from http://www.genome.jp/kegg/tool/map\_pathway.html

kegg\_plot\_all.pl kegg.path

\# This script will call three other scripts: kegg\_plot\_step\_1.pl, kegg\_plot\_step\_2.pl, kegg\_plot\_step\_3.pl

\# It will output words like:

\# Seperating original files according to Level1: kegg\_plot\_step\_1.pl....

\# Further processing output from former step: kegg\_plot\_step\_2.pl....

\# Calculating: kegg\_plot\_step\_3.pl...

\# Sorting results...

\# Drawing...

\# Done!

\# There will be three outputs: final\_kegg\_for\_ggplot (TSV output), plot3.pdf (final pdf outpuf), test.csv (CSV output)

**4. Plot the radar (spider) chart**

install.packages("fmsb")

library(fmsb)

data <- read.csv("transposon.csv",header=T)

data <- data.frame(A=c(data[1]),B=c(data[2]))

radarchart(data,

`           `cglty = 1,       # Grid line type

`           `cglcol = "gray", # Grid line color

`           `cglwd = 1,       # Line width of the grid

`           `pcol = 4,        # Color of the line

`           `plwd = 2,        # Width of the line

`           `plty = 1)        # Line type of the line

**5. Genome completeness analysis using Benchmarking Universal Single-Copy Orthologs (BUSCO)**

**5.1 Dependencies of BUSCO**

\# Make sure you have installed HMMER (HMMER > 3.1b2) , BLAST+, Augustus (>3.2.1 or 3.2.2) before using BUSCO

**5.2 run BUSCO**

nohup python2.7 BUSCO.py -i M.honghuensis\_200\_scaf.fasta -o M.honghuensis -l /opt/biosoft/busco/metazoa\_odb9 -m geno -c 24 &

nohup /opt/biosoft/python\_3/bin/python3.6 /opt/biosoft/busco/BUSCO.py -i H\_1\_all\_6\_frame -o H\_1\_all\_6\_frame -l /opt/biosoft/busco/metazoa\_odb9/ -m proteins -c 64 &

python2.7 BUSCO.py -i SEQUENCE\_FILE -o OUTPUT\_NAME -l LINEAGE -m geno

python2.7 BUSCO.py -i SEQUENCE\_FILE -o OUTPUT\_NAME -l LINEAGE -m prot

python2.7 BUSCO.py -i SEQUENCE\_FILE -o OUTPUT\_NAME -l LINEAGE -m tran

\# Main result file: short\_summary\_XXXX.txt

**6. Genome completeness analysis using Core Eukaryotic Genes Mapping Approach (CEGMA)**

**6.1 Install dependencies for CEGMA**

\# Install Genewise 2.4.1

\# Install geneid v1.4

wget ftp://genome.crg.es/pub/software/geneid/geneid\_v1.4.4.Jan\_13\_2011.tar.gz

tar zxf geneid\_v1.4.4.Jan\_13\_2011.tar.gz

rm geneid\_v1.4.4.Jan\_13\_2011.tar.gz

cd geneid/

make

cd bin/

echo "PATH=$PATH:/opt/biosoft/geneid/bin" >> ~/.bashrc

source ~/.bashrc

#Install CEGMA

wget http://korflab.ucdavis.edu/datasets/cegma/CEGMA\_v2.5.tar.gz

tar zxf CEGMA\_v2.5.tar.gz

cd CEGMA\_v2.5/

make

\# If it report the error, “Can't locate Cegma.pm in @INC (@INC contains: /usr/local/lib64/perl5 /usr/local/share/perl5 /usr/lib64/perl5/vendor\_perl /usr/share/perl5/vendor\_perl /usr/lib64/perl5 /usr/share/perl5 .) at ./cegma line 34”

\# To solve this error:

sudo cp Cegma.pm /usr/local/lib64/perl5/

sudo cp FAlite.pm /usr/local/lib64/perl5/

sudo cp geneid.pm /usr/local/lib64/perl5/

sudo cp HMMstar.pm /usr/local/lib64/perl5/

echo "PATH=$PATH:/opt/biosoft/CEGMA\_v2.5/bin" >> ~/.bashrc

export CEGMA="/opt/biosoft/CEGMA\_v2.5"

export CEGMATMP="/opt/biosoft/CEGMA\_v2.5"

export PERL5LIB="$PERL5LIB:$CEGMA/lib"

**6.2 Run CEGMA for genomes and transcriptomes**

nohup cegma -g WL\_all\_filter\_Unigene.fasta -T 8 &

\# CEGMA actually has two databases, 458 eukaryotic conserved sequences, and 248 eukaryotic ultraconserved sequences, the latter is a subset of the former.

\# All output results, except “output.completeness\_ report” are based on 248 sequences and others are based on 458 sequences.

**6.3 Run CEGMA for proteomes**

\# CEGMA was not developed for proteomes, but there is a tricky way to use it for evaluating the completeness of proteomes.

cd /opt/biosoft/CEGMA\_v2.5/data/hmm\_profiles

cat ./\*.hmm > KOG.hmm

nohup hmmpfam –cpu 32 KOG.hmm ../../MH\_MCPID.pep > cegma.output &

parse\_hmmpfam\_output.pl cegma.output ../../H\_MCPID.pep 1e-05 458 /opt/biosoft/CEGMA\_v2.5/ > result

**7. Use ggplot2 to make the top hit species distribution plot of genome sequences**

library(ggplot2)

library(ggthemes)

data <- read.csv("Species\_distribution.csv",header=T)

data <- data[seq(30),]   # only use the former 30 lines 

data <- data.frame(A=c(data[1]),B=c(data[2]),C=c(data[3]))  # Make it into a data.frame

data$Species <- factor(data$Species, levels = data$Species[order(data$Count)]) # Reorder the factors

p <- ggplot(data,aes(y=Count, x=Species))

p+geom\_bar(stat="identity", width=0.8, fill= "dodgerblue3",colour="black", cex=0.4)+coord\_flip()+ theme\_set(theme\_bw())+labs(title='The top hit species distribution')+theme(line=element\_line(size=0.5))+ ylab("BLAST Hit")+geom\_text(data=data, aes(x=Species, y=Count, label=paste0(round((Count/829)\*100, digits=2), "%")), hjust=-0.2 ,size=3)+ scale\_y\_continuous(limits=c(0,250))

ggsave("plot2.pdf", width=8, height=8)

**8. Use ggplot2 to visualize the Hisat mapping results of genome analysis (polar plot)**

data<-read.csv("map.csv",header=T) 

data <- data.frame(A=c(data[1]),B=c(data[2]))  # convert it into dataframe

p <- ggplot(data,aes(x=Type,y=Number,fill= Type)) 

p+ geom\_bar(stat="identity")

data$Type <- factor(data$Type, levels = c("Total reads","Total mapped reads","Unique mapped reads","Reads mapped in pair","Spliced mapped reads","Multiple mapped reads"))    # reorder the factors 

p <- ggplot(data,aes(x=reorder(Type,Number), y=Number, fill= Type)) + theme\_set(theme\_bw())  # use the white theme

p+geom\_bar(stat="identity",alpha=0.9,width=0.5,colour='black',cex=0.4)+coord\_polar(theta="y")+labs(title = "Mapping summary")+ theme(axis.text.x = element\_blank(),axis.text.y = element\_blank(),axis.title.x=element\_blank(), axis.title.y=element\_blank(), legend.title = element\_blank(),axis.ticks=element\_blank()) +guides(fill =guide\_legend(keywidth = 1.1, keyheight = 1.1))

**9. Use ggplot2 to make five-way Venn diagram**

library(ggplot2)

library(VennDiagram)

\# Read sample data

A <- sample(1:1000, 400, replace = FALSE);

B <- sample(1:1000, 600, replace = FALSE);

C <- sample(1:1000, 350, replace = FALSE);

D <- sample(1:1000, 550, replace = FALSE);

E <- sample(1:1000, 375, replace = FALSE);

\# Drawing

venn.plot <- venn.diagram(x = list(A = A, B = B, C = C, D = D, E = E), filename = "Venn\_5set\_pretty.png", imagetype = "png", col = "black", fill = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"), alpha = 0.50, cex = c(1.5, 1.5, 1.5, 1.5, 1.5, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 1, 1, 1, 1, 1.5), cat.col = c("dodgerblue", "goldenrod1", "darkorange1", "seagreen3", "orchid3"), cat.cex = 1.5, cat.fontface = "bold", margin = 0.05);

**10. Get the statistics of assembled genome**

\# Please check the genome\_statistic.pl in attached Scripts

**Please cite this paper as**

Guo, Q., Atkinson, S. D., Xiao, B., Zhai, Y., Bartholomew, J. L., & Gu, Z. (2022). A myxozoan genome reveals mosaic evolution in a parasitic cnidarian. BMC biology, 20(1), 1-19. 

**License**

All source code, i.e. scripts/\*.pl, scripts/\*.sh or scripts/\*.py are under the MIT license.
