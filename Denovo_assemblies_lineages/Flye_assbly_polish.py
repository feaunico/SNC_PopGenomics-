import os

os.system('makedir racon_polishing')
os.system('makedir pilon_polishing')

ON =      #path to fastq.gz with Oxford Nanopore reads

os.system('flye --nano-raw ' + ON + ' --genome-size 35 --threads 24 -o flye_genome --iterations 1')

genome = 'DSEP_Walh_flye_genome/assembly.fasta'
os.system('minimap2 -a -t 14 ' + genome + ' ' + ON + ' > racon_polishing/mappings/map1.sam')
os.system('racon/build/bin/racon -t 14 ' + ON + ' racon_polishing/mappings/map1.sam ' + genome + ' > racon_polishing/assembly_pol1.fasta')

os.system('minimap2 -a -t 14 racon_polishing/assembly_pol1.fasta ' + ON + ' > racon_polishing/mappings/map2.sam')
os.system('racon/build/bin/racon -t 14 ' + ON + ' racon_polishing/mappings/map2.sam racon_polishing/assembly_pol1.fasta > racon_polishing/assembly_pol2.fasta')

os.system('minimap2 -a -t 14 racon_polishing/assembly_pol2.fasta ' + ON + ' > racon_polishing/mappings/map3.sam')
os.system('racon/build/bin/racon -t 14 ' + ON + ' racon_polishing/mappings/map3.sam racon_polishing/assembly_pol2.fasta > racon_polishing/assembly_pol3.fasta')

os.system('minimap2 -a -t 14 racon_polishing/assembly_pol3.fasta ' + ON + ' > racon_polishing/mappings/map4.sam')
os.system('racon/build/bin/racon -t 14 ' + ON + ' racon_polishing/mappings/map4.sam racon_polishing/assembly_pol3.fasta > racon_polishing/assembly_pol4.fasta')

lib1 = #path to fastq.gz with Illumina paired end reads - mate 1
lib2 = #path to fastq.gz with Illumina paired end reads - mate 2
#1
os.system('samtools faidx racon_polishing/assembly_pol4.fasta')
os.system('bwa index racon_polishing/assembly_pol4.fasta')
os.system("bwa mem -t $(nproc) -M racon_polishing/assembly_pol4.fasta " + lib1 + " " +  lib2 + " > pilon_polishing/mappings/map5.sam")
os.system("samtools view -S -b pilon_polishing/mappings/map5.sam > pilon_polishing/mappings/map5.bam")
os.system('samtools sort pilon_polishing/mappings/map5.bam -o pilon_polishing/mappings/map5s.bam')
os.system('samtools index pilon_polishing/mappings/map5s.bam')
os.system('java -jar ../../pilon-1.24.jar --genome racon_polishing/assembly_pol4.fasta --frags pilon_polishing/mappings/map5s.bam --output assembly_pol5 --outdir pilon_polishing')

#2
os.system('samtools faidx pilon_polishing/assembly_pol5.fasta')
os.system('bwa index pilon_polishing/assembly_pol5.fasta')
os.system("bwa mem -t $(nproc) -M pilon_polishing/assembly_pol5.fasta " + lib1 + " " +  lib2 + " > pilon_polishing/mappings/map6.sam")
os.system("samtools view -S -b pilon_polishing/mappings/map6.sam > pilon_polishing/mappings/map6.bam")
os.system('samtools sort pilon_polishing/mappings/map6.bam -o pilon_polishing/mappings/map6s.bam')
os.system('samtools index pilon_polishing/mappings/map6s.bam')
os.system('java -jar ../../pilon-1.24.jar --genome pilon_polishing/assembly_pol5.fasta --frags pilon_polishing/mappings/map6s.bam --output assembly_pol6 --outdir pilon_polishing')

#3
os.system('samtools faidx pilon_polishing/assembly_pol6.fasta')
os.system('bwa index pilon_polishing/assembly_pol6.fasta')
os.system("bwa mem -t $(nproc) -M pilon_polishing/assembly_pol6.fasta " + lib1 + " " +  lib2 + " > pilon_polishing/mappings/map7.sam")
os.system("samtools view -S -b pilon_polishing/mappings/map7.sam > pilon_polishing/mappings/map7.bam")
os.system('samtools sort pilon_polishing/mappings/map7.bam -o pilon_polishing/mappings/map7s.bam')
os.system('samtools index pilon_polishing/mappings/map7s.bam')
os.system('java -jar ../../pilon-1.24.jar --genome pilon_polishing/assembly_pol6.fasta --frags pilon_polishing/mappings/map7s.bam --output assembly_pol7 --outdir pilon_polishing')

#4
os.system('samtools faidx pilon_polishing/assembly_pol7.fasta')
os.system('bwa index pilon_polishing/assembly_pol7.fasta')
os.system("bwa mem -t $(nproc) -M pilon_polishing/assembly_pol7.fasta " + lib1 + " " +  lib2 + " > pilon_polishing/mappings/map8.sam")
os.system("samtools view -S -b pilon_polishing/mappings/map8.sam > pilon_polishing/mappings/map8.bam")
os.system('samtools sort pilon_polishing/mappings/map8.bam -o pilon_polishing/mappings/map8s.bam')
os.system('samtools index pilon_polishing/mappings/map8s.bam')
os.system('java -jar ../../pilon-1.24.jar --genome pilon_polishing/assembly_pol7.fasta --frags pilon_polishing/mappings/map8s.bam --output assembly_pol8 --outdir pilon_polishing')


#5
os.system('samtools faidx pilon_polishing/assembly_pol8.fasta')
os.system('bwa index pilon_polishing/assembly_pol8.fasta')
os.system("bwa mem -t $(nproc) -M pilon_polishing/assembly_pol8.fasta " + lib1 + " " +  lib2 + " > pilon_polishing/mappings/map9.sam")
os.system("samtools view -S -b pilon_polishing/mappings/map9.sam > pilon_polishing/mappings/map9.bam")
os.system('samtools sort pilon_polishing/mappings/map9.bam -o pilon_polishing/mappings/map9s.bam')
os.system('samtools index pilon_polishing/mappings/map9s.bam')
os.system('java -jar ../../pilon-1.24.jar --genome pilon_polishing/assembly_pol8.fasta --frags pilon_polishing/mappings/map9s.bam --output assembly_pol9 --outdir pilon_polishing')

#6
os.system('samtools faidx pilon_polishing/assembly_pol9.fasta')
os.system('bwa index pilon_polishing/assembly_pol9.fasta')
os.system("bwa mem -t $(nproc) -M pilon_polishing/assembly_pol9.fasta " + lib1 + " " +  lib2 + " > pilon_polishing/mappings/map10.sam")
os.system("samtools view -S -b pilon_polishing/mappings/map10.sam > pilon_polishing/mappings/map10.bam")
os.system('samtools sort pilon_polishing/mappings/map10.bam -o pilon_polishing/mappings/map10s.bam')
os.system('samtools index pilon_polishing/mappings/map10s.bam')
os.system('java -jar ../../pilon-1.24.jar --genome pilon_polishing/assembly_pol9.fasta --frags pilon_polishing/mappings/map10s.bam --output assembly_pol10 --outdir pilon_polishing')

#7
os.system('samtools faidx pilon_polishing/assembly_pol10.fasta')
os.system('bwa index pilon_polishing/assembly_pol10.fasta')
os.system("bwa mem -t $(nproc) -M pilon_polishing/assembly_pol10.fasta " + lib1 + " " +  lib2 + " > pilon_polishing/mappings/map11.sam")
os.system("samtools view -S -b pilon_polishing/mappings/map11.sam > pilon_polishing/mappings/map11.bam")
os.system('samtools sort pilon_polishing/mappings/map11.bam -o pilon_polishing/mappings/map11s.bam')
os.system('samtools index pilon_polishing/mappings/map11s.bam')
os.system('java -jar ../../pilon-1.24.jar --genome pilon_polishing/assembly_pol10.fasta --frags pilon_polishing/mappings/map11s.bam --output assembly_pol11 --outdir pilon_polishing')


#8
os.system('samtools faidx pilon_polishing/assembly_pol11.fasta')
os.system('bwa index pilon_polishing/assembly_pol11.fasta')
os.system("bwa mem -t $(nproc) -M pilon_polishing/assembly_pol11.fasta " + lib1 + " " +  lib2 + " > pilon_polishing/mappings/map12.sam")
os.system("samtools view -S -b pilon_polishing/mappings/map12.sam > pilon_polishing/mappings/map12.bam")
os.system('samtools sort pilon_polishing/mappings/map12.bam -o pilon_polishing/mappings/map12s.bam')
os.system('samtools index pilon_polishing/mappings/map12s.bam')
os.system('java -jar ../../pilon-1.24.jar --genome pilon_polishing/assembly_pol11.fasta --frags pilon_polishing/mappings/map12s.bam --output assembly_pol12 --outdir pilon_polishing')
