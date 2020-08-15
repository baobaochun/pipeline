> 这个流程的数据是三个健康人类的粪便宏基因组ONT测序数据，我自己用于初步熟悉三代ONT数据的分析流程。数据下载地址：https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA508395 ，这个流程还有许多没有完善的地方，需要后面进行更多的尝试来补充。

- #### 从NCBI上下载数据
要先从NCBI上下载数据，推荐的方式是使用NCBI官方的下载工具sratoolkit，要下载对应的服务器系统的版本
```
#查看服务器系统
cat /proc/version
#进入软件安装目录
cd /public/home/lichunhui/software/
#下载相应版本的软件
wget -c https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.10.8/sratoolkit.2.10.8-centos_linux64.tar.gz
#解压
tar xzvf sratoolkit.2.10.8-centos_linux64.tar.gz
```
解压之后可以直接进行使用，也可以把软件添加进环境变量，我没有进行添加。NCBI里面的数据以.sra的格式进行存储，每一个数据都有一个对应的编号，下载的时候直接输入相应的编号即可
```
#将数据下载到指定文件夹下，-p显示进度
/public/home/lichunhui/software/sratoolkit.2.10.8-centos_linux64/bin/prefetch -p -O /data2/home/lichunhui/human_stool SRR8427257

#下载下来的数据是.sra格式，要对其进行拆包，转换为.fastq格式
cd /data2/home/lichunhui/human_stool/SRR8427257
/public/home/lichunhui/software/sratoolkit.2.10.8-centos_linux64/bin/fastq-dump SRR8427257.sra
rm SRR8427257.sra

#然后在该目录下就得到对应的.fastq文件,可以进行查看
less /data2/home/lichunhui/human_stool/SRR8427257/SRR8427257.fastq
```

- #### 数据质量查看
还是先用fastqc对原始数据进行质量的查看与评估。fastqc生成两个文件，查看其中的.html文件，发现SRR8427257.fastq的前几个碱基是随机引入的序列，需要进行切除，序列的Q值不高。
```
fastqc -o /data2/home/lichunhui/human_stool/00_fastqc/ -t 10 /data2/home/lichunhui/human_stool/SRR8427257/SRR8427257.fastq
```

- #### 数据初步质控
使用nanofilt软件进行初步的质量控制，初步的质控后可以再使用fastqc查看序列
```
#针对3代数据处理，创建新的conda环境
conda create -n 'nanopore' python='3.7'

#安装nanofilt并使用其来进行质控，这里筛选掉了长度<1000bp和Q值<8的序列，并切除了前40bp
NanoFilt /data2/home/lichunhui/human_stool/SRR8427257/SRR8427257.fastq -l 1000 -q 8 --headcrop 40 > /data2/home/lichunhui/human_stool/01_nanofilt/SRR8427257_nanofilt.fastq
```

- #### 去除宿主污染
minimap2是针对长读段数据的比对软件，将数据与人类参考基因组进行比对，去除宿主污染的影响
```
#minimap.sh
minimap2 -ax map-ont /public/home/renqingmiao2018/project/yak.rumen.metagenome/00.ref/human.genome/GCF_000001405.38_GRCh38.p12_genomic.fna /data2/home/lichunhui/human_stool/01_nanofilt/SRR8427257_nanofilt.fastq -t 12 | samtools view -bS > /data2/home/lichunhui/human_stool/02_minimap2/SRR8427257.bam

nohup time sh /data2/home/lichunhui/human_stool/02_minimap2/minimap.sh > /data2/home/lichunhui/human_stool/02_minimap2/minimap2.log 2>&1 &
```

- #### 比对后bam文件处理
```
#查看统计信息,主要看reads的比对率
cd /data2/home/lichunhui/human_stool/02_minimap2
samtools flagstat SRR8427257.bam

#在二代数据中我们要提取的数据的flag值是12，表示的是提取出没有比对上的序列及其另一端序列，但是在三代数据中，由于没有双端数据，因此要修改提取的FLAG值为4，否则不能提取出信息
samtools view -bu -f 4 SRR8427257.bam > SRR8427257_unmap.bam

#sort排序
samtools sort SRR8427257_unmap.bam -o SRR8427257_unmap_sort.bam

#得到fastq文件
bedtools bamtofastq -i SRR8427257_unmap_sort.bam -fq SRR8427257_meta1.fastq
```

- #### 组装
根据文献中使用的方法，使用canu v2.0和flye v2.8两个组装软件对序列进行组装。组装结果的好坏就我组装出来的结果来看受参数genomesize的影响较大。  
**canu组装**
```
#canu组装genomesize=50m
#canu.sh
canu -p SRR8427257 -d SRR8427257_canu_50m genomeSize=50m -nanopore /data2/home/lichunhui/human_stool/02_minimap2/SRR8427257_meta1.fastq

cd /data2/home/lichunhui/human_stool/04_assemble
nohup time sh canu.sh 2>&1 &

#质量评估
python /public/home/lichunhui/software/quast/quast.py -o /data2/home/lichunhui/human_stool/04_assemble/canu50m_quast /data2/home/lichunhui/human_stool/04_assemble/SRR8427257_canu_50m/SRR8427257.contigs.fasta


#canu组装genosize=100m
#canu.sh
canu -p SRR8427257 -d SRR8427257_canu_100m genomeSize=100m -nanopore /data2/home/lichunhui/human_stool/02_minimap2/SRR8427257_meta1.fastq

cd /data2/home/lichunhui/human_stool/04_assemble
nohup time sh canu.sh 2>&1 &

#质量评估
python /public/home/lichunhui/software/quast/quast.py -o /data2/home/lichunhui/human_stool/04_assemble/canu100m_quast /data2/home/lichunhui/human_stool/04_assemble/SRR8427257_canu_100m/SRR8427257.contigs.fasta
```
**flye组装**
```
#flye组装genosize=100m
flye --nano-raw /data2/home/lichunhui/human_stool/02_minimap2/SRR8427257_meta1.fastq \
--out-dir /data2/home/lichunhui/human_stool/04_assemble/SRR8427257_flye_100m \
--genome-size 100m -t 12 --meta

cd /data2/home/lichunhui/human_stool/04_assemble
nohup time sh flye.sh > flye_100m.log 2>&1 &

#质量评估
python /public/home/lichunhui/software/quast/quast.py -o /data2/home/lichunhui/human_stool/04_assemble/flye100m_quast /data2/home/lichunhui/human_stool/04_assemble/SRR8427257_flye_100m/assembly.fasta


#flye组装genomesize=250m
flye --nano-raw /data2/home/lichunhui/human_stool/02_minimap2/SRR8427257_meta1.fastq \
--out-dir /data2/home/lichunhui/human_stool/04_assemble/SRR8427257_flye_250m \
--genome-size 250m -t 12 --meta

cd /data2/home/lichunhui/human_stool/04_assemble
nohup time sh flye.sh > flye_250m.log 2>&1 &

#质量评估
python /public/home/lichunhui/software/quast/quast.py -o /data2/home/lichunhui/human_stool/04_assemble/flye250m_quast /data2/home/lichunhui/human_stool/04_assemble/SRR8427257_flye_250m/assembly.fasta
```

- #### 对reads和contig进行物种注释，对contig进行功能注释
因为三代数据长读长的特点，可以直接对reads进行注释
```
#使用kraken2对质控后的reads进行注释，kraken2极其耗内存！！！
# read_kraken2.sh
/public/home/lichunhui/software/kraken2/kraken2 --db /data2/home/yanxiaoting/biosoft/kraken2/KrakenDB/KrakenDB_201809/ --threads 10 --report-zero-counts --output SRR8427257.readinfo --report SRR8427257.report /data2/home/lichunhui/human_stool/02_minimap2/SRR8427257_meta1.fastq

cd /data2/home/lichunhui/human_stool/03_read_kraken2
nohup time sh read_kraken2.sh > read_kraken2.log 2>&1 &


#使用kraken2对组装后的contigs进行注释
#contig_kraken2.sh
/public/home/lichunhui/software/kraken2/kraken2 --db /data2/home/yanxiaoting/biosoft/kraken2/KrakenDB/KrakenDB_201809/ --threads 10 --report-zero-counts --output SRR8427257_contig.readinfo --report SRR8427257_contig.report /data2/home/lichunhui/human_stool/04_assemble/SRR8427257_canu_100m/SRR8427257.contigs.fasta

cd /data2/home/lichunhui/human_stool/05_congtig_kraken2
nohup time sh contig_kraken2.sh > contig_kraken2.log 2>&1 &


#使用eggnog-mapper对contig进行同源蛋白注释
#eggnog.sh
python /data2/home/lichunhui/database/eggnog/eggnog-mapper/emapper.py -m diamond -i /data2/home/lichunhui/human_stool/04_assemble/SRR8427257_canu_100m/SRR8427257.contigs.fasta -o eggNOG_diamond  --cpu 24 --translate > log_eggnog.log 2>&1 &

cd /data2/home/lichunhui/human_stool/06_eggnog
nohup eggnog.sh > eggnog.log 2>&1 &
```
