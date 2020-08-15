### 数据处理
- #### fastqc查看序列质量分布情况
##### 使用格式
`fastqc -o 输出文件路径 -t 线程数 fastq文件路径`  
  
fastqc生成的报告的文件名是根据输入来定的，运行结束后会有两个文件，一个.html网页文件，一个是.zip压缩文件。
```
#FastQC v0.11.9
nohup fastqc -o ~/project/lungfish_meta -t 2 \
~/data/Lungfish-Mete/Dormancy/00.original_data/03.Meta_dormancy/X101SC19101272-Z01-J007/1.rawdata/1-1_BDME202032413-1a/1-1_BDME202032413-1a_1.fq.gz >~/project/lungfish_meta/fastqc.log 2>&1 &
```


- #### 使用trimmomatic去除接头和低质量序列
##### 双端测序数据使用方法
```
java -jar trimmomatic-0.32.jar PE \
	[-threads <threads>] \
	[-phred33|-phred64] \
	[-trimlog <trimLogFile>] \
	[<inputFile1> <inputFile2>] \
	[<outputFile1P> <outputFile1U> <outputFile2P> <outputFile2U>] \
	[ILLUMINACLIP:<fastaWithAdaptersEtc>:<seed mismatches>:<palindrome clip threshold>:<simple clip threshold>:<minAdapterLength>:<keepBothReads>] \
	[SLIDINGWINDOW:<windowSize>:<requiredQuality>] \
    [LEADING:<quality>] \
	[TRAILING:<quality>] \
	[MINLEN:<length>]
```
##### 参数说明
`ILLUMINACLIP:TruSeq2-PE.fa:2:40:15`指定包含接头和引物序列的fasta文件的路径，一般情况，从测序公司拿到的数据都是去除了接头和引物序列的。  
`SLIDINGWINDOW:4:25`滑窗剪切，统计滑窗口中所有碱基的平均质量值，如果低于设定的阈值，则切掉窗口。例子：设置4bp窗口，碱基平均质量值阈值为25  
`LEADING:25`从reads的起始端开始切除质量值低于设定阈值的碱基，直到有一个碱基其质量值达到阈值  
`TRAILING:25`从reads的末端开始切除质量值低于设定阈值的碱基，直到有一个碱基质量值达到阈值  
`MINLEN:25`设定一个最短reads长度，当reads经过前面的过滤之后，如果保留下来的长度低于这个阈值，整条 reads 被丢弃。  
输出结果有4个文件，正向和反向双端序列都保留的paired文件和一段序列被丢弃的unpaired文件

##### 例子
```
trimmomatic PE -threads 10 \
     ~/data/Lungfish-Mete/Dormancy/00.original_data/03.Meta_dormancy/X101SC19101272-Z01-J007/1.rawdata/1-1_BDME202032413-1a/1-1_BDME202032413-1a_1.fq.gz \
     ~/data/Lungfish-Mete/Dormancy/00.original_data/03.Meta_dormancy/X101SC19101272-Z01-J007/1.rawdata/1-1_BDME202032413-1a/1-1_BDME202032413-1a_2.fq.gz \
     ${base}_1.qc.fq.gz ${base}_s1_se \
     ${base}_2.qc.fq.gz ${base}_s2_se \
     ILLUMINACLIP:TruSeq2-PE.fa:2:40:15 \
     LEADING:25 \
     TRAILING:25 \
     SLIDINGWINDOW:4:25 \
     MINLEN:100 
```


- #### 使用BWA建立宿主参考基因组索引
```
bwa index -a bwtsw ref.fa
```
##### 参数说明
`-a bwtsw` 用于比对的算法，有两种，一种是`is`,适用于2GB以下的参考基因组的建立，2GB以上的要使用`bwtsw`


- #### 使用BWA进行比对，并使用samtools输出为bam文件
`bwa mem -t 线程数 -M -R reads标头 参考序列文件 fq1 fq2 > | samtools view -bS 输出bam文件路径`
##### 参数说明
`bwa mem`是BWA比对算法中的一种，更快更准确  
`-t` 设置线程数  
`-M` Mark shorter split hits as secondary(for Picard compatibility)将较短的split hits标记为secondary，与picard兼容  
`-R` 设置reads标头
```
#组装单个样本
bwa mem -t 24 -M -R '@RG\tID:1-1_BDME202032413-1a\tLB:1-1_BDME202032413-1a\tPL:Illumina\tSM:1-1_BDME202032413-1a' \
/public/home/wangkun/projects/lungfish/06.bionano_scaffold/lungfish_pilon_bionano.fasta \
~/data/Lungfish-Mete/Dormancy/02.cleandata/1-1_BDME202032413-1a_1.clean.trim.fq.gz \
~/data/Lungfish-Mete/Dormancy/02.cleandata/1-1_BDME202032413-1a_2.clean.trim.fq.gz \
| samtools view -bS > ~/data2/lunfish_meta/Dormancy/bwa_map/1-1_BDME202032413-1a.bam
#Real time: 8349.134 sec; CPU: 197076.041 sec

#批量读取
#!/bin/bash
for i in ~/data/Lungfish-Mete/Dormancy/02.cleandata/*_1.clean.trim.fq.gz
do
echo $i
base=$(basename $i _1.clean.trim.fq.gz)
echo $base
bwa mem -t 10 -M -R '@RG\tID:{base}\tLB:{base}\tPL:Illumina\tSM:{base}' \
/public/home/wangkun/projects/lungfish/06.bionano_scaffold/lungfish_pilon_bionano.fasta \
~/data/Lungfish-Mete/Dormancy/02.cleandata/${base}_1.clean.trim.fq.gz \
~/data/Lungfish-Mete/Dormancy/02.cleandata/${base}_2.clean.trim.fq.gz \
| samtools view -bS > ~/data2/lunfish_meta/Dormancy/03_bwa_map/${base}.bam
done
```
  

- #### 使用samtools将sam文件转换为bam文件  
`samtools view -bS 输入sam文件 > 输出bam文件`
##### 参数说明
`-b`该参数设置输出 BAM 格式，默认下输出是 SAM 格式文件  
`-S`默认下输入是 BAM 文件，若是输入是 SAM文件，则最好加该参数，否则有时候会报错
```
samtools view -bS 1-1_BDME202032413-1a.sam > 1-1_BDME202032413-1a.bam
```


- #### 使用samtools提取出没有比对到参考序列上的比对结果
`samtools view -bu -f 12 -F 256 输入bam文件 > 输出bam文件`
##### 参数说明
`-u`以未压缩的bam格式输出  
`-f`只输出在比对flag中包含该整数的序列  
`-F`跳过比对flag中含有该整数的序列  
`-f`和`-F`后面的数字代表了不同的FLAG值，可以查一下不同FLAG值代表的含义
```
#!/bin/bash
for i in ~/data/Lungfish-Mete/Dormancy/03.mapping/*-1a.bam
do
echo $i
base=$(basename $i .bam)
echo $base
samtools view -bu -f 12 -F 256 ~/data/Lungfish-Mete/Dormancy/03.mapping/${base}.bam > \
~/data2/lunfish_meta/Dormancy/04_select_unmap/${base}_unmap.bam
done
```



- #### 使用samtools对未比对上参考序列的结果进行排序
`samtools sort [option] <in.bam> -o <out.bam>`
```
#!/bin/bash
for i in ~/data2/lunfish_meta/Dormancy/04_select_unmap/*.bam
do
echo $i
base=$(basename $i _unmap.bam)
samtools sort ${base}_unmap.bam -o ${base}_sort.bam
done
```


- #### 使用bedtools将bam文件转换为fastq文件，这一步操作完成之后得到的就是清洗好的宏基因组数据
`bedtools bamtofastq -i 输入bam文件 -fq 输出fq文件`
##### 参数说明
`-i`输入的bam文件  
`-fq`输出的fq文件  
`-fq2`如果是双端测序，就要使用这个参数
```
#!/bin/bash
for i in ~/data2/lunfish_meta/Dormancy/04_select_unmap/*_sort.bam
do
base=$(basename $i _sort.bam)
echo $base
bedtools bamtofastq -i ~/data2/lunfish_meta/Dormancy/04_select_unmap/${base}_sort.bam \
-fq ~/data2/lunfish_meta/Dormancy/05_meta_data/${base}_1.fq \
-fq2 ~/data2/lunfish_meta/Dormancy/05_meta_data/${base}_2.fq
done
```


- #### 使用MegaHit进行单个样本的组装
```
#!/bin/bash
for i in ~/data2/lunfish_meta/Dormancy/05_meta_data/*_1.fq
do
echo $i
base=$(basename $i _1.fq)
echo $base
megahit -1 ~/data2/lunfish_meta/Dormancy/05_meta_data/${base}_1.fq \
-2 ~/data2/lunfish_meta/Dormancy/05_meta_data/${base}_2.fq \
-o ~/data2/lunfish_meta/Dormancy/06_megahit_zz/${base} &
wait
cd ~/data2/lunfish_meta/Dormancy/06_megahit_zz/${base} 
mv final.contigs.fa ../${base}.megahit_contigs.fa
done
```
- #### 使用IDBA进行组装
```
#安装
cd ~/software
git clone https://github.com/loneknightpy/idba.git
cd ~/software/idba
./build.sh
./configure
make

#使用
#IDBA的组装,要将两个paired-end的fq文件保存在一个fq文件中
for i in ~/data2/lunfish_meta/Dormancy/05_meta_data/*_1.fq
do
echo $i
base=$(basename $i _1.fq)
echo $base
~/software/idba/bin/fq2fa --merge ~/data2/lunfish_meta/Dormancy/05_meta_data/${base}_1.fq \
~/data2/lunfish_meta/Dormancy/05_meta_data/${base}_2.fq \
~/data2/lunfish_meta/Dormancy/idba_zz/${base}_12.fq
done

wait
#组装，参数设置有点讲究,根据自己的实际数据来
for i in ~/data2/lunfish_meta/Dormancy/idba_zz/*_12.fq
do
echo $i
base=$(basename $i _12.fq)
echo $base
~/software/idba/bin/idba_ud -r ~/data2/lunfish_meta/Dormancy/idba_zz/${base}_12.fq \
-o ~/data2/lunfish_meta/Dormancy/idba_zz/idba_result \
--pre_correction  \
--maxk 124 \
--min_contig 200 \
--num_threads 4
done
```


- #### 使用prodigal预测ORF
##### 参数说明  
`-a`输出氨基酸文件  
`-c`不允许基因断开，也就是要求完整的ORF，有起始和终止结构  
`-d`输出预测基因的序列文件
`-f`选择输出文件的格式，有gbk，gff和sco格式可以选择  
`-g`指定密码子表，默认11  
`-i`输入文件，即需要预测的基因组的序列文件  
`-m`屏蔽基因组中的N碱基  
`-o`输出文件，默认为屏幕输出  
`-p`选择方式，是single还是meta，默认是single  
`-q`不输出错误信息到屏幕
```
#!/bin/bash
for i in ~/data2/lunfish_meta/Dormancy/06_megahit_zz/*.megahit_contigs.fa
do
echo $i
base=$(basename $i .megahit_contigs.fa)
prodigal -a ${base}.orf.faa \
-i ~/data2/lunfish_meta/Dormancy/06_megahit_zz/${base}.megahit_contigs.fa \
-f gff -o ${base}.gff -p single -q -d ${base}.orf.ffn
done
```

- #### 使用cd-hit聚类去冗余
##### 参数说明
`cd-hit-est`对核苷酸序列进行聚类  
`-i`  fasta格式的输入文件  
`-o`输出文件名  
`-n`  word_length  
`-c`相似性  
`-G`设置为0代表使用本地序列一致性  
`-M`设置为0，不限制内存的使用  
`-d`设置为0表示使用fasta标题中第一个空格前的字段作为序列名字  
`-aS`控制短序列比对严格程度的参数，默认为0，若设为0.9则表示比对区间要占到短序列的90%  
`-T`设置线程数
```
cat *.ffn > dormancy.ffn
cd-hit-est -i dormancy.ffn -o dormancy.geneSet.ffn -n 9 -c 0.95 -G 0 -M 0 -d 0 -aS 0.9 -T 10
```

- #### 分箱
MetaWRAP是一款整合了质控、拼接、分箱、提纯、评估、物种注释、丰度估计、功能注释和可视化的分析流程，纳入超140个工具软件，可一键安装；流程整合了CONCOCT、MaxBin、 metaBAT等三款分箱工具以及提纯和重组装算法。我只使用其中的分箱流程
```
#安装metawrap进行分箱
#新建conda环境
conda create -n metawrap python=2.7

#添加新频道
conda config --add channels ursky

#安装
conda install -y -c ursky metawrap-mg

#要在合并分箱后使用checkM进行分箱评估，需要提前下载好checkM数据库并添加进metawrap的搜索路径
mkdir checkm_folder
checkm data setRoot
# CheckM will prompt to to chose your storage location... Give it the path to the folder you just made.
# Now manually download the database:
cd checkm_folder
wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz
tar -xvf *.tar.gz

#分箱
nohup metawrap binning -t 24 --metabat2 --maxbin2 --concoct \
-a /data2/home/yanxiaoting/project/Golden_Snub-nosed_Monkey/Meta/Faece-8/03.assembly-IDBA/11.assembly.fa \
-o /data2/home/lichunhui/test/metawrap_bin \
/data2/home/lichunhui/test/reads_1.fastq \
/data2/home/lichunhui/test/reads_2.fastq > bin.log 2>&1 &

#合并分箱
nohup metawrap bin_refinement -o bin_refinement -t 20 -c 50 -x 10 \
-A /data2/home/lichunhui/test/metawrap_bin/metabat2_bins/ \
-B /data2/home/lichunhui/test/metawrap_bin/maxbin2_bins/ \
-C /data2/home/lichunhui/test/metawrap_bin/concoct_bins/ > bin_refinement.log 2>&1 &
```


### 注释
- #### 使用kraken2进行注释
- #### 提取kraken2的注释结果进行汇总统计
```
'''
#该脚本编写时的python版本
Python 3.7.6 | packaged by conda-forge | (default, Jun  1 2020, 18:57:50) 
[GCC 7.5.0] on linux
#该脚本实现的内容，将kraken2注释得到的.report文件作为输入，输出结果为域、门、属三个分类水平下，各个样本的微生物的相对比例丰度，输出结果为3个文件，domain_abun.csv,phylum_abun,genus_abun.csv

#python kraken2_annotation.py [输入文件路径] [输出文件路径]

#示例：python kraken2_annotation.py ~/project/lungfish_meta/kraken2/ ~/project/lungfish_meta/kraken2/result/

#注意：该脚本目前只处理过少样本的数据，对于大样本量数据的处理速度尚且未知。
(2020/07/25,李春晖)
'''
from pandas import DataFrame
import pandas as pd
import sys
import os,glob,re


def cal_abun(infile_path,samplename):
    #提取第2列，第4列和第6列的信息
    list1=[]
    list3=[]
    list5=[]

    with open(infile_path,'r') as infile_data:
        for line in infile_data:
            an_info=line.split()
            list1.append(an_info[1])
            list3.append(an_info[3])
            if len(an_info) > 3:
                list5.append('_'.join(an_info[5:len(an_info)]))
                continue
            if len(an_info) == 3:
                list5.append(an_info[5])
    #这里偷了一下懒，作为总数的root上的reads数的位置是固定的，所以直接取来使用            
    list1=[int(i)/int(list1[1]) for i in list1]
    kraken2_annotation=DataFrame({samplename:list1,'rank_code':list3,'name':list5})
    return kraken2_annotation[kraken2_annotation[samplename] != 0]

def main():
    path = sys.argv[1]#读取路径
    outpath = sys.argv[2]#读取输出路径
    infile_path=glob.glob(os.path.join(path,'*.report'))#获取路径下全部的.report文件的路径

    file_name=['rank_code','name']#准备列名
    #得到所有.report文件的文件名
    for i in infile_path:
        mode='[\w\-]*.report'
        file_name.extend(re.findall(mode, i))


    result_df=cal_abun(infile_path[0],file_name[2])#先用第一个样本的数据生成一个dataframe作为总的dataframe

    #再分别将其余样本的数据加入到总的dataframe中
    for i in range(1,len(infile_path)):
        df=cal_abun(infile_path[i],file_name[i+2])
        result_df=pd.merge(result_df,df,on=['name','rank_code'],how='outer')
    
    #按照域、门、属三类进行输出
    result_df[result_df['rank_code'].isin(['D'])][file_name].to_csv(outpath+'domain_abun.csv',index=False)
    result_df[result_df['rank_code'].isin(['P'])][file_name].to_csv(outpath+'phylum_abun.csv',index=False)
    result_df[result_df['rank_code'].isin(['G'])][file_name].to_csv(outpath+'genus_abun.csv',index=False)


if __name__ == '__main__':
    main()
```

- #### 使用diamond在数据库中进行比对注释
##### 下载NCBI-NR数据库，并构建diamond参考
```
#NCBI-NR数据库比较大，下载解压后有140G的大小（2020/8/5）
nohup wget -c ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz >log.log 2>&1 &
#多线程解压，提高速度，-k:解压后不删除压缩文件
unpigz -k -p 16 nr.gz
#构建diamond参考
nohup diamond makedb --in nr -d nr -p 12 > log.log 2>&1 &
```


##### 在NR数据库中进行比对注释
在蛋白质数据库中进行比对，最好是要把核苷酸序列提前转换为蛋白质序列，但是diamond也提供了blastx用于将核苷酸序列比对到蛋白质序列上，来实现序列的注释，更多使用细节请参考diamond的官方文档
http://www.diamondsearch.org/index.php
```
#diamond的输出文件可以有两种选择，分别代表两种后续的处理方式，方法一是输出.daa文件，该格式是DIAMOND专有的二进制格式，另外一个软件MEGAN也支持它，并允许快速导出结果。
nohup time diamond blastx --db /data2/home/lichunhui/database/nr/nr.dmnd \
-t /data2/home/lichunhui/lunfish_meta/Dormancy/08_diamond_nr/ \
-p 20 \
-q /data2/home/lichunhui/lunfish_meta/Dormancy/06_megahit_zz/dormancy.fa \
--daa /data2/home/lichunhui/lunfish_meta/Dormancy/08_diamond_nr/dormancy.daa \
> /data2/home/lichunhui/lunfish_meta/Dormancy/08_diamond_nr/log.log 2>&1 &

#方法二是直接生成.m8格式的文件
nohup time diamond blastx --db /data2/home/lichunhui/database/nr/nr.dmnd \
-t /data2/home/lichunhui/lunfish_meta/Dormancy/08_diamond_nr/ \
-p 20 \
-q /data2/home/lichunhui/lunfish_meta/Dormancy/06_megahit_zz/dormancy.fa \
-o /data2/home/lichunhui/lunfish_meta/Dormancy/08_diamond_nr/dormancy.nr.m8 \
> /data2/home/lichunhui/lunfish_meta/Dormancy/08_diamond_nr/log.log 2>&1 &
```


##### 在KEGG数据库中进行比对注释
因为KEGG数据库要下到本地是收费的（大概），所以在网上没有找到相关的下载方式，这里还是用到了KEGG数据库来进行diamond建库
```
#建库
nohup time diamond makedb --in /public/home/yanxiaoting/database/kegg_all_clean/kegg_all_clean.fa -d /data2/home/lichunhui/database/kegg -p 10 > log.log 2>&1 &

#比对注释
nohup time diamond blastx --db /data2/home/lichunhui/database/kegg/kegg.dmnd -t /data2/home/lichunhui/lunfish_meta/Dormancy/08_diamond_nr/ -p 20 -q /data2/home/lichunhui/lunfish_meta/Dormancy/06_megahit_zz/dormancy.fa -o /data2/home/lichunhui/lunfish_meta/Dormancy/08_diamond_nr/dormancy.kegg.m8 > /data2/home/lichunhui/lunfish_meta/Dormancy/08_diamond_nr/log.log 2>&1 &
```


- #### 使用eggnog-mapper进行注释
使用eggnog-mapper这个软件的初衷是为了进行eggNOG的注释，但是在查看官方文档的过程中发现了这个软件还能提供一些其他常用数据库的注释，例如nr，kegg，nog等，该软件除了有v1.0还有v2.0，且v2.0提供了更多数据库的注释。
```
#eggnog-mapper依赖python2.7
#创建新的conda环境
conda create -n eggnog python=2.7

#本地下载
wget https://github.com/eggnogdb/eggnog-mapper/archive/1.0.3.tar.gz -O ~/software/eggnog-mapper-1.0.3.tar.gz

#解压
tar zxf ~/software/eggnog-mapper-1.0.3.tar.gz -C /data2/home/lichunhui/database/eggnog/

#进入软件目录并下载数据库到软件默认的data目录下（也可以自己指定下载目录）
cd /data2/home/lichunhui/database/eggnog/eggnog-mapper-1.0.3
nohup ./download_eggnog_data.py euk bact arch viruses > log.log 2>&1 &

#注释，基于hmmer
nohup python /data2/home/lichunhui/database/eggnog/eggnog-mapper-1.0.3/emapper.py -d bact -i /data2/home/lichunhui/lunfish_meta/Dormancy/06_megahit_zz/dormancy.fa --data_dir /data2/home/yanxiaoting/biosoft/eggnog-mapper/eggnog-mapper-1.0.3/data -o eggNOG_hmmer --cpu 20 --translate > log_hmmer.log 2>&1 &
```
```
#安装v2.0版本的eggnog-mapper
git clone https://github.com/jhcepas/eggnog-mapper.git

#下载数据库，使用软件目录下的这个.py脚本自动下载eggnog数据库
download_eggnog_data.py

#注释，基于diamond
nohup python /data2/home/lichunhui/database/eggnog/eggnog-mapper/emapper.py -m diamond -i /data2/home/lichunhui/lunfish_meta/Dormancy/06_megahit_zz/dormancy.fa -o eggNOG_diamond  --cpu 24 --translate > log_diamond.log 2>&1 &
```