# 自定义StrainPanDA数据库

## 建泛基因组

### Step1: 获取所有菌株的信息

可以从Genome List in NCBI下载相关的信息

例如, 对 *E. coli*, 可以从此下载： [Genome List - Genome - NCBI (nih.gov)](https://www.ncbi.nlm.nih.gov/genome/browse#!/prokaryotes/Escherichia coli)

### Step2:

在SimStr的docker里 (see [SimStr](https://github.com/xbiome/StrainPanDA/blob/main/SimStr/README.md), 自动选择代表序列并产生泛基因组，例如:

```
fileout=Bifidobacterium-longum-202009
sh /script/Ref_str_select_ANI_mash.sh Bifidobacterium-longum-20200915.csv $fileout /script/ 99
```

其中，

Bifidobacterium-longum-20200915.csv 是step1里下载的 *B. longum* 文件; 

Bifidobacterium-longum-202009 是泛基因组数据库的名字，他有固定的格式: genus-species-yearmonth(in number) ;

99 是最大可接受的两个基因组之间的ANI，如果两个基因组之间ANI大于这个值，将会被合并到一个cluster。目前建议使用该值。如果增加ANI的相似性，StrainPanDA的结果可靠性可能反而下降。

/script/ 是SimStr在docker里面的路径。

注：如果不能使用docker，可以使用singularity进行类似的操作

### 输出文件:

只有$fileout 文件夹 (例如上述例子：Bifidobacterium-longum-202009) 是泛基因组. 他的格式和PanPhlAn是一样的 (see [Pangenome generation (old) · SegataLab/panphlan Wiki (github.com)](https://github.com/SegataLab/panphlan/wiki/Pangenome-generation-(old)))

所有其他文件都是中间文件，如果后续还要在想通路径建泛基因组，需要先移除。

> **注意**
> 1. 如果同时建多个库，每个建库需要在分开的文件夹，否则互相之间会影响.
> 2. 如果panphlan报类似如下错误 "zlib.error: Error -3 while decompressing: invalid block type", 这是因为wget下载的基因组不完整。此时，需要手动修正这些文件。用户可以通过下述命令`check_gz_integrity.sh`,检查文件的完整性。然后利用 `broken_ID.txt` 去重新下载(例如, `cd strain-pan;grep -f broken_ID.txt $fileout"_ref_dl_list.txt_sel.txt" > broken_ref_dl_list.txt; wget -c broken_ref_dl_list.txt`). 有时候，失败文件是由于`wget -c`造成的,此时你可以直接用wget重新下载，例如上述例子改成 `cd strain-pan;grep -f broken_ID.txt $fileout"_ref_dl_list.txt_sel.txt" > broken_ref_dl_list.txt; wget -i broken_ref_dl_list.txt`.
> 3. 此外，如果panphlan失败了，需要先删除所有和panphlan相关的输出文件夹再重新跑流程, 否则panphlan会报错说文件夹已存在.  (一般来说, 主要有下载不完整引起,此时用户可以下载完了以后，直接单独跑panphlan: 如: `panphlan_pangenome_generation.py -c $fileout --i_fna strain-pan/ --i_gff strain-pan/ -o $fileout`. 在本案例里, file_out 是 `Bifidobacterium-longum-202009`
> 4. 构建时使用的源于step1的 `.csv` 文件必须是最新的，否则有些下载链接可能因为NCBI更新而已经失效。

## 注释泛基因组

所有的注释过程也必须在SimStr的docker里运行，因为涉及比较多依赖包 (see [SimStr](https://github.com/xbiome/StrainPanDA/blob/main/SimStr/README.md))



### CAZy注释

```
python /script/cazy_anno.py -i . -o cazy -d path_to_CAZy_DB
```

-i 输入目录, 可以是当期目录，或者是建泛基因组时用的$fileout 输出目录.可以支持同时对多个泛基因组进行注释

-d path_to_CAZy_DB 指的时CAZy (dbCAN2) 数据库. 可以从以下方式获取: (see [Index of /dbCAN2/download (unl.edu)](https://bcb.unl.edu/dbCAN2/download/) for updates and instructions)

```
wget http://bcb.unl.edu/dbCAN2/download/CAZyDB.07312019.fa.nr && diamond makedb --in CAZyDB.07312019.fa.nr -d CAZy \
    && wget http://bcb.unl.edu/dbCAN2/download/Databases/dbCAN-HMMdb-V8.txt && mv dbCAN-HMMdb-V8.txt dbCAN.txt && hmmpress dbCAN.txt \
    && wget http://bcb.unl.edu/dbCAN2/download/Databases/tcdb.fa && diamond makedb --in tcdb.fa -d tcdb \
    && wget http://bcb.unl.edu/dbCAN2/download/Databases/tf-1.hmm && hmmpress tf-1.hmm \
    && wget http://bcb.unl.edu/dbCAN2/download/Databases/tf-2.hmm && hmmpress tf-2.hmm \
    && wget http://bcb.unl.edu/dbCAN2/download/Databases/stp.hmm && hmmpress stp.hmm 
```

#### 注释的输出

最终的输出是 cazy/cazy 目录. 每个泛基因组会有一个文件. 每个文件里, 第一列是基因家族ID, 第二列是CAZy catalog ID.



### VFDB annotation

```
python /script/vfdb_anno.py -i . -o vfdb_out -db $path_to_VFDB -pi 50 -cov 50
```

-i 输入目录, 可以是当期目录，或者是建泛基因组时用的$fileout 输出目录.可以支持同时对多个泛基因组进行注释

-d $path_to_VFDB 指的时CAZy (dbCAN2) 数据库. 可以从以下方式获取:

```
diamond makedb --in VFDB_setA_pro.fas -d $path_to_VFDB
```

The VFDB_setA_pro.fas 从此处下载 [VFDB: Virulence Factors of Bacterial Pathogens (mgc.ac.cn)](http://www.mgc.ac.cn/cgi-bin/VFs/v5/main.cgi)，下载时选择protein sequences of core dataset链接。(当前版本于2021_04_13下载)

$path_to_VFDB 可以是任何用户指定路径.

-pi and -cov 参数是用于指定blast 的 identity 和 coverage percentage的阈值. 

#### 注释的输出

最终的输出是 vfdb_out/final_out folder. 目录. 每个泛基因组会有一个文件. 每个文件里, 第一列是基因家族ID, 第二列是 VFDB ID.

### EggNOG mapper 注释

一般的功能注释，包括（KEGG KO等）可以通过以下方式进行 [EggNOG mapper](http://eggnog-mapper.embl.de/).

```
for i in ${pangenome_db_path}/*/*.ffn; do bs=`basename $i`;emapper.py -i $i --itype CDS --cpu 16 -o ${bs%%.ffn};done
for i in  ${pangenome_db_path}/* ; do bs=`basename $i` ;  cp panphlan_${bs}_centroids.emapper.annotations $i/${bs}_emapper_anno.tsv;done
```

#### 注释的输出

最终的输出会在每个泛基因组的文件夹. 每列的注释，在 `{species_version}_emapper_anno.tsv` 文件里. 例如:

```
#query	seed_ortholog	evalue	score	eggNOG_OGs	max_annot_lvl	COG_category	Description	Preferred_name	GOs	EC	KEGG_ko	KEGG_Pathway	KEGG_Module	KEGG_Reaction	KEGG_rclass	BRITE	KEGG_TC	CAZy	BiGG_Reaction	PFAMs
```

