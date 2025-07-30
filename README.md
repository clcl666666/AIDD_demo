# AIDD_course
## AI辅助酶设计
### Discovery of deaminase functions by structure-based protein clustering. 2023. Cell文章复现
#### pdb文件下载
1. https://www.ebi.ac.uk/interpro/ 网页上，使用DddA deaminase检索
2. 批量下载pfam家族PF14428中的alphafold预测的结构。首先将界面中的表格复制粘贴出来获得AlphaFold DB的accession id。得到文件pfam14428_f100_id.csv
3. 运行download_alphafold_pdb.py文件获得100个pdb文件，下载至pdb_files文件夹中

#### 使用TM-score计算两两结构之间的相似性

先安装软件：
```
wget https://zhanggroup.org/TM-score/TMscore.cpp # 下载源码
g++ -O3 -o TMscore TMscore.cpp # 编译 
./TMscore model.pdb native.pdb 
```
