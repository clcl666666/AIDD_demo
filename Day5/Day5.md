# AIDD_domo
## AI辅助酶设计
### Learning protein fitness landscapes with deep mutational scanning data from multiple sources. 2023. Cell systems文章复现
#### 软件安装
```
git clone https://github.com/cl666666/GVP-MSA.git


conda create -n gvpmsa python==3.9
conda activate gvpmsa
pip3 install torch torchvision --index-url https://download.pytorch.org/whl/cu126

pip install fair-esm 或者pip install git+https://github.com/facebookresearch/esm.git 
pip install numpy==1.24.0
pip install pandas matplotlib scipy omegaconf scikit-learn biotite biopython
```
```
# torch_geometric安装
>>> import torch
>>> print(torch.__version__)
2.7.1+cu126
>>> print(torch.cuda.is_available())
True
>>> print(torch.version.cuda)
12.6

pip install torch-geometric

pip install torch-scatter -f https://data.pyg.org/whl/torch-2.7.0+cu126.html
pip install torch-sparse -f https://data.pyg.org/whl/torch-2.7.0+cu126.html
pip install torch-cluster -f https://data.pyg.org/whl/torch-2.7.0+cu126.html
pip install torch-spline-conv -f https://data.pyg.org/whl/torch-2.7.0+cu126.html
```

```
#测试模型，在单个数据集上随机划分数据，训练和验证。
python train_single_protein_randomly.py --train_dataset_names TEM1 --n_ensembles 3  

```


#### 运行报错解决
```
python ./simple_models/addition.py --dataset_name TEM1
FileNotFoundError: [Errno 2] No such file or directory: '/home/chenlin/directed_evolution/gvp/input_data/TEM1/TEM1_single.csv'
# 文件路径需要换成相对路径。
```
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
#### Unweighted Pair Group Method with Arithmetic mean (UPGMA)聚类，Figtree美化.

UPGMA算法步骤：

1. 计算样本间的欧几里得距离矩阵。
2. 找到距离最小的两个簇，合并它们，并记录合并高度。
3. 更新距离矩阵，使用新簇与其他簇的平均距离。
4. 重复直到只剩一个簇

输出保存成Newick文件，放入figtree(http://tree.bio.ed.ac.uk/software/figtree/)中美化。


# 训练一个简单的transformer的神经网络。
1)	获得序列embedding以构建下游模型（Cell systmes文章举例），从文章github仓库中提炼序列embedding的代码并学习使用。https://github.com/fhalab/MLDE?tab=readme-ov-file#generating-encodings-with-generate_encoding.py，看懂代码中EncodingGenerator的类，将这个类方法用在我们自己的代码上，实现蛋白质序列的不同方式encoding，包括"onehot", "georgiev", “esm”系列模型。
路径是在：MLDE/generate_encoding.py这个代码里面找到MLDE/code/encode/encoding_generator.py，

2） 构建一个简单的transformer模型训练

3） Rosetta, fast relax, ddG的计算。

4）加一个ligandMPNN
