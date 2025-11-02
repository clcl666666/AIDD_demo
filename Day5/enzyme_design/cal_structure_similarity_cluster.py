from scipy.cluster.hierarchy import linkage, dendrogram,to_tree

# pip install py-upgma
from pathlib import Path
import subprocess
import numpy as np
import matplotlib.pyplot as plt
# cmd = './TMscore ../AIDD_demo/pdb_files/AF-A0A1M4Z5B2-F1.pdb ../AIDD_demo/pdb_files/AF-A0A1M4ZGM6-F1.pdb'

def linkage_to_newick(linkage_matrix, labels,output_file="tree.nwk"):
    """
    Convert a linkage matrix to Newick format.
    
    Parameters:
    linkage_matrix : numpy.ndarray
        The hierarchical clustering encoded as a linkage matrix
    labels : list
        List of labels for the samples
    
    Returns:
    newick_str : str
        The tree in Newick format
    """
    # Convert linkage matrix to a tree structure
    tree = to_tree(linkage_matrix, rd=False)
    
    def _to_newick(node, labels):
        if node.is_leaf():
            return f"{labels[node.id]}:{node.dist/2:.4f}"
        else:
            left = _to_newick(node.get_left(), labels)
            right = _to_newick(node.get_right(), labels)
            return f"({left},{right}):{node.dist/2:.4f}"
    
    newick_str = _to_newick(tree, labels) + ";"
    with open(output_file, 'w') as f:
        f.write(newick_str)
    return 

def get_tmscore(pdb1, pdb2, tmscore_exe="./TMscore"):
    try:
        # 运行 TM-score 命令，捕获输出
        result = subprocess.run(
            [tmscore_exe, pdb1, pdb2],
            capture_output=True,
            text=True,
            check=True
        )
        
        # 解析 stdout 中的 TM-score 值
        for line in result.stdout.splitlines():
            if "TM-score    =" in line:
                tm_score = float(line.split("=")[1].split()[0])
                return tm_score
        
        raise ValueError("TM-score not found in output")
    
    except subprocess.CalledProcessError as e:
        print(f"Error: {e.stderr}")
        raise
    except Exception as e:
        print(f"Error processing output: {str(e)}")
        raise

def get_distance_matrix(pdb_files):
    n = len(pdb_files)
    tmscore_matrix = np.zeros((n, n))
    labels = [f.stem for f in pdb_files]  # 使用文件名（不含路径和扩展名）作为标签

    # 只计算上三角（对称矩阵）
    for i in range(n):
        tmscore_matrix[i, i] = 1.0  # 自身比较的 TM-score 为 1
        for j in range(i + 1, n):
            tm_score = get_tmscore(pdb_files[i], pdb_files[j])
            tmscore_matrix[i, j] = tmscore_matrix[j, i] = tm_score

    return tmscore_matrix, labels

pdb_dir = Path('pdb_files_new')
pdb_files = list(pdb_dir.rglob('AF*'))
# 只取前50个进行计算
pdb_files = pdb_files[:20]

tmscore_matrix, labels = get_distance_matrix(pdb_files)

# 调用UPGMA（method='average'）
# [0,1] 1最像，0
Z = linkage(1-tmscore_matrix, method='average')
linkage_to_newick(Z, labels,output_file="tree.nwk")

# 绘制树状图
fig, ax = plt.subplots(figsize=(12, 7))
dn = dendrogram(Z, labels=labels, ax=ax)

# 设置标签倾斜与对齐
ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right', fontsize=9)
plt.title('UPGMA Clustering Dendrogram')
plt.ylabel('Distance')
plt.tight_layout()
plt.savefig('dendrogram_new.png', dpi=300)
plt.close()
