# 运行dssp计算蛋白质的二级结构（α-helix, β-sheet 等），除了需要安装biopython的包，还需要安装dssp,他其实调用外部命令行工具 mkdssp 来做计算的。
# conda install -c salilab dssp

from Bio.PDB import PDBParser, DSSP

pdb_file = "protein.pdb"
parser = PDBParser(QUIET=True)
structure = parser.get_structure("protein", pdb_file)

model = structure[0]  # 第一模型
dssp = DSSP(model, pdb_file)

for key in list(dssp.keys())[:5]:
    res_id = key[1][1]
    ss = dssp[key][2]  # 二级结构类型
    print(f"残基 {res_id} 二级结构: {ss}")