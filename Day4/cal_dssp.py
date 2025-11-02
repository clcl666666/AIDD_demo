from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP
from Bio.Data.PDBData import residue_sasa_scales
import pandas as pd
from collections import defaultdict

protein_letters_1to3 = {
    "A": "Ala",
    "C": "Cys",
    "D": "Asp",
    "E": "Glu",
    "F": "Phe",
    "G": "Gly",
    "H": "His",
    "I": "Ile",
    "K": "Lys",
    "L": "Leu",
    "M": "Met",
    "N": "Asn",
    "P": "Pro",
    "Q": "Gln",
    "R": "Arg",
    "S": "Ser",
    "T": "Thr",
    "V": "Val",
    "W": "Trp",
    "Y": "Tyr",
}
# Used for relative accessibility

residue_max_sasa = residue_sasa_scales["Sander"]
# rel_acc = acc / self.residue_max_acc[resname]
pdb_f = '/fastone/users/shiny/flask/results/Humanlization/20250106/CM000.pdb'

# 解析PDB文件
p = PDBParser()
structure = p.get_structure("protein", pdb_f)  # 替换为您的PDB文件路径
model = structure[0]  # 假设我们只处理第一个模型

# 使用DSSP分析二级结构和SASA
dssp = DSSP(model, pdb_f, dssp="mkdssp")  # 确保mkdssp在您的系统路径中

# 创建一个空的DataFrame来存储结果
sasa = defaultdict(list)

# 遍历DSSP结果并提取SASA值
for key in dssp.keys():
    value = dssp[key]
    chain_id, res_id = key
    dssp_index, amino_acid, sec_struct, rel_asa, phi, psi = value[:6]  # 获取DSSP索引、氨基酸、二级结构、相对ASA、PHI、PSI
    sasa['Chain'].append(chain_id)
    sasa['Residue_ID'].append(res_id[1])
    sasa['Residue_Name'].append(amino_acid)
    sasa['SASA_relative'].append(rel_asa)
    max_ss = residue_max_sasa[protein_letters_1to3[amino_acid].upper()]
    sasa['SASA'].append(rel_asa*max_ss)
sasa_df = pd.DataFrame(sasa)

# 将DataFrame输出为CSV文件
sasa_df.to_csv('sasa_results.csv', index=False)