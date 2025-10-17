import pandas as pd

def read_rosetta_scorefile_per_res(path):
    with open(path) as f:
        lines = [line for line in f if line.startswith("SCORE:")]
    
    # 以空格分隔，跳过 "SCORE:" 列
    data = [line.strip().split()[1:] for line in lines]
    columns = data[0]
    df = pd.DataFrame(data[1:], columns=columns)

    # 将数值列转为 float
    for col in df.columns[2:-1]:
        df[col] = df[col].astype(float)
    
    # df['hydrophobic'] = df['fa_atr']+df['fa_rep']*0.55 + df['fa_sol']*0.9375+df['fa_intra_rep']*0.005
    # df['electrostatic'] = df['fa_elec'] * 0.875
    # df['proline ring'] = df['pro_close'] * 1.25
    # df['hbond'] = 1.17*(df['hbond_sr_bb']+df['hbond_lr_bb']+df['hbond_bb_sc'])+df['hbond_sc']*1.1
    # df['disulfide'] = df['dslf_fa13']
    # df['torsion_conflicts'] = df['omega']*0.625+df['fa_dun']*0.7 +df['p_aa_pp']*0.4+df['yhh_planarity']*0.625+df['rama_prepro']*0.25

    df['hydrophobic'] = df['fa_atr']+df['fa_rep'] + df['fa_sol']+df['fa_intra_rep']
    df['electrostatic'] = df['fa_elec']
    df['proline ring'] = df['pro_close']
    df['hbond'] = df['hbond_sr_bb']+df['hbond_lr_bb']+df['hbond_bb_sc']+df['hbond_sc']
    df['disulfide'] = df['dslf_fa13']
    df['torsion_conflicts'] = df['omega']+df['fa_dun']+df['p_aa_pp']+df['yhh_planarity']+df['rama_prepro']

    return df[['pdb_id','hydrophobic','electrostatic','proline ring','hbond','disulfide','torsion_conflicts','score']]

def read_rosetta_scorefile_energy_breakdown(path):
    with open(path) as f:
        lines = [line for line in f if line.startswith("SCORE:")]
    
    # 以空格分隔，跳过 "SCORE:" 列
    data = [line.strip().split()[1:] for line in lines]
    columns = data[0]
    df = pd.DataFrame(data[1:], columns=columns)
    df = df.query("restype2 != 'onebody'")
    for col in df.columns[7:-1]:
        df[col] = df[col].astype(float)

    # df['hydrophobic'] = df['fa_atr']+df['fa_rep']*0.55 + df['fa_sol']*0.9375+df['fa_intra_rep']*0.005
    # df['electrostatic'] = df['fa_elec'] * 0.875
    # df['proline ring'] = df['pro_close'] * 1.25
    # df['hbond'] = 1.17*(df['hbond_sr_bb']+df['hbond_lr_bb']+df['hbond_bb_sc'])+df['hbond_sc']*1.1
    # df['disulfide'] = df['dslf_fa13']
    # df['torsion_conflicts'] = df['omega']*0.625+df['fa_dun']*0.7 +df['p_aa_pp']*0.4+df['yhh_planarity']*0.625+df['rama_prepro']*0.25

    df['hydrophobic'] = df['fa_atr']+df['fa_rep'] + df['fa_sol']+df['fa_intra_rep']
    df['electrostatic'] = df['fa_elec']
    df['proline ring'] = df['pro_close']
    df['hbond'] = df['hbond_sr_bb']+df['hbond_lr_bb']+df['hbond_bb_sc']+df['hbond_sc']
    df['disulfide'] = df['dslf_fa13']
    df['torsion_conflicts'] = df['omega']+df['fa_dun']+df['p_aa_pp']+df['yhh_planarity']+df['rama_prepro']

    return df[['pdbid1','restype1','pdbid2','restype2','hydrophobic','electrostatic','proline ring','hbond','disulfide','torsion_conflicts','total']]



# 每一个氨基酸的总共的
file_path = '/fastone/users/linchen/src/mycode/rosetta_InterfaceAnalyzer/per_residue_energy/per_res.sc'
residue_df = read_rosetta_scorefile_per_res(file_path)
# 氨基酸和氨基酸之间的相互作用
file_path = '/fastone/users/linchen/src/mycode/rosetta_InterfaceAnalyzer/per_residue_energy/energy_breakdown.sc'
break_residue_df = read_rosetta_scorefile_energy_breakdown(file_path)
