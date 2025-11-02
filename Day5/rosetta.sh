/fastone/users/linchen/src/rosetta_bin_linux_2018.09.60072_bundle/main/source/bin/score_jd2.static.linuxgccrelease -in:file:s PGF-CH1.pdb -ignore_waters true -out:pdb -ignore_unrecognized_res  

/fastone/users/linchen/src/rosetta_bin_linux_2018.09.60072_bundle/main/source/bin/relax.static.linuxgccrelease -s PGF-CH1_0001.pdb -out:suffix _relaxed -relax:default_repeats 5

# 获得每一个氨基酸的能量
/fastone/users/linchen/src/rosetta_bin_linux_2018.09.60072_bundle/main/source/bin/per_residue_energies.static.linuxgccrelease \
-in:file:s PGF-CH1_0001_relaxed_0001.pdb \
-out:file:silent per_res.sc

# 拆分一个氨基酸的能量为，onebody和与其他氨基酸之间的相互作用
/fastone/users/linchen/src/rosetta_bin_linux_2018.09.60072_bundle/main/source/bin/residue_energy_breakdown.static.linuxgccrelease \
-in:file:s PGF-CH1_0001_relaxed_0001.pdb \
-out:file:silent energy_breakdown.sc

# 运行flex_ddG,会自动运行inputs文件里面的所有，需要把Rosetta的位置改了。
cd flex_ddG_tutorial
python run_example_1.py
# 运行分析脚本,非常好，这一行是老版本 pandas 的写法。新的版本已经没有append这个命令了。建立一个新的环境，python=3.9,安装一个老版本的pandas。pip install pandas==1.5.3 numpy==1.23.5

python analyze_flex_ddG.py output
# 运行脚本后，终端中会输出几个表格，包括：
# 野生型界面结合能 ΔG（wt_dG）
# 突变体界面结合能 ΔG（mut_dG）
# 突变后结合能变化 ΔΔG（mutant ΔΔG）
# 这些结果也会被写入到
# analysis_output 目录下的 .csv 文件中。
