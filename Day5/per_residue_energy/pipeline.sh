
/fastone/users/linchen/src/rosetta_bin_linux_2018.09.60072_bundle/main/source/bin/score_jd2.static.linuxgccrelease \
-in:file:s PGF-CH1.pdb \
-ignore_waters true -out:pdb -ignore_unrecognized_res  

/fastone/users/linchen/src/rosetta_bin_linux_2018.09.60072_bundle/main/source/bin/relax.static.linuxgccrelease -s PGF-CH1_0001.pdb \
-out:suffix _relaxed -relax:default_repeats 5

# 获得每一个氨基酸的能量
/fastone/users/linchen/src/rosetta_bin_linux_2018.09.60072_bundle/main/source/bin/per_residue_energies.static.linuxgccrelease \
-in:file:s PGF-CH1_0001_relaxed_0001.pdb \
-out:file:silent per_res.sc

# 拆分一个氨基酸的能量为，onebody和与其他氨基酸之间的相互作用
/fastone/users/linchen/src/rosetta_bin_linux_2018.09.60072_bundle/main/source/bin/residue_energy_breakdown.static.linuxgccrelease \
-in:file:s PGF-CH1_0001_relaxed_0001.pdb \
-out:file:silent energy_breakdown.sc
