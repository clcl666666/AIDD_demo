
==== Quantifying target hydrophobicity ====

Start by using a molecular viewing tool to extract just your target protein of interest. Take care to remove solvent molecules and binding partners as these will affect the calculation. 
Save this as target.pdb. And run the following command:

$ROSETTA/bin/rosetta_scripts -parser:protocol per_res_sap.xml -beta_nov16 -renumber_pdb -s target.pdb
/fastone/users/linchen/src/rosetta.binary.ubuntu.release-408/main/source/bin/rosetta_scripts.static.linuxgccrelease -parser:protocol per_res_sap.xml -beta_nov16 -renumber_pdb -s 7xwo_part.pdb

This will output a target_0001.pdb with per-residue SAP score information. Next, open PyMOL in the same folder as your target_0001.pdb (you may need to cd within pymol),
 and open target_0001.pdb. Set your favorite viewing settings (the authors would use show cartoon; show sticks; hide (hydro)) and then run the following command from within pymol.

run per_res_sap.py

That command will label each CA with the SAP score of the sidechain and will color them more red if they have higher values. 
(If the pymol command failed, you can manually open the pdb with a text editor and look at the numbers.) 
In general, you want to select a binding site with enough possible ΔSAP. The way to think about this is that if you can block 100% of the SASA of a residue,
 you accumulate the SAP score. However, if you can only block half of a residue, you only get half. 
Cavities behind sidechains can artificially inflate the numbers here because you won't be able to cover the backside of the residue. Be careful around low-density sites.

You can use the Extended Data figure with x-axis = Target Delta SAP and y-axis = Target Success rate as a guide to these numbers. 
We've been a bit more ambitious since these early experiments and as such, these are the cutoffs that seems to make sense:

SAP Score:
 0 -  6 -- No chance of binding
 6 -  9 -- Probably won't work, but you might get lucky
 9 - 12 -- The limit. Definitely order a lot of designs
12 - 15 -- Success is likely. You should still order a lot
15 - 18 -- Easier target. Getting a binder by ordering 100 designs might be possible
18+     -- Easy target. Be mindful you don't create overly non-polar designs that aggregate


Once you are done, delete the target_0001.pdb produced here to avoid confusion in later steps.
