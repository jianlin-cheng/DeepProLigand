"""
#################################################################
# Runs fast relax
#Input: averaged structure in .pdb format
#Output: refined structure in .pdb format
#Author: Nabin Giri (ngzvh@missouri.edu)
#################################################################

"""

from pyrosetta import *
init()

def fastrelax(input_pdb):
	pose = pose_from_pdb(input_pdb)
	scorefxn = rosetta.core.scoring.ScoreFunctionFactory.create_score_function("ref2015_cst.wts")
	repeat = 5
	frelax = rosetta.protocols.relax.FastRelax(scorefxn,repeat)
	frelax.apply(pose)
	output_name = input_pdb.split('_')[0]
	out = f'{output_name}_relax.pdb'
	pose.dump_pdb(out)


if __name__ =="__main__":
	dir = '/data/'
	final_7770 = f"{dir}/7770_final.pdb"
	fastrelax(final_7770)
	final_22898 = f"{dir}/22898_final.pdb"
	fastrelax(final_22898)
	final_30210 = f"{dir}/30210_final.pdb"
	fastrelax(final_30210)
	print("############ NG: Fast Relax is Complete ############")
	
	