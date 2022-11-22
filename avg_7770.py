"""
#################################################################
#Algorithm 1 of the paper. Averages predicted structure and reference structure for cryo-EM map 7770.
#Input: superimposed protein structure (predicted structure superimposed relative to reference structure)
#Output: average structure saved in .pdb format
#Author: Nabin Giri (ngzvh@missouri.edu)
#################################################################

"""
import os
from pyrosetta import *

init()


def get_residues(pdb_file, output_dir):
    threshold = 0.5
    alt_threshold = 0.3
    pose = pyrosetta.pose_from_pdb(pdb_file)
    print(f"The number of residues are {pose.total_residue()}")
    print(f"The number of chain are {pose.num_chains()}")

    pdb_name = os.path.basename(pdb_file)
    target_name = pdb_name.split('_')[0]

    file_path = f"{output_dir}/avg_{target_name}.pdb"
    with open(file_path, 'w') as fi:
        fi.write('HEADER')
        fi.write('\t')
        fi.write('PROTEIN STRUCTURE')
        fi.write('\n')
        fi.write('TITLE')
        fi.write('\t')
        fi.write('AVERAGE OF SUPERIMPOSED 7770')
        fi.write('\n')
        fi.write('AUTHOR')
        fi.write('\t')
        fi.write('UM BML')
        fi.write('\n')

    start_A1 = pose.conformation().chain_begin(1)
    end_A1 = pose.conformation().chain_end(1)

    start_A2 = pose.conformation().chain_begin(2)
    end_A2 = pose.conformation().chain_end(2)

    start_B1 = pose.conformation().chain_begin(3)
    end_B1 = pose.conformation().chain_end(3)

    start_B2 = pose.conformation().chain_begin(5)
    end_B2 = pose.conformation().chain_end(5)

    start_C1 = pose.conformation().chain_begin(6)
    end_C1 = pose.conformation().chain_end(6)

    start_C2 = pose.conformation().chain_begin(8)
    end_C2 = pose.conformation().chain_end(8)

    start_D1 = pose.conformation().chain_begin(9)
    end_D1 = pose.conformation().chain_end(9)

    print("############ Entering A Chain ############")
    atoms = ['N', 'CA', 'C', 'O']
    count = 0
    count_avg_A = 0
    chain = 'A'

    for i in range(start_A1, end_A1 + 1):
        for j in range(start_A2, end_A2 + 1):
            if (pose.residue(j).is_ligand()):
                continue
            d = distance(pose, i, j)
            if d < threshold:
                count_avg_A = count_avg_A + 1
                for atm in range(len(atoms)):
                    count = count + 1
                    avg = average(pose, i, j, atoms[atm])
                    save_avg(pose, file_path, i, count, avg, atoms[atm], chain)
                break
            elif j == end_A2:
                for k in range(start_B2, end_B2 + 1):
                    d = distance(pose, i, k)
                    if d < threshold:
                        count_avg_A = count_avg_A + 1
                        for atm in range(len(atoms)):
                            count = count + 1
                            avg = average(pose, i, k, atoms[atm])
                            save_avg(pose, file_path, i, count, avg, atoms[atm], chain)
                        break
                    elif k == end_B2:
                        for l in range(start_C2, end_C2 + 1):
                            d = distance(pose, i, l)
                            if d < threshold:
                                count_avg_A = count_avg_A + 1
                                for atm in range(len(atoms)):
                                    count = count + 1
                                    avg = average(pose, i, l, atoms[atm])
                                    save_avg(pose, file_path, i, count, avg, atoms[atm], chain)
                                break
                            elif l == end_C2:
                                for atm in range(len(atoms)):
                                    count = count + 1
                                    save(pose, file_path, i, count, atoms[atm], chain)

    print(f"Number of residues averaged for Chain A are {count_avg_A}")
    print(f"Model-2: Reference structure residue difference of A is {end_A2 - start_A2}")

    with open(file_path, 'a') as fi:
        fi.write('TER')
        fi.write('\n')

    print("############ Entering B Chain ############")
    chain = 'B'
    count_avg_B = 0
    for i in range(start_B1, end_B1 + 1):
        for j in range(start_A2, end_A2 + 1):
            if (pose.residue(j).is_ligand()):
                continue
            d = distance(pose, i, j)
            if d < threshold:
                count_avg_B = count_avg_B + 1
                for atm in range(len(atoms)):
                    count = count + 1
                    avg = average(pose, i, j, atoms[atm])
                    save_avg(pose, file_path, i, count, avg, atoms[atm], chain)
                break
            elif j == end_A2:
                for k in range(start_B2, end_B2 + 1):
                    d = distance(pose, i, k)
                    if d < threshold:
                        count_avg_B = count_avg_B + 1
                        for atm in range(len(atoms)):
                            count = count + 1
                            avg = average(pose, i, k, atoms[atm])
                            save_avg(pose, file_path, i, count, avg, atoms[atm], chain)
                        break
                    elif k == end_B2:
                        for l in range(start_C2, end_C2 + 1):
                            d = distance(pose, i, l)
                            if d < threshold:
                                count_avg_B = count_avg_B + 1
                                for atm in range(len(atoms)):
                                    count = count + 1
                                    avg = average(pose, i, l, atoms[atm])
                                    save_avg(pose, file_path, i, count, avg, atoms[atm], chain)
                                break
                            elif l == end_C2:
                                for atm in range(len(atoms)):
                                    count = count + 1
                                    save(pose, file_path, i, count, atoms[atm], chain)

    print(f"Number of residues averaged for Chain B is {count_avg_B}")
    print(f"Model-2: Reference structure residue difference of B is {end_B2 - start_B2}")

    with open(file_path, 'a') as fi:
        fi.write('TER')
        fi.write('\n')

    print("############ Entering C Chain ############")
    chain = 'C'
    count_avg_C = 0
    for i in range(start_C1, end_C1 + 1):
        for j in range(start_A2, end_A2 + 1):
            if (pose.residue(j).is_ligand()):
                continue
            d = distance(pose, i, j)
            if d < alt_threshold:
                count_avg_C = count_avg_C + 1
                for atm in range(len(atoms)):
                    count = count + 1
                    avg = average(pose, i, j, atoms[atm])
                    save_avg(pose, file_path, i, count, avg, atoms[atm], chain)
                break
            elif j == end_A2:
                for k in range(start_B2, end_B2 + 1):
                    d = distance(pose, i, k)
                    if d < alt_threshold:
                        count_avg_C = count_avg_C + 1
                        for atm in range(len(atoms)):
                            count = count + 1
                            avg = average(pose, i, k, atoms[atm])
                            save_avg(pose, file_path, i, count, avg, atoms[atm], chain)
                        break
                    elif k == end_B2:
                        for l in range(start_C2, end_C2 + 1):
                            d = distance(pose, i, l)
                            if d < alt_threshold:
                                count_avg_C = count_avg_C + 1
                                for atm in range(len(atoms)):
                                    count = count + 1
                                    avg = average(pose, i, l, atoms[atm])
                                    save_avg(pose, file_path, i, count, avg, atoms[atm], chain)
                                break
                            elif l == end_C2:
                                for atm in range(len(atoms)):
                                    count = count + 1
                                    save(pose, file_path, i, count, atoms[atm], chain)

    print(f"Number of residues averaged for Chain C is {count_avg_C}")
    print(f"Model-2: Reference structure residue difference of C is {end_C2-start_C2}")

    with open(file_path, 'a') as fi:
        fi.write('TER')
        fi.write('\n')

    print("############ Entering D Chain ############")
    chain = 'D'
    count_avg_D = 0
    for i in range(start_D1, end_D1 + 1):
        for j in range(start_A2, end_A2 + 1):
            if (pose.residue(j).is_ligand()):
                continue
            d = distance(pose, i, j)
            if d < threshold:
                count_avg_D = count_avg_D + 1
                for atm in range(len(atoms)):
                    count = count + 1
                    avg = average(pose, i, j, atoms[atm])
                    save_avg(pose, file_path, i, count, avg, atoms[atm], chain)
                break
            elif j == end_A2:
                for k in range(start_B2, end_B2 + 1):
                    d = distance(pose, i, k)
                    if d < threshold:
                        count_avg_D = count_avg_D + 1
                        for atm in range(len(atoms)):
                            count = count + 1
                            avg = average(pose, i, k, atoms[atm])
                            save_avg(pose, file_path, i, count, avg, atoms[atm], chain)
                        break
                    elif k == end_B2:
                        for l in range(start_C2, end_C2 + 1):
                            d = distance(pose, i, l)
                            if d < threshold:
                                count_avg_D = count_avg_D + 1
                                for atm in range(len(atoms)):
                                    count = count + 1
                                    avg = average(pose, i, l, atoms[atm])
                                    save_avg(pose, file_path, i, count, avg, atoms[atm], chain)
                                break
                            elif l == end_C2:
                                for atm in range(len(atoms)):
                                    count = count + 1
                                    save(pose, file_path, i, count, atoms[atm], chain)

    print(f"Number of residues averaged for Chain D is {count_avg_D}")

    with open(file_path, 'a') as fi:
        fi.write('END')


def distance(pose, res_i, res_j):
    xyz_i = pose.residue(res_i).xyz('CA')
    xyz_j = pose.residue(res_j).xyz('CA')
    dist = (xyz_i - xyz_j).norm()
    return dist


def average(pose, res_i, res_j, atm):
    xyz_i = pose.residue(res_i).xyz(atm)
    xyz_j = pose.residue(res_j).xyz(atm)
    avg = list()
    for i in range(len(xyz_i)):
        avg.append((xyz_i[i] + xyz_j[i]) / 2)
    return avg


def save(pose, file_path, res_i, count, atm, chain):
    with open(file_path, 'a') as fi:
        fi.write('ATOM')
        fi.write('  ')
        fi.write(str(count).rjust(5))
        fi.write('  ')
        fi.write(atm.ljust(4))
        fi.write(pose.residue(res_i).name()[0:3].rjust(3))
        fi.write(' ')
        fi.write(chain)
        fi.write(str(pose.pdb_info().number(res_i)).rjust(4))
        fi.write('    ')
        xyz_atm = pose.residue(res_i).xyz(atm)
        for i in range(len(xyz_atm)):
            fi.write(str(round(xyz_atm[i], 3)).rjust(8))
        fi.write(str(1.00).rjust(5))
        fi.write(str(0.00).rjust(5))
        fi.write('           ')
        fi.write(atm[0:1].rjust(1))
        fi.write('  ')
        fi.write('\n')


def save_avg(pose, file_path, res_i, count, avg, atm, chain):
    with open(file_path, 'a') as fi:
        fi.write('ATOM')
        fi.write('  ')
        fi.write(str(count).rjust(5))
        fi.write('  ')
        fi.write(atm.ljust(4))
        fi.write(pose.residue(res_i).name()[0:3].rjust(3))
        fi.write(' ')
        fi.write(chain)
        fi.write(str(pose.pdb_info().number(res_i)).rjust(4))
        fi.write('    ')
        for i in range(len(avg)):
            fi.write(str(round(avg[i], 3)).rjust(8))
        fi.write(str(1.00).rjust(5))
        fi.write(str(0.00).rjust(5))
        fi.write('           ')
        fi.write(atm[0:1].rjust(1))
        fi.write('  ')
        fi.write('\n')


if __name__ == "__main__":
    pdb_file = '/data/7770_superimposed.pdb'
    output_dir = '/data/averaged'
    get_residues(pdb_file, output_dir)
    print("############ NG: Complete ############")
