import argparse
import glob
import os
from copy import deepcopy
from dataclasses import dataclass
from typing import List, Union

import numpy as np


@dataclass
class GroMol:
    title: str
    num_atoms: int
    resi_num: Union[List[int], np.ndarray]
    resi_name: List[str]
    atom_name: List[str]
    atom_num: Union[List[int], np.ndarray]
    xyz: np.ndarray  # nm
    vel: np.ndarray  # nm/ps
    box_size: np.ndarray = np.zeros(3)

    def __str__(self) -> str:
        gro_contents = []
        gro_contents.append(f"{self.title}")
        gro_contents.append(f"{self.num_atoms}")

        for i in range(self.num_atoms):
            gro_contents.append(
                f"{self.resi_num[i]:5d}"
                f"{self.resi_name[i]:5s}"
                f"{self.atom_name[i]:5s}"
                f"{self.atom_num[i]:5d}"
                f"{self.xyz[i][0]:8.3f}"
                f"{self.xyz[i][1]:8.3f}"
                f"{self.xyz[i][2]:8.3f}"
                f"{self.vel[i][0]:8.4f}"
                f"{self.vel[i][1]:8.4f}"
                f"{self.vel[i][2]:8.4f}"
            )

        gro_contents.append(
            f"{self.box_size[0]:10.6f}"
            f"{self.box_size[1]:10.6f}"
            f"{self.box_size[2]:10.6f}"
        )

        return "\n".join(gro_contents)


def _parse_args():
    # prepare argument parser
    parser = argparse.ArgumentParser(
        description="Convert trajectory in VENUS96 output to gro format "
        "with coordinates and velocities included."
    )
    parser.add_argument(
        "-v",
        "--venus",
        type=str,
        help="VENUS96 output file",
        required=True,
    )
    parser.add_argument(
        "-g",
        "--gro",
        type=str,
        help="Template gro file (optional)",
        default="template.gro",
        required=True
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        help="Output gro file (optional, defaults to traj_idx.gro)",
        default="traj.gro"
    )
    parser.add_argument(
        "-r",
        "--reorder",
        type=str,
        help="Reorder map (optional)",
        required=False,
    )
    parser.add_argument(
        "--no-split",
        help="If specified, write all trajectories in one file. (optional)",
        action="store_true",
    )
    args = parser.parse_args()
    return args


def _check_input_files(args):
    # check VENUS96 output file
    if not os.path.exists(args.venus):
        raise FileNotFoundError(
            f"VENUS96 output file {args.venus} does not exist!"
        )

    if not os.path.exists(args.gro):
        raise FileNotFoundError(
            f"Template gro file {args.gro} does not exist!"
        )

    if args.reorder:
        reorder_map = args.reorder
        if not os.path.exists(reorder_map):
            raise FileNotFoundError(
                f"Reorder map file {reorder_map} does not exist!"
            )


def venus2gro(
    venus_out: str,
    template_gro: str,
    out_gro: str,
    reorder: Union[str, None],
    split: bool = True,
) -> List[str]:
    """Convert VENUS96 output file to gro file.

    Args:
        venus_out (str): Path to VENUS96 output file.
        template_gro (str): Path to template gro file.
        out_gro (str): Path to output gro file(s).
        reorder (str): Path to reorder map text file.
            Defaults to None, i.e. atoms in VENUS96 output file shares the same
            order in template gro file.
        split (bool): Whether to split different trajectories.
            Defaults to True.

    Returns:
        traj_gro_list (List[str]): List of paths to all trajectories extracted
            from VENUS96 output file.
    """

    """Section 0. Load reorder map"""
    if reorder is not None:
        reorder_map = np.loadtxt(reorder, dtype=np.int32) - 1

    """Section 1. Parse template gro"""
    temp_contents = open(template_gro, "r").readlines()

    temp_title = temp_contents.pop(0).rstrip()
    temp_num_atoms = int(temp_contents.pop(0))

    temp_resi_num = []
    temp_resi_name = []
    temp_atom_name = []
    temp_atom_num = []

    for _ in range(temp_num_atoms):
        line = temp_contents.pop(0)
        temp_resi_num.append(int(line[0:5]))
        temp_resi_name.append(line[5:10].rstrip())
        temp_atom_name.append(line[10:15].rstrip())
        temp_atom_num.append(int(line[15:20]))

    temp_box_size = [float(x) for x in temp_contents.pop(0).split()]

    """Section 2. Parse VENUS96 output file"""
    contents = open(venus_out, "r").readlines()
    num_lines = len(contents)
    cur_line_idx = 0
    cur_traj_idx = 0

    anchor_num_atoms = "NUMBER OF ATOMS"
    anchor_masses = "MASSES OF ATOMS"
    anchor_traj = "TRAJECTORY NUMBER"
    anchor_cycle = "THE CYCLE COUNT IS"

    num_atoms = 0
    masses = []
    traj_list = []

    while cur_line_idx < num_lines:
        line = contents[cur_line_idx]

        # parse num_atoms
        if anchor_num_atoms in line:
            num_atoms = int(line.split("=")[-1])
            assert num_atoms == temp_num_atoms, \
                "Number of atoms in VENUS96 output and template gro file" \
                " mismatch, please check these files."

        # parse masses
        if anchor_masses in line:
            cur_line_idx += 2
            line = contents[cur_line_idx]
            masses = [float(x) for x in line.split()]

        # parse traj
        if anchor_traj in line:
            cur_traj_idx = int(line.split()[3])
            if cur_traj_idx > len(traj_list):
                print(f"Trajectory: {cur_traj_idx}")
                traj_list.append([])

        # parse current traj
        if anchor_cycle in line:
            cur_cycle = int(line.split()[4])
            cur_time = float(line.split()[6]) * 1e-02  # 10.0e-14 s -> 1.0 ps
            print(f"Cycle: {cur_cycle} Time: {cur_time:.3f} ps")

            cur_line_idx += 4

            xyz = []
            vel = []

            for i in range(num_atoms):
                cur_line_idx += 1
                line = contents[cur_line_idx]
                arr = line.split()

                # coordinates
                xyz_atom = arr[:3]  # angstrom
                xyz_atom = [0.1 * float(x) for x in xyz_atom]  # angstrom -> nm
                xyz.append(xyz_atom)

                # momenta
                p_atom = arr[3:]  # 0.1 (g/mol) * (nm/ps)
                v_atom = [10.0 * float(p) / masses[i]
                          for p in p_atom]  # (nm/ps)
                vel.append(v_atom)

            if reorder is not None:
                _xyz = deepcopy(xyz)
                _vel = deepcopy(vel)

                for (old_idx, new_idx) in enumerate(reorder_map):
                    xyz[old_idx] = _xyz[new_idx]
                    vel[old_idx] = _vel[new_idx]

            # create GroMol
            title = f"{temp_title} Traj {cur_traj_idx}, t= {cur_time:.4f} ps"
            gro_mol = GroMol(
                title=title,
                num_atoms=num_atoms,
                resi_num=temp_resi_num,
                resi_name=temp_resi_name,
                atom_name=temp_atom_name,
                atom_num=temp_atom_num,
                xyz=np.array(xyz),
                vel=np.array(vel),
                box_size=np.array(temp_box_size),
            )

            traj_list[-1].append(gro_mol)

        cur_line_idx += 1

    """Section 3. Write trajectories to gro file(s)"""
    def _auto_backup_file(file):
        if not os.path.exists(file):
            return

        file_dir = os.path.dirname(file)
        file_name = os.path.basename(file)
        old_file_pattern = f"#{file_name.replace('.gro', '')}.*.gro#"
        old_files = sorted(glob.glob(os.path.join(
            file_dir,
            old_file_pattern,
        )))
        num_old_files = len(old_files)
        new_backup_file_name = \
            f"#{file_name.replace('.gro', '')}.{num_old_files + 1}.gro#"
        os.rename(file, new_backup_file_name)

    traj_gro_list = []

    if not split:
        traj_list = np.concatenate(traj_list).reshape(1, -1).tolist()

    for (traj_idx, traj) in enumerate(traj_list):
        traj_gro_name = out_gro.replace(".gro", f"_{traj_idx}.gro")
        _auto_backup_file(traj_gro_name)
        with open(traj_gro_name, "w") as f:
            for frame in traj:
                f.writelines(f"{str(frame)}\n")

        traj_gro_list.append(traj_gro_name)

    return traj_gro_list


if __name__ == "__main__":
    # parse args
    args = _parse_args()

    # check input and output files
    _check_input_files(args)

    print(args)

    # parse output
    traj_list = venus2gro(
        venus_out=args.venus,
        template_gro=args.gro,
        out_gro=args.output,
        reorder=args.reorder,
        split=(not args.no_split),
    )
