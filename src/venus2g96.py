import argparse
import glob
import os
from copy import deepcopy
from dataclasses import dataclass
from typing import List, Tuple, Union

import numpy as np


@dataclass
class G96Mol:
    title: str
    timestep: Tuple[int, float]
    position: np.array
    velocity: np.array
    box: np.array = np.zeros(3)

    def __str__(self) -> str:
        g96_contents = []

        # title
        g96_contents.append("TITLE")
        g96_contents.append(self.title)
        g96_contents.append("END")

        # timestep
        g96_contents.append("TIMESTEP")
        g96_contents.append(f"{self.timestep[0]:15d}{self.timestep[1]:15.6f}")
        g96_contents.append("END")

        # position
        g96_contents.append("POSITIONRED")

        for pos in self.position:
            g96_contents.append(f"{pos[0]:15.9f}{pos[1]:15.9f}{pos[2]:15.9f}")

        g96_contents.append("END")

        # velocity
        g96_contents.append("VELOCITYRED")

        for vel in self.velocity:
            g96_contents.append(f"{vel[0]:15.9f}{vel[1]:15.9f}{vel[2]:15.9f}")

        g96_contents.append("END")

        # box
        g96_contents.append("BOX")
        g96_contents.append(
            f"{self.box[0]:15.9f}{self.box[1]:15.9f}{self.box[2]:15.9f}"
        )
        g96_contents.append("END")

        return "\n".join(g96_contents)


def _parse_args():
    # prepare argument parser
    parser = argparse.ArgumentParser(
        description="Convert trajectory in VENUS96 output to g96 format "
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
        "-o",
        "--output",
        type=str,
        help="Output g96 file (optional, defaults to traj_idx.g96)",
        default="traj.g96"
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

    if args.reorder:
        reorder_map = args.reorder
        if not os.path.exists(reorder_map):
            raise FileNotFoundError(
                f"Reorder map file {reorder_map} does not exist!"
            )


def venus2g96(
    venus_out: str,
    out_g96: str,
    reorder: Union[str, None],
    split: bool = True,
) -> List[str]:
    """Convert VENUS96 output file to g96 file.

    Args:
        venus_out (str): Path to VENUS96 output file.
        out_g96 (str): Path to output g96 file(s).
        reorder (str): Path to reorder map text file.
        split (bool): Whether to split different trajectories.
            Defaults to True.

    Returns:
        traj_g96_list (List[str]): List of paths to all trajectories extracted
            from VENUS96 output file.
    """

    """Section 0. Load reorder map"""
    if reorder is not None:
        reorder_map = np.loadtxt(reorder, dtype=np.int32) - 1

    """Section 1. Parse VENUS96 output file"""
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

        # parse masses
        if anchor_masses in line:
            cur_line_idx += 2
            line = contents[cur_line_idx]
            masses = np.array([float(x) for x in line.split()])

        # parse traj
        if anchor_traj in line:
            cur_traj_idx = int(line.split()[3])
            if cur_traj_idx > len(traj_list):
                print(f"Trajectory: {cur_traj_idx}")
                traj_list.append([])

        # parse current traj
        if anchor_cycle in line:
            cur_cycle = int(line.split()[4])
            # integration stepsize in units of 10-14 sec
            # 1.0e-14 s -> 1.0e-02 ps
            cur_time = float(line.split()[6]) * 1.0e-02
            print(f"Cycle: {cur_cycle} Time: {cur_time:.3f} ps")

            cur_line_idx += 4

            xyz = []
            vel = []

            for i in range(num_atoms):
                cur_line_idx += 1
                line = contents[cur_line_idx]
                arr = line.split()

                # coordinates in units of Angstrom
                # 1.0 Angstrom -> 0.1 nm
                xyz_atom = 0.1 * np.array([float(x) for x in arr[:3]])
                xyz.append(xyz_atom)

                # momenta in units of amu * Angstrom / (1.0e-14 s)
                # 1.0 amu * Angstrom / (1.0e-14 s) -> 10.0 (g/mol) * (nm/ps)
                p_atom = np.array([float(x) for x in arr[3:]])
                v_atom = 10.0 * p_atom / masses[i]
                vel.append(v_atom)

            if reorder is not None:
                _xyz = deepcopy(xyz)
                _vel = deepcopy(vel)

                for (old_idx, new_idx) in enumerate(reorder_map):
                    xyz[old_idx] = _xyz[new_idx]
                    vel[old_idx] = _vel[new_idx]

            # create G96Mol
            g96_mol = G96Mol(
                title="",
                timestep=[cur_cycle, cur_time],
                position=xyz,
                velocity=vel,
            )

            traj_list[-1].append(g96_mol)

        cur_line_idx += 1

    """Section 2. Write trajectories to g96 file(s)"""
    def _auto_backup_file(file):
        if not os.path.exists(file):
            return

        file_dir = os.path.dirname(file)
        file_name = os.path.basename(file)
        old_file_pattern = f"#{file_name.replace('.g96', '')}.*.g96#"
        old_files = sorted(glob.glob(os.path.join(
            file_dir,
            old_file_pattern,
        )))
        num_old_files = len(old_files)
        new_backup_file_name = \
            f"#{file_name.replace('.g96', '')}.{num_old_files + 1}.g96#"
        os.rename(file, new_backup_file_name)

    traj_g96_list = []

    if not split:
        traj_list = np.concatenate(traj_list).reshape(1, -1).tolist()

    for (traj_idx, traj) in enumerate(traj_list):
        traj_g96_name = out_g96.replace(".g96", f"_{traj_idx + 1}.g96")
        _auto_backup_file(traj_g96_name)
        with open(traj_g96_name, "w") as f:
            for frame in traj:
                f.writelines(f"{str(frame)}\n")

            traj_g96_list.append(traj_g96_name)

    return traj_g96_list


if __name__ == "__main__":
    # parse args
    args = _parse_args()

    # check input and output files
    _check_input_files(args)

    # parse output
    traj_list = venus2g96(
        venus_out=args.venus,
        out_g96=args.output,
        reorder=args.reorder,
        split=(not args.no_split),
    )
