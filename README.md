# venus2gro

Extract coordinates and momenta from VENUS96 output and write trajectories in GROMACS gro format.

Author: mizu-bai

## Requirements

- python 3.7 and above
- numpy

## Usage

```bash
$ python3 venus2gro.py -h
usage: venus2gro.py [-h] -v VENUS -g GRO [-o OUTPUT] [-r REORDER] [--no-split]

Convert trajectory in VENUS96 output to gro format with coordinates and velocities included.

optional arguments:
  -h, --help            show this help message and exit
  -v VENUS, --venus VENUS
                        VENUS96 output file
  -g GRO, --gro GRO     Template gro file (optional)
  -o OUTPUT, --output OUTPUT
                        Output gro file (optional, defaults to traj_idx.gro)
  -r REORDER, --reorder REORDER
                        Reorder map (optional)
  --no-split            If specified, write all trajectories in one file. (optional)
```

## Example

In folder `example/`, there is a QCT trajectory of methane molecule.

- `ch4.dt5`: VENUS96 input file.
- `ch4.out`: VENUS96 output file, containing 5 trajectories.
- `template.gro`: Template gro file, `venus2gro.py` will read title, number of atoms, residue numbers, residue names, atom names, atom numbers and box size stored in it.
- `reorder.txt`: In VENUS96 output, the order of atoms is `H H H H C`, while we want to use order in `template.gro`, which is `C H H H H`. Thus, a reorder file should be used, the first line is `5`, indicating that the 5-th atom (`C`) in VENUS96 output should be the 1-st atom in gro file. If the atom order is same in VENUS96 and `template.gro`, this file and `-r/--reorder` option are not necessary.

Run this command, then 5 trajectories in gro format will be created.

```bash
$ python3 ../src/venus2gro.py -v ch4.out -g template.gro -o ch4_traj.gro -r reorder.txt
```

## Reference

(1) Hase, W. L.; Duchovic, R. J.; Hu, X.; Komornicki, A.; Lim, K. F.; Lu, D.-H.; Peslherbe, G. H.; Swamy, K. N.; Vande Linde, S. R.; Varandas, A., Wang, H.; Wolf, R. J. VENUS96: A general chemical dynamics computer program, _Quantum Chemical Program Exchange (QCPE) Bulletin_, **1996**, _16_ (4), 671. https://www.depts.ttu.edu/chemistry/Venus/index.php

(2) Abraham, M. J.; Murtola, T.; Schulz, R.; Páll, S.; Smith, J. C.; Hess, B.; Lindahl, E. GROMACS: High Performance Molecular Simulations through Multi-Level Parallelism from Laptops to Supercomputers. _SoftwareX_ **2015**, _1–2_, 19–25. https://doi.org/10.1016/j.softx.2015.06.001.
