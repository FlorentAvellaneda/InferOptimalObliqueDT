# Learning Optimal Oblique Decision Trees with (Max)SAT

The provided source code is designed for compilation and execution on GNU/Linux systems.


## File organization

- `datasets/` : contains the datasets used in the experiments
- `nuwls/` : contains the binary of the NuWLS solver
- `main/` : contains the source code of the project
- `lib/` : contains the source code of the libraries used in the project


## Dependencies

### Not included
- g++ (`sudo apt install g++` on Ubuntu)
- cmake (`sudo apt install cmake` on Ubuntu)
- Z3 (`sudo apt install libz3-dev` on Ubuntu)
- libsvm (`sudo apt install libsvm-dev` on Ubuntu)
- libgmp (`sudo apt install libgmp-dev` on Ubuntu)
- libgmpxx (`sudo apt install libgmpxx-dev` on Ubuntu)
- xdot (`sudo apt install xdot` on Ubuntu) : to see the learned tree

### Included
- Cadical : https://github.com/arminbiere/cadical
- EvalMaxSAT : https://github.com/FlorentAvellaneda/EvalMaxSAT
- NuWLS : https://github.com/filyouzicha/NuWLS
- CLI11 : https://github.com/CLIUtils/CLI11
- fmt::fmt : https://github.com/fmtlib/fmt


## Build

```bash
mkdir build
cd build
cmake ..
make
```


## Usage


**Usage:** `./build/main/LearnODT [OPTIONS] CSV_file`

**Options:**

  `-h,--help`                   Print this help message and exit

  `--exp UINT`                  Experiment for Fig2 and Fig3 (1, 2, 3 or 4)

  `--seed UINT`                 seed (default:0, means random seed)

  `-v INT`                      vebosity (default: 1)

  `-d UINT`                     depth (default: 2)

  `--exact`                     Exact classification

  `--Z3`                        Use Z3 encodage

  `--kcross INT`                k cross validation (default: 0 (no cross validation), -1 for one out cross validation)

  `--shati`                     Use Shati et al. encodage

  `--noshow`                    Don't show the tree

  `--solver TEXT`               External solver cmd (internal solver by default)

## Examples

From the root directory of the project, you can run the following examples (and do not forget to install xdot to see the learned tree).

Learn an optimal oblique decision tree on the `datasets/mouse.csv` dataset with a depth of 2 and our MaxSAT encoding:
```bash
./build/main/LearnODT datasets/mouse.csv -d 2
```

Learn an optimal oblique decision tree on the `datasets/mouse.csv` dataset with a depth of 2 and the Z3 encoding:
```bash
./build/main/LearnODT datasets/mouse.csv -d 2 --Z3
```

Learn an optimal oblique decision tree on the `datasets/mouse.csv` dataset with a depth of 2 and the Shati et al. encoding:
```bash
./build/main/LearnODT datasets/mouse.csv -d 2 --shati
```

Learn an optimal oblique decision tree on the `datasets/mouse.csv` dataset with a depth of 2 and our MaxSAT encoding using the external solver `nuwls` with timeout of 3 seconds:
```bash
./build/main/LearnODT datasets/mouse.csv -d 2 --solver "timeout 3 ./nuwls/nuwls"
```

Run experiment 1, 2, 3 and 4 used to build Fig2:
```bash
./build/main/LearnODT datasets/mouse.csv --exp 1
./build/main/LearnODT datasets/mouse.csv --exp 2
./build/main/LearnODT datasets/mouse.csv --exp 3
./build/main/LearnODT datasets/mouse.csv --exp 4
```
