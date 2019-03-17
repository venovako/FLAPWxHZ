# FLAPWxHZ
The Hari-Zimmermann complex generalized SVD and EVD.

## Building

### Prerequisites

A recent 64-bit Linux (e.g., CentOS 7.6) or macOS (e.g., Mojave) is needed.

Have the Intel MKL (Math Kernel Library) installed.

Then, clone and build [JACSD](https://github.com/venovako/JACSD) in a directory parallel to this one.

### Make options

Run ``make`` as follows:
```bash
cd src
make [CPU=x64|x200] [NDEBUG=0|1|2|3|4|5] [all|clean|help]
```
where ``CPU`` should be set for the Intel C++ and Fortran compilers to ``x64`` for Xeons, or to ``x200`` for Xeon Phi KNLs, respectively.
If ``CPU`` is not set, GNU C and Fortran compilers will be used instead.

Here, ``NDEBUG`` should be set to the desired optimization level (``3`` is a sensible choice).
If unset, the predefined debug-mode build options will be used.

For example, ``make CPU=x200 NDEBUG=3 clean all`` will trigger a full, release-mode rebuild for the KNLs.

## Execution

### Command line

In the examples below, ``TPC`` stands for threads-per-core.
If the hyperthreading is not desired, it should be set to ``1``.

``FN`` is the input and output file name prefix (without an extension).

#### Phase 0

```bash
/path/to/phase0.exe input.bin output.bin
```

Phase 0 is a data conversion phase from a custom data format to a set of plain binary files.

#### Phase 1

```bash
OMP_NUM_THREADS=T OMP_PLACES=CORES OMP_PROC_BIND=SPREAD,CLOSE /path/to/phase1.exe FN L a G TPC
```

``L``, ``a``, and ``G`` are the problem-specific parameters.

#### Phase 2

```bash
OMP_NUM_THREADS=T OMP_PLACES=CORES OMP_PROC_BIND=SPREAD,CLOSE /path/to/phase2.exe FN M N TPC
```

#### Phase 3

```bash
OMP_NUM_THREADS=T OMP_PLACES=CORES OMP_PROC_BIND=SPREAD,CLOSE /path/to/phase3.exe FN M N TPC JSTRAT1 NSWP1 JSTRAT2 NSWP2
```

``JSTRAT1`` is the inner, and ``JSTRAT2`` the outer Jacobi strategy.
``JSTRAT1`` can be ``2`` for ``cycwor`` or ``4`` for ``mmstep``.
``JSTRAT2`` can be ``3`` for ``cycwor`` or ``5`` for ``mmstep``.

``NSWP1`` and ``NSWP2`` are the maximal numbers of the inner and of the outer sweeps allowed, respectively.

#### Phase 4

```bash
OMP_NUM_THREADS=T OMP_PLACES=CORES OMP_PROC_BIND=SPREAD,CLOSE /path/to/phase4.exe FN N TPC
```

### Data format

To be explained in more detail...

#### Phase 0

Outputs ``FN.X``, ``FN.T``, ``FN.U``.

#### Phase 1

Input: ``FN.X``, ``FN.T``, ``FN.U``.

Output: ``FN.YY``, ``FN.WW``, ``FN.JJ``.

#### Phase 2

Input: ``FN.YY``, ``FN.WW``, ``FN.JJ``.

Output: ``FN.Y``, ``FN.W``, ``FN.J``.

#### Phase 3

Input: ``FN.Y``, ``FN.W``, ``FN.J``.

Output: ``FN.YU``, ``FN.WV``, ``FN.Z``; ``FN.EY``, ``FN.EW``, ``FN.E``; ``FN.SY``, ``FN.SW``, ``FN.SS``.

#### Phase 4

Input: ``FN.Z``.

Output: ``FN.ZZ``.

This work has been supported in part by Croatian Science Foundation under the project IP-2014-09-3670 ([MFBDA](https://web.math.pmf.unizg.hr/mfbda/)).
