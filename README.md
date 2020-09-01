# FLAPWxHZ
The Hari-Zimmermann complex generalized hyperbolic SVD and EVD.

A part of the supplementary material for the paper arXiv:[1907.08560](https://arxiv.org/abs/1907.08560 "The LAPW method with eigendecomposition based on the Hari–Zimmermann generalized hyperbolic SVD") \[math.NA\].

## Building

### Prerequisites

A recent 64-bit Linux (e.g., CentOS 7.7 with devtoolset-8) or macOS (e.g., Catalina) is needed.

Have the Intel MKL (Math Kernel Library) installed.

Then, clone and build [JACSD](https://github.com/venovako/JACSD) in a directory parallel to this one.

### Make options

Run ``make`` as follows:
```bash
cd src
make [COMPILER=gnu|x64|x200|nvidia] [NDEBUG=0|1|2|3|4|5] [all|clean|help]
```
where ``COMPILER`` should be set for the Intel C/C++ and Fortran compilers (version 19.1+/2020+ recommended) to ``x64`` for Xeons, or to ``x200`` for Xeon Phi KNLs, respectively.
If ``COMPILER`` is not set, GNU C/C++/Fortran compilers will be used instead.

GNU Fortran 9 and 10 are *not* supported!
Please take a look [here](https://gcc.gnu.org/gcc-9/changes.html) for the explanation regarding the MAX and MIN intrinsics.
Currently, only GPU Fortran *8* is fully supported.
On RHEL/CentOS it is provided by, e.g., devtoolset-8.

Here, ``NDEBUG`` should be set to the desired optimization level (``3`` is a sensible choice).
If unset, the predefined debug-mode build options will be used.

For example, ``make COMPILER=x200 NDEBUG=3 clean all`` will trigger a full, release-mode rebuild for the KNLs.

## Execution

### Command line

In the examples below, ``TPC`` stands for threads-per-core.
If the hyperthreading is not desired, it should be set to ``1``.

``FN`` is the input and output file name prefix (without an extension).

#### Phase 0

```bash
/path/to/phase0.exe input.bin FN
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

``JSTRAT1`` can be ``2`` for ``cycwor`` or ``4`` for ``mmstep`` (recommended).

``JSTRAT2`` can be ``3`` for ``cycwor`` (recommended if a particular number of threads is supported) or ``5`` for ``mmstep``.

``NSWP1`` (``1`` for block-oriented) and ``NSWP2`` (``30`` should suffice in most cases) are the maximal numbers of the inner and of the outer sweeps allowed, respectively.

#### Phase 4

```bash
OMP_NUM_THREADS=T OMP_PLACES=CORES OMP_PROC_BIND=SPREAD,CLOSE /path/to/phase4.exe FN N TPC
```

### Data format

All data is stored in the Fortran array order.

The testing dataset is available for [download](http://euridika.math.hr:1846/Jacobi/FLAPW-data/) (please, conserve the bandwidth by downloading only what is of interest to you).

An example of data format of the test cases:

| file name |      data type |      rows | columns |
| --------- | -------------- | --------- | ------- |
| ``FN.X``  | ``COMPLEX(8)`` | ``2*L*a`` | ``G``   |
| ``FN.T``  | ``COMPLEX(8)`` | ``2*L``   | ``2*L`` |
| ``FN.U``  |    ``REAL(8)`` | ``L*a``   | ``1``   |
| ``FN.YY`` | ``COMPLEX(8)`` | ``2*L*a`` | ``G``   |
| ``FN.WW`` | ``COMPLEX(8)`` | ``2*L*a`` | ``G``   |
| ``FN.JJ`` | ``INTEGER(8)`` | ``2*L*a`` | ``1``   |
| ``FN.Y``  | ``COMPLEX(8)`` | ``G``     | ``G``   |
| ``FN.W``  | ``COMPLEX(8)`` | ``G``     | ``G``   |
| ``FN.J``  | ``INTEGER(8)`` | ``G``     | ``1``   |
| ``FN.P``  | ``INTEGER(8)`` | ``G``     | ``1``   |
| ``FN.O``  | ``INTEGER(8)`` | ``G``     | ``1``   |
| ``FN.YU`` | ``COMPLEX(8)`` | ``G``     | ``G``   |
| ``FN.WV`` | ``COMPLEX(8)`` | ``G``     | ``G``   |
| ``FN.Z``  | ``COMPLEX(8)`` | ``G``     | ``G``   |
| ``FN.EY`` |    ``REAL(8)`` | ``G``     | ``1``   |
| ``FN.EW`` |    ``REAL(8)`` | ``G``     | ``1``   |
| ``FN.E``  |    ``REAL(8)`` | ``G``     | ``1``   |
| ``FN.SY`` |    ``REAL(8)`` | ``G``     | ``1``   |
| ``FN.SW`` |    ``REAL(8)`` | ``G``     | ``1``   |
| ``FN.SS`` |    ``REAL(8)`` | ``G``     | ``1``   |
| ``FN.ZZ`` | ``COMPLEX(8)`` | ``G``     | ``G``   |

#### Phase 0

Outputs ``FN.X``, ``FN.T``, ``FN.U``.

#### Phase 1

Input: ``FN.X``, ``FN.T``, ``FN.U``.

Output: ``FN.YY``, ``FN.WW``, ``FN.JJ``.

#### Phase 2

Input: ``FN.YY``, ``FN.WW``, ``FN.JJ``.

Output: ``FN.Y``, ``FN.W``, ``FN.J``, ``FN.P``, ``FN.O``.

#### Phase 3

Input: ``FN.Y``, ``FN.W``, ``FN.J``.

Output: ``FN.YU``, ``FN.WV``, ``FN.Z``; ``FN.EY``, ``FN.EW``, ``FN.E``; ``FN.SY``, ``FN.SW``, ``FN.SS``.

#### Phase 4

Input: ``FN.Z``.

Output: ``FN.ZZ``.

This work has been supported in part by Croatian Science Foundation under the project IP-2014-09-3670 ([MFBDA](https://web.math.pmf.unizg.hr/mfbda/)).
