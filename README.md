# mon_projet_codex

This repository now includes a MATLAB helper `compute_uio_gains.m` that
solves the LMIs of Theorem 1 and computes the observer gains K1.

If you prefer a ready-to-run script, execute `matlab/compute_uio_gains_script.m`
in MATLAB or Octave. It loads `ETAPE1.mat`, solves the LMIs subsystem by
subsystem, prints the matrices `P` and `K1`, and saves them to `ETAPE2.mat`.
