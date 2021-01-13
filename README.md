# XPA Algorithm in Continuous Time

This repository includes programs for solving the Krusell-Smith model in continuous time by the KS, XPA, and REITER algorithms as described in

* Masakazu Emoto and Takeki Sunakawa (2021) "Applying the Explicit Aggregation Algorithm to Heterogeneous Agent Models in Continuous Time," which is avaialble at https://tkksnk.github.io/wps/

* To replicate the results in the paper,

  1. make sure that
  ```
  diagnose = 0; % display diagnoses and graphs
  loadtemp = 1; % load temp.mat and parameter sigma; turn on when executing run.m
  UpwindKZ = 1; % upward scheme for K and Z
  KFEnoKZ  = 1; % transition matrix without terms wrt K and Z
  ```
  in `main_XPA_v1.m` and `main_KS_v1.m`.

  2. run `run.m` in each XPA, KS, and REITER folders, which solves model for different values of the standard deviation of aggregate productivity by calling `main_XPA_v1.m`, `main_KS_v1.m`, or `main_REITER_v1.m`.

  3. Then run `Figure1and2.m` `Figure3.m` `Table3.m` to replicate figures and tables in the paper.

* For robustness checks in Appendix, you may want to change the above meta parameters in `main_XPA_v1.m` and `main_KS_v1.m` and rerun `run.m`

For questions about the programs, please send an email to: Masakazu Emoto <masakazu.emoto@gmail.com>

The source code is licensed MIT. The website content is licensed CC BY 4.0, see LICENSE.
