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

  2. run `run.m` in each XPA, KS, and REITER folders, which solves model for different values of the standard deviation of aggregate productivity by calling `main_XPA_v1.m`, `main_KS_v1.m`, or `main_REITER_v1.m`. It may take time for a while especially for KS. Also note that `UpwindKZ = 0` does not work for XPA.

  3. Then run `Figure1and2.m` and `Figure3.m` to replicate the figures and tables in the paper.

* You can also directly run each `main_XPA_v1.m`, `main_KS_v1.m`, or `main_REITER_v1.m` by setting `diagnose = 1` to see how each algorithm works.

* For the robustness checks in Appendix, you may want to change the parameter $\mu$ and the above meta parameters in `main_XPA_v1.m` and `main_KS_v1.m` and rerun `run.m`

For questions about the programs, please send an email to: Masakazu Emoto <masakazu.emoto@gmail.com>

The source code is licensed MIT, see LICENSE.
