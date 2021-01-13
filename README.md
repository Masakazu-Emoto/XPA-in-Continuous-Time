# XPA in Continuous Time

These are programs for solving Krusell-Smith model in continuous time by Krusell-Smith, XPA, and REITER as described in

* Masakazu Emoto and Takeki Sunakawa (2021) "Applying the Explicit Aggregation Algorithm to Heterogeneous Agent Models in Continuous Time," which is avaialble at https://tkksnk.github.io/wps/

* To replicate the results in the paper,
  1. make sure that
  ```
  diagnose = 0;
  loadtemp = 1;
  UpwindKZ = 1;
  KFEnoKZ  = 1;
  ```
  in `main_XPA_v1.m` and `main_KS_v1.m`.
  2. run `run.m` in each XPA, KS, and REITER folders, which solves model for different values of the standard deviation of aggregate productivity by calling `main_XPA_v1.m`, `main_KS_v1.m`, and `main_REITER_v1.m`.
  3. Then run `Figure1and2.m` `Figure3.m` `Table3.m` to replicate figures and tables in the paper.

* For robustness checks in Appendix, you may want to change the above meta parameters in `main_XPA_v1.m` and `main_KS_v1.m` and rerun `run.m`

For questions about the programs, please send an email to: Masakazu Emoto <masakazu.emoto@gmail.com>
