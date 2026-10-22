# Interventional Bayesian Causal Discovery (IBCD)

IBCD is an empirical Bayesian causal discovery method using interventional data to infer causal graphs with user-selectable scale-free or Erdős–Rényi structural priors. You can read the [IBCD paper](https://arxiv.org/abs/2510.01562) on arXiv.

### Installation ###

Create a new environment and install the required packages:

```
conda create -n ibcd python=3.10 -y
conda activate ibcd
pip install -r requirements.txt
```


### Quick start ### 

To run IBCD use the following command

```
python ibcd.py --data /PATH/TO/data.csv \
               --prior sf \
               --output_dir OUTPUT_FOLDER
```

Use `--prior sf` (scale-free) or `--prior er` (Erdős–Rényi) depending on the expected graph structure of your data. If the structure is unknown, we recommend using the `sf` prior as a default. For details on all available configuration options, see the **Arguments** section below.


### Input file ###

IBCD takes a single CSV file where each row is a sample, columns `V1, V2, …, Vp` contain expression values, and the final column `target` indicates whether the sample is `control` or belongs to a specific interventional condition (e.g., `V1`, `V2`, `V3`, …). See the input data example CSV [here](https://github.com/bcbg-bio/IBCD/blob/main/data/input/data.csv).

#### Input format example ####

Control samples appear first; interventional samples follow with a label indicating which experimental condition they belong to.

```
V1,   V2,   V3,   V4,   V5,   target
4.92, 5.13, 6.02, 6.44, 6.91, control
4.81, 5.05, 5.87, 6.29, 6.65, control
4.76, 5.01, 5.92, 6.41, 6.78, control
4.88, 5.10, 5.95, 6.38, 6.70, control
...
3.01, 3.36, 4.67, 4.93, 5.44, V5
3.82, 2.64, 3.93, 4.64, 4.90, V5
2.27, 2.94, 4.14, 4.97, 5.48, V5
2.45, 3.24, 4.36, 4.67, 5.13, V5
```

### Output files ###

IBCD produces four output files. See the output files example [here](https://github.com/bcbg-bio/IBCD/tree/main/data/output).

- **G.csv**: Inferred causal graph given as the posterior-mean weighted adjacency matrix.  
- **G_draws.npy**: Posterior samples of the weighted adjacency matrix $G$, saved as a NumPy array across MCMC draws.
- **pip.csv**: Posterior inclusion probability for each edge, which measures how strongly the posterior supports the existence of an edge. 
- **lfsr.csv**: Local false sign rate, the posterior probability that the inferred sign of an edge is incorrect.

`G_draws.npy` is produced by every run but is not included in the example
directory, as it scales with the number of draws (`num_chains` × `num_samples`
× `p` × `p` float32 values) and is too large to keep in the repository.

The example outputs in `data/output` were generated from `data/input/data.csv`
(p = 50) with:

```
python src/ibcd.py --data data/input/data.csv \
                   --prior sf \
                   --output_dir data/output \
                   --num_warmup 300 \
                   --num_samples 200 \
                   --num_chains 3
```

`num_warmup` and `num_samples` are reduced from the defaults of 1000 to keep
the example inexpensive to reproduce on a CPU; the run above takes roughly 50
minutes on a 16-core laptop. Use the defaults for real analyses.


### Arguments ###
```
usage: ibcd.py [-h] --data DATA --prior {sf,er} --output_dir OUTPUT_DIR [--sf_anchor {em,none}]
               [--pi0_floor PI0_FLOOR] [--alpha_er ALPHA_ER] [--num_warmup NUM_WARMUP]
               [--num_samples NUM_SAMPLES] [--num_chains NUM_CHAINS]
               [--target_accept_prob TARGET_ACCEPT_PROB] [--max_tree_depth MAX_TREE_DEPTH]
               [--chain_method {parallel,sequential,vectorized}] [--save_diagnostics] [--seed SEED]
               [--rho_penalty {none,gaussian}] [--rho_sigma RHO_SIGMA] [--init_opt_steps INIT_OPT_STEPS]
               [--slab_width {none,em}] [--init_opt_lr INIT_OPT_LR] [--init_jitter INIT_JITTER]
               [--epsilon EPSILON] [--oracle_G ORACLE_G] [--global_prior] [--likelihood {mn,mvn}]
               [--mvn_rank MVN_RANK] [--mvn_jitter MVN_JITTER]

IBCD pipeline. 1) Load data.csv (observation + intervention) 2) Run 2SLS 3) Choose SF (scale-free) or ER
(Erdős–Rényi) empirical prior 4) Fit empirical Bayesian spike-and-slab prior on matrix normal model 5)
Output G, PIP, and LFSR.

options:
  -h, --help            show this help message and exit
  --data DATA           Path to input data CSV.
  --prior {sf,er}       Choice of empirical prior: 'sf' = scale-free, 'er' = Erdős–Rényi.
  --output_dir OUTPUT_DIR
                        Directory to save all outputs.
  --sf_anchor {em,none}
                        SF prior only: where the overall sparsity level comes from. 'em' takes the global
                        spike weight from the EM of equation 12 and lets the degree profile share it out
                        across rows. 'none' is the published max-normalisation, pi0_i = 1 - theta_i /
                        max(theta, phi), which has no level of its own. Default em.
  --pi0_floor PI0_FLOOR
                        SF prior only: smallest per-node spike proportion, so no row is left entirely
                        unshrunk. The row budget is capped at (1 - pi0_floor) * (D - 1) and the excess
                        redistributed.
  --alpha_er ALPHA_ER   Alpha for EM in ER prior. Controls shrinkage strength. Default=2.0.
  --num_warmup NUM_WARMUP
                        Number of NUTS warm-up iterations. Default = 1000.
  --num_samples NUM_SAMPLES
                        Number of posterior samples per chain after warm-up. Default = 1000.
  --num_chains NUM_CHAINS
                        Number of parallel MCMC chains. Default = 3.
  --target_accept_prob TARGET_ACCEPT_PROB
                        NUTS target acceptance probability. Raising it shrinks the adapted step size, which
                        reduces divergences at the cost of longer trajectories. Default = 0.9.
  --max_tree_depth MAX_TREE_DEPTH
                        Maximum NUTS tree depth; each iteration costs at most 2^depth - 1 leapfrog steps.
                        Default = 12.
  --chain_method {parallel,sequential,vectorized}
                        How to draw multiple chains. 'vectorized' maps them onto one device, which is the
                        only form of within-process parallelism available when CUDA exposes a single device,
                        and avoids the post-hoc stack that 'sequential' pays for. 'parallel' needs one
                        visible device per chain and falls back to sequential otherwise. Default =
                        vectorized.
  --save_diagnostics    Write diagnostics.json to the output directory: divergences, leapfrog steps, split
                        R-hat, ESS, spectral-radius quantiles, rejected-draw counts, runtime and the run
                        configuration. Off by default because R-hat and ESS are O(D^2) over the entries of
                        G.
  --seed SEED           PRNG seed for MCMC. Default = 42.
  --rho_penalty {none,gaussian}
                        Soft acyclicity constraint on G: 'gaussian' is the Appendix H prior N(0,
                        rho_sigma^2) on the spectral radius of G, estimated by 50 steps of power iteration.
                        Default none.
  --rho_sigma RHO_SIGMA
                        Scale of the Gaussian rho penalty. Default 0.5, the Appendix H value.
  --init_opt_steps INIT_OPT_STEPS
                        NUTS starts from init_to_median(num_samples=50), an essentially empty G, after this
                        many Adam steps on the log posterior, which leave the posterior unchanged; 0 starts
                        NUTS at the (jittered) init_to_median draws themselves. Default 2000.
  --slab_width {none,em}
                        'none' for the horseshoe as published; 'em' for the regularised horseshoe with its
                        slab width set to the RMS scale of the EM's fitted slab. The regularised slab's tail
                        beyond the width is Gaussian rather than Cauchy. Default em.
  --init_opt_lr INIT_OPT_LR
                        Adam learning rate for the descent before NUTS. Default 0.01.
  --init_jitter INIT_JITTER
                        Noise added to log(lam) and eps after optimising, so the chains do not start at one
                        point. Default 0.1.
  --epsilon EPSILON     Threshold for computing PIP: edges with |G| > epsilon are counted as active. Default
                        = 0.05.
  --oracle_G ORACLE_G   Ablation (Table 7, Appendix E): CSV of the true G; the ER prior's spike weight is
                        set to its share of zero off-diagonal entries instead of the EM estimate. ER prior
                        only.
  --global_prior        Ablation (Table 7): give every edge the EM's global spike weight instead of the
                        edge-specific weights. ER prior only.
  --likelihood {mn,mvn}
                        'mn' is the matrix-normal likelihood. 'mvn' is the Table 7 ablation: a multivariate
                        normal on vec(R_hat) whose covariance is the full covariance S of equation 9
                        truncated to its top --mvn_rank components. S is D^2 x D^2, so small D only. Default
                        mn.
  --mvn_rank MVN_RANK   Components of S kept for --likelihood mvn. Default 10, as in the paper.
  --mvn_jitter MVN_JITTER
                        Diagonal added to the truncated S for --likelihood mvn, which is otherwise singular.
                        Default 1e-5.
```

