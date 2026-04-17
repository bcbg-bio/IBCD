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
               --causal_order false \
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

IBCD produces three output files. See the output files example [here](https://github.com/bcbg-bio/IBCD/tree/main/data/output).

- **G.csv**: Inferred causal graph given as the posterior-mean weighted adjacency matrix.  
- **G_draws.npy**: Posterior samples of the weighted adjacency matrix $G$, saved as a NumPy array across MCMC draws.
- **pip.csv**: Posterior inclusion probability for each edge, which measures how strongly the posterior supports the existence of an edge. 
- **lfsr.csv**: Local false sign rate, the posterior probability that the inferred sign of an edge is incorrect.


### Arguments ###
```
usage: ibcd.py [-h] --data DATA --prior {sf,er} --causal_order {true,false} --output_dir OUTPUT_DIR
               [--alpha_sf ALPHA_SF] [--alpha_er ALPHA_ER] [--num_warmup NUM_WARMUP]
               [--num_samples NUM_SAMPLES] [--num_chains NUM_CHAINS] [--epsilon EPSILON]

IBCD pipeline. 1) Load data.csv (Y_matrix + target) 2) Run 2SLS 3) Choose SF (scale-free) or ER
(Erdős–Rényi) empirical prior 4) Fit empricial Bayesian spike-and-slab on matrix normal model 5) Output G
draws, posterior mean G, PIP, and LFSR.

options:
  -h, --help            show this help message and exit
  --data DATA           Path to input data CSV.
  --prior {sf,er}       Choice of empirical prior: 'sf' = scale-free, 'er' = Erdős–Rényi.
  --causal_order {true,false}    'false' = input variables are unordered, as in most real data (e.g. perturb-seq); 'true' = input variables are already in causal/topological order.
  --output_dir OUTPUT_DIR
                        Directory to save all outputs.
  --alpha_sf ALPHA_SF   Penalty parameter for SF prior (scale-free row-wise optimization). Default=1.0.
  --alpha_er ALPHA_ER   Alpha for EM in ER prior. Controls shrinkage strength. Default=2.0.
  --num_warmup NUM_WARMUP
                        Number of NUTS warm-up iterations. Default = 300.
  --num_samples NUM_SAMPLES
                        Number of posterior samples per chain after warm-up. Default = 1000.
  --num_chains NUM_CHAINS
                        Number of parallel MCMC chains. Default = 3.
  --epsilon EPSILON     Threshold for computing PIP: edges with |G| > epsilon are counted as active.
                        Default = 0.05.
```

