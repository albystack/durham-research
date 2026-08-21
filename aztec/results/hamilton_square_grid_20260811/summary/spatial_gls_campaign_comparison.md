# Spatial block-GLS campaign comparison

The table below selects the disorder covariance from each analysis. 
Intervals and model comparisons come from joint environment bootstraps.

| Campaign | min-order tag | c | 95% interval | P(c>0) | delta BIC | held-out RMSE log / log2 |
|---|---:|---:|---:|---:|---:|---:|
| directed_gamma_k05_min128 | 128 | -0.002513 | [-0.007742, 0.002546] | 0.1657 | -2.977 | 0.096 / 0.092 |
| directed_gamma_k05_min512 | 512 | -0.0248 | [-0.04424, -0.006096] | 0.0050 | 2.851 | 0.095 / 0.086 |
| directed_gamma_k05_min1024 | 1024 | -0.004031 | [-0.0553, 0.04446] | 0.4162 | -3.441 | 0.079 / 0.088 |
| directed_gamma_k1_min128 | 128 | 0.001144 | [-0.003937, 0.005915] | 0.6590 | -3.670 | 0.096 / 0.096 |
| directed_gamma_k1_min512 | 512 | 0.01189 | [-0.00676, 0.03014] | 0.8850 | -2.091 | 0.103 / 0.098 |
| directed_gamma_k1_min1024 | 1024 | 0.02968 | [-0.02061, 0.07804] | 0.8738 | -2.118 | 0.104 / 0.102 |
| directed_lognormal_var2_min128 | 128 | 0.001218 | [-0.003915, 0.006297] | 0.6655 | -3.649 | 0.148 / 0.147 |
| directed_lognormal_var2_min512 | 512 | -0.007155 | [-0.02596, 0.01052] | 0.2054 | -3.109 | 0.155 / 0.155 |
| directed_lognormal_var2_min1024 | 1024 | -0.03973 | [-0.09073, 0.01054] | 0.0581 | -1.074 | 0.197 / 0.208 |
| directed_uniform_0_2_min128 | 128 | 0.00337 | [-0.001555, 0.008088] | 0.9090 | -1.992 | 0.165 / 0.161 |
| directed_uniform_0_2_min512 | 512 | -0.004514 | [-0.02134, 0.01212] | 0.2896 | -3.413 | 0.152 / 0.167 |
| directed_uniform_0_2_min1024 | 1024 | -0.001267 | [-0.04652, 0.04222] | 0.4529 | -3.463 | 0.163 / 0.200 |
| undirected_gamma_k05_min128 | 128 | 0.003652 | [-0.001268, 0.008608] | 0.9253 | -1.778 | 0.117 / 0.115 |
| undirected_gamma_k05_min512 | 512 | 0.006979 | [-0.01159, 0.02509] | 0.7703 | -3.114 | 0.109 / 0.110 |
| undirected_gamma_k05_min1024 | 1024 | 0.02425 | [-0.02417, 0.07264] | 0.8337 | -2.502 | 0.108 / 0.122 |

These are finite-size numerical comparisons, not asymptotic proofs.
