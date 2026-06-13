# Draft Update To Prof. Chhita

Dear Prof. Chhita,

I have started running simulations for the fixed site-dependent transition probability model. I implemented the edge-weight version you suggested, with

```text
w_N = 1,  w_E = 1,  w_S = u,  w_W = v,
```

where `u` and `v` are sampled from positive mean-one Gamma distributions and then normalised into transition probabilities.

One issue has appeared quite clearly in the numerical diagnostics. Although the raw weights have mean one, the normalised transition probabilities seem to produce a strong effective north-east drift, especially for moderate or strong disorder. For example, with `Gamma(k, 1/k)` weights and `L = 256`, the mean boundary hit location was approximately:

```text
k = 1:   (0.894 L, 0.889 L)
k = 0.5: (0.916 L, 0.918 L)
```

The winding variance then appears quite flat, which seems consistent with the drifted-walk behaviour you mentioned rather than with the intended anomalous winding question.

As a diagnostic comparison, I also tried a balanced model where all four outgoing weights are independent mean-one Gamma variables. That model stays centred in the boundary-hit diagnostics, but up to `L = 256` its winding variance looks closer to the usual logarithmic behaviour than to clear `(log L)^2` growth.

Before I run larger simulations, I wanted to check whether the two-weight model above is definitely the intended model, or whether I should instead use a locally more balanced random environment, for example by randomising all four outgoing edge weights symmetrically.

Best wishes,

Alberto

