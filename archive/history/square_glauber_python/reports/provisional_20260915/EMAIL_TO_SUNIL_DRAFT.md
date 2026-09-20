Subject: Glauber dynamics experiment: method, results and sampling difficulty

Dear Sunil,

I wanted to summarise my understanding of Glauber dynamics, explain how I am using it, and share the results so far.

My understanding is that Glauber dynamics samples configurations by making small random local updates. For the dimer model, I choose a square uniformly at random. If it contains two parallel dimers, I choose between the horizontal and vertical arrangements with probabilities proportional to the products of their edge weights. Otherwise, nothing changes. This update preserves the desired weighted equilibrium distribution. The practical difficulty is knowing how long to run before the starting configuration no longer matters; an apparently stable trace alone does not establish this.

To address that, I am using coupling from the past. I start from the lowest and highest possible height configurations and evolve both using exactly the same random updates, from a chosen time in the past to time zero. Because these updates preserve height order, if the two configurations agree at time zero, every possible starting configuration gives the same answer. This provides an equilibrium sample without choosing an arbitrary burn-in. If they do not agree, I double how far back I start, keeping the previously assigned random history unchanged, and try again.

For the experiment, I generate random edge weights from a gamma distribution with shape 1 and mean 1, then keep those weights fixed throughout sampling. I also run a control with every edge weight equal to 1. The grid sizes are L=8, 12, 16, 20, 24, 28 and 32.

For each environment, I run the whole sampling procedure twice with independent random streams. Each run has its own coupled lower and upper configurations, so the two final samples are independent conditional on the shared weights. I measure the height at the centre and height differences between separated points. Comparing the two samples across environments lets me measure the variation caused by the random weights, separately from fluctuations within one fixed environment.

The analysis is now complete for 23,799 of the 23,800 planned pairs, using 10,000 bootstrap resamples of whole environments. The main findings are:

- The gamma model has positive centre disorder covariance, increasing overall from about 0.71 at L=8 to 1.41 at L=32. The uniform controls fluctuate around zero, although the L=24 centre control has a nominal 95% interval just below zero.
- A simple logarithmic fit describes the centre data reasonably well over these sizes. Adding a (log L)^2 term gives a coefficient of 0.045, with a 95% bootstrap interval from -0.344 to 0.425. The spatial analysis also does not establish a positive quadratic term.
- The fitted curvature changes considerably when I change the size range. I therefore do not think these results distinguish the asymptotic growth laws reliably yet.

The remaining difficulty is one gamma environment at L=32. Both sampling runs still fail to couple after a history of roughly one billion sweeps. At this size, one sweep means 961 attempted local updates for each bounding chain. Doubling the history again would therefore mean roughly another 30 hours, without a guarantee of success. I have paused further sampling.

The current results are explicitly provisional: the unresolved environment is omitted, and the bootstrap intervals do not account for possible bias from that omission. It is only one environment, but its unusually slow sampling means I cannot assume it is irrelevant. The long coupling time is also not, by itself, evidence for a particular covariance growth law.

Does this approach seem appropriate to you, and would you recommend trying to finish the last environment, or reconsidering the sampling strategy and size range before spending more computation?

Best,
Alberto
