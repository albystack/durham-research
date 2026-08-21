# Sunil Chhita — Research Correspondence Archive

> Purpose: source-of-truth archive of the supervisor messages supplied for this research project.
> The archive is intentionally **not** a default Codex context file. Use `context/RESEARCH_BRIEF.md`,
> `context/RESULTS.md`, and `context/NEXT_STEPS.md` first; consult this file only when exact supervisor wording
> or historical context matters.
>
> Scope: Sunil's messages only. Formatting/blank lines have been lightly normalised; wording is
> otherwise preserved from the supplied email exports.

## Thursday, 2 April 2026 at 16:47 — Re: Uniform Spanning Tree Project Outcomes and Results

Dear Alberto,

Thank you for the messages.

The implementation looks good - well done!

Let me think about next steps/extensions and get back to you after the Easter weekend.

Best wishes

Sunil

---

## Monday, 13 April 2026 at 11:28 — Re: Uniform Spanning Tree Project Outcomes and Results

Dear Alberto,

Thanks for the message. Apologies in the delay getting back to you.

The outputs that you sent before look great.  I have had a chance to think about possible extensions and have listed them below.  I should stress that you should prioritise your degree work and revision over this project, so there is no pressure from my side to make any progress.

Extension 1: Try to put implement spanning trees on the triangular lattice ([https://en.wikipedia.org/wiki/Hexagonal\_lattice](https://eur01.safelinks.protection.outlook.com/?url=https%3A%2F%2Fen.wikipedia.org%2Fwiki%2FHexagonal_lattice\&data=05%7C02%7Calberto.rescigno%40durham.ac.uk%7C28b88802cd634440c05a08decc4295ec%7C7250d88b4b684529be44d59a2d8a6f94%7C0%7C0%7C639172783515506623%7CUnknown%7CTWFpbGZsb3d8eyJFbXB0eU1hcGkiOnRydWUsIlYiOiIwLjAuMDAwMCIsIlAiOiJXaW4zMiIsIkFOIjoiTWFpbCIsIldUIjoyfQ%3D%3D%7C0%7C%7C%7C\&sdata=xsXMZu0vI6%2BInZ63Qzv5lPvg5mZz07k%2Bnub%2FcEKPvPo%3D\&reserved=0 "https://eur01.safelinks.protection.outlook.com/?url=https%3A%2F%2Fen.wikipedia.org%2Fwiki%2FHexagonal_lattice\&data=05%7C02%7Calberto.rescigno%40durham.ac.uk%7C28b88802cd634440c05a08decc4295ec%7C7250d88b4b684529be44d59a2d8a6f94%7C0%7C0%7C639172783515506623%7CUnknown%7CTWFpbGZsb3d8eyJFbXB0eU1hcGkiOnRydWUsIlYiOiIwLjAuMDAwMCIsIlAiOiJXaW4zMiIsIkFOIjoiTWFpbCIsIldUIjoyfQ%3D%3D%7C0%7C%7C%7C\&sdata=xsXMZu0vI6%2BInZ63Qzv5lPvg5mZz07k%2Bnub%2FcEKPvPo%3D\&reserved=0")).  The dual graph of the triangular lattice is the honeycomb graph.  The random walk step has 6 possible directions, so a simple symmetric random walk on the triangular lattice chooses each direction with probability 1/6.

For your original code and this extension, Temperley’s bijection applies which means that the configurations that you have drawn correspond to a random tiling.  The picture you sent previously then corresponds to a random domino tiling, whereas Extension 1 corresponds to a random lozenge tiling.

Extension 2: Change the probabilities on the random walk behind the spanning trees so that they are no longer 1/4 (or 1/6 for the triangular lattice).  These will then correspond to simulations of drifted random walks.

Extension 3:  Change the probabilities of the random walk so that at each site, you have random probabilities.  There are a couple of cases here

1. The random walk steps are fixed in time (i.e. around each site, they always have the same probabilities)
2. The random walk steps change at each time step.

I’m guessing that case 1 will be the most interesting so let’s stick with this setup.  You can select any distribution for the probabilities, but I’m thinking having a normal distribution would be sufficient.  For instance, let v and w be two samples from the normal distribution with mean 0 and standard deviation s.  Then the probability of going left, right, up or down is equal to (1/2+u)(1/2+v),(1/2-u)(1/2+v),(1/2+u)(1/2-v),(1/2-u)(1/2-v) respectively.  We would want to see what happens when s changes.

I’m primarily interested in the winding of a single branch of the spanning tree, with the single branch started at the centre of the spanning tree.  Here, the windings of a branch are the number of left turns minus by the number of right turns that the branch makes until it hits the root.  For instance, a branch could wrap around itself twice counter clockwise, so the winding would be 8 (or 9 depending on how it hits the root).

I’m interested in the variance of (many samples) of the windings for large trees with different values of the standard deviation s. The mean should hopefully be zero. I’m primarily interested in large sample sizes.  For the case you simulated before (i.e. each random walk step has probability 1/4), the variance will be C log L, where L is the size of the box and C is some constant.

It is an open problem in the case outlined above - we believe it should be C (log L)^2  (where C is a constant) but there isn’t much computational evidence and there is certainly no mathematical proof.  In order to see anything of interest, it would be good to plot a range of large sizes.

Let me know how you get on.  Hopefully, extension 2 and perhaps part of extension 3 are implementable quite easily from the code that you already have.

Best wishes

Sunil

---

## Monday, 8 June 2026 at 09:40 — Re: Uniform Spanning Tree Project Outcomes and Results

Dear Alberto,

Thanks for the message and good to hear about that are interested in continuing the project.  Hope that the exams all went smoothly.

The plan that you mentioned sounds good.  You are welcome to choose any distribution for the random walk jumps that you would like, e.g. choosing exponential or gamma will get rid of having negative probabilities.

For the probabilities themselves, I think what I wrote below wasn’t quite correct.  It’s probably easiest to have edge weights and then get the probabilities there. So if we label the edges from each site N, E, S, W (for the usual compass directions), we could select the N and E edges to have weight 1, while the S and W edges have weight u and v, where u and v are selected from some distribution with mean 1.  Then, the probabilities would be 1/(2+u+v), 1/(2+u+v), u/(2+u+v) and v/(2+u+v) for transitioning N, E, S and W.

With regards to implementation, what you wrote sounds good. It might be simplest to just focus on drawing the loop erased path from the origin as opposed to drawing the whole tree (to reduce simulation time).  Finally, if the random variables u and v do not have mean zero (so the random walk has a drift), then you should get that the variance of the winding of the loop-erased path starting from the origin should be finite.

Best wishes

Sunil

---

## Monday, 15 June 2026 at 10:30 — Re: Uniform Spanning Tree Project Outcomes and Results

Dear Alberto,

Thanks for the message. Good to see about the progress and the updates.

From what you mentioned, I think you should do as you propose - choose a more locally balanced random environment.  Just to confirm, each site should have a different sample.

It would also be good to run larger simulations and also experiment with different choices of distributions if the gamma distribution does not indicate (log L)^2 behaviour for larger simulations.

There is an interest in the community about the conjectured (log L)^2 behaviour, so if we can nail down the distribution, I think this should give some interest.

Incidentally, I will be on campus on Wednesday and Thursday, so we can set up an appointment to discuss if you are available.

Best wishes

Sunil

---

## Wednesday, 17 June 2026 08:32:28 — Re: Uniform Spanning Tree Project Outcomes and Results

Dear Alberto,

This sounds good. 

How about 2pm on Thursday in my office? 

Best wishes

Sunil

---

## Friday, 10 July 2026 at 09:43 — Re: New Julia results

Dear Alberto,

Hope you are having a good summer.

Many thanks for the update.  It’s good to see that you are up and running with Julia.  The graphs look great, especially page 3! It does seem that Julia is able to handle larger domains too, which is good.

I agree, it does look like we are not seeing any (log L)^2 behaviour.

Could I check one aspect?  Could I check that the annealed version means that you are running a different background for each loop erased random walk?  (i.e. so you do not reuse the same environment twice).

Best wishes

Sunil

---

## Monday, 13 July 2026 at 11:31 — Re: New Julia results

Dear Alberto,

Thanks for the clarification and explanation of the code.

Yes, if you could rerun with one walk per environment, that would be super helpful. I’m probably being overly cautious but we want to remove any type of additional averaging over the environment. It will probably lead to the same results, but it is good to check nonetheless.

Best wishes

Sunil

---

## Thursday, 16 July 2026 at 17:06 — Re: New Julia results

Hi Alberto,

Thanks for the message.

I agree, it’s looking like C(Log L)^p with p=1.   The plots and analysis look good.

I exchanged some messages with a colleague who is working on the super-rough fluctuations (i.e. C(log L)^2) in a slightly different model and they mentioned that they have a proof.  They use the height function from the double-dimer model, which takes two independent copies of the dimer model and subtracts the difference between the height functions.  The model that they consider is a domino tiling of an Aztec diamond with random weights.

I think we should try to see if we can see anything in the double-dimer version of our model.  What this means in terms of spanning trees and our environments is the following:

1. For each sample of the random environment (i.e. the random weights), run two independent copies of the loop erased walks from the origin until they hit the boundary.
2. Computing the number of left turns minus the number of right turns for each independent copy of the loop erased walks, i.e. the windings that you were computing previously
3. Take the difference between the two windings.

It turns out by doing this, you will get the height function of the double dimer model.  Again, I think we probably want to look at the annealed picture.

Have a go with this and see whether you get anything different.  I’m expecting the same behaviour that we observed before, since it is just under copies, but it is probably worth checking to see if there is anything.

Best wishes

Sunil

---

## Thursday, 23 July 2026 at 16:05 — Re: New Julia results

Hi Alberto,

All fine here, I hope you are doing well too.

Thanks for the messages and the plots.

Sorry about the delay getting back to you.  I have been away this week in the US at a conference and have just returned.

From your messages, I think it looks like we are converging on (C log L)^p with p=1 for this model.  That’s interesting that there is some correlation between the two walks on the same environment; I wasn’t expecting this.

Aside from checking the code (which will need to get done at some stage) or optimizing it to be able to run larger simulations, if you agree, I think the next step is to compare the behaviour in the Aztec diamond numerically.  As mentioned last time, a couple of colleagues (Duits and van Peski) have mentioned in passing that they have a proof of the (C log L)^2 behaviour. If we could see that there is a difference between the numerics from the Aztec diamond and the windings of loop erased random walks, then that would be interesting.   Here, we would be wanting to focus on the height of the face at one point of the Aztec diamond, such as the face in the center.

The downside with Aztec diamond, is that the maximum size you may be able to sample is probably around 600 to 800 for random weights.  I do have another colleague (Leonid Petrov) that has implemented the sampling algorithm and optimized part of it using Claude; see [https://lpetrov.cc/domino/](https://eur01.safelinks.protection.outlook.com/?url=https%3A%2F%2Flpetrov.cc%2Fdomino%2F\&data=05%7C02%7Calberto.rescigno%40durham.ac.uk%7C68830ab986c44a748e9a08defa12e7a4%7C7250d88b4b684529be44d59a2d8a6f94%7C0%7C0%7C639223156330320814%7CUnknown%7CTWFpbGZsb3d8eyJFbXB0eU1hcGkiOnRydWUsIlYiOiIwLjAuMDAwMCIsIlAiOiJXaW4zMiIsIkFOIjoiTWFpbCIsIldUIjoyfQ%3D%3D%7C0%7C%7C%7C\&sdata=sF%2BXVXc5%2B9ya4THRCfxTgLH0EL9lyL2%2FhD3Yk2CugGg%3D\&reserved=0 "https://eur01.safelinks.protection.outlook.com/?url=https%3A%2F%2Flpetrov.cc%2Fdomino%2F\&data=05%7C02%7Calberto.rescigno%40durham.ac.uk%7C68830ab986c44a748e9a08defa12e7a4%7C7250d88b4b684529be44d59a2d8a6f94%7C0%7C0%7C639223156330320814%7CUnknown%7CTWFpbGZsb3d8eyJFbXB0eU1hcGkiOnRydWUsIlYiOiIwLjAuMDAwMCIsIlAiOiJXaW4zMiIsIkFOIjoiTWFpbCIsIldUIjoyfQ%3D%3D%7C0%7C%7C%7C\&sdata=sF%2BXVXc5%2B9ya4THRCfxTgLH0EL9lyL2%2FhD3Yk2CugGg%3D\&reserved=0") for an example of tilings of the Aztec diamond.   Let me think if there is anything else we could do to reduce the complexity as we only want a specific statistic.

Best wishes

Sunil

---

## Friday, 24 July 2026 at 18:38 — Re: New Julia results

Hi Alberto,

Thanks for the message.

Yes, the plan is next to explore tilings of the Aztec diamond with random weights and compare with your other results.

That’s good to point about the double-dimer representation. I believe that this should be the eventual goal for the Aztec diamond model.  For Duits and van Peski’s result, I believe that they are using Gamma distributed edge weights, so it probably makes sense to use these.

I’m not sure how accessible Leo Petrov’s code is for the Aztec diamond. I have some old Julia code that I could share next week if that is helpful, but it would probably need to be optimized.  Since we are only requiring the height function of one face, there may be some simplifications that could be done - let me get back to you early next week along with some pointers on the domino shuffling algorithm.

Best wishes

Sunil

---

## Sunday, 26 July 2026 at 11:17 — Re: New Julia results

Hi Alberto,

Apologies about the weekend email (I have a few admin duties on Monday morning, so I doubt I will get to my email until the afternoon).

I think this link does look like Leo’s code. If you go back one file directory (to domino tilings), then I think you can find different packages there. I think a lot of it was coded with agents (from Claude).

I was thinking that perhaps we could have one last try with the loop erased walk case before starting on the Aztec diamond. This time, instead of generating the random environment first, sample the steps directly using the relevant distributions.  That is, at each point that the random walk visits, you will need to generate the probabilities that the walk takes steps N,E, W and S each time.

This will mean that the transition probabilities will not be fixed with time. This interpretation fits better with random walk in random environment.  I’m hoping that this should be simpler to implement than the having the fixed environment and may give something different.

I will get back to you on the Aztec diamond case later on Monday, if not, then on Tuesday.

Best wishes

Sunil

---

## Monday, 27 July 2026 at 15:08 — Re: New Julia results

Hi Alberto,

Thanks for the message and also running these.

Yes, these should hopefully be much simpler to implement, as you do not need to keep track of the whole environment, just sampling each step as needed.

Would it be possible to run one of the distributions, say the Gamma distributed random variable with relatively large variance, to size 5000?  You might be able to implement larger simulations than 5000 because you shouldn’t need as much memory to run these (i.e. there isn’t a 4\*5000^2 block of numbers to store).

Another observable that you could look at is the length of the LERW.  This will be given by the length of the list that you obtain that tracks the loop erased path from (0,0) to the boundary box.  If the observable works like the base case (i.e. where each step has probability 1/4), then one should get the length to be L^{5/4}.  That is, you will get for the base case

log (Length of LERW)/ log L

to be equal to 5/4.

Let me know if you see anything different happening here or whether you get the same features. I haven’t found anything simpler with the Aztec diamond just yet, but will confirm tomorrow.

Best wishes

Sunil

---

## Tuesday, 28 July 2026 at 17:09 — Re: New Julia results

Hi Alberto,

Thanks for the message and the update.

This looks good - so it looks like even the loop erased random walk in random environment is behaving similar to the usual random walk case as well.

I think the next step would be to focus on the Aztec diamond.  Do you have a copy of Mathematica - I know that Durham has a site license and you might be able to view a copy from AppsAnywhere.

I have attached an old copy of some Julia code that randomly simulates the Aztec diamond with chosen weights.  The output is to give a matrix and then in the past, I have used mathematica to draw the resulting picture. The input for the simulations is a table of entries that gives the weights (there are some random weights that I have given into this matrix). Do you want to have a go at trying to simulate an Aztec diamond with random weights?  Once this is done, we can try to get the height function and hopefully see numerically the (log L)^2 behavior.

Best wishes

Sunil

---

## Tuesday, 28 July 2026 at 18:20 — Re: New Julia results

Hi Alberto,

Thanks for the message. It looks like you have a feel on what to do next.

You are welcome to run it with any distribution you want. I originally coded it with uniformly random weights. It is probably good to get a picture with this and then try choosing gamma distributed weights.

Eventually, we probably want to use Gamma distributed weights as there are some simplifications that we can use later with these weights that Duits and van Peski found.

Best wishes

Sunil

---

## Wednesday, 29 July 2026 at 14:05 — Re: New Julia results

Hi Alberto,

Thanks for the email.  I think the next step is then to get height function and then compute the height in the center of the Aztec diamond for different realization.  The heights are defined on the faces of the graph (where the graph is the one with dimers as opposed to dominoes).

The height function in this case is equal to +/-3 if you cross a dimer between a shared edge and -/+1 if do not cross a dimer between the shared edge. Here the sign depends on the parity, that is, if you go between two adjacent faces, the white vertex can be on either the left or the right and the sign is different for each case.  The construction is such if you go around each white vertex clockwise, the height function increases by 1 each time, unless you cross a dimer where the height function drops by 3. For black vertices, the convention is the opposite.

Here is the Julia code for getting the height function:

\#

function heightfunction(x0)

      n=size(x0)[1]

      m=Int(floor(n/2))

      A=zeros(n+1,m+1)

      for j in 1\:m+1

            A[1,j]=2\*j-2

      end

      for i in 2\:n+1

            for j in 1\:m

                  if Int(x0[i-1,2\*j-1])==1

                        A[i,j]=A[i-1,j]+3

                  else

                        A[i,j]=A[i-1,j]-1

                  end

            end

      end

      return A

end

\#

The even rows are offset from the odd rows so that we can store the data in a table.

It’s probably worth trying to get statistics of the height function in the center of the Aztec diamond. That is, compute the height at the center of the Aztec diamond with random weights in a similar fashion to what you did before.  I’m guessing that you will be able to get to size 600.  If we want to go larger, we will need to do some optimizing of the code.  Hopefully, we see the (Log L)^2 behavior here.  If not, we may have to look at the double dimer version of this model, but let’s see how this part goes first.

Best wishes

Sunil

---

## Wednesday, 29 July 2026 at 15:06 — Re: New Julia results

Correction: When I meant in a similar fashion to before, I meant computed many samples.  For computing the height at the center, you just need an element of the array.

Best wishes

Sunil

---

## Thursday, 30 July 2026 at 17:40 — Re: New Julia results!

Hi Alberto,

Thanks for the update and sending these through. It’s looking promising.  Just a clarification when you mention:

\>>For the random weights, I used the Duits-Van Peski Gamme model with independent a\_ij distributed as Gamma(0.2,1) and b\_ij distributed as Gamma(0.25,1) with the other two edge weights fixed to 1. I applied their weight reduction recurrence and then sampled the tiling by deletion, then sliding and creation. To make the computation smoother for larger sizes, I drew the independent creation choices while reducing the weights and storing them as bits. Choices corresponding to sites that are not holes during shuffling are not used, so this has the same distribution as drawing a choice only when a new pair of dominoes is required.

I’m not exactly sure what you mean by this.  Are you saying that you choose the weights that Duits and Van Peski use? This should hopefully give a simplification so that you do not need to simulate using the whole table of weights.  If you are using this, then hopefully you should be able to get larger samples (as you are not needing to keep track of all the weights at each step).

For the double dimer model, you are going to have to sample two configurations with the same weights and then take the height difference.  But hopefully, it looks like we will see the (Log L)^2 for the height function variance without needing to do this.

Best wishes

Sunil

---

## Thursday, 30 July 2026 at 18:12 — Re: New Julia results!

Hi Alberto,

No problem and thanks for the clarification, I understand what you are saying now.

I will have a look at their paper closer.  I’m hoping that there is a simplification of the weights so that we do not have to keep track of the reduction recurrence.  If so, this will need a lot less computational power and allow larger simulations much faster.

Best wishes

Sunil

---

## Monday, 3 August 2026 at 18:36 — Re: New Julia results!

Hi Alberto,

Thanks for the email.  This is really great progress! I haven’t had a chance to digest your email fully yet (as have been busy with a new starter today), but it certainly looks like it is heading in a good direction.

From what I can tell, you have managed to find the (log L)^2 behavior. This is what I had in mind, but it is more subtle than I realized.  For the results that you have, are these for the Aztec diamond? Do you think you can get similar results for the dimer model on the square grid (i.e. the spanning tree model that you were studying earlier).  It might be good to try a few different weights as well with both models.

Yes, I think it would be good to see the code. If we do want to write up the results into a numerics paper, the code will have to be available along with some commentary of the results/analysis.

Best wishes

Sunil

---

## Tuesday, 4 August 2026 at 10:16 — Re: New Julia results!

Hi Alberto,

Thanks for the message.

Yes, this is what I had in mind.

\>>For the square grid case, would you like me to reproduce the same paired analysis for the Aztec diamond - namely sampling two independent configs in the same random environment and separating the conditional component from the shared environment covariance using the spatial height increments?

Yes, this would be great since we want to see if we get the (Log L)^2 behaviour in this case too.

Best wishes

Sunil

---

## Friday, 7 August 2026 at 16:48 — Re: New Julia results!

Dear Alberto,

Thanks for the update.

Looking forward to seeing the results when they come back in.

If the results come back in and they are as expected, we can then start to think about the write-up.  This stage might take some time (journals are quite strict here so we won’t be able to use any AI systems to help with this part).

Best wishes

Sunil

---

## Undated latest supplied reply — after the 11 August 2026 update

> The export did not include a visible Date/Subject header for this topmost message, so no date has been inferred.

Hi Alberto,

Thanks for running this and sorry for the delay getting back to you.   

That’s a pity that we are not seeing the (Log L)^2 behaviour here in the square grid. 

I would like to check that we are actually doing the right simulation here, that is whether the correspondence with the spanning trees is does work in the additional randomness setting. There are a couple of checks that we could do here:

1. We could run with random walk with random environment with random weights as we did before, calling a different random weight at each step (as opposed to fixing the environment). This is the same simulation that you did before, but you will be getting different results due to using the different observable that you found.  
2. Use a different simulation approach using dimer spin flips (Glauber dynamics).  Here, take a dimer configuration on the grid and apply flips. Here a flip is around a square face, and only works if there are dimer covering the opposite edges around the square, and the flip is rotating the edges. That is, the edges of the square are labelled A,B,C and D in clockwise order, and dimers cover A and C, then a flip will mean that dimers cover B and D.  The flips occur with a local probability that depends on the weights.  

Let’s see if 1 comes back with anything. If not, then we can try 2, but let me know before starting as I should give you more information on the spin flips and provide some pointers. 

Best wishes

Sunil
