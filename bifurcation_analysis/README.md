# Bifurcation analysis of TGFβ1, THBS1, and FMOD network



<br>

## Figure S5A
Fitting our model to the experimental data using Gnuplot (version 5.4),
we inferred the model parameters.

Below is our bash command:

```sh
cd 001_fitting
gnuplot gscript_fitting.txt
```

The results are stored in `001_fitting/fit.log`.



<br>

## Figures 4B and 4C
We computed a bifurcation diagram of [TGFβ1] with respect to parameters for
THBS1 and FMOD.

Below is our bash command:

```sh
cd 002_control_pq
bash aout_001.sh
bash aout_002.sh
gnuplot gscript_pq.txt
```

Also, we computed bifurcation diagrams of [THBS1], [FMOD], and [TGFβ1]
with respect to a parameter for THBS1 by the following command:

```sh
cd 003_control_p_THBS1
bash aout_001.sh
bash aout_002.sh
gnuplot gscript_p.txt
```

Similarly, we computed bifurcation diagrams of [THBS1], [FMOD], and [TGFβ1]
with respect to a parameter for FMOD by the following command:

```sh
cd 004_control_q_FMOD
bash aout_001.sh
bash aout_002.sh
gnuplot gscript_q.txt
```
