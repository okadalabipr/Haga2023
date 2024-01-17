# Bifurcation analysis of TGF-β1, THBS1, and FMOD network



<br>

## Figure S8A and S8C
Fitting our model to the experimental data using Gnuplot (version 5.4),
we inferred the model parameters.

Below is our bash command:

```sh
cd 001_fitting
gnuplot gscript_fitting_Fig3A.txt
mv 001_fitting/fit.log 001_fitting/fit_gscript_fitting_Fig3A.log

gnuplot gscript_fitting_Fig2AB.txt
mv 001_fitting/fit.log 001_fitting/fit_gscript_fitting_Fig2AB.log
```

The results are stored in `001_fitting/fit_gscript_fitting_Fig3A.log`
and `001_fitting/fit_gscript_fitting_Fig2AB.log`.



<br>

## Figures 4B, S8B, and 4D
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

Likewise, we computed bifurcation diagrams of endogenous [THBS1], [FMOD],
and [TGFβ1] with respect to a parameter gamma for TGF-β1 using the
following command:

```sh
cd 010_control_r_PDL
bash aout_001.sh
bash aout_002.sh
gnuplot gscript_r.txt
```

## Figures 4C and S10B
We computed distributions of [THBS1] and [FMOD] with respect to a parameters for THBS1.

Below is our bash command:

```sh
cd 020_stochastic_simulation_forK1
bash aout.sh
```

Similarly, we computed distributions of [THBS1] and [FMOD] with respect to a parameters for FMOD.

Below is our bash command:

```sh
cd 021_stochastic_simulation_forK2
bash aout.sh
```

Likewise, we computed distributions of endogenous [THBS1] and [FMOD] with respect to a parameter gamma for TGF-β1 using the
following command:

```sh
cd 030_stochastic_simulation_forGamma
bash aout.sh
```
