# Bifurcation analysis of TGF-β1, THBS1, and FMOD network



<br>

## TGF-β1, THBS1, and FMOD network model
Let $x_1$, $x_2$, and $x_3$ be [THBS1], [FMOD], and [TGF-β1], respectively.
Based on the experimental results, we constructed the following mathematical model.

$$
\begin{align}
\frac{\displaystyle {\rm d} x_1}{\displaystyle {\rm d} t}
&= a \left( \frac{\displaystyle x_3^{n_a}}{\displaystyle x_3^{n_a} + K_a} \right)
- d_1 x_1, \tag{1} \\[3pt]
\frac{\displaystyle {\rm d} x_2}{\displaystyle {\rm d} t}
&= b \left( \frac{\displaystyle K_b}{\displaystyle x_3^{n_b} + K_b} \right)
- d_2 x_2, \tag{2} \\[3pt]
\frac{\displaystyle {\rm d} x_3}{\displaystyle {\rm d} t}
&= c \left( \frac{\displaystyle x_1}{\displaystyle x_1 + K_1} \right)
\left( \frac{\displaystyle K_2}{\displaystyle x_2 + K_2} \right) - d_3 x_3, \tag{3}
\end{align}
$$

where $a/d_1$, $b/d_2$, and $c/d_3$ are maximal values of [THBS1], [FMOD], and [TGF-β1], respectively;
$K_a$, $K_b$, $K_1$, and $K_2$ are the half saturation constants; and
$n_a$ and $n_b$ are Hill coefficients.



<br>

## Steady state problem
Setting the left-hand sides of Eqs. (1)–(3) to zero, the steady state solutions of $x_1$, $x_2$, and $x_3$ can be obtained as follows:

$$
\begin{align}
\bar{x_1} &= A
\left( \frac{\displaystyle \bar{x_3}^{n_a}}{\displaystyle \bar{x_3}^{n_a} + K_a} \right),
\tag{4} \\[3pt]
\bar{x_2} &= B
\left( \frac{\displaystyle K_b}{\displaystyle \bar{x_3}^{n_b} + K_b} \right),
\tag{5} \\[3pt]
\bar{x_3} &= C
\left( \frac{\displaystyle \bar{x_1}}{\displaystyle \bar{x_1} + K_1} \right)
\left( \frac{\displaystyle K_2}{\displaystyle \bar{x_2} + K_2} \right), \tag{6}
\end{align}
$$

where $A = a / d_1$, $B = b / d_2$, and $C = c / d_3$.

By fitting (4) and (5) to the experimental data, the parameters can be inferred using a nonlinear least-squares method with Gnuplot (version 5.4).
Below is our bash command:

```sh
cd 001_fitting
gnuplot gscript_fitting.txt
```

We obtained $A = 2.36$, $B = 0.33$, $K_a = 0.016$, $K_b = 0.002$, $n_a = 1.6$, and $n_b = 1.7$ (see `001_fitting/fit.log`).



<br>

Substituting (4) and (5) into (6), and denoting

$$
f(x) = C \left[ \frac{\displaystyle x^{n_a}}{\displaystyle x^{n_a} + K_A (x^{n_a} + K_a)} \right]
\left[ \frac{\displaystyle K_B (x^{n_b} + K_b)}{\displaystyle K_b + K_B (x^{n_b} + K_b)} \right] 
- x, \tag{7}
$$

where $K_A = K_1 / A$ and $K_B = K_2 / B$, the following equation holds:

$$f(\bar{x}_3) = 0. \tag{8}$$

Hereinafter, let $x = \bar{x}_3$.

Assuming $n_a \ge 1$ (this is validated later), one can see that $x = 0$ is a trivial solution of (8).
Divided by $x$ ($\ne 0$), (8) is reduced to the following equation:

$$
\begin{align}
C K_B x^{n_a - 1} (x^{n_b} + K_b) - \left[ x^{n_a} + K_A (x^{n_a} + K_a) \right]
\left[ K_b + K_B (x^{n_b} + K_b) \right] = 0. \tag{9}
\end{align}
$$

Here, let us denote

$$
\begin{align}
g(x) &= C K_B x^{n_a - 1} (x^{n_b} + K_b) \\[3pt]
&- \left[ x^{n_a} + K_A (x^{n_a} + K_a) \right]
\left[ K_b + K_B (x^{n_b} + K_b) \right], \tag{10} \\[3pt]
\frac{\displaystyle {\rm d}g}{\displaystyle {\rm d}x}
&= (n_a - 1) C K_B x^{n_a - 2} (x^{n_b} + K_b) + n_b C K_B x^{2(n_a - 1)} \\
&- \left[ n_a (1 + K_A) x^{n_a - 1} \right] \left[ K_b + K_B (x^{n_b} + K_b) \right] \\
&- \left[ x^{n_a} + K_A (x^{n_a} + K_a) \right] \cdot n_b K_B x^{n_b - 1}. \tag{11}
\end{align}
$$

Now (9) can be numerically solved using Newton's method as follows:

$$
x^{(n + 1)} = x^{(n)} - \frac{g(x^{(n)})}{\frac{{\rm d}g}{{\rm d}x}(x^{(n)})},
$$

where $x^{(0)}$ is an initial guess for a root of $g(x)$ and $\{ x^{(n)} \}$ is expected to converge to the root.



<br>

## Stability analysis
Stability of a steady state solution of $f(x) = 0$ can be determined by the sign of $\frac{{\rm d} f}{{\rm d} x}$, which is calculated as follows:

$$
\begin{align}
\frac{{\rm d} f}{{\rm d} x}
&= C \left[ \frac{\displaystyle n_a K_A K_a x^{n_a - 1}}{\displaystyle \left\{ x^{n_a} + K_A (x^{n_a} + K_a) \right\}^2} \right]
\left[ \frac{\displaystyle K_B (x^{n_b} + K_b)}{\displaystyle K_b + K_B (x^{n_b} + K_b)} \right]
\\[3pt]
&+ C \left[ \frac{\displaystyle x^{n_a}}{\displaystyle x^{n_a} + K_A (x^{n_a} + K_a)} \right]
\left[ \frac{\displaystyle n_b K_B K_b x^{n_b - 1}}{\displaystyle \left\{ K_b + K_B (x^{n_b} + K_b) \right\}^2} \right] - 1.
\tag{12}
\end{align}
$$

If $\frac{{\rm d} f}{{\rm d} x}(\bar{x}) > 0$ (resp. $\frac{{\rm d} f}{{\rm d} x}(\bar{x}) < 0$), the steady state solution $\bar{x}$ is unstable (resp. stable).

For instance, $\left. \frac{{\rm d} f}{{\rm d} x} \right|_{x = 0} = -1 < 0$.
Hence, the trivial solution $\bar{x} = 0$ is always stable.



<br>

## Bifurcation analysis
By arbitrarily setting $C = 1$, one can compute a bifurcation diagram of $\bar{x}_3$ with respect to $K_1$ and $K_2$.

Below is our bash command:

```sh
cd 002_control_pq
bash aout_001.sh
bash aout_002.sh
gnuplot gscript_pq.txt
```

Also, bifurcation diagrams of $\bar{x}_1$, $\bar{x}_2$, and $\bar{x}_3$ with respect to $K_1$ can be computed by the following command:

```sh
cd 003_control_p_THBS1
bash aout_001.sh
bash aout_002.sh
gnuplot gscript_p.txt
```

Similarly, bifurcation diagrams with respect to $K_2$ can be computed by the following command:

```sh
cd 004_control_q_FMOD
bash aout_001.sh
bash aout_002.sh
gnuplot gscript_q.txt
```
