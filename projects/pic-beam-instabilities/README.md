This file can be compiled to pdf with redered equations: `pandoc README.md -o README.pdf`

# Derivations

Here are some derivations on how some quantitiesa are calculated from the numerical values.

These assume that units $m_0$ and $q_0$ are chosen such that for electrons
$|\hat{q}| = \hat{m}$.

Following relation is used: $\Delta x = \frac{4\pi q_0^2 \hat{c}^2}{m_0 c^2}$

## Normalized plasma energy density

\begin{align*}
\frac{B^2}{8\pi n_0 \gamma_b m_e c^2}
= \frac{\hat{B}^2 B_0^2 \Delta x^3}{8\pi \hat{n}_0 \gamma_b m_e c^2}
= \frac{\hat{B}^2} {2 \hat{n}_0 \gamma_b} \frac{B_0^2\Delta x^3} {4\pi m_ec^2}
\end{align*}

\begin{align*}
\frac{B_0^2\Delta x^3} {4\pi mc^2}
&= \frac{m_0^2 c^2} {q_0^2 \hat{c}^2 \Delta t^2} \frac{\Delta x^3} {4\pi m_0 \hat{m}_e c^2} \\
&= \frac{m_0 c^2} {4\pi q_0^2 \hat{c}^2 } \frac{\Delta x^3} {\Delta t^2 \hat{m}_e c^2} \\
&= \frac{\Delta x^2} {\Delta t^2 \hat{m}_e c^2} \\
&= \frac{1} {\hat{m}_e \hat{c}^2}
\end{align*}

\begin{align*}
\frac{B^2}{8\pi n_0 \gamma_b m_e c^2}
= \frac{\hat{B}^2}{2 \hat{n}_0 \gamma_b \hat{m}_e \hat{c}^2}
\end{align*}


## Numerical classic plasma frequency

\begin{align*}
\omega_p^2 &= \frac{4\pi n q^2}{m} \\
&= \frac{4\pi \hat{n} q_0^2 \hat{q}^2}{\Delta x^3 m_0 \hat{m}} \\
&= \frac{4\pi q_0^2 \hat{c}^2}{m_0 c^2}
\frac{\hat{n} \hat{q} c^2}{\Delta x^3 \hat{c}^2} \\
&= \frac{\hat{n} \hat{q} c^2}{\Delta x^2 \hat{c}^2}
\end{align*}

Using $\omega^2_p = \frac{\hat{\omega}_p^2}{\Delta t^2}$, we get:

\begin{align*}
\hat{\omega}_p^2 = \frac{\hat{n} \hat{q} c^2 \Delta t^2}{\Delta x^2 \hat{c}^2} = \hat{n} \hat{q}
\end{align*}

