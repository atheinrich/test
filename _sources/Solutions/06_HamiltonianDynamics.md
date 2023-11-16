---
jupytext:
  formats: ipynb,md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.13.0
kernelspec:
  display_name: Python 3 (phys-521)
  language: python
  name: phys-521
---

```{code-cell}
:tags: [hide-input]

import mmf_setup; mmf_setup.nbinit(quiet=True)
import logging;logging.getLogger("matplotlib").setLevel(logging.CRITICAL)
```


# Assignment 6: Hamiltonian Dynamics

## Arbitrary Dispersion: Negative Mass and Special Relativity

Consider a particle moving in 1D with  kinetic energy $E(p)$ moving under a constant
force $F$.  Use Hamilton's equations to find the general solution for the position
$x(t)$ of the particle.  Check your answer with the familiar solution for $E(p) =
p^2/2m$.  Discuss the physical meaning of $E'(p)$ and $E''(p)$ in terms of Newton's law
and the particle motion.

:::{note}
  This approach also works for $E(p) = \sqrt{p^2c^2 + m^2c^4}$ where $c$ is the speed of
  light. This gives the motion of a particle under constant force in special relativity
  The corresponding coordinate transformation into a co-moving constantly accelerating
  frame gives rise to [Rindler
  coordinates](https://en.wikipedia.org/wiki/Rindler_coordinates), which are applicable
  close to the surface of the earth.  These have interesting properties associated with
  general relativity, including time-dilation at different heights, and an event horizon
  at a distance $d=mc^2/F$ below the observer.  For the earth, this distance is about
  $d\approx0.3$[pc](https://en.wikipedia.org/wiki/Parsec), well past the limit where the
  approximation of a constant gravitational field breaks down. 
:::

### Solution

::::{solution}

Hamilton's equations follow from the Hamltonian:

\begin{gather*}
  H(x, p) = E(p) + V(x), \\
  \dot{p} = -\pdiff{H}{x} = F(x) = -V'(x), \qquad
  \dot{x} = E'(p).
\end{gather*}

From the last equation, we see that $E'(p)$ is the velocity of the particle.  In fluid
dynamics, E&M ect. this corresponds to the **group velocity**, not the phase velocity.

Newton's law is obtained by differentiating this equation once more to obtain the acceleration:

\begin{gather*}
  \ddot{x} = E''(p)\dot{p} = E''(p)F(x) = -E''(p)V'(x).
\end{gather*}

Comparing this to Newton's law $F = ma$, we have $\ddot{x} = m^{-1}F$, hence $E''(p)
\equiv m^{-1}$ should be interpreted as the **inverse effective mass**.  Note that in
higher dimensions, this is formally a tensor:

\begin{gather*}
  [\mat{M}^{-1}]_{ij} = \frac{\partial^2 E(\vect{p})}{\partial p_i\partial p_j}.
\end{gather*}

Now for constant force $V(x) = -Fx$, we can integrate Hamilton's equations.

:::{warning}
  Several people tried to "integrate" these equations doing something like this: since
  $\dot{x} = E'(p)$, 
  
  \begin{gather*}
    x(t) = \int \dot{x} \d{t} = \int E'(p) \d{t} \neq E'(p) t + x_0.
  \end{gather*}
  
  This is only correct if $E'(p)$ is a constant, but is generally incorrect since $p =
  p(t)$ depends on time.  Specifically in this case, $p(t) = p_0 + Ft$.  Please be
  careful to make sure that your integration variables make sense and that all
  dependencies are made explicit before integrating.
:::

The first equation of motion is straightforward: $\dot{p}(t) = F$, therefore $p(t) =
p_0 + Ft$.  The second equation uses a trick that, for a constant force, we can change
variables from $\d{p} = F\d{t}$, hence we can change the integration variable and use
the fact that the integrand $E'(p)$ is a total derivative:

\begin{gather*}
  E'(p) = \diff{E(p)}{p} = \diff{E}{t}\diff{t}{p} = \diff{E}{t}\dot{p}^{-1} = \diff{E}{t}F^{-1}.
\end{gather*}

Let $t\in[0, t]$, with $p(0) = p_0$, $E_0 = E(p_0)$, etc.  Then:

\begin{gather*}
  \dot{p} = F, \qquad
  \int_{p_0}^{p}\d{p} = \int_{0}^{t} F \d{t},\qquad
  p(t) = p_0 + F t,\\
  \dot{x} = E'(p), \qquad
  \int_{x_0}^{x}\d{x} = \int_{0}^{t} E'(p)\d{t}
                      = \frac{1}{F}\int_{0}^{t} \diff{E(p)}{t}\d{t}
                      = \frac{1}{F}\int_{E_0}^{E}\d{E}\\
  x(t) - x_0 = \frac{E\bigl(p(t)\bigr) - E_0}{F}.
\end{gather*}

Hence,

\begin{gather*}
  x(t) = x_0 + \frac{E(p_0 + F t) - E(p_0)}{F}.
\end{gather*}

For regular Newtonian mechanics, $E(p) = p^2/2m$, and we have the familiar parabola:

\begin{gather*}
  x(t) = x_0 + \frac{(p_0 + F t)^2 - p_0^2}{2mF}
       = x_0 + \frac{2p_0 F t + F^2 t^2}{2mF}\\
       = x_0 + \frac{p_0}{m} t + \frac{F}{m}\frac{t^2}{2}
       = x_0 + v_0 t + a\frac{t^2}{2}.
\end{gather*}

For the relativistic dispersion, we have

\begin{gather*}
  x(t) = x_0 + \frac{\sqrt{(p_0 + F t)^2c^2 + m^2c^4} - \sqrt{p_0^2c^2 + m^2c^4}}{F}.
\end{gather*}

Choosing coordinates so that $x_0 = -\sqrt{p_0^2c^2 + m^2c^4}/F$ and starting our clock
at $t =  -p_0/F$, this simplifies:

\begin{gather*}
  x(t) = \sqrt{t^2c^2 + \frac{m^2c^4}{F^2}} = \sqrt{t^2c^2 + d^2}.
\end{gather*}

From the following figure, we see that light starting out a distance $d$ behind the
accelerating observer will never reach the observer.  This corresponds to the Rindler
horizon.  An object dropped by the accelerating observer (green circle) will continue on
a straight path, crossing the horizon (green x).  From the point of view of the
accelerating observer (not obvious from this plot -- you must transform [Rindler
coordinates](https://en.wikipedia.org/wiki/Rindler_coordinates) to see this), the
dropped object will fall slower and slower, never crossing the horizon.  To the
accelerating observer, the time for the dropped object will slow down as it approaches
the horizon. This is the equivalent of gravitational redshift.

::::

```{code-cell} ipython3
:tags: [hide-input]

%matplotlib inline
import numpy as np, matplotlib.pyplot as plt

c = d = 1.0

t_max = 3.0
t_drop = 0.3

t = np.linspace(0, t_max)
x = np.sqrt((c*t)**2 + d**2)
x_drop = np.sqrt((c*t_drop)**2 + d**2)
v_drop = c*t_drop / x_drop
t_drop_horizon = (x_drop - v_drop * t_drop) / (c-v_drop)
t_o = np.linspace(t_drop, t_max)
x_o = x_drop + v_drop * (t_o - t_drop)
fig, ax = plt.subplots()
ax.plot(x/d, c*t, '-C0', label="Accelerating observer")
ax.plot(c*t/d, c*t, ':C1', label="Light")
ax.plot(x_o/d, c*t_o, '--C2', label="Dropped object")
ax.plot([x_drop / d], [c*t_drop], 'oC2')
ax.plot([c*t_drop_horizon/d], [c*t_drop_horizon], 'xC2')
ax.set(xlabel="x/d", ylabel="$ct$", aspect=1)
ax.grid(True)
ax.legend();
```








## Hamilton Jacobi Equation

Your assignment is to analyze the motion of a harmonic oscillator using the Hamiltonian
formalism.  Please follow the outline given below for analyzing a free particle and
complete the same type of analysis for the harmonic oscillator.

### Free Particle (Sample Analysis)

#### Lagrangian Analysis

*Analyze the problem using the Lagrangian formalism.*

1. We use the generalized coordinate $q = x$ and velocity $\dot{q} = v$ so that the Lagrangian is:

   $$
     \mathcal{L}(q, \dot{q}, t) = \frac{m}{2}\dot{q}^2.
   $$

2. The canonical momentum is:

   $$
     p = \pdiff{\mathcal{L}}{\dot{q}} = m\dot{q}.
   $$

3. The Euler-Lagrange equation is:

   $$
     \dot{p} = m\ddot{q} = \pdiff{\mathcal{L}}{q} = 0
   $$

4. The general solution is

   $$
      \dot{q}(t) = v, \qquad
      q(t) = q_0 + v_0t.
   $$

+++

#### Hamiltonian Analysis

+++

1. From above, we can form the Hamiltonian using the Legendre transform.  First we invert $p = m\dot{q}$ to find $\dot{q} = p/m$, then we perform the Legendre transformation:

   $$
     H(q, p, t) = p \dot{q}(p) - \mathcal{L}(q, \dot{q}(p), t) = \frac{p^2}{2m}.
   $$

2. Now we can express the problem in terms of Hamilton's equation of motion:

   \begin{align*}
     \begin{pmatrix}
       \dot{q}\\
       \dot{p}
     \end{pmatrix}
     &=
     \begin{pmatrix}
       0 & 1\\
       -1 & 0
     \end{pmatrix}
     \cdot
     \begin{pmatrix}
       \pdiff{H}{q}\\
       \pdiff{H}{p}
     \end{pmatrix}
     =
     \begin{pmatrix}
       0 & 1\\
       -1 & 0
     \end{pmatrix}
     \cdot
     \begin{pmatrix}
       0\\
       p/m
     \end{pmatrix}\\
     &=
     \frac{1}{m}
     \begin{pmatrix}
       0 & 1\\
       0 & 0
     \end{pmatrix}
     \cdot
     \begin{pmatrix}
       q\\
       p
     \end{pmatrix}.
   \end{align*}

3. The solution can be expressed as:

   \begin{align*}
    \begin{pmatrix}
       q\\
       p
     \end{pmatrix}
     &=
     e^{\left(\begin{smallmatrix} 0 & 1\\ 
                                  0 & 0\end{smallmatrix}\right)\frac{t}{m}}
    \cdot
    \begin{pmatrix}
       q_0\\
       p_0
     \end{pmatrix}                                 
     =
     \left(
     \mat{1} 
     + \begin{pmatrix} 
         0 & 1\\ 
         0 & 0
       \end{pmatrix}\frac{t}{m}
     \right)
    \cdot
    \begin{pmatrix}
       q_0\\
       p_0
     \end{pmatrix}\\
     &=
     \begin{pmatrix} 
       1 & \frac{t}{m}\\ 
       0 & 1
     \end{pmatrix}
    \cdot
    \begin{pmatrix}
       q_0\\
       p_0
     \end{pmatrix}
     =
    \begin{pmatrix}
       q_0 + p_0t/m\\
       p_0
     \end{pmatrix}                                 
   \end{align*}
   
   where we have used the fact that $\mat{A}^n = \mat{0}$ for $n>1$ where
   $\mat{A}=\bigl(\begin{smallmatrix} 0 & 1\\ 0 & 0\end{smallmatrix}\bigr)$.
   
   
4. The Hamilton-Jacobi equation is:

   $$
     H\left(q, \pdiff{S}{q}, t\right) + \pdiff{S}{t}
     = \frac{1}{2m}\left(\pdiff{S}{q}\right)^2 + \pdiff{S}{t} = 0\\
     = \frac{1}{2m}S_{,q}^2 + S_{,t} = 0.
   $$
   
5. This equation is separable, and we may place the $q$'s on one side, and the $t$'s on
   the other to obtain:

   $$
     \frac{1}{2m}S_{,q}^2 = E = -S_{,t}.
   $$
   
   Integrating each side, we obtain:
   
   $$
      S(q, t) = \sqrt{2mE} q + f(t), \qquad
      S(q, t) = W(q) - Et,
   $$
   
   where $f(t)$ and $W(q)$ are the integration constants of each of the pieces.  The solution is thus:
   
   $$
     S(q,t) = \sqrt{2mE}q - Et + S_0.
   $$
   
   Note that $W(q) = \sqrt{wmE}q + S_0$ is called "*Hamilton's characteristic function*" (i.e. in Fetter and Walecka) or sometimes the "*abbreiviated action*" (Landau and Lifshitz) and the form $S(q,t) = W(q) - Et$ is always valid when $H(q, p, t) = H(q, p)$ is independent of time.

6. We are free to choose any new coordinate $Q$ as long as the invertability requirement still holds:

   $$
     \left(\frac{\partial^2 S}{\partial q \partial Q}\right) = \sqrt{\frac{m}{2E(Q)}}E'(Q) \neq 0.
   $$
   
   Since the new Hamiltonian $H'(Q, P, t) = H(q, p, t) + S_{,t} = 0$ by construction, the equations of motion are $\dot{Q} = \dot{P} = 0$ and $P$ and $Q=E$ are constants of motion.
   
   A convenient choice is $E(Q) = Q$: i.e. introducing the energy $E$ as the new coordinate.
   
7. After this choice is made, we have the *generating function* for the canonical transformation:

   $$
     S(q, Q, t) = \sqrt{2mQ}q - Qt + S_0.
   $$   
     
   The canonical momentum follows from the following, which may be inverted to express $q(Q, P, t)$:
   
   $$
     P = - \pdiff{S(q,Q, t)}{Q} = t - \sqrt{\frac{m}{2E}} q\\
     q(Q, P, t) = \sqrt{\frac{2E}{m}}(t - P) = v_0(t - t_0) = q_0 + v_0t, \\
     p(Q, P, t) = S_{,q} = \sqrt{2mE} = mv_0.
   $$
   
   Thus, we see that the canonical momentum $P$ to the energy $Q=E$ is the initial time $t_0$.
   
8. Inverting these, we have the explicit canonical transformation:

   $$
     Q(q, p, t) = E = \frac{p^2}{2m}, \qquad
     P(q, p, t) = t - \sqrt{\frac{m}{2E}}q = t - \frac{mq}{p}.
   $$
   
   We can now explicitly check that the Poisson bracket (I am using the convention of Landau and Lifshitz here) satisfies the canonical commutation reationships:
   
   $$
     [P, Q] = \pdiff{P}{p}\pdiff{Q}{q} - \pdiff{P}{q}\pdiff{Q}{p}
            = \frac{mq}{p^2}\cdot 0 - \left(-\frac{m}{p}\right)\cdot \frac{p}{m}
            = 1.
   $$

9. *(This analysis is a little harder to do for the oscillator, so do not feel you have
   to do it.)* Armed with the solutions, we may construct the action function for the
   path connecting $(q_0, t_0)$ to $(q_1, t_1)$:

   $$
     q(t) = q_0 + \frac{q_1 - q_0}{t_1 - t_0}(t-t_0), \qquad
     \dot{q}(t) = \frac{q_1 - q_0}{t_1 - t_0}.     
   $$
   
   Hence, we have the action:
   
   $$
     \bar{S}(q_0, t_0; q, t) = \int_{t_0}^{t} \mathcal{L}\d{t}
     = \int_{t_0}^{t} \frac{m}{2}\frac{(q_1 - q_0)^2}{(t_1 - t_0)^2}\d{t}
     = \frac{m}{2}\frac{(q_1 - q_0)^2}{t - t_0}.
   $$
   
   This allows us to construct the general solution to any initial-value problem for the
   Hamilton-Jacobi equation:
   
   $$
     H(q, S_{,q}(q, t), t) + S_{,t}(q,t) = 0, \qquad S(q, t_0) = S_0(q)
   $$
   
   as
   
   $$
     S(q, t) = S_0(q_0) + \int_{t_0}^{t} L(q(t), \dot{q}(t), t)\;\d{t}
   $$
   
   where the action is computed over the trajectory starting from $q(t_0) = q_0$ with
   initial momentum $p_0 = S_0'(q_0)$ and ending at $q(t) = q$ at time $t$.  For this
   problem $p_0 = mv_0$ so we have 
   
   $$
     q(t) = q_0 + v_0(t - t_0) = q_0 + \frac{S'_0(q_0)}{m}(t - t_0)
   $$
   
   which must be inverted to find $q_0 = q_0(q, t, t_0)$.  The explicit solution here is
   expressed in terms of this:
   
   $$
     S(q, t) = S_0(q_0) + \frac{1}{2m}[S'_0(q_0)]^2 (t - t_0).
   $$
   
   Since $H$ is independent of time, we can take $t_0 = 0$ without loss of generality.
   Now, consider an example problem from Arnold where he asks for the solution to this
   problem with initial conditions $S_0(q) = \frac{mq^2}{2T}$ (though he chooses units
   where $m=T=1$).  This can be explicity constructed:
   
   $$
     S'_0(q_0) = \frac{m q_0}{T} = mv_0, \qquad
     q = q_0 + \frac{q_0}{T}t, \qquad
     q_0 = \frac{q}{1 + \frac{t}{T}} = \frac{q T}{T + t}.
   $$
   
   The explicit solution is thus
   
   $$
     S(q, t) = \frac{m q^2 T}{2(T + t)^2} + \frac{1}{2m}\left(\frac{mq}{T + t}\right)^2 t
     =  \frac{mq^2}{2(T + t)}.
   $$
   
   Note that this does *not* have the same form as the separable solution we constructed
   above.  This is due to the choice of initial condistions $S_0(q)$.  In particular,
   our separable solution corresponds with the initial conditions $S_0(q) = m v_0 q$
   instead.  Mathematically, however, the Hamilton-Jacobi equations can be solved with
   any arbitrary initial conditions.

### Harmonic Oscillator (Solution)

Here is a similar analysis for the harmonic oscillator.


#### Lagrangian Analysis

*Analyze the problem using the Lagrangian formalism.*

1. Use the generalized coordinate $q = x$ and velocity $\dot{q} = v$ so that the
   Lagrangian is:

   $$
     \mathcal{L}(q, \dot{q}, t) = \frac{m}{2}\dot{q}^2 - \frac{k}{2}q^2.
   $$

2. Find the the canonical momentum?

   $$
     p = m\dot{q}.
   $$

3. Write the Euler-Lagrange equation:

   $$
     \dot{p} = -kq
   $$

4. Write down the general solution:

   $$
      q(t) = q_0\cos(\omega t) + \frac{\dot{q}_0}{\omega}\sin(\omega t).
   $$


#### Hamiltonian Analysis


*Analyze the problem using the Hamiltonian formalism*.

1. Use the Legendre transform to write the Hamiltonian:

   $$
     H(q, p, t) = p\dot{q} - L = \frac{p^2}{2m} + \frac{k q^2}{2}.
   $$

2. Now we can express the problem in terms of Hamilton's equation of motion:

   $$
     \begin{pmatrix}
       \dot{q}\\
       \dot{p}
     \end{pmatrix}
     = 
     \begin{pmatrix}
       H_{,p} = p/m\\
       -H_{,q} = -kq
     \end{pmatrix}
     = \overbrace{
       \begin{pmatrix}
         0 & m^{-1}\\
         -k & 0
       \end{pmatrix}
     }^{\mat{A}}
     \cdot
    \begin{pmatrix}
       q\\
       p
     \end{pmatrix}.
   $$

3. To compute the matrix exponential, we note that:
   \begin{gather*}
     \mat{A}^2 = 
     -\frac{k}{m}\mat{1}, \qquad
     \mat{A}^3 = -\frac{k}{m}\mat{A},
   \end{gather*}

   etc.  Thus, letting $\omega^2 = k/m$, we have
   
   \begin{gather*}
     e^{\mat{A}t} = 
     \overbrace{\left(
       1 - \frac{\omega^2 t^2}{2} + \frac{\omega^4 t^4}{4!} \cdots\right)
     }^{\cos(\omega t)}\mat{1}\\
     +
     \overbrace{\left(
       t - \frac{\omega^2 t^3}{3!} - \frac{\omega^4 t^5}{5!} \cdots\right)
     }^{\omega^{-1}\sin(\omega t)}\mat{A}\\
     = \cos(\omega t)\mat{1} + \omega^{-1} \sin(\omega t)\mat{A}.
   \end{gather*}

   The solution can be expressed as:

   $$
    \begin{pmatrix}
       q\\
       p
     \end{pmatrix}
     =
     e^{\mat{A}t}
    \begin{pmatrix}
       q_0\\
       p_0
     \end{pmatrix}
    =
    \begin{pmatrix}
      \cos(\omega t) &  \sin(\omega t)/\sqrt{km}\\
      -\sqrt{km}\sin (\omega t) & \cos(\omega t)
    \end{pmatrix}
    \begin{pmatrix}
       q_0\\
       p_0
     \end{pmatrix}.
   $$
   
   ::::{warning}

   Note that we recover the **oscillations** implicit in the name "Harmonic Oscillator"
   and found in from the Lagrangian analysis.  There is **no justification** for
   truncating the exponential.  All terms must be kept, and the series summed.
   
   If you
   are not comfortable with the "trick" presented here, there are many other ways of
   doing this:
   
   :::{admonition} Details
   :class: dropdown
   
   **Diagonalize the matrix:**
   
   This is the most work, but straightforward.  The eigenvalues satisfy $\det(\mat{A} -
   \lambda \mat{1}) = \lambda^2 + k/m = 0$, hence $\lambda_{\pm} = \pm \I \omega$ where
   $\omega^2 = k/m$.  Solving for the eigenvectors, we have
   \begin{gather*}
     \mat{A} =
     \overbrace{
     \begin{pmatrix}
       m^{-1} & \lambda_-\\
       \lambda_+ & -k
     \end{pmatrix}}^{\mat{S}}
     \overbrace{
     \begin{pmatrix}
       \lambda_+ & 0\\
       0 & \lambda_-
     \end{pmatrix}}^{\mat{D}}
     \overbrace{
     \begin{pmatrix}
       m^{-1} & \lambda_-\\
       \lambda_+ & -k
     \end{pmatrix}^{-1}}^{\mat{S}^{-1}}\\
     =
     \begin{pmatrix}
       m^{-1} & -\I\omega\\
       \I\omega & -k
     \end{pmatrix}
     \begin{pmatrix}
       \I\omega & 0\\
       0 & -\I\omega
     \end{pmatrix}
     \frac{
       \begin{pmatrix}
       -k & \I\omega\\
       -\I\omega & m^{-1}
     \end{pmatrix}}{\underbrace{\frac{-k}{m}-\omega^2}_{-2\omega^2}}.
   \end{gather*}
   Note: the matrix is not Hermitian, and the eigenvectors are not orthonormal, so we
   must use the matrix inverse on the right.
   
   The matrix exponential is thus
   \begin{gather*}
     e^{\mat{A}t} = \mat{S}e^{\mat{D}t}\mat{S}^{-1}
     e^{\mat{A}t} = \mat{S}
     \begin{pmatrix}
       e^{\I\omega t} & 0\\
       0 & e^{-\I\omega t}
     \end{pmatrix}
     \mat{S}^{-1}\\
     =
     \frac{1}{\frac{-k}{m}-\omega^2}
     \begin{pmatrix}
       m^{-1} & -\I\omega\\
       \I\omega & -k
     \end{pmatrix}
     \begin{pmatrix}
       e^{\I\omega t} & 0\\
       0 & e^{-\I\omega t}
     \end{pmatrix}
     \begin{pmatrix}
       -k & \I\omega\\
       -\I\omega & m^{-1}
     \end{pmatrix}\\
     = 
     \begin{pmatrix}
       \cos(\omega t) & 
       \frac{1}{m\omega}\sin(\omega t)\\
       \frac{-k}{\omega}\sin(\omega t) &
       \cos(\omega t)
     \end{pmatrix},
   \end{gather*}
   as claimed above.
   
   **Convert the system to a second-order equation and solve conventionally:**
   
   \begin{gather*}
      \dot{q} = \frac{p}{m}, \qquad \dot{p} = -kq\\
       \ddot{q} = \frac{\dot{p}}{m} = -\frac{k}{m}q\\
       q(t) = a\cos(\omega t) + b\sin(\omega t), \qquad \omega^2 = \frac{k}{m}.
   \end{gather*}
   Now solve for $a$ and $b$ in terms of the initial conditions.
   
   
   :::


      
   
   
   
   ::::

   
4. The Hamilton-Jacobi equation is:

   $$
     H\left(q, S_{,q}, t\right) + S_{,t} = \frac{S_{,q}^2}{2m} + \frac{kq^2}{2} + S_{,t} = 0.
   $$
   
5. This equation is separable, and therefore each side is a constant:

   \begin{gather*}
     \frac{S_{,q}^2}{2m} + \frac{kq^2}{2} = -S_{,t} = E.
   \end{gather*}

   Integrating, we get the "*abbreviated action*" $W(q)$ and the generating function $S(q,t)$:

   \begin{align*}
     S(q, t) &= -Et + f(q), \qquad
     \diff{S}{q} = \sqrt{2mE - mkq^2}, \\
     S(q, t) &= \int\d{q} \sqrt{2mE - mkq^2} + g(t)\\
     &= \underbrace{\frac{q}{2}\sqrt{2mE - kmq^2} 
        + \frac{mE}{\sqrt{km}}\tan^{-1}\frac{\sqrt{mk}q}{\sqrt{2mE - kmq^2}}}_{W(q)} - Et.
   \end{align*}
   
6. You are free to choose any new coordinate $Q$ as long as the invertability
   requirement still holds.  Here we let $Q = E$:

   $$
     \left(\frac{\partial^2 S}{\partial q \partial Q}\right) = 
     \frac{m}{\sqrt{2mQ - mkq^2}} \neq 0.
   $$
   
7. Now express the *generating function* $S(q,Q,t)$ for the canonical transformation in
   terms of your chosen coordinate:

   $$
     S(q, Q, t) = \frac{q}{2}\sqrt{2mQ - kmq^2} 
     + \frac{mQ}{\sqrt{km}}\tan^{-1}\frac{\sqrt{mk}q}{\sqrt{2mQ - kmq^2}} - Qt.
   $$
     
   Use this to determine the canonical momenta:
   
   \begin{align*}
     P = -S_{,Q} &= t - \overbrace{\sqrt{\frac{m}{k}}}^{1/\omega}\tan^{-1}\frac{\sqrt{mk}q}{\sqrt{2mQ - kmq^2}},\\
     p(Q, P, t) = S_{,q} &= \sqrt{2mQ - mkq^2},\\
     q(Q, P, t) &= \sqrt{\frac{2Q/k}{1 + \cot^2\bigl(\omega(t - P)\bigr)}}.
   \end{align*}

8. Invert these to get the explicit canonical transformation:

   \begin{align*}
     Q(q, p, t) &= \frac{p^2}{2m} + \frac{m\omega^2 q^2}{2},\\
     P(q, p, t) &= t - \frac{1}{\omega}\tan^{-1}\frac{m\omega q}{p}.
   \end{align*}
   
   Explicitly check that the Poisson bracket (I am using the convention of Landau and
   Lifshitz here) satisfies the canonical commutation relationships: 
   
   $$
     [P, Q] = \pdiff{P}{p}\pdiff{Q}{q} - \pdiff{P}{q}\pdiff{Q}{p}
            = \frac{m^2q^2\omega^2}{m^2q^2\omega^2 + p^2} - 
              \frac{-p^2}{m^2q^2\omega^2 + p^2}
            = 1.
   $$

9. Armed with the generating function designed so that $H'(Q, P, t) = 0$, argue that $Q$
   and $P$ are constant, and hence use your result from step 7. to write down the
   complete solution and compare with your previous results. 

   \begin{align*}
     q(t) &= \sqrt{\frac{2Q/k}{1 + \cot^2\bigl(\omega(t - P)\bigr)}}\\
          &=\sqrt{\frac{2Q}{k}}\sin\bigl(\omega(t - P)\bigr),\\
     p(t) &= \sqrt{2mQ - mkq^2} = \sqrt{\frac{2mQ}{1+\tan^2\bigl(\omega(t - P)\bigr)}}\\
          &= \sqrt{2mQ}\cos\bigl(\omega(t - P)\bigr).
   \end{align*}
   
   These are the same after identifying $P$ with an appropriate phase, and $Q$ with the
   energy.  We may now construct the action function for the path connecting $(q_0, t_0)$ to
   $(q_1, t_1)$:
   
   \begin{gather*}
     q(t) = A\cos\bigl(\omega (t-t_0)\bigr) + B\cos\bigl(\omega (t-t_1)\bigr),\\
     \begin{pmatrix}
       q_0\\
       q_1
     \end{pmatrix}
     =
     \begin{pmatrix}
       1 & \cos\bigl(\omega(t_0-t_1)\bigr)\\
       \cos\bigl(\omega(t_1-t_0)\bigr) & 1
     \end{pmatrix}
     \cdot
     \begin{pmatrix}
       A\\
       B
     \end{pmatrix}.
   \end{gather*}
   
   Let $c = \cos\bigl(\omega(t_1-t_0)\bigr)$ and $s = \sin\bigl(\omega(t_1-t_0)\bigr)$,
   then $A = (q_0-cq_1)/s^2$, $B=(q_1-cq_0)/s^2$ so
   
   \begin{align*}
     q(t) &= \frac{
       (q_0 - cq_1)\cos\bigl(\omega (t-t_0)\bigr)
       +           
       (q_1 - cq_0)\cos\bigl(\omega (t-t_1)\bigr)}{s^2},\\
     p(t) &= m\dot{q}.
   \end{align*}
   
   After some tedious algebra, we can compute the action:
   
   \begin{align*}
     S(t_0, q_0; t_1, q_1) 
     &= \int_{t_0}^{t_1}\left(\frac{p^2}{2m} - \frac{m\omega q^2}{2}\right)\d{t}\\
     &= \frac{m \omega}{2}
        \frac{(q_0^2 + q_1^2)\cos\bigl(\omega (t_1 - t_0)\bigr) - 2q_0q_1}
             {\sin\bigl(\omega (t_1 - t_0)\bigr)}.
    \end{align*}
    
    As expected from the time-independence of $H$, this is a function of $t_1 - t_0$.
    
10. Since this problem is both separable and periodic, we can express the solution in
    terms of *action-angle variables*.  The action variable is simply:
    
    \begin{gather*}
      J = \oint p \d{q} = \int_{0}^{2\pi} \frac{p^2}{m\omega} \d{\theta}
        = \frac{2 m E}{m\omega} \int_{0}^{2\pi} \cos^2\theta \d{\theta}
        = \frac{2 \pi E}{\omega}, \\
      E(J) = \frac{\omega J}{2\pi},
    \end{gather*}
    
    where we let $\theta = \omega (t - P)$ so that $m\omega \d{q} = p \d{\theta}$.  The
    abbreviated action $W(J)$ and corresponding angle variable $w$ are
    
    \begin{gather*}
      W(J) = Et = \frac{\omega J t}{2\pi}, \qquad
      w = \diff{W(J)}{J} = \frac{\omega}{2\pi}t.
    \end{gather*}
    
    The linear frequency is
    
    \begin{gather*}
      \nu = \diff{E(J)}{J} = \frac{\omega}{2\pi}
    \end{gather*}
    
    as expected.  In this case, $Q = w$ and $P = J$ are canonically conjugate:
    
    \begin{gather*}
      J(p, q) = \frac{2\pi E(p, q)}{\omega} 
              = \frac{2\pi}{\omega}\left(\frac{p^2}{2m} + \frac{m\omega^2 q^2}{2}\right),\\
      w(p, q) = \frac{\omega t}{2\pi} = 
                \frac{1}{2\pi} \cos^{-1}\frac{p}{\sqrt{2mE}},
    \end{gather*}
    
    thus, the Poisson braket is:
    
    \begin{gather*}
      [P, Q] = [J, w] = \pdiff{J}{p}\pdiff{w}{q} - \pdiff{w}{p}\pdiff{J}{q} = 1.
    \end{gather*}


