---
jupytext:
  formats: ipynb,md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.13.8
kernelspec:
  display_name: Python 3 (phys-521)
  language: python
  name: phys-521
---

```{code-cell}
:tags: [hide-cell]

import mmf_setup;mmf_setup.nbinit()
%matplotlib inline
import numpy as np, matplotlib.pyplot as plt
```

# Assignment 4: Small Oscillations

**Due: 11:59pm Friday 28 October 2022**

+++

## Masses and Springs 1


Consider an equiliateral triangle of three equal masses $m$ and springs with spring
constant $k$ and equilibrium length $l$ as discussed in class.  Consider this free to
move in three dimensions (i.e. on the International Space Station).

1. How many normal modes do you expect for this system?
2. How many different frequencies do you expect?
3. What would be the degeneracy of these frequencies?  Describe the modes with a sketch.
4. Compute at least one non-zero frequency.

+++

### Solution

:::::{solution}
:show:

1. There should be $9 = 3\times 3$ modes since there are three particles in three
   dimensions.  This should be compared with the $6=3\times 2$ modes in 2D.  Note: these
   6 modes should remain since planar motion is still allowed in 3D. Thus, we are
   looking for an addition 3 modes that take us out of the plane.  A little thought
   should convince you that none of these new out-of-the plane modes should have a
   finite frequency since this motion is orthogonal to the spring stretching, and hence
   would not cause a linear increase in the spring force.  Clearly there is 1 additional
   translation zero-mode, so the other two must be rotations.
   
   :::{margin}
   
   It may not be obvious that there are 3 rotational modes -- in particular, one might
   worry that there are only 2 linearly independent modes.  The counting here can help
   resolve this.  We should certainly have all 6 modes from the 2D case plus an
   additional translational zero mode.  This leaves $2 = 9-6-1$ additional modes which
   are the two rotational modes.  The only other possibility would be one extra
   rotational mode and another single non-zero frequency, but the only non-zero
   frequency modes must be planar oscillations.
   :::
2. Thus, there should be 3 translation and 3 rotational modes with zero frequencies.  The
   remaining 3 modes are the same planar modes discussed in class: one breathing mode,
   and then two linearly independent modes with a different frequency.
3. In order words, there should be 3 different frequencies, $\omega = 0$ ($\times 6$),
   $\omega_1$ ($\times 1$), and $\omega_2$ ($\times 2$). 
4. The equilibrium state has the particles sitting in an equilateral triangle with sides
   of length $l$.  We can represent this in the complex plane as having the particles at
   the points $x_{n} = re^{2\pi\I n/3}$.  The length of the springs is: 

   \begin{gather*}
     r\abs{e^{2\pi \I/3} - e^{-2\pi \I/3}} = 2r\sin\frac{2\pi}{3} 
     = r\sqrt{3}.
   \end{gather*}

   The mode can be described in terms a dilation where each particle moves away-from or
   towards the origin in synchrony.  The Lagrangian can thus be expressed as:
   
   \begin{gather*}
     L[r,\dot{r}] = \frac{3}{2}m\dot{r}^2 - \frac{3k}{2}\left(r\sqrt{3} - l\right)^2
     = \frac{3}{2}m\dot{r}^2 - \frac{9k}{2}(r - r_0)^2,
   \end{gather*}
   
   where $r_0$ is the displacement at the equilibrium length of the springs.  This is
   quadratic in the displacement $\delta = r-r_0$ so we immediately have the frequency:
   
   \begin{gather*}
     \omega_1^2 = \frac{k_{\mathrm{eff}}}{m_{\mathrm{eff}}}
     = \frac{9k}{3m}
     = \frac{3k}{m}.
   \end{gather*}
   
   The other mode is a little more difficult to compute but has frequency
   
   \begin{gather*}
     \omega_2^2 = \frac{3k}{2m}.
   \end{gather*}

:::{admonition} Comment 1
Some people made a mistake computing the potential energy.  The potential energy
of a spring is:

\begin{gather*}
  \frac{k}{2}(r-l)^2
\end{gather*}

where $r$ is the length of the spring and $l$ is the equilibrium length.  For two points
in 2D we have $r = \sqrt{(x_1-x_2)^2 + (y_1-y_2)^2}$, hence:

\begin{gather*}
  \frac{k}{2}(r-l)^2 = \frac{k(r^2-2rl+l^2)}{2}\\
  = \frac{k}{2}\left((x_1-x_2)^2 + (y_1-y_2)^2 - l^2 + 2l\sqrt{(x_1-x_2)^2 + (y_1-y_2)^2}\right).
\end{gather*}

Note in particular that, even under small oscillations about the equilibrium position,
this depends on the equilibrium length $l$.
:::

:::{admonition} Comment 2
Some people had the idea of introducing three coordinates $r_{1}$, $r_{2}$, and $r_{3}$
corresponding to the lengths of each spring.  This makes the computation of the
potential energy very simple:

\begin{gather*}
  V = \frac{k}{2}\bigl((r_1-l)^2 + (r_2-l)^2 + (r_3 - l)^2\bigr).
\end{gather*}

Unfortunately, computing the kinetic energy is not so simple with this approach.  In
particular, the following problems remain:

1. We have only introduced 3 coordinates, whereas we need a total of 9.  If you argue
   carefully about the zero modes, then you might think that these three coordinates are
   enough to describe the three non-zero modes.  They are, but this argument is quite
   subtle.
2. To get the kinetic energy, one needs to compute the actual motion of the three
   particles.  This requires some complicated trigonometry, as well as ensuring that the
   center of mass does not move so that the previous argument removing the 6 zero modes
   remains valid.
   
This approach can be used for the breathing mode where one can argue that $r=r_1=r_2=r_3$
and that the distance from the center $d = r/\sqrt{3}$, hence the
Lagrangian in terms of the deviation $\delta = r-l$ is

\begin{gather*}
  L = 3\frac{m}{2}\dot{d}^2 - 3\frac{k}{2}(r-l)^2
    = \frac{m}{2}\dot{r}^2 - 3\frac{k}{2}(r-l)^2
    = \frac{m}{2}\dot{\delta}^2 - 3\frac{k}{2}\delta^2.
\end{gather*}

Hence, $\omega^2 = (3k/2) / (m/2) = 3k/m$ as we found above.

:::

:::{admonition} Comment 3

Several people describe the modes using a diagram like the following:

```{glue:} rotation_mode_fig
```

While this might have a germ of a correct argument which would follow along the lines
of: "The modes must form a representation of the symmetry group of the object, and there
is a 3-dimensional representation that characterizes this mode..." but this argument is
quite tricky to make.

As it stands, however, it has a fundamental problem: the equilibrium state is
**different** in each picture.  In the first, the equilibrium triangle lies in the
$x$-$y$ plane, while in the second, it lies in the $x$-$z$ plane.  This is not allowed:
when considering normal modes, you must use the **same** equilibrium state, and look at
fluctuations about this state.

Consider the same argument made about the breathing mode: one might be tempted to
conclude that there are 3 breathing modes, but there is only one.
:::
:::::

```{code-cell} ipython3
:tags: [remove-cell]

from myst_nb import glue

fig, axs = plt.subplots(1, 3, figsize=(3*1.3+1, 1.3))

for ax in axs:
    r = 0+0j
    dr = 1.0 + 0j
    R = np.exp(2j*np.pi/3)
    v = -0.2*R
    for n in range(3):
        ax.plot([r.real], [r.imag], 'oC0')
        ax.arrow(r.real, r.imag, v.real, v.imag, width=0.02)
        r += dr; dr *= R; v *= R

axs[0].set(xlabel='x', ylabel='y', aspect=1, xticks=[], yticks=[])
axs[1].set(xlabel='x', ylabel='z', aspect=1, xticks=[], yticks=[])
axs[2].set(xlabel='y', ylabel='z', aspect=1, xticks=[], yticks=[]);

glue("rotation_mode_fig", fig, display=False);
```
+++

## Masses and Springs 2

Repeat the process for a tetrahedron with 4 equal masses and 6 equal springs.
You may want to check your results with the code on CoCalc.

+++

### Solution

A tetrahedron in 3D has 4 masses, so we should now have $12 = 4\times 3$ modes.  The
tetrahedron spontaneously breaks translation and rotational symmetries, implying that
there should be 3 zero modes ($\omega = 0$) associated with pure translation, and 3 zero
modes corresponding to rotations (see discussion above about why there are 3 -- not 2 --
rotational modes).  This leaves 6 non-zero modes.

One of these should be a spherically symmetric breathing mode which has frequency
$\omega^2 = 4k/m$ as we calculate below.

Another mode corresponds to one face breathing symmetrically while oscillating towards
and away from the other vertex.  Note that this must be orthogonal to the breathing
mode, so that while the face is contract, it must be moving *away* from the other vertex
(unlike the breathing mode, which moves towards the other node while contracting).
There are 4 apparent modes here corresponding to the four faces, but it is possible that
they are linearly dependent.  Counting supports this hypothesis and there are 3 of
these.  This mode has frequency $\omega^2 = k/m$ but I do not have a trick for
calculating it.

There are also 3 modes where two adjacent nodes contract while the other two expand.
Again these might be linearly dependent.  Counting supports this hypothesis and there
are 2 of these.  A little thought shows that only the two springs discussed here
compress and expand in this mode -- the others just rotate.  Thus, the frequency is just
that of two masses connected by a single spring which is easily computed by noting that
the center of the spring does not move, so we have a single mass with half a spring, or
$\omega^2 = 2k/m$. *(Recall that half a spring is stiffer than a whole spring)*

Since there are only 5 remaining mode, we can conclude that the two sets discussed above
must actually be degenerate as discussed.

:::{margin}
A tetrahedron has 4 vertices and 6 edges.
:::
The breathing mode is fairly easy to calculate.  If the equilibrium length of the
springs is $l$, then, after a bit of geometry, we find that the distance from the center
of the tetrahedron to the vertex is $d = l \sqrt{3/8}$ *(we will check this below)*.
The Lagrangian is thus:

\begin{gather*}
  L = 4\frac{m}{2}\frac{3}{8}\dot{\delta}^2 - 6\frac{k}{2}\delta^2, \qquad
  \omega_{\text{breathing}}^2 = \frac{6k/2}{12m/16} = \frac{4k}{m}.
\end{gather*}

$$
  \newcommand{\vect}[1]{\mathbf{#1}}
  \newcommand{\uvect}[1]{\mathbf{\hat{#1}}}
  \newcommand{\norm}[1]{\lVert#1\rVert}
$$

Here I attempt a brute force solution to analysis the tetrahedron normal-mode problem.
To start, I will simply write the force law for the 4 points in terms of a spring force
where $l$ is the equilibrium length:

$$ \vect{F}_{i} = \sum_{j\neq i} -k (\norm{\vect{r}_{ji}} - l)\uvect{r}_{ji} =
  \sum_{j\neq i} -k (\norm{\vect{x}_i - \vect{x}_j} - l)\frac{\vect{x}_i -
  \vect{x}_j}{\norm{\vect{x}_i - \vect{x}_j}}.
$$

Our strategy is to evolve these with some damping $\Gamma$ to allow the system to relax
into the ground state.  Then we kick the system and see what modes we excite.


$$
  \ddot{\vect{x}}_{i} = \frac{\vect{F}_{i}}{m} - \Gamma \dot{\vect{x}}_{i}.
$$

[YouTube Discussion](https://youtu.be/ahJLD6_DAJc)

```{code-cell} ipython3
:tags: [hide-input]

%matplotlib inline
import numpy as np, matplotlib.pyplot as plt
from IPython.display import clear_output
from mpl_toolkits.mplot3d import Axes3D

from scipy.integrate import solve_ivp

m = 1
k = 1
l = 1


def F(x):
    """Return the array of forces.

    x = [[x_0, y_0, z_0],
         [x_1, y_1, z_1]
         [x_2, y_2, z_2],
         ...
        ]

    N = 4
    3 dimensions
    x.shape = (4, 3)
    """
    F = []
    for i in range(4):
        F_i = 0
        for j in range(4):
            if j == i:
                continue
            xij_ = x[i, :] - x[j, :]
            xij = np.linalg.norm(xij_)
            F_i += -k * (xij - l) * xij_ / xij
        F.append(F_i)
    return np.array(F)


def energy(q, dq):
    """Return the total energy of the system."""
    T = m / 2 * (dq**2).sum()
    V = 0
    for i in range(4):
        for j in range(4):
            if j == i:
                continue
            qij_ = q[i, :] - q[j, :]
            qij = np.linalg.norm(qij_)
            V += k * (qij - l)**2 / 2.0
    return T + V


def rhs(t, x, gamma=0):
    """RHS of the set of ODEs.

    Arguments
    ---------
    gamma : float
       Damping coefficient.
    """
    q, dq = x.reshape(2, 4, 3)
    ddq = F(q) / m - gamma * dq
    return np.array([dq, ddq]).ravel()


def plot(x, fig=None, ax=None):
    """Create a 3D plot of the molecule."""
    q, dq = x.reshape(2, 4, 3)
    xs, ys, zs = q.T
    if ax is None:
        if fig is None:
            fig = plt.gcf()
        ax = fig.add_subplot(111, projection='3d')

    ax.scatter3D(xs, ys, zs)
    order = [0, 1, 2, 0, 3, 1, 2, 3]
    ax.plot(xs[order], ys[order], zs[order])
    return fig, ax


rng = np.random.default_rng(seed=2)
x0 = rng.random(2 * 4 * 3)

ax = None
fig = None

dt = 1.0
x = 1 * x0
E = 1
args = (1.0, )  # Damping to find ground state

t0 = 0
while abs(E) > 1e-8:
    t1 = t0 + dt
    res = solve_ivp(rhs, t_span=(t0, t1), y0=x, args=args)
    x = res.y[:, -1]
    t0 = t1
    global ax, fig
    if ax:
        ax.clear()
    fig, ax = plot(x, fig=fig, ax=ax)
    q, dq = x.reshape(2, 4, 3)
    E = energy(q, dq)
    plt.title("E={}".format(E))
    display(fig)
    clear_output(wait=True)

x0 = x
```

Note that the final energy is $E=0$ as we would expect since the springs reach their
equilibrium length and there is no kinetic energy.

Next we compute the distance from the center of mass to a vertex and verify that $d/l = \sqrt{3/8}$

```{code-cell} ipython3
qs = x0.reshape(2, 4, 3)[0]
q_cm = qs.sum(axis=0)/4
ds = np.linalg.norm(qs - q_cm[np.newaxis, :], axis=-1)
assert np.allclose((ds/l)**2, 3/8, rtol=1e-4)
```

Now we give the system a random
displacement and perform a Fourier analysis of the resulting motion.

```{code-cell} ipython3
rng = np.random.default_rng(seed=2)
dx = rng.random(x.shape)
dx[4 * 3:] = 0
x1 = x0 + 0.001 * dx
#x1 = 1.01*x0  # Exite dilation mode mostly
T = 500.0
Nt = 2**9
ts = np.linspace(0, T, Nt)
res = solve_ivp(rhs, t_span=(0, T), y0=x1, t_eval=ts)
xs = res.y.T
plt.plot(ts, xs)
plt.ylim(0.714, 0.718)
plt.xlim(0, 100);
```

```{code-cell} ipython3
for x_ in xs.T:
    plt.psd(
        x_,
        NFFT=Nt,  # Use all the points for the FFT
        Fs=Nt / T,  # This is the sample rate in samples/s
        detrend='linear',  # Remove the linear trends (f=0 modes)
    )

# Here are the analytic normal mode frequencies:
ws = np.sqrt(np.array([k / m, 2 * k / m, 4 * k / m]))
fs = ws / (2 * np.pi)  # Linear frequencies
for f in fs:
    plt.axvline(f, color='y')
```

See also [Normal Modes.ipynb].

[Normal Modes.ipynb]: <https://cocalc.com/projects/31c120c9-2956-4420-9d6f-374a6ee32df3/files/Normal%20Modes.ipynb>

## 2 Masses

Provide a complete solution describing the motion of two coupled oscillators as shown
below that move in 1D (left-right) without friction etc.

![Two masses and springs](figSprings2m.png)

Be sure to check various limiting cases to make sure your answer is
correct, such as equal masses and spring constants, setting some spring constants to
zero, etc.

:::{margin}
Hint: Your final answer will probably be simpler if your replace your masses and spring
constants with the frequencies

\begin{align*}
  \omega_{j}^2 &= \frac{k_{j} + k_{12}}{m_{j}}, \\
  \omega_{12}^2 &= \frac{k_{12}}{\sqrt{m_1m_2}}.
\end{align*}

:::

Use your answer to characterize how energy is transferred between two harmonic
oscillators.  I.e. consider the case where the middle spring has a very small spring
constant $k_{12} \ll k_{1, 2}$.  This corresponds to two weakly coupled harmonic
oscillators.  Suppose you excite the left oscillator with angular frequency $\omega_1
\approx \sqrt{k_1/m_1}$: what conditions must be true in order that a significant
portion of the energy can be transferred to the second oscillator with angular frequency
$\omega_2 \approx \sqrt{k_2/m_2}$?

Be sure to apply the checklist described at the start of the course when you consider
the completeness of your solution: 

1. Is the problem well defined and well formulated?  (Does it have a definite answer?
   Could you tell a computer how to solve this problem?)
2. How will you describe the solution?  (What variables will you use?  What parameters
   are required?)  Note: the solution will often appear much simpler if you describe the
   solution with appropriate variables.
3. What is your intuition for the problem? What do you think the system will do? Record
   your predictions before you start solving the equations.
4. Formulate the equations that will distinguish the physical solution from all other
   possible behaviors admitted by your description.  What physics do you include?  What
   effects do you ignore?  Justify briefly.  (If you were to setup an experiment to test
   your results, would you need to consider the effect of Jupiter's gravitational field?
   Why or why not?)
5. Solve the equations to get the solution.  Consider different strategies here as some
   approaches will be much cleaner and more efficient than others.  I tend to start
   working quickly to see how things progress, then once I see my way through to a
   portion of the solution, I often start over again using a more streamlined approach
   that will ensure I make fewer mistakes with my algebra etc. Some suggestions:
   1. Identify all dimensionless parameters.  The qualitative behavior of the system will
      depend (only) on these.
   2. Are there any simplifying limits where you can quickly solve the problem?  If you
      get stuck, these might help you find the more general solution.
   3. If finding an analytic solution is challenging, maybe try quickly implementing a
      numerical solution so you can build a better intuition for what should happen.
      *This introduces the added complication of making sure your numerical
      implementation is correct, however, I find I often learn most about a problem
      while trying to making the numerical implementation agree with my analytic
      calculations, so this extra effort is usually not wasted.*

6. Check your solution:
   1. Does it agree with your intuition?  If not, you need to either find your mistake,
      or possibly change your mental model of how this system works (once you verify
      that your solution is indeed correct).
   2. Is it dimensionally correct? 
   3. Are there limiting case you can check easily? 
   4. Can you compare with a colleague who independently solved the problem?
   5. Can you solve it a different way to compare? (Maybe with different variables or in
      a different coordinate system?)
   6. Can you do an experiment or relate this to a physical system where you know what
      happens?

7. Generalize: Sometimes once you have solved a problem, you can generalize your
   solution.  In this case, do you see how you could solve a generalized problem with 3,
   4, or more masses? Does your method generalize to the cased of different masses or
   different spring constants?  Even if you do not actually solve the generalized
   problem, thinking about this can help you find more efficient techniques for solve
   similar problems in the future.

### Solution

:::{admonition} Strategy

1. Get the equations of motion using a Lagrangian.

   *I may need to be careful about working with springs because the equilibrium length
   can play an important role, but in 1D this should be able to be neglected.*
2. Find the equilibrium points.

   *In this case, I can probably choose my coordinates so the equilibrium corresponds to
   $x_1=x_2=0$.*
3. Expand the equations of motion about the equilibrium points.

   *As these are springs, the equations will probably already be linear, so there will
   be no need to expand.*
4. Form a matrix equation, and diagonalize to get the normal modes.

   *Since we are interested in small $k_{12}$, these might simplify if we expand in
   $k_{12}$.  Will have to see.*

To figure out how efficiently energy is transferred to the second spring, I plan to
start from rest with $x_1(0) = \epsilon$ and $x_{2}(0) = 0$.  This system will have
total energy $E = (k_1 + k_{12})x_1^2(0)/2 \approx k_1\epsilon^2/2$. (The approximation
is neglecting the energy of the middle spring which is small in the weakly coupled
case.)  I then plan to look at the maximum displacement of $x_2(t)$ at some later time,
so I can estimate the maximum amount of energy transferred:

$$
  \frac{E_2}{E} \approx \frac{k_2(\max x_2)^2}{k_1\epsilon^2}.
$$

To do this, we need to express the initial state as a sum of the normal modes.
Formally, we can solve the problem with matrix as follows:

$$
  \diff[2]{\vect{x}}{t} = -\mat{M}\cdot\vect{x}, \qquad
  \vect{x}(t) = \cos(\mat{M} t)\cdot \begin{pmatrix}
    \epsilon\\
    0
  \end{pmatrix}
$$

but we need to compute the cosine of a matrix which can be done by diagonalizing.  It is
probably simple to just compute the eigenvectors and the solve the equations.

I expect that, when $k_{12}$ is small, the situation will be analogous to that of the
one harmonic oscillator with frequency $\approx \omega_1$ acting as an external source
for the other oscillator with frequency $\approx \omega_2$.  From the solution to a
driven harmonic oscillator, we know that one only has significant amplification when
$\omega_2 \approx \omega_1$, so I expect efficient transfer only close to this resonance.

:::

As discussed in class, and as is easily checked, for motion in 1D, we need not worry
about the equilibrium length of the springs.  We can thus introduce the two coordinates
$x_1$ and $x_2$ as the displacements from equilibrium, leading to the following
Lagranian and equations of motion:

:::{margin}
As an example of how numerical checks helped me, I originally forgot to include the
diagonal contribution of $k_{12}$, writing the matrix as:

$$
  \begin{pmatrix}
    -k_1 & k_{12}\\
    k_{12} & -k_2
  \end{pmatrix}.
$$

This gave an answer that did not make sense. I quickly found the mistake by running the
26 lines of code below (not including plotting) to find that energy was not conserved,
indicating an error in my formulations.  For many problems, it is very fast to write a
simple code to check your algebra and this can be an efficient use of your time. 
:::
\begin{gather*}
  L = \frac{m_1\dot{x}_1^2 + m_2\dot{x}_2^2}{2} - \frac{k_1x_1^2 + k_2x_2^2 +
  k_{12}(x_2-x_1)^2}{2},\\
  \begin{aligned}
    m_1 \ddot{x}_1 &= -k_1 x_1 - k_{12}(x_1 - x_2),\\
    m_2 \ddot{x}_2 &= -k_2 x_2 - k_{12}(x_2 - x_1).
  \end{aligned}\\
  \begin{pmatrix}
    m_1\\
    & m_2
  \end{pmatrix}
  \cdot
  \begin{pmatrix}
    \ddot{x}_1\\
    \ddot{x}_2
  \end{pmatrix}
  =
  \begin{pmatrix}
    \overbrace{-k_1 - k_{12}}^{-\tilde{k}_1} & k_{12}\\
    k_{12} & \underbrace{-k_2 - k_{12}}_{-\tilde{k}_2}
  \end{pmatrix}
  \cdot
  \overbrace{
  \begin{pmatrix}
    x_1\\
    x_2
  \end{pmatrix}}^{\vect{x}}.
\end{gather*}

At this point, we implement a simple numerical solution to make sure that these make
sense.  We use the variable:

$$
  \vect{y} = \begin{pmatrix}
    x_1\\
    x_2\\
    \dot{x}_1\\
    \dot{x}_2
  \end{pmatrix}.
$$

```{code-cell} ipython3
:tags: [hide-input]

from scipy.integrate import solve_ivp

rng = np.random.default_rng(seed=2)
k1, k2, k12, m1, m2 = 1.0 + rng.random(5)
k12 *= 0.1 # make coupling weak.

k12 *= 3
m2 /= 2

def rhs(t, y):
    x1, x2, dx1, dx2 = y
    # ddx1 = (-k1 * x1 + k12 * x2) / m1  # Bug!
    # ddx2 = (-k2 * x2 + k12 * x1) / m2  # Bug!
    ddx1 = (-k1 * x1 - k12 * (x1 - x2)) / m1
    ddx2 = (-k2 * x2 - k12 * (x2 - x1)) / m2
    return (dx1, dx2, ddx1, ddx2)

y0 = (1.0, 0.0, 0.0, 0.0)

T = 40.0

res = solve_ivp(rhs, t_span=(0, T), y0=y0, atol=1e-8, rtol=1e-8)
t = res.t
x1, x2, v1, v2 = res.y
E1 = m1*v1**2/2 + k1*x1**2/2
E2 = m2*v2**2/2 + k2*x2**2/2
E12 = k12 * (x2-x1)**2/2
E = E1 + E2 + E12

# Check energy conservation.
assert np.allclose(E, E[0])

fig, ax = plt.subplots(figsize=(8, 3))
l1, = ax.plot(t, x1, label="$x_1$")
l2, = ax.plot(t, x2, label="$x_2$")
ax.set(xlabel="$t$", ylabel="$x$")
ax1 = ax.twinx()
ax1.plot(t, E1/E[0], "--", c=l1.get_c(), label="$E_1$")
ax1.plot(t, E2/E[0], "--", c=l2.get_c(), label="$E_2$")
ax1.plot(t, E/E[0], "--k", label="$E$")
ax1.plot(t, E12/E[0], ":C4", label="$E_{12}$")
ax1.set(ylabel="$E_{i}/E$")
ax.legend()
ax1.legend();
```

Let $\tilde{k}_{1,2} = k_{1,2} + k_{12}$, then the frequencies of the normal modes have the form:

\begin{gather*}
  \DeclareMathOperator{\Re}{Re}
  \vect{x}(t) = \Re(\vect{A} e^{\I \omega t}), \qquad
  \det\left\vert
    \begin{matrix}
      -\omega^2 m_1 + \tilde{k}_1 & -k_{12}\\
      -k_{12} & -\omega^2 m_2 + \tilde{k}_2
    \end{matrix}
  \right\vert = 0,\\
  \omega^4 m_1 m_2 - (m_2\tilde{k}_1 + m_1\tilde{k}_2) \omega^2 + \tilde{k}_1\tilde{k}_2 - k_{12}^2 = 0,\\
  \omega^2 = \frac{m_2\tilde{k}_1 + m_1\tilde{k}_2}{2m_1m_2} \pm 
      \sqrt{\frac{(m_2\tilde{k}_1 - m_1\tilde{k}_2)^2 + 4m_1m_2k_{12}^2)}{4m_1^2m_2^2}}.
\end{gather*}

This can be simplified by expanding and introducing some angular frequency variables:

\begin{gather*}
    \omega_{1} = \sqrt{\frac{k_1+k_{12}}{m_1}}, \qquad
    \omega_{2} = \sqrt{\frac{k_2+k_{12}}{m_2}}, \qquad
    \omega_{12}= \sqrt{\frac{k_{12}}{\sqrt{m_1m_2}}},\\
  \omega^2_{\pm} = \underbrace{\frac{\omega_1^2 + \omega_2^2}{2}}_{\omega_0^2} \pm 
    \sqrt{\left(\frac{\omega_1^2 - \omega_2^2}{2}\right)^2 + \omega_{12}^4}, \qquad
    \omega_0^2 = \frac{\omega_1^2 + \omega_2^2}{2}.
\end{gather*}

As a check, we can consider some limits:

1. $k_{12} = 0$.  Here we should have two independent springs with frequencies
   $\omega_{i}^2 = k_i/m_i$, which we do.
2. $k_1=k_2=0$ and $m_1=m_2=m$.  This is two masses connected by a spring.  There should
   be one zero mode corresponding to translation, and one mode with $\omega^2 = 2k/m =
   2\omega_0^2$ as discussed above.  In this cases $\omega_1 = \omega_2 = \omega_{12} =
   \omega_0$ and $\omega_{\pm}^2 = \omega_0^2 \pm \sqrt{\omega_0^4}$, as expected.
3. If $m_1 \rightarrow \infty$, we should have one mode $\omega^2 \rightarrow 0$
   and the other as a single mass $m_2$ oscillating with two
   springs: $\omega^2 = (k_2 + k_{12})/m_2 = \omega_2^2$.  In this case $\omega_1, \omega_{12}
   \rightarrow 0$, so $\omega_-^2 \rightarrow 0$ and $\omega_+^2 \rightarrow \omega_2^2$
   as expected.
4. If $k_1 = k_2 = k_{12}$ and $m_1 = m_2$ we have $\omega_1 = \omega_2 = \sqrt{2}
   \omega_{12} = \sqrt{k/m}$ and so $\omega_{\pm} = 2\omega_{12}^2 \pm \omega_{12}^2$.
   This corresponds to one mode where the masses move in synchrony so that the middle
   spring is never compressed $\omega_{-}^2 = k/m$, and another with $\omega_{+}^2 = 3k/m$.

:::{margin}
There are many ways to make dimensionless parameters.  You need to play around a bit to
see which ones will be most useful.  Here we want to consider a weak coupling, so
something proportional to $k_{12}$ is reasonable.  The second parameter appears
naturally after solving the problem.
:::
To characterize the problem, we note that the following dimensional parameters can be
formed:

\begin{gather*}
  \delta_{12} = \frac{\omega_{12}^2}{\omega_0^2} 
              = \frac{k_{12}/\sqrt{m_1m_2}}{\omega_0^2},\qquad
  \eta = \frac{\omega_{1}^2 - \omega_{2}^2}{\omega_1^2 + \omega_2^2},\\
  \frac{\omega_2^2}{\omega_0^2} = 1 - \eta,\qquad
  \frac{\omega_1^2}{\omega_0^2} = 1 + \eta,
\end{gather*}

In terms of these, we have

\begin{gather*}
  \omega^2_{\pm} = \omega_0^2\Bigl(1 \pm \sqrt{\eta^2 + \delta_{12}^2}\Bigr).
\end{gather*}

To answer the question about how much energy can be transferred from one mode to
another, however, we need to know what the eigenvectors look like:

:::{margin}
Once we insert the expression for the frequencies $\omega_{\pm}$, the bottom two
equations are equivalent.  We will use the first one for convenience.
:::
\begin{gather*}
  \begin{pmatrix}
    -\omega^2_{\pm} m_1 + \tilde{k}_1 & -k_{12}\\
    -k_{12} & -\omega^2_{\pm} m_2 + \tilde{k}_2
  \end{pmatrix}
  \cdot
  \begin{pmatrix}
    x^{\pm}_1\\
    x^{\pm}_2
  \end{pmatrix} = 0, \\
  x^{\pm}_2 = \frac{\tilde{k}_1 - \omega^2_{\pm} m_1}{k_{12}}x^{\pm}_1, \qquad
  x^{\pm}_1 = \frac{\tilde{k}_2 - \omega^2_{\pm} m_2}{k_{12}}x^{\pm}_2.
\end{gather*}

Let's assume there are no degeneracies, so we can write the (un-normalized) modes as:

\begin{gather*}
  \ket{\pm} = 
  \begin{pmatrix}
    x^{\pm}_1\\
    x^{\pm}_2
  \end{pmatrix}
  =
  \begin{pmatrix}
    1\\
    \alpha_{\pm}
  \end{pmatrix}, \\
  \alpha_{\pm} = \frac{\tilde{k}_1 - \omega^2_{\pm} m_1}{k_{12}}
  = \sqrt{\frac{m_1}{m_2}}
    \frac{\overbrace{\eta \mp \sqrt{\eta^2 + \delta_{12}^2}}^{\frac{\omega_1^2 - \omega_{\pm}^2}{\omega_0^2}}}
         {\delta_{12}}.
\end{gather*}

Starting at rest with oscillations of amplitude $\epsilon$ only in $x_1$, we look for coefficients
$a_{\pm}$ such that: 

\begin{gather*}
  \ket{+}a_{+} + \ket{-}a_{-} =
  \begin{pmatrix}
    1\\
    0
  \end{pmatrix}, \qquad
  a_{+} = \frac{-\alpha_{-}}{\alpha_{+} - \alpha_{-}}, \qquad
  a_{-} = \frac{\alpha_{+}}{\alpha_{+} - \alpha_{-}}.
\end{gather*}

The time-evolution is thus:

\begin{align*}
  \frac{x_1(t)}{\epsilon}
  &= a_{-}\cos(\omega_- t) + a_{+}\cos(\omega_+ t)\\
  &=\frac{\alpha_{+}\cos(\omega_- t) - \alpha_{-}\cos(\omega_+ t)}{\alpha_{+} - \alpha_{-}},\\
  \frac{x_2(t)}{\epsilon} 
  &= a_{-}\alpha_{-}\cos(\omega_- t) + a_{+}\alpha_{+}\cos(\omega_+ t), \\
  &= \frac{\alpha_{+}\alpha_{-}}{\alpha_{+} - \alpha_{-}}
     \Bigl(\cos(\omega_- t) - \cos(\omega_+ t)\Bigr).
\end{align*}

Unless $\omega_{+} = \omega_{-}$, $x_2$ will have a maximum amplitude of

\begin{gather*}
  \max \frac{x_2}{\epsilon} =
  \frac{2\alpha_{+}\alpha_{-}}{\alpha_{+} - \alpha_{-}}
  =
  \sqrt{\frac{m_1}{m_2}}\frac{\delta_{12}}{\sqrt{\eta^2+\delta_{12}^2}}
\end{gather*}

We now look at the energy transfer.  If the coupling is small, then the potential energy
will be $E \approx kx^2/2$ where $x$ is the maximum extension.  Thus, the maximum amount
of energy transferred to the second mass will be:

$$
  \frac{E_{2}}{E} 
   \approx \frac{k_2 (\max x_2)^2}{\tilde{k}_1 \epsilon^2}
   \approx \frac{\omega_2^2}{\omega_1^2}\frac{\delta_{12}^2}{\eta^2 + \delta_{12}^2}
   = \frac{1-\eta}{1+\eta}\frac{\delta_{12}^2}{\eta^2 + \delta_{12}^2}.
$$

Here is what this looks like for a variety of values of $\eta \in (-1, 1)$:

```{code-cell} ipython3
:tags: [hide-input]

eta = np.linspace(-1, 1, 200)[1:-1]
fig, ax = plt.subplots()
for d12 in [0.01, 0.1, 0.2, 0.3]:
    l, = ax.plot(eta, (1-eta)/(1+eta)*d12**2/(eta**2+d12**2), 
                 label=f"$\delta_{{12}}^2={d12}$")
    #kw = dict(c=l.get_c(), ls=":")
    #ax.plot(eta, (1-eta-d12)/(1+eta)*d12**2/(eta**2+d12**2), **kw)

ax.legend()
ax.set(xlabel=r"$\eta$", ylabel=r"$\max E_2/E$", ylim=(0, 1));
```

```{code-cell} ipython3
:tags: [hide-cell]

# Here is some code that explores the troubling regime.  We find that
# there in this regime, there is no good way to partition the energy.
from scipy.integrate import solve_ivp

xi = -0.9
d12 = np.sqrt(0.3)
eta = xi*np.sqrt(1-d12**2)
w0 = 1.0
w1, w2 = w0*np.sqrt([1+eta, 1-eta])
m1 = 10.0
m2 = 2.0
k12 = d12*w0**2 * np.sqrt(m1*m2)
k1, k2 = w1**2*m1 - k12, w2**2*m2 - k12
print(k1, k2, k12)

def rhs(t, y):
    x1, x2, dx1, dx2 = y
    # ddx1 = (-k1 * x1 + k12 * x2) / m1  # Bug!
    # ddx2 = (-k2 * x2 + k12 * x1) / m2  # Bug!
    ddx1 = (-k1 * x1 - k12 * (x1 - x2)) / m1
    ddx2 = (-k2 * x2 - k12 * (x2 - x1)) / m2
    return (dx1, dx2, ddx1, ddx2)

y0 = (1.0, 0.0, 0.0, 0.0)

T = 40.0

res = solve_ivp(rhs, t_span=(0, T), y0=y0, atol=1e-8, rtol=1e-8)
t = res.t
x1, x2, v1, v2 = res.y
E1 = m1*v1**2/2 + (k1+k12)*x1**2/2
E2 = k2*x2**2/2
E12 = m2*v2**2/2 + k12*x2**2/2 - k12 * x1 *x2
E = E1 + E2 + E12

# Check energy conservation.
assert np.allclose(E, E[0])

fig, ax = plt.subplots(figsize=(8, 3))
l1, = ax.plot(t, x1, label="$x_1$")
l2, = ax.plot(t, x2, label="$x_2$")
ax.set(xlabel="$t$", ylabel="$x$")
ax1 = ax.twinx()
ax1.plot(t, E1/E[0], "--", c=l1.get_c(), label="$E_1$")
ax1.plot(t, E2/E[0], "--", c=l2.get_c(), label="$E_2$")
ax1.plot(t, E/E[0], "--k", label="$E$")
ax1.plot(t, E12/E[0], ":C4", label="$E_{12}$")
ax1.set(ylabel="$E_{i}/E$")
ax.legend()
ax1.legend();
```

The maximum near $\eta = 0$ is reasonable -- I expect maximal transfer when the
frequencies are resonant $\omega_2 \approx \omega_1$.  It is a little disturbing that
the value exceeds $1$, however, I suspect this is because of the approximation $k_2
\rightarrow \tilde{k}_2$ in the numerator.  The behaviour for small $\eta \approx -1$,
however, is puzzling.  This comes from the denominator $1+\eta$.  At this point, I was
not completely certain that I did not make some sort of error, so I compared a few
expressions with the following code.


```{code-cell} ipython3
#:tags: [hide-cell]

tk1, tk2 = k1 + k12, k2 + k12
w1, w2, w12 = np.sqrt([tk1/m1, tk2/m2, k12/np.sqrt(m1*m2)])
w0 = np.sqrt((w1**2+w2**2)/2)
delta12 = 2*w12**2/w1/w2
delta12 = w12**2/w0**2
eta = (w1**2 - w2**2)/(w1**2 + w2**2)
print(eta, delta12)
pm = np.array([+1, -1])
wpm = wp, wm = w0 * np.sqrt(1 + pm * np.sqrt(eta**2 + delta12**2))
alpha_pm = alpha_p, alpha_m = (tk1 - wpm**2 * m1) / k12
ap, am = a_pm = np.divide([-alpha_m, alpha_p], alpha_p - alpha_m)

Ap = np.array(
    [[-wp**2*m1 + tk1, -k12],
     [-k12, -wp**2*m2 + tk2]])
Am = np.array(
    [[-wm**2*m1 + tk1, -k12],
     [-k12, -wm**2*m2 + tk2]])
assert np.allclose(0, np.linalg.det(Ap))
assert np.allclose(0, np.linalg.det(Am))
assert np.allclose(delta12, k12/np.sqrt(m1*m2)/w0**2)
# plt.plot(t, x1)
# plt.plot(t, am*np.cos(wm*t) + ap*np.cos(wp*t))

assert np.allclose(
    (w1**2-wpm**2)/w0**2,
    eta - pm*np.sqrt(eta**2 + delta12**2))
assert np.allclose(
    alpha_pm,
    np.sqrt(m1/m2)*(eta - pm*np.sqrt(eta**2 + delta12**2)) / delta12)

assert np.allclose(
    2*alpha_p*alpha_m/(alpha_p - alpha_m),
    np.sqrt(m1/m2)*delta12/np.sqrt(eta**2+delta12**2))

assert np.allclose([(w1/w0)**2, (w2/w0)**2], 
                   [1+eta, 1-eta])

print("The following are my numerical vales vs my formula")
print("x2 max:", x2.max(), 2*alpha_p*alpha_m/(alpha_p - alpha_m))
print("E2/E:", E2.max()/E.max(), (1-eta)/(1+eta) * delta12**2 / (eta**2 + delta12**2))
print("E2/E:", k2*(x2.max())**2/k1/(x1.max())**2,
              (1-eta)/(1+eta) * delta12**2 / (eta**2 + delta12**2))
eps = 1.00
print(E.max(), tk1*eps**2/2, E2.max(), tk2*alpha_p*alpha_m/(alpha_p - alpha_m))
print(E.max(), tk1*eps**2/2, E2.max(), k2*alpha_p*alpha_m/(alpha_p - alpha_m))
```

I found several bugs, for example
($\delta_{12}$ vs $\delta_{12}^2$ in some places, missing factors of 2 etc.  Finally, I
use the numerical solution to find the maximum extent $x_2$ and print out the
numbers. This is close to the value given above, so I think I have done things correctly.

So what does this mean?  Take a moment to see if you can reason about what is happening
here?

:::{admonition} Spoiler
:class: dropdown

While it appears that $\eta$ can take values in $(-1, 1)$, note that we actually have a
more restrictive criterion that $\eta^2 + \delta_{12}^2 < 1$, otherwise one of the
frequencies $\omega_-$ becomes imaginary indicating an instability.

Keeping $k_2$ rather than $\tilde{k}_2$ we can obtain the bound:

$$
  \frac{E_{2}}{E} 
   \lessapprox \frac{k_2 (\max x_2)^2}{\tilde{k}_1 \epsilon^2}
   = \frac{(\omega_2^2 - k_{12}/m_2)}{\omega_1^2}\frac{\delta_{12}^2}{\eta^2 + \delta_{12}^2}\\
   = \frac{(1 - \eta - \delta_{12}\sqrt{m_1/m_2})}{1+\eta}\frac{\delta_{12}^2}{\eta^2 +
  \delta_{12}^2}
$$

It might seem that you can get back to the original problem by making $m_1/m_2$ small,
but there are non-trivial constraints on the parameters such that you can't
independently adjust everything while holding $\eta$ and $\delta_{12}$ fixed.  This
needs a bit more exploration, and perhaps a different choice of dimensionless
parameters, but after playing for a bit, I am pretty convinced that this is correct.
:::

## (Optional) Double Pendulum

Consider a double pendulum as shown below:

```{code-cell}
:tags: [hide-input]

import sympy
sympy.init_printing()
from phys_521 import double_pendulum
double_pendulum.draw_figure()
```

The formal task here is to characterize and quantify the normal modes of this system,
but I would like you to take this further and perform a thorough analysis of this
system, investigating both analytic features (such as the normal modes) and features
that are only apparent numerically.  Upon request, I will provide code examples to help
you analyze this system, but would like you to first try to solve it from scratch
yourself.  (Please use the notation given above so that we can compare answers.) 

Here are some recommendations for analyzing this system.  Note, I recommend first
finding out as much information as possible using techniques that take very little of
**your** time and effort.  These quick solutions might be approximate, or very slow to
run on a computer but get you something you can use to check your more refined solutions
later.  Once you have a set of quick and dirty results that you know are correct, you
can start to implement a more elegant or efficient solution, checking it against the
other solutions.  *As an example, consider performing a nasty integral - the integral
might admit a nice analytic form after performing some clever substitutions, but you
might first try simply approximating it numerically to see if you can get an idea of the
answer.  Then as you implement the more sophisticated transformations, you can keep
checking to make sure that you have not make any mistakes.* 

**Brief Exploration:** Spend a few minutes thinking about the problem, maybe writing a
Lagrangian, getting familiar with the coordinates and structure.  The equations for this
system are moderately complicated, so correctly deriving them and implementing them can
be prone to errors, thus it is crucial that you think about how you can check your
answers.  I will suggest some general approaches here: not all will be relevant or
useful for this problem.

**Dimensionless Parameters:** Identify the dimensionless parameters of the problem.
What will set the scale for periods?  When will the motion be simple, when will it be
complicated?  In your normal mode analysis, what will constitute "small"?  These
parameters will characterize the different qualitative behaviors of the system, allowing
you to classify them.

**Limiting Cases:** Identify simple limiting cases.  Can you reduce the system to a
simpler system with appropriate limits?  Write these limits down explicitly, and if
possible, determine how the leading order corrections will scale as you move away from
the ideal limits.  This will provide a very useful way of quantitatively testing your
code.  *(One issue is that some useful limits might be impossible to implement
numerically because of divergence or singularities.  Thus your numerical implementation
will have to simulate a system close to, but not at, the limit.  Since the system cannot
be exactly situated at the limit in these cases, understanding how the errors will scale
is crucial for testing.)*

**Conservation Laws and Symmetry:** Identify conserved quantities, or relationships that
follow from symmetry.

**Brute Force Solution:** Is there a simple way to implement brute force solution?
Sometimes a problem can be formulated in a dead simple way that takes **you** almost no
time, but perhaps is very slow for the computer, or can only give an approximate answer.
If so, get a simple or limited version of your code working.

**CAS:** If your problem is messy, can you use a computer algebra system
[CAS](https://en.wikipedia.org/wiki/Computer_algebra_system) like
[SymPy](http://www.sympy.org/en/index.html) or [Sage](http://www.sagemath.org/) (both
can be used on CoCalc) to make sure your results are correct?  Sometimes these can even
be used to generate the functions used in your numerical code.  *(Warning: Although I
really like this approach in principle, I find that often I waste far more time trying
to automatically generate code than I would spend if I just sat down and did the thing
by hand once.  However, if you need to generate many similar but different pieces of
code with changing functional forms, then code generation becomes a very viable option.
Keep in mind, however, that the final functional form generated by a CAS will likely not
be optimal in terms of performance, or roundoff error.  Often hand-manipulation will
greatly improve the result, but this is where I end up wasting lots of time - trying to
automate these simplifications.)*

**Numerical Implementation:** Once you have a numerical implementation, test it
carefully against the limiting case you found earlier.  Also check that your units
etc. are correct by running the code with different scalings of the parameters, but the
same dimensionless combinations.  The results should be identical, but with different
units.

**Testing:** Check that your code obeys all the conservation laws and symmetry
properties you derived earlier.  If you find something that seems wrong, how do the
errors behave?  If they are smooth, there is probably some systematic error (missing
factor, etc.).  If they appear to fluctuate rapidly or randomly, they might be numerical
errors due to roundoff, or an algorithm not converging.  *(With
[`solve_ivp`](https://docs.scipy.org/doc/scipy/reference/reference/generated/scipy.integrate.solve_ivp.html)
or similar integrators, you might have problems if you take to large steps.  These can
be fixed by taking smaller step sizes, or increasing the equivalent of the `max_step`
parameter to allow the algorithm to do more work for each time step.)*

**Exploration:** Armed with a validated numerical solution, start exploring various
aspects of the physics.  Check your intuition, and try to test analytic results and
approximations.  If you find something surprising, can you quantify it or understand
what is happening qualitatively from the analytics?

```{margin}
[<img style="float:right; padding-left:1em" width="150"
src="https://upload.wikimedia.org/wikipedia/commons/8/87/Double_pendulum_flips_graph.png">](https://commons.wikimedia.org/wiki/File:Double_pendulum_flips_graph.png)
```

**Generalization**:
Consider generalizing the problem: are there any ways you can generalize your solution
with a minimal amount of work to extend its range of validity?  Can you test other
aspects of physics or apply your solution to related problems?  Once you have a working
numerical solution to a problem, and have put some time into understanding the behavior
of the system, it is a good idea to see if you can get some mileage out of it.  *(For
example, the system under consideration here exhibits chaotic dynamics.   It is thus
provides a great opportunity for exploring and learning about
[chaos](https://en.wikipedia.org/wiki/Chaos_theory), [Poincaré
sections](https://en.wikipedia.org/wiki/Poincar%C3%A9_map), [Lyapunov
exponents](https://en.wikipedia.org/wiki/Lyapunov_exponent), etc.:*
