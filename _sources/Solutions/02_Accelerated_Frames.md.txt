---
jupytext:
  formats: ipynb,md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.11.1
kernelspec:
  display_name: Python 3 (phys-521)
  language: python
  name: phys-521
---

# Assignment 2: Accelerated Frames

+++

## Merry-go-round

Consider a free particle moving in 2D ($x$-$y$ plane) with a constant velocity.  Find
the trajectory for this "free" particle in an accelerated frame rotating with $\omega(t)
= \alpha t$ (i.e. with angle $\theta(t) = \alpha t^2/2$).  In principle you could solve
the problem in the rotating frame, but there is an easier way to find the trajectory.
Show that your trajectory explicitly satisfies Newton's law in the accelerating frame
including all three corrections.

+++

### Solution

::::{solution}

:::{admonition} Comments

The key here is to find the transformation between the rotating frame and the inertial
frame.  Solve the problem in the inertial frame (easy $\vect{x} = \vect{x}_0 +
\vect{v}_0 t$) and then use this to transform to rotating frame where the transformed
solution should satisfy the modified force law (with coriolis, centrifugal etc. terms).
:::

The final solution is to note that the rotating coordinates are related to the inertial
coordinates by the following transformation
\begin{gather*}
  \begin{pmatrix}
    x_r(t)\\
    y_r(t)
  \end{pmatrix}
  =
  \begin{pmatrix}
    \cos \theta & \sin \theta\\
    -\sin \theta & \cos \theta
  \end{pmatrix}
  \begin{pmatrix}
    x_i(t)\\
    y_i(t)
  \end{pmatrix}.
\end{gather*}
Then one simply needs to solve the problem in the inertial frame:
\begin{gather*}
  \vect{x}_i(t) = \vect{x}_0 + \vect{v}_0 t.
\end{gather*}

To derive the transformations, recall that transformations to a rotating frame satisfy:

\begin{gather*}
  \left.\pdiff{}{t}\right|_{\text{body}} = 
  \left.\pdiff{}{t}\right|_{\text{inertial}} + \vec{\omega}\times
\end{gather*}

with the caveat that the coordinates be aligned $\vect{X}=\vect{x}$ at the time
considered. Let $\vect{x}$ be the coordinate in the inertial frame and $\vect{X}$ in the
rotating frame, then we have the following equations of motion:

\begin{gather*}
  m\ddot{\vect{x}} = \vect{F}_{\text{ext}} = \vect{0}, \qquad
  \vect{x}(t) = \vect{x}_0 + \vect{v}_0t,\\
  \begin{aligned}
    m\left(\left.\pdiff{}{t}\right|_{\text{inertial}}\right)^2\vect{x}
    &= m\left(\left.\pdiff{}{t}\right|_{\text{body}} + \vec{\omega}\times\right)^2\vect{x}\\
    &= m\left(\left.\pdiff{}{t}\right|_{\text{body}} + \vec{\omega}\times\right)
     (\dot{\vect{X}} + \vec{\omega}\times\vect{X})\\
    &= m(\ddot{\vect{X}} + \dot{\vec{\omega}}\times\vect{X} + 2\vec{\omega}\times\dot{\vect{X}}
    + \vec{\omega}\times\vec{\omega}\times \vect{X})
  \end{aligned}\\
  0 = \ddot{\vect{X}}(t) + \dot{\vec{\omega}}\times\vect{X}
                     + 2\vec{\omega}\times\dot{\vect{X}}
                     + \vec{\omega}\times(\vec{\omega}\times\vect{X}),\\
  \ddot{\vect{X}}(t) = - \dot{\omega}\uvect{z}\times\vect{X}
                       - 2\omega \uvect{z}\times\dot{\vect{X}}
                       + \omega^2\vect{X}.
\end{gather*}

I can never remember the sign in this transformation, however, so let's derive it properly:

\begin{gather*}
  \vect{X} = \mat{R}_{\theta}\cdot\vect{x} = e^{-\vec{\theta}(t)\times}\cdot\vect{x}
  = \begin{pmatrix}
      \cos\theta & \sin\theta\\
      -\sin\theta & \cos\theta
  \end{pmatrix}
  \cdot
  \vect{x},\\
  [\vec{\theta}\times]_{kj} = \theta\epsilon_{2jk} = 
  \begin{pmatrix}
   0&-\theta(t)\\
   \theta(t)&0
  \end{pmatrix}.
\end{gather*}

The definition $\vec{\theta}\times$ follows from $\vect{A}\times\vect{B} = \vect{C}$
where $C_k = \sum_i\epsilon_{ijk}A_iB_j$ and $\epsilon_{ijk}$ is the totally
antisymmetric [Levi-Civita symbol](https://en.wikipedia.org/wiki/Levi-Civita_symbol)
with $\epsilon_{012} = 1$.  To check the signs, consider the $x$ axis in the rotating
frame $\vect{X} = (1,0)$, in the inertial frame this should rotate as $\vect{x} =
(\cos\theta, \sin\theta)$ which is what we have here.  Thus, the trajectory is:

:::{math}
:class: full-width

\begin{align*}
  \vect{X}(t) &= \begin{pmatrix}
    \cos\theta & \sin\theta\\
    -\sin\theta & \cos\theta
  \end{pmatrix}\cdot(\vect{x}_0 + \vect{v}_0 t), 
  \\
  \dot{\vect{X}}(t) &= 
  \dot{\theta}\begin{pmatrix}
    -\sin\theta & \cos\theta\\
    -\cos\theta & -\sin\theta
  \end{pmatrix}\cdot(\vect{x}_0 + \vect{v}_0 t) 
  + 
  \begin{pmatrix}
    \cos\theta & \sin\theta\\
    -\sin\theta & \cos\theta
  \end{pmatrix}\cdot\vect{v}_0, 
  \\
  \ddot{\vect{X}}(t) &= 
  \begin{pmatrix}
    -\sin\theta & \cos\theta\\
    -\cos\theta & -\sin\theta
  \end{pmatrix}\cdot[\ddot{\theta}(\vect{x}_0 + \vect{v}_0 t) + 2\vect{v}_0]
  -
  \dot{\theta}^2\begin{pmatrix}
    \cos\theta & \sin\theta\\
    -\sin\theta & \cos\theta
  \end{pmatrix}\cdot(\vect{x}_0 + \vect{v}_0 t) 
\end{align*}
:::

Now we can simply differentiate the transformation to obtain the equations of motion.
*(Note: since we want the final expressions in terms of $\vect{X}$, it is fastest to
differentiate the inverse transform.  This is not necessary, but if you work with the
original expression, you will have one more step of substituting previous relationships
to replace $\vect{x}$ and $\dot{\vect{x}}$ with the corresponding expressions in terms
of $\vect{X}$)*:

Thus we obtain the following equations of motion in the rotating frame 

\begin{gather*}
  \vect{x} = \mat{R}^{-1}_{\theta}\cdot\vect{X}, \\
  \dot{\vect{x}} = \dot{\mat{R}}_{\theta}^{-1}\cdot\vect{X} 
                 + \mat{R}^{-1}_{\theta}\cdot\dot{\vect{X}},\\
  \ddot{\vect{x}} = \ddot{\mat{R}}_{\theta}^{-1}\cdot\vect{X} 
                  + 2\dot{\mat{R}}_{\theta}^{-1}\cdot\dot{\vect{X}}
                 + \mat{R}^{-1}_{\theta}\cdot\ddot{\vect{X}}.
\end{gather*}
 
Now, since the original particle is free, $\ddot{\vect{x}} = \vect{0}$, so we can multiply the last equation by $\mat{R}_{\theta}$ on the left:

\begin{gather*}
  \mat{R}_{\theta}\cdot\ddot{\vect{x}} = \mat{0} \qquad\Longrightarrow \qquad
     \ddot{\vect{X}} = -\mat{R}_{\theta}\ddot{\mat{R}}_{\theta}^{-1}\cdot\vect{X} 
     -2\mat{R}_{\theta}\dot{\mat{R}}_{\theta}^{-1}\cdot\dot{\vect{X}}.
\end{gather*}

Finally, we need to perform some matrix algebra to simplify the expressions
$\dot{\mat{R}}^{-1}$ and $\ddot{\mat{R}}_\theta^{-1}$.  These follow from
$\mat{R}_{\theta}\cdot\mat{R}^{-1}_{\theta} = \mat{1}$ and $\dot{\mat{1}} = \mat{0}$.
Since the rotation axis does not change direction, all matrices commute and we can
simply differentiate the exponential without any concern about commutation.

\begin{align*}
  \mat{R_\theta}^{-1} 
  &= \mat{R_{-\theta}},
  \\
  \dot{\mat{R}_\theta}^{-1} 
  &= e^{\mat{\vec{\theta}\times}}\cdot\mat{\vec{\dot{\theta}}\times}
  = e^{\mat{\vec{\theta}\times}}\cdot\mat{\vec{\omega}\times},
  \\
  \ddot{\mat{R}_\theta}^{-1}
  &= e^{\mat{\vec{\theta}\times}}\cdot\mat{\vec{\omega}\times} \cdot \mat{\vec{\omega}\times}
  +
  e^{\mat{\vec{\theta}\times}}\cdot\mat{\vec{\dot{\omega}}\times}
  \\
  \mat{R_\theta}\cdot\dot{\mat{R}_\theta}^{-1}
  &= \mat{\vec{\omega}\times}
  \\
  \mat{R_\theta}\cdot\ddot{\mat{R}_\theta}^{-1}
  &= \mat{\vec{\omega}\times} \cdot \mat{\vec{\omega}\times} + \mat{\vec{\dot{\omega}}\times}, 
  \qquad
  \mat{\vec{\omega}\times} = \begin{pmatrix}
    0 & -\omega\\
    \omega & 0
  \end{pmatrix}, \qquad
  \mat{\vec{\omega}\times} \cdot \mat{\vec{\omega}\times} = -\omega^2\mat{1}
\end{align*}

Collecting the results, we have:

\begin{align*}
   \ddot{\vect{X}} &= -\vec{\omega}\times\vec{\omega}\times\vect{X}
   - \dot{\vec{\omega}}\times\vect{X} - 2\vec{\omega}\times\dot{\vect{X}}
   = \omega^2\cdot\vect{X} - \dot{\vec{\omega}}\times\vect{X} - 2\vec{\omega}\times\dot{\vect{X}}
   \\
   &= \omega^2\vect{X} -
   \begin{pmatrix}
     0 & -1\\
     1 & 0
   \end{pmatrix}\cdot\left(
     \dot{\omega}\vect{X}
     +2\omega\dot{\vect{X}}\right),\\
   \begin{pmatrix}
    \ddot{x}\\
    \ddot{y}
  \end{pmatrix} &=
  \begin{pmatrix}
    \omega^2 x + \dot{\omega}y + 2\omega\dot{y}\\
    \omega^2 y - \dot{\omega}x + 2\omega\dot{x}
  \end{pmatrix} =
  \begin{pmatrix}
    \alpha^2 t^2 x + \alpha y + 2\alpha t\dot{y}\\
    \alpha^2 t^2 y - \alpha x + 2\alpha t\dot{x}
  \end{pmatrix}.
\end{align*}

*(Quick check: dimensions are correct and this is equivalent to (10.2) in the book and
to what we remembered above.)*

Now we can substitute the original solution into the right hand side which quickly
allows us to verify this as the solution.

:::{admonition} Explicit solution in rotating frame.
:class: dropdown

If you want to actually solve these equation, one can express them as complex numbers,

\begin{gather*}
  Z = X + \I Y,\\
  \ddot{Z} = \omega^2 Z - \I(\dot{\omega}Z + 2\omega\dot{Z})
           = (\omega^2 - \I\dot{\omega})Z - 2\I\omega\dot{Z}\\
  \ddot{Z} + q(t)Z + p(t)\dot{Z} = 0,\\
  p(t) = \I\dot{\omega}(t) - \omega^2(t), \qquad
  q(t) = 2\I\omega(t).
\end{gather*}

This is a homogeneous second-order ODE in standard form.  Note that one immediately sees
a trivial solution $Z=0$.  This corresponds to the particle sitting at rest at the
origin which is the same in both frames.  Solving such equations can be tricky, but
there exist nice theorems such as [Abel's
identity](https://en.wikipedia.org/wiki/Abel%27s_identity) that allow you to construct
the general solution if you can find a specific solution.  One could use our results
here to solve a general class of such problems where, for arbitrary $q(t) =
2\I\omega(t)$,

\begin{gather*}
  \ddot{Z} + q(t)Z + p(t)\dot{Z} = 0, \qquad
  p(t) = \frac{q^2(t)}{4} + \frac{\dot{q}(t)}{2},\\
  Z(t) = e^{\I\theta(t)}z(t) = e^{\I\theta(t)}(z_0 + v_0 t),\qquad
  \theta(t) = \frac{1}{2\I}\int_0^{t}\d{t}\;q(t),
\end{gather*}

These equations can also be expressed these as a set of coupled first-order equations
with time-dependent coefficents:

\begin{gather*}
  \vect{q}(t) = \begin{pmatrix}
    Z(t)\\
    \dot{Z}(t)
  \end{pmatrix}, \qquad
  \dot{\vect{q}}(t) = \mat{A}(t)\cdot\vect{q}(t), \\
  \mat{A}(t) = \begin{pmatrix}
    0 & 1\\
    \omega^2(t) - \I\dot{\omega}(t) & -2\I\omega(t)
  \end{pmatrix}
  = \begin{pmatrix}
    0 & 1\\
    \alpha^2t^2 - \I\alpha & -2\I\alpha t
  \end{pmatrix}.
\end{gather*}

The challenge now is that the matrix $\mat{A}(t)$ does not commute at different times.
These types of equations can be solved formally as:

\begin{gather*}
  \begin{pmatrix}
    Z(t)\\
    \dot{Z}(t)
  \end{pmatrix}
  = T\left[e^{\int_0^{t}\mat{A}(t)\d{t}}
  \right]
  \begin{pmatrix}
    Z(0)\\
    \dot{Z}(0)
  \end{pmatrix}
\end{gather*}

where $T[]$ is the time-ordering operator which places all terms in the expansion of the
exponential with later times appearing on the left.  If $\mat{A}(t)$ were to commute
with itself at all times, this would be the complete solution, but since it does not, we
face a similar problem to that of solving the time-dependent Schrödinger equation in
quantum mechanics.  Of course, such a solution can be easily integrated numerically.
:::

::::

## Larmor's Theorem: Problem 2.1

+++

Do problem 2.1 from {cite:p}`Fetter:2003`:

> **2.1** **Larmor's theorem**
>
> 1. The Lorentz force implies the equation of motion $m\ddot{\vect{r}} = e(\vect{E} +
>     c^{-1}\vect{v}\times\vect{B})$.  Prove that the effect of a weak uniform magnetic
>     field $\vect{B}$ on the motion of a charged particle in a central electric field
>     $\vect{E} = E(\norm{\vect{r}})\uvect{r}$ can be removed by transforming to a
>     coordinate system rotating with an angular frequency $\omega_L =
>     -(e/2mc)\vect{B}$. (State precisely what "weak" means.)
> 2. Extend this result to a system of particles of given ratio $e/m$ interacting
>     through potentials $V_{ij}(\norm{\vect{r}_i - \vect{r}_j})$.

In class you have seen how trajectories and dynamics in a rotating frame can look much
more complicated than in an inertial frame.  Here you will explore how the motion of a
charged particle in a weak electric field looks like motion in a rotating frame.  Hence,
by removing the "rotation" you can simplify the problem to permit an analysis that is
more akin to mechanics in a non-rotating frame.

Notes:
1. When justifying what makes something small, be sure to compare quantities of the same
   dimension.  *(E.g., it make no sense to say that the mass of an electron $m$ is
   small.  Compared to your mass, it is indeed small, but compared to the mass of an
   electron neutrino, it is huge!)*
2. The point of the second part is to explain why the charge-to-mass ratio must be the
   same for all particles, and why the potentials must only depend on their separation.
   How would the analysis break down if these conditions were not satisfied?

+++

### Solution

::::{solution}

:::{admonition} Comments

This is the reverse of the previous problem.  Here we can get a simpler problem by
transforming the system to a rotating frame where the coriolis term cancels part of our
original problem.  Solving in the rotating frame is mathematically simpler, and then we
can transform back to get the physical solution.

Note: Some people noted that a free particle moving in a constant magnetic field moves
in a circle or helix, and suggested that you could transform to a rotating frame so that
this motion is simple for any field, not just weak fields.  This fails for two reasons:

1. The transformation here requires choosing the rotation axis.  This depends on the
   position and velocity of the particle, so it does not offer a general solution.  In
   the problem discussed, the rotation axis goes through the center of the fixed
   electrostatic potential (nucleus in an atom for example) but applies for all
   electrons.
2. This transformation only works for free particles – which are easy to solve – so
   does not help simplify the problem when you really need help.
3. The angular frequency here is the cyclotron frequency which differs from the
   solution discussed here by a factor of 2.

As noted in the problem, one must make a dimensionless quantity small compared to one.
Note also that one cannot compare vectors with inequality signs: $\vect{E} < \vect{B}$
makes no sense.  One needs to take the norm.

Another strange thing people did was to write:

\begin{gather*}
  \vect{E} = E \norm{\vect{r}} \uvect{r}.
\end{gather*}

This makes no sense dimensionally, and leads to a missing dependence on $r$ in the
small-field condition.  Presumably it is a misinterpretation of the notation in the
question where the electric field is written **as a function of $r = \norm{\vect{r}}$**
but pointing in the $\uvect{r}$ direction:
\begin{gather*}
  \vect{E} = E(r) \uvect{r}.
\end{gather*}
$E(\norm{\vect{r}}) = E(r)$ means that the magnitude of the electric field **is a
function of the radius $r$**, not that the dimensions of $E$ suddenly changed to have an
extra factor of $1/r$.
:::

The key here is that the Lorentz force equation looks like a transformation to a
rotating frame.  Using the previous notation where $\vect{X}$ are coordinates in a frame
rotating with angular velocity $\vect{\omega}$ and $\vect{x}$ are the inertial
coordinates.  Compare (10.2) with the Lorentz force law:

\begin{gather*}
  \ddot{\vect{r}} = \frac{\vect{F}}{m} 
     = \frac{e}{m}\left(\vect{E} + \frac{\vect{v}\times\vect{B}}{c}\right)
     = \frac{e}{m}\vect{E} - \frac{e\vect{B}}{mc}\times\vect{v}\\
  \ddot{\vect{X}} = \frac{\vect{F}_{b}}{m} - 2\vec{\omega}\times \dot{\vect{X}} - (\vec{\omega}\,\times)^2 \vect{X}.
\end{gather*}

Thus, the effect of a magnetic field looks like the Coriolis force.  We may thus cancel
this term by boosting to an appropriately rotating frame where the Coriolis term will
cancel the effects of the magnetic field.  Effecting this transformation $\vect{R} =
e^{-t\vec{\omega}\times}\vect{r}$:

\begin{gather*}
  \vect{r} = e^{t\vec{\omega}\times}\vect{R} , \qquad
  \vect{v} = e^{t\vec{\omega}\times}(\dot{\vect{R}} + \vec{\omega}\times\vect{R}), \\
  \ddot{\vect{r}} = e^{t\vec{\omega}\times}\left(
    (\vec{\omega}\times)^2\vect{R} 
    + 2\vec{\omega}\times\dot{\vect{R}} 
    + \ddot{\vect{R}}\right)
    = \frac{e}{m}\vect{E} - \frac{e\vect{B}}{cm}\times e^{t\vec{\omega}\times}(\dot{\vect{R}} + \vec{\omega}\times\vect{R}).
\end{gather*}

Thus, we can cancel the $\dot{\vect{R}}$ terms if we set 

\begin{gather*}
  \vect{\omega} = -\frac{e\vect{B}}{2mc}
\end{gather*}

to obtain

\begin{gather*}
    \ddot{\vect{R}}
    = \frac{e}{m}e^{-t\vec{\omega}\times}\vect{E} 
    + \left(\frac{e\vect{B}}{2cm}\times\right)^2\vect{R}.
\end{gather*}

:::{note} 
Note the sign on this term: if you have the wrong sign, you might have mistakenly used
$\vect{v} = \dot{\vect{R}}$ instead of converting $\dot{\vect{R}}$ back to the
inertial frame.

Note also that there is no restriction here for the motion to be planar.  There is no
reason for $\vect{B}$ to be perpendicular to $\vect{R}$.  The potential is central, but
the magnetic field breaks rotational invariance.  In other words, once we move to the
rotating frame, we have a central potential (for weak $\vect{B}$) and the orbit can be
in any plane, including one that is not orthogonal to $\vect{B}$.
:::

Now, if $\vect{E} = -\vect{\nabla}\phi(r)$ is centrally symmetric, then
$e^{-t\vec{\omega}\times}\vect{E} = \vect{E}$ is constant.  Finally, we can negect the
last term if

\begin{gather*}
  \frac{e^2\abs{B}^2}{4m^2c^2} \ll \frac{\abs{\ddot{\vect{R}}}}{\abs{\vect{R}}}, 
  \qquad
  \text{or}
  \qquad
  \frac{e\abs{B}^2 \abs{\vect{R}}}{4mc^2 \abs{\vect{E}}} \ll 1.
\end{gather*}

The small quantity on the left is the square of the classical Larmor frequency.

:::{note} 
The radial extent of the system is relevant, not just the ratio of the electric and
magnetic fields.  This make sense, because the term we are dropping is analogous to the
centrifugal force which gets larger as you move away from the center.  In particular, if
$\vect{E} = \vect{0}$, then neglecting the second term implies the motion will be a
straight line, but in fact, it will be a circle.  Eventually, this condition will be
violated and we cannot neglect the second term.  However, if the electric field is
strong enough, the electron will be tightly bound, and the condition will be valid for
all times.
:::

Thus, for a central electric potential with a small magnetic field we have the
simplified equation:

\begin{gather*}
  \ddot{\vect{R}} \approx \frac{e}{m}\vect{E}.
\end{gather*}

*For some reason, several students invoked an argument about the motion of the center of
mass.  This is not relevant for the problem here.*

The key last part of the problem is to note that if the interactions in a system are
pairwise, then since $\abs{\vect{r}_i - \vect{r}_j} = \abs{\vect{R}_i - \vect{R}_j}$ is
invariant under rotations, the equations can be applied to more complex systems.
However, for the cancellation to work for all particles, we must have the same $\omega$
for each particle, which requires that the charge-to-mass ratio $e/m$ be constant.


Although this is not part of the problem, note that if instead, we choose a frame with
the [**cyclotron
frequency**](https://en.wikipedia.org/wiki/Cyclotron_resonance#Cyclotron_Resonance_frequency)
(which differs from the Larmor frequency by a factor of 2)

\begin{gather*}
  \vect{\omega} = -\frac{e\vect{B}}{mc}
\end{gather*}

we obtain

\begin{gather*}
  \ddot{\vect{R}}
  = \frac{e}{m}e^{-t\vec{\omega}\times}\vect{E} - \vec{\omega}\times\dot{\vect{R}}.
\end{gather*}

Unlike the previous equations, this is not so useful for general solutions, but if
$\vect{E}$, $\vect{B}$ and $\dot{\vect{R}}$ are all parallel, then we have a description
of motion in 1D $\ddot{R} = eE/m$ which describes a helical motion of a charged particle
in a cyclotron.

Thus, we have two cases where transforming to a rotating frame actually simplifies the
equations of motion.

::::

+++

## Motion on the Earth: Problem 2.2

Do problem 2.2 from {cite:p}`Fetter:2003`:

> **2.2** Assume that over the time interval of interest, the center of mass of the
>    earth moves with approximately constant velocity with respect to the fixed stars
>    and that $\vect{\omega}$, the angular velocity of the earth, is a
>    constant. Rederive the terrestrial equations of particle motion (11.8) and (11.6)
>    by writing Newton's law in a body-fixed frame with origin at the surface of the
>    earth (Fig. 11.2).
> \begin{gather*}
    \vect{g} = -(GM_eR_e^{-2} - \omega^2 R_e\sin^2\theta) \hat{\vect{r}}
             + \tfrac{1}{2}\sin 2\theta\; \omega^2 R_e \hat{\vect{\theta}} \tag{11.6}\\
  m\ddot{\vect{r}} = m \vect{g} - 2m \vect{\omega}\times \dot{\vect{r}} \tag{11.8}
\end{gather*}

Perhaps the most important daily experience we have with accelerating frames is the
motion due to the rotation of the eath.  Here you are asked to rederive two expressions
in the book: one for the local gravitational acceleration $\vect{g}$, and the other for
Newton's laws.  The $\vect{g}$ here is what you would measure if you were to drop an
object in the laboratory.  The corrections due to the motion of the earth would be
critical if you were, for example, using such a measurement to determine the
distribution of mass in the earth since that is only one component of the total
acceleration.

Notes:
* The point of this problem is not to do a bunch of algebra: it is to figure out what
    the difference between Newton's laws expressed relative to an origin at the center
    of the earth, and those expressed in a terrestrial frame.  Hint: When considering
    the origin at the center of the earth, $\vect{r}$ is very long (radius of the earth
    plus any deviations in the lab).  This means that the centrifugal force is quite
    strong.  However, when derived in a terrestrial frame with the origin attached to
    the surface of the earth, $\vect{r}$ is short (few meters), so the centrifugal force
    is much smaller.  Obviously the motion of a particle in a terrestrial lab does not
    depend on where you place your origin.  How do you reconcile this apparent paradox?

The corrections to the dynamical laws affect the motion of projectiles as you will use
below.

+++

### Solution

::::{solution}

:::{admonition} Comments

The key here is that when measuring displacements from the surface, one is missing the
distance from the center of the earth.  Hence the centrifugal term appears much smaller
than it should be.  The missing piece comes from the fact that the origin of the surface
reference frame is accelerating, hence we need to include this acceleration as an
additional fictitious force, which exactly corrects
for the missing centrifugal piece obtained in the center-of-the-earth frame.

Specifically, one should use (10.2) (see below).  In the earth-centered frame, the $\d^2
\vect{a}/\d{t}^2 = 0$ term vanish since the center of the earth is not accelerating, but
$\vect{r}$ is huge, including the radius of the earth.  In contrast, a frame on the
surface of the earth is accelerating because it is moving in a circle due to the earth's
rotation.  Thus  $\d^2 \vect{a}/\d{t}^2 \neq 0$, and exactly compensates for the fact
that now $\vect{r}$ (which we call $\vect{q}$ below) is now short and lacking the radius
of the earth.

:::

The question here is to rederive equations (11.6) and (11.8) by considering a frame
moving at the surface of the earth:

\begin{gather}
  \vect{g} = -\left(\frac{GM_e}{R_e^2} - \omega^2R_e\sin^2\theta\right)\uvect{r}
    + \frac{1}{2}\sin 2\theta\; \omega^2R_e\uvect{\theta},\tag{11.6}\\
  m\ddot{\vect{r}} = \vect{F}' - 2m\vect{\omega}\times \vect{v} + m\vect{g}.\tag{11.8}
\end{gather}

In the book, these were derived from (11.3) and (11.4) by considering a frame rotating
at the center of the earth where the coordinate $\vect{r}$ is the position of a particle
relative to the center of the earth:

\begin{gather}
  m\left(\diff[2]{\vect{r}}{t}\right)_e 
  = \vect{F}_{g} + \vect{F}' 
    - 2m\vect{\omega}\times\left(\diff{\vect{r}}{t}\right)_e 
    - m \vect{\omega}\times(\vect{\omega} \times \vect{r}), \tag{11.3}\\
    \vect{F}_{g} = - \frac{GM_e m \vect{r}}{r^3},
    \tag{11.4}
\end{gather}

which in turn were derived from (10.2):

:::{math}
:class: full-width

\begin{gather}
  m\left(\diff[2]{\vect{r}}{t}\right)_{\rlap{\mathrm{body}}} = 
  \vect{F}^{(e)} - m \left(\diff[2]{\vect{a}}{t}\right)_{\rlap{\mathrm{inertial}}}
  - 2m \vect{\omega}\times\left(\diff{\vect{r}}{t}\right)_{\rlap{\mathrm{body}}}
  - m\vect{\omega}\times(\vect{\omega} \times \vect{r})
  - m \diff{\vect{\omega}}{t}\times \vect{r},
  \tag{10.2}
\end{gather}

:::

*For the earth, $\dot{\omega} < 0$ due to the tidal effects of the moon, but this is
still very small (~1-2ms per century).*

with $\ddot{\vect{a}} = \vect{0}$ (this is the comment in the problem about the center of mass
of the earth moving with approximately constant velocity) and $\dot{\vect{\omega}}
\approx 0$.

The book is now asking us to derive these equations in a new set of coordinates
$\vect{q}$ which is the position of the particle relative to the surface of the earth.
Thus:

\begin{gather*}
  \vect{q} = \vect{r} - \vect{R} 
\end{gather*}

where $\vect{R}$ is the vector from the center of the earth to the surface of the
earth. Note that $\vect{q}$ is now much smaller (meters) than $\vect{r}$ (which includes
the radius of the earth).  Thus, the centrifugal piece is much smaller.  This is exactly
compensated for by the fact that now the center of the coordinate system is
*accelerating* in addition to the rotation with angular velocity $\vect{\omega} = \omega
\uvect{Z}$.  Note that the acceleration of the vector $\vect{R}$ is

\begin{gather*}
  \vect{a} = \vect{R} = e^{\vect{\omega}t\times}\vect{R}_0, \qquad
  \dot{\vect{R}} = \vect{\omega}\times e^{\vect{\omega}t\times}\vect{R}_0
                 = \vect{\omega}\times \vect{R}, \\
  \left(\diff[2]{\vect{a}}{t}\right)_{\rlap{\mathrm{inertial}}} = \ddot{\vect{R}} 
           = \vect{\omega}\times(\vect{\omega}\times\vect{R}),
\end{gather*}

which has the same form as the centrifugal term.  The equations of motion now follow from (10.2):

:::{math}
:class: full-width

\begin{align}
  m\left(\diff[2]{\vect{q}}{t}\right)_{\rlap{\mathrm{body}}} &= 
  \vect{F}_{g} + \vect{F}' - m \left(\diff[2]{\vect{a}}{t}\right)_{\rlap{\mathrm{inertial}}}
  - 2m \vect{\omega}\times\left(\diff{\vect{q}}{t}\right)_{\rlap{\mathrm{body}}}
  - m\vect{\omega}\times(\vect{\omega} \times \vect{q}),\\
&= 
  \vect{F}_{g} + \vect{F}' - m \vect{\omega}\times(\vect{\omega} \times \vect{R})
  - 2m \vect{\omega}\times\left(\diff{\vect{q}}{t}\right)_{\rlap{\mathrm{body}}}
  - m\vect{\omega}\times(\vect{\omega} \times \vect{q}),\\  
&= 
  \vect{F}_{g} + \vect{F}'
  - 2m \vect{\omega}\times\left(\diff{\vect{q}}{t}\right)_{\rlap{\mathrm{body}}}
  - m\vect{\omega}\times\bigl(\vect{\omega} \times (\vect{q} + \vect{R})\bigr),\\
&= 
  \vect{F}_{g} + \vect{F}'
  - 2m \vect{\omega}\times\left(\diff{\vect{q}}{t}\right)_{\rlap{\mathrm{body}}}
  - m\vect{\omega}\times(\vect{\omega} \times \vect{r}).\tag{11.3}  
\end{align}
:::

Thus, we have rederived (11.3) in the new frame.

:::{important}

The key point is that, whereas the last term had a purely centrifugal
origin in the center-of-mass frame, it comes from a combination of the frame
acceleration term and the much weaker centrifugal term in the surface frame.  

:::

From this point the derivation of the (11.6) and (11.8) is the same.

::::

+++

## Projectile Motion: Problem 2.5

> **2.5** A cannon is placed on the surface of the earth at colatitude (polar angle)
>    $\theta$ and pointed due east.
>
> 1. If the cannon barrel makes an angle $\alpha$ with the horizontal, show that the
>     lateral deflection of a projectile when it strikes the earth is
>     $(4V_0^3/g^2)\omega \cos\theta\;\sin^2\alpha \cos\alpha$, where $V_0$ is the
>     initial speed of the projectile and $\omega$ is the earth's angular-rotation
>     speed.  What is the direction of this deflection?
> 2. If $R$ is the range of the projectile for the case $\omega=0$, show that the change
>     in the range is given by $(2R^3/g)^{1/2} \omega \sin \theta \bigl[(\cot
>     \alpha)^{1/2} - \tfrac{1}{3}(\tan \alpha)^{3/2}\bigr]$. Neglect terms of order
>     $\omega^2$ throughout.

Here you can use the equations of motion you derived in problem 2.2 to analyze the
motion of a projectile.

I suggest (but don't require) that you try this problem two ways:

1. Use the equations of motion you derived above, expanding for small $\omega$ as
   suggested in the problem.
2. Treat the canon ball as an orbiting body and compare its elliptical motion in the
   central potential of the earth with the solid body rotation of the earth.
   (I.e. solve the problem in the inertial frame.)

If part 2. is too messy for you to do in a reasonable time, try problem 2.3 instead.  I
would like you to see how working in the accelerated frame is the same as working in the
inertial frame so you can choose whichever method is simpler when you come across a
real problem of this type.

> **2.3** An observer that rotates with the earth drops a particle from a height $h$
> above the earth's surface.  Analyze the motion from an inertial frame of reference and
> rederive the net eastward deflection of $\tfrac{1}{3}\omega g \sin \theta\;
> (2h/g)^{3/2}$, where $\theta$ is the observer's polar angle.

Notes:
* Don't use the provided answers to "guide" you.  Work through the problem first on your
    own, then use these to check your work.  If you do not get the correct answer, then
    you probably made a mistake.  There are several places where it is easy to forget an
    important piece of physics in these problems that will spoil your calculation.  To
    check your answer with the form of the answer give, it might be best to plot the two
    results or look at them numerically to avoid doing unnecessary algebra.
* If you can't figure out what is missing, perhaps try a numerical solution to make sure
    that the formula presented here are indeed correct.

### Solution

:::::{solution}

:::{admonition} Comments

When first working on this problem, it seems quite tricky, but after some inspection and
a few false starts it should be clear how to proceed.  Trying to correctly reproduce the
final answer will ensure that you realize all of the physical points.  In particular,
the solution requires a careful accounting for the corrections to leading order
$O(\omega)$.  In particular:

* The rotation causes both lateral and vertical deflections.  The vertical deflection
  can be expressed in terms of the modified effective gravitational constant
  $g_{\mathrm{eff}} = g(1-\delta)$ where $\delta \sim O(\omega)$ arises from the
  vertical component of the coriolis force.  This effect extends the time of flight $t_f
  = t_f^{0} + O(\omega)$.
* For part a), the lateral deflection is entirely due to the rotation, so the leading term
  $O(\omega)$ suffices and one can use so the leading-order expression for the time of flight.
* For part b), however, the leading effect is $O(\omega^0)$, so the correction includes
  both contributions from the coriolis force **and** from the extended time of flight.
* For some reason – probably to get something like the final result given in the book –
  some people expressed the extended range as $r_x + r_y$ or $r_x+r_z$.  Adding
  components like this makes **no sense**!  They would have to be added in quadrature
  $\sqrt{r_x^2+r_y^2}$: the correction to the range from this hypotenuse is
  $O(\omega^2)$.  The missing linear correction comes from the extended time of flight.

:::

Here we work in the surface frame with axis $\uvect{x}$ pointing east, $\uvect{y}$
pointing north and $\uvect{z}$ pointing up.  Due to the earth being a geoid, the
gravitational acceleration (11.6) points exactly in the $-\uvect{z}$ direction:
$\vect{g} = -g\uvect{z}$.  This completely accounts for the centrifugal piece, so it
only remains to determine the Coriolis acceleration $-2m\vect{\omega}\times\vect{v}$.

::::{note}

One might worry that "up" as defined by $\uvect{z}$ differs from the the $\uvect{r}$
which points radially from the center of the earth.  The "polar" angle $\theta$
represents the angle between $\vect{\omega}$ and $\uvect{r}$, so in principle, we must
include the effects of the geoid through an angle $\beta$ such that $\theta - \beta$ is
the angle between $\vect{\omega}$ and $\uvect{z}$.  Being careful, the angular velocity
can be expressed as:

\begin{gather*}
  \vect{\omega} = \omega\cos\theta\; \uvect{r} - \omega\sin\theta\; \uvect{\theta}
                = \omega\cos(\theta-\beta)\; \uvect{z} + \omega\sin(\theta-\beta)\;\uvect{y}.
\end{gather*}

However, for the given problem we can neglect these corrections as they are very small:

\begin{gather*}
  \tan\beta = \frac{g_\theta}{g_r} 
  = \frac{\sin 2\theta\; \omega^2 R_e/2}
         {g_0 - \omega ^2 R_e \sin^2\theta}
  = \epsilon\frac{\tfrac{1}{2}\sin 2\theta}
            {1 - \epsilon\sin^2\theta}\\
  \epsilon = \frac{\omega^2 R_e}{g_0} \approx 0.003.
\end{gather*}

Thus, we can safely ignore this effect and simply consider $\theta$ to be the angle
between $\vect{\omega}$ and $\uvect{z}$.

:::{important}
Notice here that I am computing a dimensionless quantity $\epsilon$ here: you must do
this to clearly say that something is "small" compared to $1$.
:::

::::

Thus, in the surface frame, the acceleration $\vect{a} = \vect{g} -
2\vect{\omega}\times\vect{v} + \order(\omega^2)$ and initial velocity $\vect{v}_0$ are:

\begin{gather*}
  \vect{\omega} =
  \begin{pmatrix}
    \omega_x\\
    \omega_y\\
    \omega_z\\
  \end{pmatrix}
  =
  \begin{pmatrix}
    0\\
    \omega \cos\theta\\
    \omega \sin\theta\\
  \end{pmatrix}, \qquad
  \vect{v_0} = 
  \begin{pmatrix}
    V_0\cos\alpha\\
    0 \\
    V_0\sin\alpha
  \end{pmatrix},\\
  \vect{a} = 
  \begin{pmatrix}
    2(v_y\omega_z-v_z\omega_y)\\
    2v_x\omega_z \\
    - (g - 2v_x\omega_y)
  \end{pmatrix}
  \approx
  \begin{pmatrix}
    2\omega(v_y\cos\theta - v_z\sin\theta)\\
    2v_x\omega\cos\theta \\
    - (g - 2\omega v_x\sin\theta)
  \end{pmatrix}.
\end{gather*}

Now if $\omega = 0$ we have the solution:

\begin{gather*}
  \vect{r} = \begin{pmatrix}
    tV_0\cos\alpha\\
    0 \\
    tV_0\sin\alpha - \frac{g t^2}{2}
  \end{pmatrix}, \qquad
  \vect{v} = \begin{pmatrix}
    V_0\cos\alpha\\
    0 \\
    V_0\sin\alpha - g t
  \end{pmatrix}, \\
  t_f = \frac{2V_0\sin\alpha}{g}, \qquad
  R = t_fV_0\cos\alpha = \frac{2V_0^2\sin\alpha\cos\alpha}{g}
  = \frac{V_0^2\sin(2\alpha)}{g},
\end{gather*}

where $t_f$ is the time to impact and $R$ is the range.

Physically, the rotation has three effects:

1. It reduces $g$, thereby increasing $t_f$ and $R_0$.
2. It adds a lateral deflection due to the acceleration $2v_x\omega_z \uvect{y}$.
3. It alters the velocity in the $\uvect{x}$ direction.

To fully account for these is quite complicated because they all affect each other: for
example, once we account for the Coriolis force changing $v_x$, we must go back and see
how this changes the other velocities... and so forth.  However, the base effect is
order $\order(\omega)$, so if we only want the answer to leading order, we can neglect
the higher order corrections, each step of which will bring in another factor of
$\omega$.  Thus, note that $v_y \propto \omega$ so $a_x \approx 2v_z\omega_y +
\order(\omega^2)$ for example.

Hence the only leading order effects are the deflection and the reduced $g$.  We start
with the deflection:

\begin{gather*}
  v_y = \int_0^t a_y\d{t} = 2v_x\omega_z t + \order(\omega^2), \qquad  
  y = \int_0^t v_y\d{t} = v_x\omega_z t^2 + \order(\omega^2).
\end{gather*}

The total lateral deflection is thus:

:::{math}
:class: full-width

\begin{gather*}
  y = v_x \omega_z t_f^2 = V_0\cos\alpha \; \omega\cos\theta \; \frac{4V_0^2\sin^2\alpha}{g^2}
    = \frac{4V_0^3}{g^2}\omega \cos\theta \sin^2\alpha\cos\alpha + \order(\omega^2), 
\end{gather*}

:::

as stated in the book.  Note that this deflection is **towards the equator**. (It is only
south in the northern hemisphere: in Australia, the deflection would be north.)

The vertical motion is only affected by the reduced $g_{\mathrm{eff}}$:

\begin{gather*}
  g_{\mathrm{eff}} = g - 2v_x\omega_y 
    =  g - 2\omega V_0 \cos\alpha\sin\theta + \order(\omega^2)
    = g(1-\delta) + \order(\omega^2),\\
  \delta = \omega\frac{2V_0\cos\alpha\sin\theta}{g}.
\end{gather*}

:::{margin}

In the previous calculation for $y$ to linear order in $\omega$, one was able to simply
use the leading order approximation $t_f \approx 2V_0\sin\alpha / g$ since the prefactor
of $\omega_z = \omega \cos\theta$ was already small.  Here, however, the next correction
cannot be neglected and we must include the next order correction by keeping $g_{\text{eff}}$.

:::

This alters the vertical motion so that the time to impact increases:

\begin{align*}
  t_f &= \frac{2V_0\sin\alpha}{g_{\mathrm{eff}}} \approx \frac{2V_0\sin\alpha}{g}(1+\delta) + \order(\omega^2), \\
  v_z &= V_0\sin\alpha -  g_{\mathrm{eff}}t = g_{\mathrm{eff}}\left(\frac{t_f}{2} - t\right),
\end{align*}

Now consider the motion in the $\uvect{x}$ direction:

\begin{align*}
  a_x &= -2\omega_yv_z = -2\omega_y g_{\mathrm{eff}}\left(\frac{t_f}{2} - t\right), \\
  v_x &= V_0\cos\alpha - \omega_y g_{\mathrm{eff}}(t_f t - t^2), \\
  x &=  V_0\cos\alpha t - \omega_y g_{\mathrm{eff}}\left(\frac{t_f t^2}{2}- \frac{t^3}{3}\right)\\
  &=  V_0\cos\alpha t - \omega_y g\left(\frac{t_f t^2}{2}- \frac{t^3}{3}\right) + \order(\omega^2).
\end{align*}

The final displacement in the $x$ directions to order $\order(\omega^2)$ is thus:

\begin{align*}
  R_x &= V_0\cos\alpha \; t_f - \omega_y g \frac{t_f^3}{6}\\
      &= R + \omega\frac{4V_0^3\sin\theta\sin\alpha\cos^2\alpha}{g^2}
          - \omega\frac{g\sin\theta}{6}\left(\frac{2V_0\sin\alpha}{g}\right)^3\\
      &= R + \frac{4V_0^3}{g^2}\omega\sin\theta \left(
        \cos^2\alpha\sin\alpha  - \frac{\sin^3\alpha}{3}\right).
\end{align*}

Thus, the change to the range is:

\begin{align*}
   \delta R &= 
   \frac{4V_0^3}{g^2}\omega \sin\theta\left(
     \cos^2\alpha\sin\alpha - \frac{\sin^3\alpha}{3}\right)\\
  &= \sqrt{2\frac{R^3}{g}}\omega \sin\theta
  \frac{\cos^2\alpha\sin\alpha - \frac{\sin^3\alpha}{3}}{\left(\frac{\sin(2\alpha)}{2}\right)^{3/2}}.
\end{align*}

:::{note}
Note: the range is also extended by the motion in the $y$ direction, but this will
appear at order $\order(\omega^2)$, so we drop it.  To see this, note that the final
location of impact is $\vect{R} = R_x\uvect{x} + y \uvect{y} = (R + x)\uvect{x} +
y \uvect{y}$, where $x\sim \order(\omega)$ and $y \sim \order(\omega)$ are small.  The
range is thus:

\begin{gather*}
  \norm{R} = \sqrt{(R + x)^2 + y^2} = \sqrt{R^2 + 2Rx + x^2 + y^2}\\
  = R\sqrt{1 + 2x/R + \order(\omega^2)} = R + x + \order(\omega^2).
\end{gather*}

Thus, even though the lateral deflection $y \sim \order(\omega)$, its contribution to
the range is $\order(\omega^2)$ and we can neglect it.

:::

A quick check shows that the last term can be rearranged to match the form in the book.
We do this graphically below, plotting the two equations without the prefactor:

:::::

```{code-cell} ipython3
%pylab inline --no-import-all
from IPython.display import clear_output; clear_output()

fig, axs = plt.subplots(1, 2, figsize=(6.5, 3))
alpha = np.linspace(0, 1, 100)[1:-1]
ax = axs[0]
ax.plot(alpha, (np.cos(alpha)**2*np.sin(alpha) - np.sin(alpha)**3/3), 
        label="Ours")
ax.plot(alpha, (np.sqrt(1./np.tan(alpha)) 
                - np.tan(alpha)**(3/2)/3)*(np.sin(alpha)*np.cos(alpha))**(3./2), '+',
        label="Book")

ax = axs[1]
ax.plot(alpha, (np.cos(alpha)**2*np.sin(alpha) 
                - np.sin(alpha)**3/3)/(np.sin(2*alpha)/2)**(3./2.), ':',
        label="Ours")
ax.plot(alpha, np.sqrt(1./np.tan(alpha)) - np.tan(alpha)**(3./2.)/3, '+',
        label="Book")
for ax in axs:
    ax.legend()
    ax.set(xlabel=r"$\alpha$");
```
