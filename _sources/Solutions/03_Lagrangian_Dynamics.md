---
jupytext:
  formats: md:myst,ipynb
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

# Assignment 3: Lagrangian Dynamics

## Brachistochrone Problem

Consider a bead of mass $m$ sliding down a wire from the point $\vect{P} = (x_0, y_0)$.

1. Write and expression for the total time $T$ it takes for the bead to slide from
   $\vect{P}$ to the origin as a functional of the path of the wire $y(x)$.
2. Write the Euler-Lagrange equations for the solution.
3. Show that the solution is a [cycloid](https://en.wikipedia.org/wiki/Cycloid) (or
   solve the problem directly).

+++

### Solution

::::{solution}
::::

Let the trajectory be $y(x)$ with $y(0) = 0$ and $y(x_0) = y_0$ held fixed.  We shall
express everything in terms of the coordinate $x$ so that the speed is

\begin{gather*}
  v = \sqrt{\dot{x}^2 + \dot{y}^2} = \dot{x}\sqrt{1+(y')^2}.
\end{gather*}

We shall also take $g$ positive but let $y$ increase downwards so that the final
solution is at a positive $y_0 = y(x_0) > 0$ (which is below $0$, but this removes minus
signs from the equations).  Note then that the potential energy is $-mgy = mgh$ relative
to the starting point.

Once we find $v$, how do we get a time from this?  Several approaches are valid, which derive from
\begin{gather*}
  t = \int_0^{t}\d{t}
\end{gather*}
by re-expressing $\d{t}$.

1. The particle must move from $x=0$ to $x=x_0$, and the component of the velocity in
   this direction is
   \begin{gather*}
     \dot{x} = \sqrt{v^2 - \dot{y}^2} = \frac{v}{\sqrt{1+(y')^2}}.
   \end{gather*}
   Thus, we can integrate:
   \begin{gather*}
      t = \int_0^{t}\d{t} = \int_{0}^{x_0}\frac{\d{t}}{\d{x}}\d{x} 
                          = \int_{0}^{x_0} \frac{\d{x}}{\dot{x}}
                          = \int_{0}^{x_0} \frac{\sqrt{1+(y')^2}}{v}\d{x}.
   \end{gather*}
2. If the Brachistochrone ends up being monotonic, then the particle also must descend
   from $y(0) = 0$ to $y(t) = y_0$, and the vertical component of the velocity is
   \begin{gather*}
     \dot{y} = \sqrt{v^2 - \dot{x}^2} = \frac{v}{\sqrt{1+(x')^2}},
   \end{gather*}
   where we now treat $x(y)$.  This gives:
   \begin{gather*}
      t = \int_0^{t}\d{t} = \int_{0}^{y_0}\frac{\d{t}}{\d{y}}\d{y} 
                          = \int_{0}^{y_0} \frac{\d{y}}{\dot{y}}
                          = \int_{0}^{y_0} \frac{\sqrt{1+(x')^2}}{v}\d{y}.
   \end{gather*}
   
3. We can directly use $v$ and note that this is the change in the arc-length $v =
   \d{l}/\d{t}$:
   \begin{gather*}
      t = \int_0^{t}\d{t} = \int_{0}^{L}\frac{\d{t}}{\d{l}}\d{l} 
                          = \int_{0}^{L} \frac{\d{l}}{v}
                          = \int_{0}^{L} \frac{\sqrt{\d{x}^2 + \d{y}^2}}{v},
   \end{gather*}
   which leads to either form upon changing variables.

How to find $v$?  This simply comes from conservation of energy (recall that $y$
increases in the downward direction):

\begin{gather*}
  E = \frac{1}{2}mv^2 - mgy = 0, \qquad
  v = \sqrt{2gy}.
\end{gather*}

:::{admonition} Lagrangian Approach

Since we are practicing Lagrangian dynamics, it might be instructive to derive this
rather trivial result.  The Lagrangian and equations of motion are:

\begin{gather*}
  \mathcal{L}(x, \dot{x}) = \frac{1}{2}mv^2 - mgh
    = m\frac{1+(y')^2}{2}\dot{x}^2 + mgy, \\
  p = m(1+(y')^2)\dot{x}, \qquad
  \dot{p} = my'y''\dot{x}^2 + mgy'.
\end{gather*}

The Euler-Lagrange equations are thus:

\begin{gather*}
  \ddot{x} = \frac{y'y''\dot{x}^2 + gy'}{1+(y')^2}
\end{gather*}

These look rather nasty, but since the Lagrangian does not depend on $t$ explicitly, the
Hamiltonian will be conserved:

\begin{gather*}
  H = p\dot{x} - \mathcal{L} = m\frac{1+(y')^2}{2}\dot{x}^2 - mgy = E = 0.
\end{gather*}
:::

Thus, we have
\begin{gather*}
  \dot{x} = \frac{v}{\sqrt{1+(y')^2}} = \sqrt{2g}\sqrt{\frac{y}{1+(y')^2}},\\
  t = \frac{1}{\sqrt{2g}}\int_{x_0}^{0}\sqrt{\frac{1+(y')^2}{y}}\d{x} 
    = \frac{1}{\sqrt{2g}} \int_{0}^{x_0}I(y, y')\d{x} 
\end{gather*}
which we must extremize.

We treat this with calculus of variations, to obtain the following Euler-Lagrange
equations (using my index notation from class for specifying partial derivatives):

\begin{gather*}
  I_{,y'}' = I_{,y},\\
  \diff{}{x}\frac{y'}{\sqrt{y(1+(y')^2)}} = \frac{1}{2}\frac{\sqrt{1+(y')^2}}{y^{3/2}}, \\
  2yy'' = 1 + (y')^2.
\end{gather*}

Here we check our algebra using SymPy:

```{code-cell} ipython3
import sympy
from sympy import var, sqrt, simplify
var('x', )
sympy.init_printing()
x = var('x')
y, dy, ddy = var("y, y', y''", cls=sympy.Function)
I = sqrt((1+dy(x)**2)/y(x))

# These make things look nicer:
substitutions = [
    (y(x).diff(x), dy(x)),
    (dy(x).diff(x), ddy(x))]
euler_lagrange_eq = (I.diff(dy(x)).diff(x) - I.diff(y(x))).subs(substitutions)
simplify(euler_lagrange_eq)
```

We see our equation in the numerator.  This equation is easily solved by the usual trick
$y'' = y'\d{y'}/\d{y}$:

\begin{gather*}
  2yy'\diff{y'}{y} = 1 + (y')^2\\
  \int \frac{2y'}{1+y'^2}\d{y'} = \int \frac{1}{y}\d{y}\\
  \ln[1+(y')^2] = C - \ln y\\
  1+(y')^2 = \frac{2R}{y},\\
  (y')^2 = \frac{2R}{y} - 1
\end{gather*}

where $2R = e^{C}$ is the integration constant.  This is exactly the [equation for a
cycloid](https://en.wikipedia.org/wiki/Cycloid#Equations) shifted by $y_0$ and
reflected.

#### Trick 1

Noticing that the integrand for $t$ is independent of $x$, we can make the analogy with
Lagrangian mechanics, and identify that the equivalent of the Hamiltonian is conserved
(sometimes called the [Beltrami
identity](https://en.wikipedia.org/wiki/Beltrami_identity):

\begin{gather*}
  L[y, y', x] \equiv I(y,y') = \sqrt{\frac{1+(y')^2}{y}}\\
  p_y \equiv \pdiff{Y(y,y')}{y'} = y'\sqrt{\frac{1}{y[1+(y')^2]}}\\
  H = p_y y' - L \equiv p_yy' - Y 
    = (y')^2\sqrt{\frac{1}{y[1+(y')^2]}} - \sqrt{\frac{1+(y')^2}{y}} = \tilde{C},\\
 \tilde{C} = \frac{-1}{\sqrt{y[1+(y')^2]}}.
\end{gather*}

This is equivalent to the above.

#### Trick 2

There is another interesting trick here.  Instead of parametrizing the path with $y(x)$,
one can instead use $x(y)$ with $\d{x} = x'(y)\d{y}$.  This may seem problematic since
$x(y)$ is likely not going to be a single-valued function, but this is easily
incorporated into the action formulation by simply specifying that the integral should
be taken over the entire curve, meaning that there will be multiple contributions from
the multi-valued regions.  The final equations are also local, so this is not an issue:

\begin{gather*}
  \sqrt{2g}t = \int_{0}^{y_0} \sqrt{\frac{(x')^2 + 1}{y}}\d{y} = \int_{0}^{y_0} J(x, x')\d{y}
\end{gather*}

Now the integrand depends only on $x'(y)$ and explicitly on $y$, but not on $x(y)$ so we
have the conservation law:

\begin{gather*}
  \pdiff{J}{x'} = \frac{x'}{\sqrt{y((x')^2 + 1)}} = \text{const.}
\end{gather*}

Squaring and isolating $x'$, we have

\begin{gather*}
  1 = C(y)\left(\frac{1}{(x')^2} + 1\right)
\end{gather*}

which immediately gives us the solution once we note that $y'(x) = 1/x'(y)$ without the
need to integrate.

+++

## Constrained Optimization

Find the maximum of the function $f(x,y) = -x^2-y^2$ subject to the constraint that
$(x-x_0)^2 + y^2 = R^2$ using the method of Lagrange multipliers.

```{code-cell} ipython3
%pylab inline --no-import-all

x0 = 5
R = 1
x = np.linspace(-10,10,100)[:, None]
y = np.linspace(-10,10,100)[None, :]
f = -x**2 - y**2
c = (x-x0)**2 + y**2

plt.figure(figsize=(10,10))
cs = plt.contour(x.ravel(), y.ravel(), f.T, 
                 20,
                 linestyles='-', colors='k')
plt.clabel(cs, fmt='%1.0f')
cs = plt.contour(x.ravel(), y.ravel(), c.T, 
                 levels=[0,1,2,3,4,5], linestyles='-', colors='b')
plt.clabel(cs, fmt='%1.0f')
plt.gca().set_aspect(1)
```

### Solution

::::{solution}

:::{admonition} Comments

The main challenge in this problem is carefully taking care of all of the conditions.
Some hangups:

* $-2y=2y\lambda$ implies that $y=0$ **or** $\lambda = -1$, not both.
* Don't forget multiple solutions: If $y\neq 0$, then $\lambda = \pm x_0/R - 1$.
* Check all solutions to determine which ones area really maxima.  Solving for the roots
  of the gradient only guarantees extrema -- they could be maxima, minima, or saddle
  points.  Either check the Hessian, or use physical/mathematical reasoning to ensure
  you have the correct solution $y=0$, $x=x_0 - R\sgn(x_0)$. 
:::

We extremize:

\begin{gather*}
  f(x,y) - \lambda \Bigl(c(x, y) - R^2\Bigr)\\
  (-x^2 - y^2) - \lambda \Bigl((x-x_0)^2 + y^2 - R^2\Bigr).
\end{gather*}

The partials are:

\begin{gather*}
  -2x = 2\lambda (x-x_0), \qquad
  -2y = 2\lambda y.
\end{gather*}

Hence, $y = 0$, and $x = \lambda x_0/(1 + \lambda)$.  The constraint equation requires
$(x-x_0)^2 + y^2 = R^2$ so

\begin{gather*}
  \Bigl(\frac{\lambda x_0 - (1+\lambda)x_0}{1+\lambda}\Bigr)^2 = R^2,\qquad
  \lambda = \pm\frac{x_0}{R} - 1.
\end{gather*}

Hence our final solution is:

\begin{gather*}
  x = \frac{x_0}{\lambda^{-1} + 1} = x_0 \pm R, \qquad y = 0.
\end{gather*}

This is easily seen by looking at the contours – the solution must lie on the $x$ axis,
and satisfy the constraint $(x-x_0)^2 = R^2$ so $x = x_0 \pm R$.

Which is the maximum?  This needs to be checked.  The value of the function at these
points is:

\begin{gather*}
  f = -(x_0 \pm R)^2.
\end{gather*}

Clearly this is maximized if we make $\abs{x_0 \pm R}$ as small as possible, so the
solution for the maximum is:

\begin{gather*}
  x = x_0 - R\sgn(x_0), \qquad y = 0, \qquad f_{\max} = -(\abs{x_0} - \abs{R})^2.
\end{gather*}

::::

+++

## Bead on a Circular Wire

Derive the equations of motion for a bead of mass $m$ moving on a frictionless circular
wire of radius $a$ (problem 3.1 below but for now consider the wire to be stationary).  Do
this three ways:

1. Use Newtonian mechanics by projecting the forces as appropriate.
2. Using the Euler-Lagrange equations for the coordinate $\theta(t)$ shown in class or
   in Fig. 16.1 in the book.
3. Using the Euler-Lagrange equations for the coordinates $r(t)$ and $\theta(t)$ but
   introducing the constraint $r(t) = a$ with the method of Lagrange multipliers.  Show
   how to use this formulation to find an expression for force of constraint exerted by
   the wire to keep the bead in place.  (I.e. the force you would feel if you were to
   hold the wire in place while the bead is moving.)
4. Derive the equations of motion for the bead if the wire moves up and down with a
   time-dependent function $h(t)$.  (This is a driven pendulum.  We will analyze it
   later.)  Could you have guessed this form?  (The answer should make sense physically.)

+++

### Solution

::::{solution}

<!-- $10=2+3+3+2$ -->

:::{admonition} Comments

This is an example of a problem where using the Lagrangian is not obviously simpler than
simply using Newton's laws.  The Lagrangian formulation does allow you to forgo breaking
the vectors into components, but there are intermediate steps in the derivation that
make the problem look more complicated.  These eventually cancel out giving a simple
final set of equations that can be more directly obtained with good physical intuition.

Don't always jump straight to formal Lagrangian/Hamiltonian techniques... sometimes good-old Newton
will get you to the correct answer faster.

:::

1. Choosing coordinates $(x, z) = (r\sin\theta, -r\cos\theta)$, we have the following:

   \begin{align*}
     \vect{v} &= r\dot{\theta}\underbrace{(\cos\theta, \sin\theta)}_{\uvect{\theta}},\\
     \vect{a} &= -r\dot{\theta}^2\underbrace{(\sin\theta, -\cos\theta)}_{\uvect{r}}
                 +r\ddot{\theta}\underbrace{(\cos\theta, \sin\theta)}_{\uvect{\theta}},\\
     \vect{F}_g &= (0, -mg) = mg\cos\theta\; \uvect{r} - mg \sin\theta \;\uvect{\theta},\\
     \vect{F}_N &= -F_N\uvect{r}.
   \end{align*}

   Newton's second law in polar coordinates thus states:
   
   \begin{gather*}
     \overbrace{-mr\dot{\theta}^2 = mg \cos\theta - F_N}^{\uvect{r}}, \qquad
     \overbrace{mr\ddot{\theta} = -mg \sin\theta}^{\uvect{\theta}}.
   \end{gather*}
   
   The first is simply the fact that, if the particle is moving with angular velocity
   $\omega = \dot{\theta}$, then it is experiencing a centripetal acceleration of
   $r\omega^2$.  The second is the additional tangential acceleration that increases or
   decreases the angular velocity.
   
   Thus, we have the solution:
   
   \begin{gather*}
     F_N = mr\dot{\theta}^2 + mg \cos\theta, \qquad
     \ddot{\theta} = -\frac{g}{r} \sin\theta.
   \end{gather*}
   
2. Since the bead moves in a circle, the speed is $v = a\dot{\theta}$ and the potential
   is $mgh = -mga\cos\theta$ so the Lagrangian is:

   \begin{gather*}
     L = \frac{m}{2} a^2\dot{\theta}^2 + mga\cos\theta,\qquad
     p = \pdiff{L}{\dot{\theta}} = m a^2 \dot{\theta},\\
     \dot{p} = \pdiff{L}{\theta} = -mga\sin\theta, \qquad
     \ddot{\theta} = \frac{-g}{a}\sin\theta
   \end{gather*}
   
   This is the second equation of motion obtained in part 1.
3. Now we include $m\dot{r}^2/2$ in the kinetic energy and replace $a$ with $r$.  The
   constraint is $r = a$ so we add:

   \begin{gather*}
     L = \frac{m}{2} (r^2\dot{\theta}^2 + \dot{r}^2)+ mgr\cos\theta - \lambda(r-a),\\
     p_\theta = m r^2 \dot{\theta},\qquad
     p_r =  m \dot{r},\\
     \dot{p}_\theta = 2mr\dot{r} \dot{\theta} + mr^2\ddot{\theta} = -mgr\sin\theta, \\
     \dot{p}_r = m\ddot{r} = mr\dot{\theta}^2 + mg\cos\theta - \lambda, \\
      2r\dot{r} \dot{\theta} + r^2\ddot{\theta} = -gr\sin\theta, \qquad
     \ddot{r} = r\dot{\theta}^2 + g\cos\theta - \frac{\lambda}{m}.
   \end{gather*}

   To enforce the constraint, we must have $\ddot{r} = \ddot{a} = 0$ so we need
   
   \begin{gather*}
     \lambda = mr\dot{\theta}^2 + gm\cos\theta
   \end{gather*}
   
   which is the constraint force $F_N$ found in part 1.  This is the tension in the
   pendulum or central force exerted by the wire.
   
4. Simply replace $g$ with $g + \ddot{h}$ since we are now in an accelerating frame.
   Formally, we have $x=a\sin\theta$, $\dot{x} = a\dot{\theta}\cos\theta$ , $y=h(t) -
   a\cos\theta$, $\dot{y} = \dot{h}+a\dot{\theta}\sin\theta$:

  \begin{gather*}
    L = \frac{m}{2}(\dot{x}^2 + \dot{y}^2) + mg\bigl(a\cos\theta + h(t)\bigr)\\
      = \frac{m}{2}\bigl(a^2\dot{\theta}^2 + 2a\dot{h}\dot{\theta}\sin\theta + \dot{h}^2\bigr)
      + mg\bigl(a\cos\theta + h(t)\bigr),\\
    p_\theta = m(a^2\dot{\theta} + a\dot{h}\sin\theta), \\
    \dot{p}_\theta = m(a^2\ddot{\theta} + a\ddot{h}\sin\theta + a\dot{h}\dot{\theta}\cos\theta)
    = ma\dot{h}\dot{\theta}\cos\theta  - mga\sin\theta,\\
    \ddot{\theta} = - \frac{g+\ddot{h}}{a}\sin\theta.
  \end{gather*}
  
  :::{note}

  Some people worked this in Cartesian coordinates with a constraint:
  
  \begin{gather*}
    L = \frac{m}{2}(\dot{x}^2 + \dot{y}^2) - mgy - \frac{m\lambda}{2}(x^2 + (y-h)^2 - a^2)\\
    p_{x} = m\dot{x}, \qquad p_y = m\dot{y}, \\
    \ddot{x} = -\lambda x , \qquad
    \ddot{y} = -g - \lambda(y-h).
  \end{gather*}
  
  This is correct, but now you need to use the constraint to find the Lagrange
  multiplier $\lambda$.  This is not straightforward:
  
  \begin{gather*}
    a\dot{a} = x\dot{x} + (\dot{y}-\dot{h})(y-h) = 0, \\
    \dot{x}^2 + x\ddot{x} + (\dot{y}-\dot{h})^2 + (\ddot{y}-\ddot{h})(y-h) = 0,\\
    \dot{x}^2 + (\dot{y}-\dot{h})^2 - \lambda \overbrace{\bigl(x^2 + (y-h)^2\bigr)}^{a^2} - (g+\ddot{h})(y-h) = 0,\\
    \lambda = \frac{\dot{x}^2 + (\dot{y}-\dot{h})^2 - (g+\ddot{h})(y-h)}{a^2}.
  \end{gather*}
  
  Transforming to the accelerating frame by letting $q=y-h$, makes manifest the effect
  that $g \rightarrow g + \ddot{h}$:
  
  \begin{gather*}
    \ddot{x} = -\lambda x , \qquad
    \ddot{q} = -(g + \ddot{h}) - \lambda q,\\
    x\dot{x} + \dot{q}q = 0, \qquad
    \lambda = \frac{\dot{x}^2 + \dot{q}^2 - (g+\ddot{h})q}{a^2}.
  \end{gather*}
   
  One can now use these to eliminate $x$ or $q$, or to transform to polar coordinates,
  but I am not sure this is a simplification over the previous procedures.
  :::
  
  
::::

:::{margin}
This is a canonical example of a case where the Hamiltonian is not the total energy.
Think carefully about the meaning of this in part **f)**.
:::
## Problem 3.1: Bead on a Rotating Circular Wire

> **3.1** A point mass $m$ slides without friction along a wire bent into a vertical
> circle of radius $a$. The wire rotates with constant angular velocity $\Omega$ about
> the vertical diameter, and the apparatus is placed in a uniform gravitational field
> $\vect{g}$ parallel to the axis of rotation.
> * **a)** Construct the lagrangian for the point mass using the angle $\theta$ (angular
>   displacement measured from the downward vertical) as generalized coordinate. 
> * **b)** Show that the condition for an equilibrium circular orbit at angle $\theta_0$
>   is $\cos\theta_0 = g/a\Omega^2$. Explain this result by balancing forces in a
>   co-rotating coordinate system. 
> * **c)** Demonstrate that this orbit is stable against small displacements along the
>   wire and show that the angular frequency of oscillation about the equilibrium orbit is
>   given by $\omega^2 = \Omega^2\sin^2\theta_0$. *Hint:* Write $\theta = \theta_0
>   +\eta(t)$, where $\eta(t)$ is a small quantity. 
> * **d)** What happens if $a\Omega^2 < g$?

:::{margin}
Hint: the total energy is not conserved.  Why?  Sometimes I like to call this a
"blender" potential, since the wire rotating at a constant angular velocity acts like a
blender.  Think about sticking something into a blender.  Is the energy of the object
conserved?  What is conserved instead?
:::
Additionally:

* **e)** Identify the conserved quantity/quantities and explain the underlying symmetries of the
  action that gives rise to it/them.
  
* **f)** Explain the physical meaning of the conserved quantity associated with time translation
  invariance.

### Solution

::::{solution}

<!-- $10=3+2+3+2$ -->

To the previous equations we must add the tangential component of the velocity $v_\perp
= \Omega r$ where $r = a\sin\theta$:

\begin{align*}
   L &= \frac{ma^2}{2} (\dot{\theta}^2 + \Omega^2\sin^2\theta) + mga\cos\theta,\\
   p &= \pdiff{L}{\dot{\theta}} = m a^2 \dot{\theta},\\
   \dot{p} &= \pdiff{L}{\theta} = ma^2\Omega^2\sin\theta\cos\theta - mga\sin\theta, \\
   \ddot{\theta} &= \left(\Omega^2\cos\theta - \frac{g}{a}\right)\sin\theta.
\end{align*}

The equilibrium circular orbit is where there is no force on $\theta$ so $\ddot{\theta} =
0$ which gives:

\begin{gather*}
  \cos\theta_0 = \frac{g}{a\Omega^2}.
\end{gather*}

This balances the centrifugal force $ma\omega^2$ with the horizontal force due to the
wire $-mgy'(x) = -mga\cos\theta_0$.

To check if this is stable, we write $\theta = \theta_0 + \eta(t)$ and expand with
$\cos(\theta_0 + \eta) = \cos\theta_0\cos\eta - \sin\theta_0\sin\eta \approx
\cos\theta_0 - \eta\sin\theta_0$:

\begin{gather*}
   \ddot{\eta} = -\eta\Omega^2\sin^2\theta_0
\end{gather*}

which has the solution $\eta = \cos(\phi_0 + t\Omega\sin\theta_0)$ which is oscillitory
about $\theta = \theta_0$, hence it is a stable orbit.

If $a\Omega^2 < g$ then the only stable orbit is $\theta = 0$ where $\sin\theta$
vanishes.  Expanding the equations for small $\theta$ we have:

\begin{gather*}
   \ddot{\theta} = \left(\Omega^2 - \frac{g}{a}\right)\theta + \order(\theta^2),
\end{gather*}

so the particle oscillates about $\theta=0$ with frequency $\omega = \sqrt{g/a -
\Omega^2}$.  Note that this becomes imaginary when $\Omega^2>g/a$ at which point this
stable equilibrium point becomes unstable to oscillations about the aforementioned
circular orbit.

:::{admonition} Comment

Some people claim that if $a\Omega^2 < g$, then no equilibrium exists.  This is careless
wording.  As we show above, $\theta=0$ is a perfectly fine equilibrium condition as one
can expect from simple physical considerations.  It just cannot be considered simply as
a circular orbit.
:::

Since the Lagrangian does not depend on time, the Hamiltonian is conserved:
\begin{gather*}
  H = \dot{\theta}p - L = 
  \frac{ma^2}{2} (\dot{\theta}^2 - \Omega^2\sin^2\theta) - mga\cos\theta
  = E - ma^2\Omega^2\sin^2\theta.
\end{gather*}
Not, however, that this is not the total energy of the system.  Instead, the Hamiltonian
includes the work done on the bead by the wire.

To understand this, think about conservation of angular momentum.  If the wire were
*not* constrained to rotate at a constant angular velocity, but were free to rotate,
then the system would have conserved energy and conserved angular momentum.  This would
cause the system to rotate faster when the bead fell closer to the core, similar to how
a figures skater rotates faster when they pull in their arms.

To compensate for this and keep rotating at a constant $\Omega$, the wire must apply a
force in the direction (or against) the motion.  What is missing from the analysis is
the Coriolis force, which points into the page in the $\uvect{\phi} =
\uvect{r}\times\uvect{\theta}$ direction:
\begin{gather*}
  \vect{F}_{\text{coriolis}} = -2m\vect{\Omega}\times\vect{v}, \qquad
  \vect{v} = a\dot{\theta}\uvect{\theta}.
\end{gather*}
Note: $\vect{v}$ here is the velocity in the rotating frame.

We can compute the cross product noting that $\theta$ is the angle between $\vect{\Omega}$ and
$\vect{r}$ and that $\pi - \theta$ is the angle between $\vect{\Omega}$ and
$\uvect{\theta}$:
\begin{gather*}
  \vect{\Omega}\times\uvect{\theta} = \Omega\cos\theta\; \uvect{\phi}, \qquad
  \vect{\Omega}\times\uvect{r} = \Omega\sin\theta\; \uvect{\phi},\\
  \vect{F}_{\text{coriolis}} = -2ma\Omega\dot{\theta}\cos\theta\; \uvect{\phi}.
\end{gather*}

This points perpendicular to the rotating plane and adds
energy at a rate (power):
\begin{gather*}
  \vect{v}_{\perp} = \vect{\Omega}\times\vect{r} = 
  a\Omega \sin\theta \; \uvect{\phi},\\
  P = \vect{F}_{\text{coriolis}}\cdot\vect{v}_{\perp}
  = -2ma^2\Omega^2\dot{\theta}\cos\theta\sin\theta
  = \diff{}{t}\left(-ma^2\Omega^2\sin^2\theta\right).
\end{gather*}
This is the derivative the extra contribution to $H$.  Thus, the conserved Hamiltonian
$H$ includes this contribution, taking into account the work done by the rotating wire
which must be supplied by the motor in order to keep the wire moving at a constant rate
$\Omega$.

:::{comment}
One might be tempted to include the centrifugal force in the calculation, however, this
is already taken into account in the regular dynamics in the rotating frame and is
responsible (along with gravity) for converting potential energy into kinetic energy in
this frame.
:::
::::
