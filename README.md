# FortranRbSerialChain

Serial RB chain solver using articulated method with arrays.

## Problem Setup & Equations

For each body the following equations are needed for the problem to be defined. 

 - The kinematics relating the motion of the previous body to the motion 
of the current body.

    $$ \begin{aligned}\boldsymbol{v}_{i} & =\boldsymbol{v}_{i-1}+\boldsymbol{s}_{i}qp_{i}\\
\boldsymbol{a}_{i} & =\boldsymbol{a}_{i-1}+\boldsymbol{s}_{i}\ddot{q}_{i}+\boldsymbol{v}_{i}\times\boldsymbol{s}_{i}qp_{i}\\
\end{aligned}$$

 - The force balance for each body, considering the reaction force coming from the next body.
 
     $$ \boldsymbol{f}_{i}={\bf I}_{i}\boldsymbol{a}_{i}+\boldsymbol{v}_{i}\times{\bf I}_{i}v_{i}+\boldsymbol{f}_{i+1}$$
	 
 - The joint power from the actuation torque
 
     $$ Q_i = \boldsymbol{s}_i^\intercal \boldsymbol{f}_i $$
	 
## Problem solution

Because the equations for each body depend on the previoius body motion _and_ the next body forces at the same time,
the solution invloves a recursive procedure developed by Roy Featherstone (https://homes.cs.washington.edu/~todorov/courses/amath533/FeatherstoneOrin00.pdf)	 

There are three stages to the solution.

### 1. Body Kinematics

Recursion from base to tip

$$\begin{aligned}\vec{r}_{i} & =\vec{r}_{i-1}+R_{i-1}\vec{\ell}\\
\vec{R}_{i} & =\vec{R}_{i-1}\mathrm{rot}(\hat{k},q_{i})\\
\vec{cg}_{i} & =\vec{r}_{i}+R_{i}\vec{cg}_{{\rm body}}\\
{\bf I}_{i} & =\mathrm{spi}(m,I_{{\rm body}},R_{i},\vec{cg}_{i})\\
\boldsymbol{s}_{i} & =\mathrm{twist}(\hat{k},\vec{r}_{i},0)\\
\boldsymbol{v}_{i} & =\boldsymbol{v}_{i-1}+\boldsymbol{s}_{i}qp_{i}\\
\boldsymbol{\kappa}_{i} & =\boldsymbol{v}_{i}\times\boldsymbol{s}_{i}qp_{i}\\
\boldsymbol{w}_{i} & =\mathrm{wrench}(m\,\vec{g},\vec{cg},0)\\
\boldsymbol{p}_{i} & =\boldsymbol{v}_{i}\times{\bf I}_{i}\boldsymbol{v}_{i}-\boldsymbol{w}_{i}
\end{aligned}$$

### 2. Articulated Inertia

Recursion from tip to base

$$\begin{aligned}\boldsymbol{d}_{i} & =\boldsymbol{p}_{i}+\boldsymbol{T}_{i+1}Q_{i+1}+\boldsymbol{\Phi}_{i+1}\left({\bf A}_{i+1}\boldsymbol{\kappa}_{i+1}+\boldsymbol{d}_{i+1}\right)\\
{\bf A}_{i} & ={\bf I}_{i}+\boldsymbol{\Phi}_{i+1}{\bf A}_{i+1}\\
\boldsymbol{T}_{i} & ={\bf A}_{i}\boldsymbol{s}_{i}\left(\boldsymbol{s}_{i}^{\intercal}{\bf A}_{i}\boldsymbol{s}_{i}\right)^{-1}\\
\boldsymbol{\Phi}_{i} & =1-\boldsymbol{T}_{i}\boldsymbol{s}_{i}^{\intercal}
\end{aligned}$$

### 3. Joint Dynamics 

Recursion from base to tip

$$\begin{aligned}\ddot{q}_{i} & =\left(\boldsymbol{s}_{i}^{\intercal}{\bf A}_{i}\boldsymbol{s}_{i}\right)^{-1}\left(Q_{i}-\boldsymbol{s}_{i}^{\intercal}\left({\bf A}_{i}\left(\boldsymbol{a}_{i-1}+\boldsymbol{\kappa}_{i}\right)+\boldsymbol{d}_{i}\right)\right)\\
\boldsymbol{a}_{i} & =\boldsymbol{a}_{i-1}+\boldsymbol{s}_{i}\ddot{q}_{i}+\boldsymbol{\kappa}_{i}\\
\boldsymbol{f}_{i} & ={\bf A}_{i}\boldsymbol{a}_{i}+\boldsymbol{d}_{i}
\end{aligned}$$

### Check Force Balance

Iteration from base to tip

$$\boldsymbol{f}_{i}={\bf I}_{i}\boldsymbol{a}_{i}+\boldsymbol{p}_{i}+\boldsymbol{f}_{i+1}$$

### Pseudo Fortranesque Code

```fortran
function calc_acc(q,qp,Q) result(qpp)
do i=1,n
	pos(i) = pos(i-1) + R(i-1)*i_*length
	ori(i) = ori(i-1)*rot(k_,q(i))
	R(i) = rot(ori(i))
	cg(i) = pos(i) + R(i)*i_*length/2
	s(i) = twist(k_, r(i), 0)
	v(i) = v(i-1) + s(i)*qp(i)
	k(i) = v(i)×s(i)*qp(i)
	I(i) = spi(mass,mmoi_body,cg(i),ori(i))
	w(i) = wrench(mass*g, cg(i), 0)
	p(i) = v(i)×I(i)*v(i)-w(i)
end do
do i=n,1
	A(i) = I(i) + RU(i+1)*A(i+1)
	d(i) = p(i) + T(i+1)*Q(i+1) + RU(i+1)*(A(i+1)*k(i+1)+d(i+1))
	T(i) = A(i)*s(i)/(dot(s(i), A(i)*s(i)))
	RU(i) = 1 - outer(T(i), s(i))
end do
do i=1,n
	qpp(i) = (Q(i)-dot(s(i), A(i)*(a(i-1)+k(i))+d(i)))/dot(s(i), A(i)*s(i))
	a(i) = a(i-1) + s(i)*qpp(i) + k(i)
end do
```
