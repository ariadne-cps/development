/***************************************************************************
 *            linear_programming.h
 *
 *  Copyright  2004-7  Pieter Collins
 *  Pieter.Collins@cwi.nl
 ****************************************************************************/

/*
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */


/*! 

\file linear_programming.h
\brief Documentation on linear programming



\page linearprogramming Linear Programming

\section standardprimaldual Standard primal and dual problems

The standard linear programming problem is
\f[ \text{(P)} \qquad \min c^Tx \text{ s.t. } Ax=b;\ x\geq0.  \f]
Without loss of generality we can take \f$b\geq0\f$.
We let \f$x^*\f$ be an optimal point.

The dual to the standard problem is
\f[ \text{(D)} \qquad \max b^T y \text{ s.t. } A^Ty\leq c \f]
or alternatively
\f[ \text{(D)} \qquad \max b^T y \text{ s.t. } A^Ty+z=c;\ z\geq 0. \f]
The variables \f$y\f$ are called <em>dual variables</em> and the variables \f$z\f$ are <em>slack variables</em>.

<b>Theorem</b> If \f$x\f$ is feasible for (P) and \f$y\f$ is feasible for (D), then \f[c^Tx \geq b^Ty.\f]
Further, if \f$x^*\f$ is optimal for (P) and \f$(y^*,z^*)\f$ is optimal for (D), then \f[c^Tx^* = b^Ty^* \qquad \text{and} \qquad x^*\!\cdot\! z^*=0 .\f]

The second condition is called <em>complementary slackness</em>.


\section optimalbasis Optimal basic solutions 

A basic solution is given by a set of column indices \f$B\f$ such that the square matrix \f$A_B\f$ formed by the \f$B\f$ columns of \f$A\f$ is nonsingular.
Then
\f[x_B=A_B^{-1}b, \ x_N=0; \quad y=(A_B^T)^{-1}c_B; \quad z_B=0, \ z_N=c_N^T-c_B^TA_B^{-1}A_N; \qquad c^Tx = c_B^TA_B^{-1}b = b^T(A_B^T)^{-1}c_B = b^Ty. \f]

Suppose \f$x^*\f$ is optimal. Then \f$x_B\geq0\f$ and \f[c^Tx = c_B^T x_B + c_N^T x_N = c_B^TA_B^{-1}(b-A_Nx_N)+c_N^Tx_N = c_B^TA_B^{-1}b + (c_N^T-c_B^TA_B^{-1}A_N)x_N = c^Tx^* + (c_N^T-c_B^TA_B^{-1}A_N)x_N,\f] and hence \f$z_N=c_N^T-c_B^TA_B^{-1}A_N\geq0\f$.


\section robustprimaldual Robust primal and dual problems

The robust standard linear programming problem is
\f[ \text{(RP)} \qquad \min c^Tx \text{ s.t. } Ax=b;\ x>0.  \f]
and additionally find column indices \f$B\f$ such that \f$A_B\f$ is nonsingular.
The dual to the standard problem is
\f[ \text{(RD)} \qquad \max b^T y \text{ s.t. } A^Ty<c \f]
or alternatively
\f[ \text{(RD)} \qquad \max b^T y \text{ s.t. } A^Ty+z=c;\ z>0. \f]


\section feasibleprimaldual Primal and dual feasibility problems

The primal feasibility problem is
\f[ \text{(PF)} \qquad Ax=b;\ x\geq0.  \f]
with dual problems
\f[ \text{(DPF)} \qquad A^Ty \leq0;\ -b^Ty\leq-\!1 \quad \text{ or } \quad \max b^Ty \text{ s.t. } A^Ty\leq 0.  \f]
Problem (PF) is solvable iff the first form of (DPF) is unsolvable or the optimal value in the second form is strictly positive.

The dual feasibility problem is
\f[ \text{(DF)} \qquad A^Ty\leq c \f]
which has dual
\f[ \text{(DDF)} \qquad Ax=0,\ c^Tx=-\!1;\ x\geq0 \quad \text{ or } \quad \min c^Tx \text{ s.t. } Ax=0,\ x\geq 0. \f]
Problem (DF) is solvable iff the first form of (DDF) is unsolvable or the optimal value in the second form is negative.

\section robustprimaldual Robust feasibility certificates

A <em>robust feasibility certificate</em> for the primal problem is a base \f$B\f$ such that \f$A_B\f$ is nonsingular, and a vector \f$x_N>0\f$ such that \f$-A_B^{-1}A_Nx_N>0\f$, and an infeasibility certificate is a vector \f$y>0\f$ such that \f$A^Ty<c\f$ and \f$b^Ty>0\f$.

A <em>robust feasibility certificate</em> for the dual problem is a vector \f$y\f$ such that \f$A^Ty<c\f$, and an infeasibility certificate is a vector \f$x_N\f$ such that \f$-A_B^{-1}A_Nx_N>0\f$ and \f$c^Tx<0\f$.


\section dualcertificate Certificates of infeasibility / Farka's Lemma

A certificate of infeasibility of \f$Ax=b;\ x\geq0\f$ is a vector \f$y\f$ such that \f$A^Ty\leq0\f$ and \f$b^Ty>0\f$.
For then \f$0\geq y^TAx = y^Tb > 0\f$, a contradiction.

A certificate of infeasibility of \f$A^Ty\leq c\f$ is a vector \f$x\geq0\f$ such that \f$Ax=0\f$ and \f$c^Tx<0\f$.
For then \f$0=y^TAx \leq c^Tx < 0\f$, a contradiction.

\section robustcertificate Certificates of robust (in)feasibility

A robust certificate for the primal feasibility problem \f$Ax=b;\ x\geq0\f$ is a base \f$B\f$ such that \f$A_B\f$ is nonsingular, and an \f$x_N>0\f$ such that \f$-A_B^{-1}A_Nx_N>0\f$; this is equivalent to \f$Ax=b;\ x>0\f$.
A robust certificate of infeasibility is a point \f$y\f$ such that \f$A^Ty<0\f$ and \f$b^Tx>0\f$.

<b>Theorem</b>
Suppose \f$Ax=b;\ x\geq0\f$ is robust. Then either \f$A\f$ has full row rank and there exists \f$x>0\f$ such that \f$Ax=b\f$, or there exists \f$y\f$ such that \f$A^T<0\f$ and \f$b^Ty>0\f$.

A robust certificate for the dual feasibility problem \f$A^Ty\leq c\f$ is a point \f$y\f$ such that \f$A^Ty<c\f$.
A robust certificate of infeasibility is a base \f$B\f$ such that \f$A_B\f$ is nonsingular, and a vector \f$x>0\f$ such that \f$Ax=0\f$ and \f$c^Tx<0\f$.
We prove \f$Ax=0\f$ by setting \f$x_B=-A_B^{-1}A_Nx_N\f$, so that \f$Ax=A_Bx_B+A_Nx_N=-A_BA_B^{-1}A_Bx_N+A_Nx_n=0\f$.

<b>Theorem</b>
Suppose \f$A^Ty\leq c\f$ is robust. Then either there exists \f$y\f$ such that \f$A^Ty<c\f$, or \f$A\f$ has full row rank and there exists \f$x>0\f$ such that \f$Ax=0\f$ and \f$c^Tx<0\f$.

<i>Proof</i> 
If \f$A^Ty<c\f$, then this holds also for perturbations of \f$y\f$.
If \f$Ax=0\f$ and \f$A_B\f$ is nonsingular, then \f$x_B=-x_NA_NA_B^{-1}\f$.
Perturbing \f$A,b\f$ and keeping \f$x_N\f$ constant, we obtain a perturbation of \f$x_B\f$, and hence a certificate for the perturbed problem.
<br>
Conversely, suppose the problem is robustly solvable.
Then the problem \f$A^Ty\leq c-\epsilon p\f$ is solvable for \f$p>0\f$ and \f$\epsilon\f$ sufficiently small. Hence there exists \f$y\f$ such that \f$A^Ty<c\f$.
<br>
Suppose the problem is robustly unsolvable
Let \f$P\f$ be a matrix with all positive entries. Since \f$A^Ty\leq c\f$ is robustly unsolvable, \f$(I+\epsilon P)^T A^Ty\leq (I+\epsilon P)^T c\f$ is unsolvable for some \f$\epsilon>0\f$. Then there exists \f$x_\epsilon\f$ such that \f$A(I+\epsilon P)x_\epsilon=0\f$, \f$x_\epsilon\geq0\f$ and \f$(I+\epsilon P)x_\epsilon c<0\f$. Then if \f$x=(I+\epsilon P)x_\epsilon\f$, then \f$Ax=0\f$, \f$c^Tx<0\f$ and \f$x>0\f$ since \f$x_\epsilon\geq0;\ x_\epsilon\neq0\f$ and \f$(I+\epsilon P)>0\f$.


\section robustsolve Converting robust feasibility problems to feasibility problems

To solve the robust primal feasibility problem,
\f[ \text{(RPF)} \qquad Ax=b;\ x>0 \f]
we choose \f$p>0\f$ and consider the problem
\f[ \min -s \text{ s.t. } Ax + Ap\,s = b;\ x,s\geq0 .\f]
Let \f$\hat{x}^T=(x\;s)\f$, \f$\hat{A}=(A\;Ap)\f$ and \f$\hat{c}^T=(0\;\mbox{}-\!1)\f$.
Then we obtain the standard primal optimisation problem 
\f$ \min \hat{c}^T\hat{x} \text{ s.t. } \hat{A}\hat{x} = b;\ \hat{x}\geq0 . \f$
If the optimal value is negative, then we have found \f$x^*,s^*\f$ such that \f$A(x^*+ps^*)=b; \ x^*\geq0,\ s^*>0\f$, so taking \f$\tilde{x}=x^*+ps^*\f$, we have \f$A\tilde{x}=b;\ \tilde{x} = x^*+s^*p \geq s^*p > 0\f$.
If the optimal value is non-negative, then we can attempt to solve the dual robust optimisation problem
\f$ \max b^Ty \text{ s.t. } Ay<0 \f$.
Since if this problem is feasible, it has unbounded solutions, we can instead choose \f$q>0\f$ look for a positive optimal value of
\f[ \max b^Ty \text{ s.t. } Ay\leq-q . \f]

To solve the robust dual feasibility problem, 
\f[ \text{(RDF)} \qquad A^T y + z = c; \ z>0 \quad \text{or} \quad  A^Ty < c, \f]
we choose \f$q>0\f$ and consider the problem
\f[ \max t \text{ s.t. } A^T y + qt \leq c . \f]
Let \f$ \hat{A}^T = (A^T \; q)\f$, \f$\hat{b}^T = (0^T \; 1)\f$ and \f$\hat{y}^T = ( y^T\;t )\f$.
Then we obtain the standard dual optimisation problem 
\f$ \max \hat{b}^T \hat{y} \text{ s.t. } \hat{A}^T \hat{y} \leq c; \f$
If the optimal value is positive, then we have found \f$y^*,t^*\f$ such that \f$A^Ty^*\leq c-t^*q < c\f$.
If the optimal value is zero or negative, then we can attempt to solve the primal robust feasibility problem 
\f$ \min b^T x \text{ s.t. } A x = 0, \ q^T x = 1, \ x\geq 0 . \f$
However, even if the dual feasibility problem is unsolvable (i.e. \f$t^* < 0\f$), the primal may become solvable by a perturbation of \f$A\f$.
We therefore consider a robust version
\f$ \min b^T x \text{ s.t. } A x = 0,\ x>0 . \f$
and introduce \f$p>0\f$ to make a problem
\f$ \min b^T x \text{ s.t. } A x = 0,\ x-p\geq0 . \f$
Note that if this problem has negative value, then we have found \f$x^*\f$ such that \f$b^Tx^*<0\f$, \f$Ax^*=0\f$ and \f$x^*>0\f$, which implies that the original dual problem has no solution.
Taking \f$\tilde{x} = x-p\f$, we obtain
\f[ \min b^T \tilde{x} + b^T p \text{ s.t. } A\tilde{x} = -Ap,\ \tilde{x}\geq 0 . \f]



\section reducedform Reduced linear programming problem

The reduced form of the standard linear programming problem is
\f[ \min c^Tx \text{ s.t. } x_B+\tilde{A}x_N=b;\ x\geq0 . \f]
Here, \f$x_B\f$ are the basic variables, and \f$x_N\f$ the non-basic variables.

If \f$A\f$ is of full row rank, by choosing a basis such that the basis matrix \f$A_B\f$ is invertible, we can put any linear programming problem in standard form.


\section boundedoptimisation  Lower and upper bounds on variables

The constrained primal linear programming problem is
\f[ \text{(CP)} \qquad \min c^Tx \text{ s.t. } Ax=b;\ l\leq x\leq u . \f]
Given lower bounded variables \f$x_L\f$ and upper bounded variables \f$x_L\f$, the problem becomes 
\f[ \min c^Tx \text{ s.t. } Ax=b;\ x_L\geq l_L;\ x_U\leq u_U\f]
and the dual problem is
\f[ \text{(DCP)} \qquad \max\ (b^T - l_L^TA_L^T - u_U^T A_U^T)y+l_L^Tc_L + u_U^Tc_U \ \text{ s.t. }\  A_L^Ty\leq c_L;\ A_U^Ty\geq c_U . \f]
Note that the objective function can be written \f$b^T y +l_L^T (c_L - A_L^Ty) + u_U^T(c_U - A_U^Ty)\f$

<b>Theorem</b>
Suppose \f$Ax=b;\ l\leq x\leq u\f$ is robust. Then either \f$A\f$ is nonsingular, and there exists \f$x\f$ such that \f$Ax=b\f$ and \f$l<x<u\f$, or there exists column index sets \f$L,U\f$ and \f$y\f$ such that \f$ (b^T-l_L^TA_L^T-u_U^TA_U^T)y>0,\ A_L^Ty<0,\ A_U^Ty>0 . \f$

To solve the robust problem, we first find a solution to the standard problem, and then try to force saturated constraints to be positive as before.


\section simplexalgorithm The simplex algorithm

  Suppose we wish to update a basis of the standard linear programming problem.
   - The current point \f$v=\mathrm{A}_B^{-1} (b - \mathrm{A}_N x_N)\f$ (typically, \f$x_N=0\f$).
   - The reduced costs are \f$c_N-(c_B\mathrm{A}_B^{-1})\mathrm{A}_N\f$.
   - If the reduced costs are all positive, the algorithm terminates. Otherwise, select \f$j\f$ such that \f$c_j<0\f$.
   - The direction to move is \f$d=\mathrm{A}_B^{-1}a_j\f$.
   - Choose \f$t\f$ maximal so that \f$l_B \leq v-td\leq u_B\f$; if the update is being used for feasibility, constraints violated by \f$v\f$ may be violated by \f$v-td\f$. Choose \f$k\f$ such that \f$i=\pi_k\f$ corresponds to a saturated constraint.
   - Replace \f$x_j\f$ by \f$x_k\f$ in the basis and update \f$\mathrm{B}:=\mathrm{A}_B{-1}\f$.
     - We have \f$\mathrm{B} \mathrm{A}_{\pi_i}=e_i\f$ for \f$i\neq k\f$ and we want \f$\mathrm{B} \mathrm{A}_{\pi_k}=e_k\f$.
     - Let \f$\mathrm{B}\mathrm{A}_{\pi_k} = a\f$.
     - For \f$i\neq k\f$, subtract \f$\mathrm{B}_{kj}\,a_i/a_k\f$ from \f$\mathrm{B}_{ij}\f$ for all \f$j\f$.
     - Then divide \f$\mathrm{B}_{kj}\f$ by \f$a_k\f$ for all \f$j\f$.


\section simplexefficiency Efficiency of the simplex algorithm

For a linear programming problem of standard form, with \f$A\f$ an \f$m\times n\f$ matrix, the number of iterations of the simplex algorithm for practical problems grows approximately as \f$m\log n\f$.

\section feasibilityalgorithms Algorithms for feasibility 

 - Constrained feasibility problem with equalities  \f$ Ax=b;\ l\leq x\leq u\ (m\leq n)\f$.<br>
   Find a set of basic variables \f$B\f$ so that \f$ A_B\f$ is nonsingular, where \f$ A_B\f$ is the matrix formed from the columns of \f$A\f$ in \f$B\f$.
   Initialise \f$x_N\f$ to \f$l\f$ for non-basic variables, and set \f$x_B=A_B^{-1}(b-A_Nx_N)\f$. Then \f$Ax=b\f$, but possibly not \f$l_B\leq x_B\leq u_B\f$.
   <br>
   Let \f$c_i=-1\f$ if \f$x_i<l_i\f$ and \f$c_i=-1\f$ if \f$x_i>u_i\f$.
   Now minimise \f$c^Tx\f$, but relax the currently violated constraints.

   \b Remark: Since we do not assume the existence of \f$\pm\infty\f$ in our number types, we use \f$l=0,\ u=-1\f$ for the constraint \f$x\geq0\f$; this is the only unbounded constraint we allow.

   See Chvatal [Chapter 8, pp 129] for more details.
 
 - Unconstrained feasibility problem with inequalities  \f$ Ax\leq b;\ l\leq x\leq u\ (m\geq n)\f$.<br>
   Let \f$I\f$ be a set of basic indices such that \f$A_I\f$ is invertible, where \f$A_I\f$ is the matrix formed from the \em rows corresponding to the \f$I\f$.
   Set \f$v=A_I^{-1}b_I\f$, and let \f$S\f$ be the constraints satisfied by \f$v\f$.
   Find \f$i\f$ such that \f$a_iv>v\f$ where \f$a_i\f$ is the \f$i^\mathrm{th}\f$ row of \f$A\f$ and let \f$c=a_i\f$.
   Change basis until either \f$a_i v\leq b\f$ or it is impossible to reduce \f$a_i x\f$ without violating \f$A_Sx\leq b_S\f$.





\section geometricfeasibility Feasibility problems for geometric operations

In the Geometry module, we need to solve the following linear programming problems to test intersection.
\f[ \begin{array}{|l||c|c|c|c|}\hline
      &\text{Polyhedron}&\text{Polytope}&\text{Zonotope}\\\hline\hline
      \text{Point} & Ap\leq b & p=Vs;\ 1\!\cdot\!s=1;\ s\geq0 & p=c+Ge;\ -1\leq e\leq1 \\\hline
      \text{Rectangle} & Ax\leq b;\ l\leq x\leq u & x=Vs;\ 1\!\cdot\!s=1;\ l\leq x\leq u;\ s\geq0 & x=c+Ge;\ l\leq x\leq u; \ -1\leq e\leq1 \\\hline
      \text{Zonotope} & A(c+Ge)\leq b;\ -1\leq e\leq 1 & Vs=c+Ge;\ 1\!\cdot s=1;\ -1\leq e\leq1;\ s\geq0 & c_1+G_1e_1=c_2+G_2e_2;\ -1\leq e_1,e_2\leq1 \\\cline{0-3}
      \text{Polytope} & AVs\leq b;\ 1\!\cdot\!s=1;\ s\geq0 & V_1s_1=V_2s_2;\ 1\!\cdot s_1=1;\ 1\cdot s_2=1;\ s_1,s_2\geq0 \\\cline{0-2}
      \text{Polyhedron} & A_1x\leq b_1;\ A_2x\leq b_2 \\\cline{0-1}
    \end{array}
\f]
We notice that by introducing slack variables, we can convert all problems into a standard linear programming problem with constraints.
 - Standard primal feasibility problem \f$ Ax=b;\ x\geq 0\f$ 
 - Constrained primal feasibility problem \f$ Ax=b;\ l\leq x\leq u\f$ 

 - Standard dual feasibility problem \f$ Ax\leq b\f$
 - Constrained dual feasibility problem \f$ Ax\leq b;\ l\leq x\leq u\f$ 

We can convert the standard dual feasibility problem into a primal linear programming problem \f$\min b^Ty\text{ s.t. } A^Ty=0\f$, but it is not so straightforward to convert a constrained dual feasibility problem into its dual. Instead we add slack variables and solve
\f$ Ax+z=b;\ l\leq x\leq u\f$. We can use the reduced simplex algorithm to take advantage of sparseness.

\section ariadnelpsolvers Linear programming solvers provided by Ariadne.

 - lpstp() Perform one step of the standard linear programming problem.
      Input: \f$\mathrm{A},b,c\f$, InOut: \f$\pi,\mathrm{A}_B^{-1}\f$.

 - lpcstp() Perform one step of the standard linear programming problem with constraints \f$l\leq x\leq u\f$.
     Any constraints which are violated are assumed to remain violated; this allows for constraints with infinities.
     Input: \f$\mathrm{A},b,c\f$, InOut: \f$\pi,\mathrm{A}_B^{-1}\f$.

 - lpupd() Update the matrix \f$\mathrm{A}_B^{-1}\f$ so that \f$\mathrm{A}_B^{-1}a=e_i\f$ by pivoting on the \f$i^\textrm{th}\f$ row.
      Input: \f$a,i\f$; InOut\f$\mathrm{A}_B^{-1}\f$.

 - lpslv() Solve the standard linear programming problem \f$\min c^Tx \text{ s.t. }Ax=b;\ x\geq0\f$.
      Input: \f$A,b,c\f$; Output: \f$\pi,\ \mathrm{A}_B^{-1},\ x^*,\ y^*,\ z^*\f$.

 - lpcslv() Solve the standard linear programming problem with constraints \f$\min c^Tx \text{ s.t. } Ax=b;\ l\leq x\leq u\f$.
      Input: \f$A,b,c\f$; Output: \f$\pi,\ \mathrm{A}_B^{-1},\ x^*,\ z^*\f$.

 - lprslv() Solve the standard linear programming problem with constraints given in reduced form \f$\min \tilde{c}^Tx_N \text{ s.t. } \tilde{A}x_N+x_B=b;\ l\leq x\leq u\f$ \f$\tilde{A}x_N+x_B=b\f$. (This is useful if a starting basis can easily be found, and there are almost as many constraints as variables.)
      Input: \f$\tilde{A},b,c,\pi\f$; Output: \f$\pi, x^*, z^*, y^*\f$.

 - lpfeas() Solve the feasibility problem \f$Ax=b;\ x\geq0\f$.

 - lpcfeas() Solve the feasibility problem \f$Ax=b\f$ with \f$l\leq x\leq u\f$.

 - lpdfeas() Solve the dual feasibility problem \f$Ax\leq b\f$ directly.

*/
