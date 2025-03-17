# Introduction to B-Splines 

According to @deboor1978, the definition of the $i$th B-spline (basis spline) of order $k$ with knots $\{t_i\}$ is given by
$$
B_{i, k, t}(x) = (t_{i+k} - t_{i})\;[t_i,\cdots,t_{i+k}](t - x)_{+}^{k-1}
$$
Where the brackets denote the divided difference on the function $t\rightarrow (t-x)^{k-1}_+$, which is the leading coefficient of the polynomial of degree $d$ (order of the set minus 1) which agrees with the function at each of the points $\{t_i, \cdots, t_{i+k}\}$. 

The divided difference can be defined recursively as
$$
[t_i,\ldots,t_{i+k}]g(t) = \begin{cases}
\frac{[t_{i},\cdots,t_{r-1},t_{r+1}, \cdots_{i+k}]g(t)
- [t_i,\cdots, t_{s-1},t_{s+1},\cdots,t_{i+k}]g(t)
}{t_{s} - t_r} 
& t_s \ne t_r, t_s, t_r \in \{t_i\ldots t_{i+k}\} \\
g^{(k)}(t_i)/k! & t_i=\ldots = t_{i+k}, g\in C^k
\end{cases}
$$
and $[t_i]g(t) = g(t_i)$ is the case for the zeroth-degree (order 1) polynomial and all $t$ are unique and $t$ is ordered.  This definition works as if $p_n(t)$ is the k-th degree polynomial agreeing at al t less than n, then $p_{n}(t) = p_{n-1}(t) + [t_i,\cdots,t_{i+n}]g(t)\;\Pi_{i<n+1} (t - t_i)$

We note the following properties of the basis:

- $B_{i}$ is zero if $x \notin (t_i, t_{i+d+1})$. 
- $B_i$ is $C^n$ continuous at $x$ where $n = k - \nu$ where $\nu$ is the number of $t_i = x$. 
- As a result, $B_i$ is represented as a degree $d$ piecewise polynomial with breakpoints at each $t_i$.
- $B_i$ goes to zero at either end of the domain
- $\sum B_i(t_i) = 1$ 
- Each $B_i$ are linearly independent, such that $B_i$ forms the basis for the piecewise polynomial space with $d + 1 - \nu$ constraints at each breakpoints $\xi_i$

As such, any spline can be represented as a sum of basis splines, so for any $f \in \mathbb{P}_{k, \xi, \nu}$, 
$$
f(x) = \sum_i \alpha_i B_i(x)
$$
 where $\alpha_i \in \mathbb{R}$.



To compute the values of B-splines more efficiently, we can use the recursive definition:
$$
B_{i, 1}(x) = \begin{cases}
c & x_i \leq x < x_{i+1} \\
0 & {\rm otherwise}
\end{cases}
$$
where $\theta$ is the heaviside step function (1 for nonzero numbers, 0 otherwise). and $c$ is some normalization coefficient.
$$
B_{i, k}(x) = \frac{x-x_i}{x_{i + k-1} - x_i} B_{i, k-1}(x) + \frac{x_{i+k} - x}{x_{i+k} - x_{i+1}} B_{i+1, k-1}(x)
$$


Another major advantage of B-splines is that their integrals and derivatives have closed-form representations as B-splines as well.
$$
\frac{d}{dx} \sum_i \alpha_i B_{i, k}= \sum_i \alpha_i^-\;B_{i, k-1} \\
\alpha_i^- = (k-1)\,\frac{\alpha_i - \alpha_{i-1}}{t_{i+k-1} - t_i}
$$

which follows from the recursive definition of the divided difference and the derivative of $(t-x)^d_+ = -d (t-x)^{d-1}_+$ (de Boor, eq. X.10). Since the derivative of a spline is a spline, we can also solve for the antiderivative of a spline:
$$
\int_{t_1}^x \sum_i \alpha_i B_{i, k}(x) dx = \sum_i^{s-1} \alpha^+_i B_{i, k+1}(x)\\
\alpha^+_i = \sum_{j=1}^i \alpha_j \frac{t_{j+k} - t_j}{k}
$$

provided $x \leq t_s$ (eq. X.22). If $\alpha_i = \delta_{i, l}$, then $\alpha_i^+ = \frac{t_{l+k} - t_l}{k}$ and 
$$
\int B_{l, d}(x) \,dx = \alpha_l^+  \sum_s B_{i, k}(t_s) = \alpha_l^+ = \frac{t_{l + k} - t_l}{k}
$$
So that the total area of a spline is given by
$$
\int_{-\infty}^{\infty} \sum_i \alpha_i B_{i, d}(x)\,dx = \frac{1}{k}\sum_i \alpha_i (t_{i+k} - t_i)
$$
