# py_frobenius
**py_frobenius** is a package implementing the Frobenius number and associated computations in Python.

Let $a_1, a_2, \ldots, a_n$ be coprime positive integers ($a_i>0$ for all $i=1,\ldots,n$ and $\mathrm{gcd}(a_1,\ldots,a_n)=1$). The *Frobenius number* of $a_1,\ldots,a_n$, denoted $$g(a_1,\ldots,a_n),$$ is the largest integer $g$ for which the following equation has no solution $(x_1,\ldots,x_n)\in\mathbb{Z}_{\geq0}^n$:
$\begin{equation}a_1x_1+\cdots+a_nx_n=g.\end{equation}$
For $n=2$, Sylvester proved in [Syl84] that
$$g(a_1,a_2)=a_1a_2-a_1-a_2.$$
For $n>2$, Curtis proved in [Cur90] that $a_1,\ldots,a_n$ and $g(a_1,\ldots,a_n)$ are algebraically independent over $\mathbb{C}$, and thus that there is not simple closed formula for $g(a_1,\ldots,a_n)$.
<br /><br />
In [Syl84], Sylvester also showed that, for $n=2$, the number of positive integers $g$ for which $(1)$ has no non-negative solutions (denoted $N(a_1,\ldots,a_n)$ in general) is
$N(a_1,a_2)=\frac{(a_1-1)(a_2-1)}{2}.$
## Included Functionality
Within this package the following functions are available:
### Frobenius Number
The frobenius number of $a_1,\ldots,a_n$. If $n=2$, we use Sylvester's formula, otherwise we use the BFD (Breadth-First Decreasing) algorithm presented in [BHNW05].
### Finding Representations
One can find all non-negative integer solutions to $(1)$. This is performed by a recursive algorithm: given an integer $N$, all representations of the integers
$$N-a_1, N-a_2, \ldots, N-a_n,$$
are found. The resulting list is then screened for duplicates. Efficiencies can be found by using graph methods that track duplicates during the algorithm and avoid searching already completed branches or branches known to provide no representations.
### Sylvester Denumerants
The Sylvester denumerant function was first investigated by  (and since named after) Sylvester in [Syl57], it is the function
$d(m\mid a_1,\ldots,a_n)=\#\left\{(x_1,\ldots,x_n)\in\mathbb{Z}_{\geq0}^n\mid a_1x_1+\cdots+a_nx_n=m\right\},$ that is the number of ways that $m$ can be represented as a non-negative integer combination of $a_1,\ldots,a_n$. The series of denumerants has the following generating function: $\sum_{m=0}^{\infty}d(m\mid a_1,\ldots,a_n)z^m=\prod_{i=1}^{n}\frac{1}{1-z^{a_i}}.$ The denumerant of $m$, therefore, is the $m\mathrm{th}$ coefficient of the generating function. We compute this in the following way:
- We can express $\frac{1}{1-z^{a_i}}$ as $\sum_{k=0}^{\infty}z^{ka_i}$. Any representation of $m$ as a non-negative integer combination of $a_1,\ldots,a_n$ must contain no more than $\left\lfloor\frac{m}{a_i}\right\rfloor$ copies of $a_i$, and so we can restrict the support to only powers of $z$ that are at most $m$.
- This gives us $n$ polynomials $\sum_{k=0}^{p_i}z^{ka_i}$ where $p_i=\left\lfloor\frac{m}{a_i}\right\rfloor$, the product of which yields a polynomial where the first $m$ coefficients match the first $m$ coefficients of the generating function for the denumerants. Thus, the denumerant is the $m\mathrm{th}$ coefficient of this polynomial.
- To compute the product of these polynomials efficiently, we employ the Cooley-Tukey Fast Fourier Transform (FFT) algorithm [CT65].
### Unrepresentable Numbers
The number of numbers that cannot be represented is the value of $N(a_1,\ldots,a_n)$. To compute this number in general (not just the specific $n=2$ case), we first compute the Frobenius number $g(a_1,\ldots,a_n)$, and then follow the procedure for obtaining the first $g(a_1,\ldots,a_n)$ correct coefficients of the generating function for the denumerants. Then $N(a_1,\ldots,a_n)$ is the number of terms with degree at most $g(a_1,\ldots,a_n)$ with coefficient $0$, since every number greater than the Frobenius number can be represented. The exact numbers that cannot be represented are the degrees of the terms with coefficient $0$.
# References
[BHNW05] D. Beihoffer, J. Hendry, A. Nijenhuis, & S. Wagon, Faster Algorithms for Frobenius Numbers, *The Electronic Journal of Combinatorics*, *12*(1), #R27 (2005).<br />
[CT65] J. Cooley, J. Tukey, An algorithm for the machine calculation of complex Fourier series, *Math. Comp.* **19** (1965), 297-301.<br />
[Cur90] F. Curtis, On formulas for the Frobenius number of a numerical semigroup, *Math. Scand.* **67** (1990), 190-192. <br />
[Syl57] J.J. Sylvester, On the partition of numbers, *Quart. J. Pure Appl. Math.* **1** (1857), 141-152.<br />
[Syl84] J.J. Sylvester, Problem 7382, *Educational Times*, **37** (1884), 26; reprinted in: Mathematical questions with their solution, *Educational Times* (with additional papers and solutions) **41** (1884), 21. <br />
