See the [exercises definition](https://www.pbr-book.org/4ed/Monte_Carlo_Integration/Exercises#)

# 1-2
The code is tested and approaches the real solution for basic polynomials. Importance sampling is tested when sampling from a uniform or Gaussiam disribution. The .py scripts plots how the invert of the variance grows with the number of samples for the MC estimator. Since this is a straight line, we can deduce that the variance evovles in $O(1/n)$.

# 3
## Case where $a,b<0$

The goal is to sample $|f(x)|=-f(x)$, since $f$ is always negative. $-f$ is the linear between $\tilde{a}=-a>0$ and $\tilde{b}-b>0$ on [0,1]: Hence we can apply the result we got when $a$ and $b$ were positive: If $\xi$ is uniformly distributed on [0,1], then X will be distributed proportionally to $|f|$ if:
$$
X = \frac{\xi(\tilde{a}+\tilde{b})}{\tilde{a}+\sqrt{(1-\xi)\tilde{a}^2+\xi\tilde{b}^2}} \\=
\frac{\xi(a+b)}{a - \sqrt{(1-\xi)a^2+\xi b^2}}
$$

## Case where $a<0$ and $b>0$


In this case, there will be a point $x_0\in[0,1]$ where $f(x_0)=0$. Hence, $|f|=-f\times 1_{[0, x_0]} + f\times 1_{[x_0, 1]}$.

We can solve for $x_0$:

$$
(1-x_0)a+x_0b = x_0(b-a) + a= 0 \iff x_0=\frac{a}{a-b}
$$

Now let's compute $p(x)$, the pdf proportional to $|f|$:
$$
\int_0^1|f| = -\int_0^{x_0} f + \int_{x_0}^1 f = \frac{1}{2}(b(1-x_0) - ax_0) = \frac{1}{2}(b-x_0(a+b)) \\=
\frac{1}{2}(b-a\frac{a+b}{a-b}) \triangleq c^{-1}
$$

or 

$$
c = 2 \frac{a-b}{b(a-b)-a(a+b)} = 2 \frac{b-a}{a^2+b^2}
$$
hence $p(x) = c|f(x)|$

Now let's compute the CDF $P(x)$. To make the study easier let's notice that $|f|=g\times1_{[0, 1]}$ where g is symmetric wrt the $x=x_0$ axis. Let $G(x)=\int_{x_0}^{x}g$. Then $G(x_0+x)=-G(x_0-x)$

Now $c^{-1}P(x) = \int_{-\infty}^{x}|f|= \begin{cases}
0 & \text{if } x <0 \\
\int_0^xg =G(x)-G(0)& \text{if } 0 \leq x \leq 1 \\
c^{-1} & \text{if } x > 1
\end{cases}$

However, for $x<x_0$:
$$
G(x) = \int_{x_0}^xg = -\int_{x_0}^x t(b-a)+adt = -\frac{1}{2}(x^2-x_0^2)(b-a) - a(x-x_0)
$$
So 
$$
c^{-1}P(x)=\begin{cases}
0 & \text{if } x < 0 \\
G(x)-G(0)=-\frac{b-a}{2}x^2-ax & \text{if } 0 \leq x \leq x_0 \\
-G(2x_0-x) - G(0) & \text{if } x_0 \leq x \leq 1 \\
1 & \text{if } 1 < x \\
\end{cases}
$$

It is important to note at this point that $p$ is positive on [0, 1] and only cancels when $x=x_0$. So $P$ is stricly growing. Now comes the inversion : let $y=P(x)$ for $x\in[0,1]$. If $y\leq P(x_0)$, then $x\in[0, x_0]$. Hence:
$$
x=P^{-1}(y)=G^{-1} (c^{-1}y + G(0)) = G^{-1}(z)
$$
where $G^{-1}$ is the invert of $G$ on $[0, x_0]$. 
$$
\Delta=a^2-4(\frac{a-b}{2})(x_0^2\frac{b-a}{2}+ax_0-z) \\=
a^2+2(b-a)(\frac{a^2}{2(b-a)}+\frac{a^2}{a-b}-z) \\=
a^2+a^2 - 2a^2-2(b-a)z \\=
2(a-b)z \geq 0
$$

$$
\boxed{\begin{array}{l}
\text{Note: } \Delta \geq 0 \\
\text{Indeed, } c^{-1}y \leq c^{-1}P(x_0)=G(x_0)-G(0)=-G(0) \\
\text{Hence } z=c^{-1}y+G(0) \leq 0\text{, and } (a-b)z \leq 0.
\end{array}}
$$

There are 2 possible solution, but the only valid one is in $[0, x_0]=[0, \frac{a}{a-b}]$ :
$$
G^{-1}(z) = \frac{a}{a-b} - \frac{1}{b-a}\sqrt{2z(a-b)}.
$$
This defines $G^{-1}$ for $y\in[0, P(x_0)]$, and threfore on [0, 1] using the symmetry. Hence, the transform we need to apply is:

$$
X = \begin{cases}
G^{-1} (c^{-1}\xi + G(0)) & \text{if } \xi < P(x_0) \\
2x_0 - G^{-1} (-c^{-1}\xi - G(0)) & \text{{else.}}
\end{cases}
$$

Which reduces to
$$
X = \begin{cases}
\frac{1}{a-b}(a+\sqrt{a^2-(a^2+b^2)\xi}) & \text{{if }}\xi < \frac{a^2}{a^2+b^2} \\
\frac{1}{a-b}(a-\sqrt{(a^2+b^2)\xi - a^2}) & \text{{else.}}
\end{cases}
$$


## An alternative strategy
It is intuitively possible to proceed in 2 steps to avoid all the complex calculation, and get more sens out of this formula. The probability domain is bimodal (split in two): let's use that.

Let's compute the probability of sampling $x\leq x_0$:
$$
c\int_0^{x_0}|f|=-\frac{cx_0a}{2}=-\frac{a}{2}\frac{2(b-a)}{a^2+b^2}\frac{a}{a-b}=\frac{a^2}{a^2+b^2} =P(x_0)
$$
We could compare $\xi$ to $P(x_0)$ to assign it conditionally to one of the 2 mode, with the correct probability. Once this choice is made, $\xi$ is still uniformly distributed, but on a narrower interval $[0, x_0]$ or $[x_0, 1]$. To do so, we can rescale it to obtain $\eta$, uniform on $[0, 1]$, then apply a transform that maps $\eta$ into a linear distribution on $[0, x_0]$ or $[x_0, 1]$.

Let's start with the case where $\xi<P(x_0)$. The transform is:
$$
\eta=\xi\frac{a^2+b^2}{a^2}
$$
We treat $\eta$ as uniformly distributed on [0, 1], and want to transform into the linear distribution that decreases from a positive value to $0$. We know the formula:
$$
\tilde{X}=1-\sqrt{1-\eta}=1+\frac{1}{a}\sqrt{a^2 - (a^2+b^2)\xi}.
$$
To make the distribuion linear on $[0, x_0]$ instead of $[0, 1]$, we can simply scale the variable:
$$
X = x_0\tilde{X}=\frac{a}{a-b}+\frac{1}{a-b}\sqrt{a^2 - (a^2+b^2)\xi}
$$


When $\xi\geq P(x_0)$, The transform is
$$
\eta = (\xi - \frac{a^2}{a^2+b^2})\frac{a^2+b^2}{b^2}=\xi\frac{a^2+b^2}{b^2} - \frac{a^2}{b^2}
$$
and to make the distribution linear and increasing from $0$ on $[0, 1]$:
$$
\tilde{X} = \sqrt{\eta}=\frac{1}{b}\sqrt{(a^2+b^2)\xi-a^2}
$$
Lastly, we apply an affine transform to have $X$ linear of $[x_0, 1]$:
$$
X=x_0+\tilde{X}(1-x_0)=\frac{a}{a-b}-\frac{1}{a-b}\sqrt{(a^2+b^2)\xi-a^2}
$$

The final formula is
$$
X = \begin{cases}
\frac{1}{a-b}(a+\sqrt{a^2-(a^2+b^2)\xi}) & \text{{if }}\xi < \frac{a^2}{a^2+b^2} \\
\frac{1}{a-b}(a-\sqrt{(a^2+b^2)\xi - a^2}) & \text{{else.}}
\end{cases}
$$


