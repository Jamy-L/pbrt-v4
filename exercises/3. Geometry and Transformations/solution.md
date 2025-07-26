# 1

The axis-aligned bounding boxes are characterised by
- The origin $o$ (the corner with the minimal coordinates, that actually corresponds to the public member pMin of the Bounds3 object)
- 3 side measure ($l_x, l_y, l_z>0$). They corresponds the coordinates of the vector pMax-pMin

The 8 corners can be defined as;
$$
p^{i,j,k} = o + i l_x e_x + jl_y e_y + kl_z e_z,
$$
where $e_x,e_y,e_z$ are the canonic unit vectors, and $i,j,k\in\{0,1\}^3$.
Given a transformation matrix $T$, we want to find $p^{\prime}_0$ such that
$$
o^\prime = \begin{pmatrix}
\min_{i,j,k}\{Tp^{i,j,k}\cdot e_x\} \\
\min_{i,j,k}\{Tp^{i,j,k}\cdot e_y\} \\
\min_{i,j,k}\{Tp^{i,j,k}\cdot e_z\} \\
1
\end{pmatrix}
$$

Let's start to develop $o^\prime_x$:
$$
o^\prime_x = \min_{i,j,k}\{Tp^{i,j,k}\cdot e_x\} = \min_{i,j,k}\{T_{0,0}p^{i,j,k}_x + T_{0,1}p^{i,j,k}_y+T_{0,2}p^{i,j,k}_z+T_{0,3}\} \\=
\min_{i,j,k}\{T_{0,0}(o_x+il_x) + T_{0,1}(o_y+jl_y)+T_{0,2}(o_z+kl_z)\}+T_{0,3} \\=
\min_{i,j,k}\{T_{0,0}il_x + T_{0,1}jl_y+T_{0,2}kl_z\}+To\cdot e_x \\=
\min_i\{iT_{0,0}l_x\} + \min_j\{jT_{0,1}l_y\}+\min_k\{kT_{0,2}l_z\}+To\cdot e_x
$$
Since $l_x,l_y,l_z$ are positive, the min depends on the sign of the $T$ coefficient. Let's call $T^+$, a version of $T$ where each negative coefficient is forced to $0$, and $T^-$ where each negative coefficient is forced to 1. Then, developping $o^\prime_y$ and $o^\prime_z$ reveals that
$$
o^\prime = T^-\begin{pmatrix}
l_x\\
l_y\\
l_z\\
0
\end{pmatrix}
+To =
T^-(q-o)+To = T^-q + (T-T^-)o = T^-q + T^+o
$$
A similar developpment for the point $q^\prime$ (that correspond to pMax in the code, i.e the vertex with the maximum coordinates) yields the expression
$$
q^\prime = T^+\begin{pmatrix}
l_x\\
l_y\\
l_z\\
0
\end{pmatrix}
+To = 
T^+(q-o)+To = T^+q + (T-T^+)o = T^+q + T^-o .
$$

We will write the code as a for-loop, iterating on the coefficients of $T$ so let's see what happens under the hood for one couple $i,j$:
- If $T[i,j] \geq 0$, then $T[i,j]o[j]$ is added to $o^\prime[j]$ and $T[i,j]q[j]$ to $q^\prime[j]$.
- Else, then $T[i,j]o[j]$ is added to $q^\prime[j]$ and $T[i,j]q[j]$ to $o^\prime[j]$.

Both of these quantities must therefore necessarily be computed at each iteration, but they are added to a different buffer depending on the sign of $T[i,j]$. This could typically be implemented by a conditional move at each iteration, but interestingly, this sign condition can be expressed differently! Since, by definition, $q[j]\geq o[j]$ for any $j$, then $T[i,j] \geq 0$ is equivalent to $T[i,j]q[j]\geq T[i,j]o[j]$. This conditional move can cleverly be replaced by min and max operators, as is done in the algorithm proposed in the Graphics Gems I. This requires to make the same comparison twice which may sound less optimised than the conditional move, but min and max operations are very presumably better fitted for vectorisation.

It is also worth noting that the final column of $T$ is a special case: The 4th component of both $q$ and $o$ is 1, so $T[i,j]$ is added to the both buffers regardless of its sign. This is why I have put this part out of the loop in my implementation, and delegated this cost directly to the constructor.

Here is my final implementation:
```c++
PBRT_CPU_GPU Bounds3f Transform::operator()(const Bounds3f &b) const {
    Bounds3f bt;
    // Initialize with the translation component directly, instead of adding at the end.
    Point3f pMin_new {m[0][3], m[1][3], m[2][3]};
    Point3f pMax_new {m[0][3], m[1][3], m[2][3]};

    for (uint i = 0; i < 3; ++i) {
        for (uint j = 0; j < 3; ++j) {
            Float mpmin = m[i][j] * b.pMin[j];
            Float mpmax = m[i][j] * b.pMax[j];
            pMin_new[i] += std::min(mpmin, mpmax);
            pMax_new[i] += std::max(mpmin, mpmax);
        }
    }
    bt.pMin = pMin_new;
    bt.pMax = pMax_new;

    return bt;
}
```

However, it may not be worth the troubles to add this implementation in the code ; I don't think bounding box transforms take a relevant amount of time during rendering.