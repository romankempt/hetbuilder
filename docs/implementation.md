# Theoretical Background

Coincidence lattices are determined with the algorithm outlined by Schwalbe-Koda ([1]).

[1]: https://doi.org/10.1021/acs.jpcc.6b01496 ". Phys. Chem. C 2016, 120, 20, 10895-10908"

Two 2D lattice bases (lattice vectors are given as column vectors) are given by:

```math
\mathbf{A} = \begin{pmatrix} a_{11} & a_{21} \\ a_{12} & a_{22} \end{pmatrix}
```
```math
\mathbf{B} = \begin{pmatrix} b_{11} & b_{21} \\ b_{12} & b_{22} \end{pmatrix}
```

Each point in the 2D plane is given by the coefficients:

```math
P(m_1, m_2) = m_1 \vec{a}_1 + m_2 \vec{a}_2
```
```math
P(n_1, n_2) = n_1 \vec{b}_1 + n_2 \vec{b}_2
```

The two bases can be rotated with respect to each other:

```math
\mathbf{R}(\theta) = \begin{pmatrix} \cos(\theta) & -\sin(\theta) \\ \sin(\theta) & \cos(\theta) \end{pmatrix}
```

Two lattice points of the two bases coincide under the following condition:

```math
\begin{pmatrix} \vec{a}_1 & \vec{a}_2 \end{pmatrix} \begin{pmatrix} m_1 \\ m_2 \end{pmatrix}
= \mathbf{R}(\theta) \begin{pmatrix} \vec{b}_1 & \vec{b}_2 \end{pmatrix}
\begin{pmatrix} n_1 \\ n_2 \end{pmatrix} \\
```
```math
\mathbf{A} \vec{m} = \mathbf{R}(\theta) \mathbf{B} \vec{n}
```

As a tolerance criterion, coincidence is accepted if the distance between the coinciding lattice points is smaller than a threshold:

```math
| \mathbf{A} \vec{m} - \mathbf{R}(\theta) \mathbf{B} \vec{n} | \leq tolerance
```

Solving this system of linear equations yields a set of associated vectors for each angle:

```math
s(\theta) = \{ (\vec{m_1}, \vec{n_1}), (\vec{m_2}, \vec{n_2}), ..., (\vec{m_s}, \vec{n_s}), ..., (\vec{m_k}, \vec{n_k}) \} \\
```

From any pair of these associated vectors that is linearly independent, one can construct supercell matrices from the row vectors:

```math
\mathbf{M} = \begin{pmatrix} m_{s1} & m_{s2} \\ m_{k1} & m_{k2} \end{pmatrix}~~~~
\mathbf{N} = \begin{pmatrix} n_{s1} & n_{s2} \\ n_{k1} & n_{k2} \end{pmatrix}
```

This yields a set $S(\theta)=\{(\mathbf{M}_i, \mathbf{N}_i)\}$ of supercell matrices.

# Implementation Details

The four coefficients $m_{s1}$, $m_{s2}$, $m_{k1}$ and $m_{k2}$ are determined by a grid search. Therefore, one has to iterate through all possible combinations for all given angles $\theta_i$

```math
(-N_{max} \leq s \leq -N_{min} < N_{min} \leq s < N_{max} \\
 \text{ and } -N_{max} \leq k \leq -N_{min} < N_{min} \leq k < N_{max}) ~\forall~\theta_i ~,
```

where $N_{min}$ and $N_{max}$ are the minimum and maximum number of translations, respectively. This yields $((2 \cdot (N_{max} - N_{min})))^4) * N_{angles})$ grid points to search through, which is done in C++ employing OpenMP parallelism.

This results in a set $S(\theta)$ that typically is very large. Before the supercells are generated, it is reduced by practical criteria.

1. All unit cell multiples are removed by ensuring that their absolute greatest common divisor $\text{gcd}$ equals 1:

    ```math
    \text{abs}(\text{gcd}(m_{s1}, m_{s2}, m_{k1}, m_{k2}, n_{s1}, n_{s2}, n_{k1}, n_{k2})) \overset{!}{=} 1
    ```

2. Two supercell matrices yield the same area if they have the same determinant for a given angle. For all supercell matrices with the same determinant, the ones with positive and symmetric entries are preferred:

    ```math
    \mathbf{M}_i = \mathbf{M}_j ~~\text{if det}(\mathbf{M}_i) =~\text{det}(\mathbf{M}_j)
    ```

After the set $S(\theta)$ is reduced, the atomic supercell configurations are generated. This may still be a larger number, so further reduction steps are necessary.

3. The lower and upper supercells are given by the product of the supercell matrices $\mathbf{M}_i$ and $\mathbf{N}_i$ with the primitive bases $\mathbf{A}$ and $\mathbf{B}$, respectively. The common supercell $\mathbf{C}_i$ can be built by a linear combination of the two:

    ```math
    \mathbf{M}_i \mathbf{A} \approx \mathbf{N}_i \mathbf{R}(\theta_i) \mathbf{B}
    ```
    ```math
    \mathbf{C}_i = \mathbf{M}_i \mathbf{A} + w \cdot (  \mathbf{N}_i \mathbf{R}(\theta_i) \mathbf{B} -  \mathbf{M}_i \mathbf{A})
    ```

    Where the weight factor $w$ ranges from 0 to 1 and determines if the common unit cell $\mathbf{C}_i$ is either completely given by the lattice of $\mathbf{A}$ or $\mathbf{B}$ or in between. This yields a set of atomic configurations for all angles $P = \{\mathbf{(C}_i,\theta_i)\}$.

4. The set  $P$ can further be reduced by removing symmetry-equivalencies. This is achieved by standardization via the Space Group library [spglib](https://atztogo.github.io/spglib/python-spglib.html). After standardization, two atomic configurations are considered equivalent if:
   
   - their space groups match
   - they have the same number of atoms
   - they yield the same area

## Definition of the Stress Measure

Dropping the index $i$, one can define the target transformation matrices:

```math
\mathbf{T_A} \mathbf{MA} = \mathbf{C} ~~~~ \mathbf{T_B}  \mathbf{N} \mathbf{R}(\theta) \mathbf{B} = \mathbf{C}
```

These are two-dimensional and can be polarly decomposed:

```math
\mathbf{T}_A = \mathbf{U}_A(\phi)\mathbf{P}_A
```
```math
\mathbf{T}_B = \mathbf{U}_B(-\phi)\mathbf{P}_B
```

In two dimensions, this decomposition is unique. We interpret the matrix $\mathbf{P}$ as strain on the lattice vectors with small rotations being removed by the rotation $\mathbf{U}(\phi)$. This allows to define a stress tensor $\varepsilon$ by substracting unity:

```math
\mathbf{\varepsilon}_A = \mathbf{P}_A - \mathbb{I}
```
```math
\mathbf{\varepsilon}_B = \mathbf{P}_B - \mathbb{I}
```

As an average value, we calculate the stress measure $\bar{\varepsilon}$:

```math
\bar{\varepsilon} = \sqrt{\frac{\varepsilon_{xx}^2 + \varepsilon_{yy}^2 +\varepsilon_{xy}^2 + \varepsilon_{xx} \cdot \varepsilon_{yy}}{4}}
```

And a total stress measure on both unit cells:
```math
\bar{\varepsilon}_{tot} = \bar{\varepsilon}_A + \bar{\varepsilon}_B
```