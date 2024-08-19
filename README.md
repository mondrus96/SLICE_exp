# Sparse + Low Rank Inverse Covariance Estimation (SLICE) Experiments
SLICE proposes an effective algorithm for the sparse + low rank Gaussian graphical model. The objective we seek to minimize is:

```math
\underbrace{- \mathcal{L}(\boldsymbol{\hat{S}};(\boldsymbol{\tilde{\Sigma}}^{-1} - \boldsymbol{\hat{L}})^{-1}) + \rho \|\boldsymbol{\hat{S}}\|_1}_{\text{penalized negative log likelihood}} + \underbrace{\|\boldsymbol{\tilde{\Sigma}}(\boldsymbol{\hat{S}} + \boldsymbol{\hat{L}}) - \boldsymbol{I}\|_F^2}_{\text{covariance fidelity}} \\
\text{s.t. } \mathcal{R}(\boldsymbol{\hat{L}}) = r, \ \text{where $0 < r < p$ }\\
```

We propose an efficient pseudo-EM algorithm which proceeds by alternating between two steps in block coordinate descent fashion. The first step is a closed-form solution for the low rank update

```math
\boldsymbol{\hat{L}} = \textit{SVD}_r(\boldsymbol{\tilde{\Sigma}}^{-1} - \boldsymbol{\hat{S}}).
```

The second step is a sparse Gaussian graphical model problem of the sample covariance matrix conditioned on the low rank matrix. This can be solved using standard methods such as GLASSO, CLIME, among others

```math
{\boldsymbol{\hat{S}}}^{(i+1)}={\underset{\boldsymbol{\hat{S}}}{\operatorname {arg\,max}}}\ Q({\boldsymbol{\hat{S}}}\mid {\boldsymbol{\hat{S}}}^{(i)}).
```