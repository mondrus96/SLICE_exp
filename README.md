# Sparse + low rank inverse covariance estimation (SLICE) experiments
SLICE proposes an effective algorithm for the sparse + low rank Gaussian graphical model. The objective we seek to minimize is:

\[
        \underbrace{- \mathcal{L}(\bs{\hat{S}};(\bs{\tilde{\Sigma}}^{-1} - \bs{\hat{L}})^{-1}) + \rho \|\bs{\hat{S}}\|_1}_{\text{penalized negative log likelihood}} + \underbrace{\|\bs{\tilde{\Sigma}}(\bs{\hat{S}} + \bs{\hat{L}}) - \bs{I}\|_F^2}_{\text{covariance fidelity}} \\
        \text{s.t. } \mathcal{R}(\bs{\hat{L}}) = r, \ \text{where $0 < r < p$ }\\
\]

We propose an efficient pseudo-EM algorithm which proceeds by alternating between two steps in block coordinate descent fashion. The first step is a closed-form solution for the low rank update

\[
\boldsymbol{\hat{L}} = \textit{SVD}_r(\bs{\tilde{\Sigma}}^{-1} - \boldsymbol{\hat{S}}).
\]

The second step is a sparse Gaussian graphical model problem of the sample covariance matrix conditioned on the low rank matrix. This can be solved using standard methods such as GLASSO, CLIME, among others

\[
{\boldsymbol{\hat{S}}}^{(i+1)}={\underset{\boldsymbol{\hat{S}}}{\operatorname {arg\,max}}}\ Q({\boldsymbol{\hat{S}}}\mid {\boldsymbol{\hat{S}}}^{(i)}).
\]

