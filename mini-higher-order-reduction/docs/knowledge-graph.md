# Knowledge Graph - mini-higher-order-reduction

## L1: Definitions
| Concept | Definition | Implementation |
|---------|-----------|----------------|
| State-space model | dx/dt = A x + B u, y = C x + D u | StateSpace struct |
| Transfer function | G(s) = num(s)/den(s) | TransferFunction struct |
| Reduced-order model | Approximate low-order system | ReducedModel struct |
| Pole | Eigenvalue s = sigma + j*omega | PoleInfo struct |
| Eigenvalue decomposition | A = V * Lambda * V^{-1} | EigenDecomp struct |
| Gramians | Controllability/observability Gramians | Gramians struct |
| Routh array | Stability analysis table | RouthArray struct |
| Model order | Dimension n of state vector | int n |
| Hankel singular values | sigma_i = sqrt(lambda_i(Wc*Wo)) | hsv[] array |
| Dominant pole | Pole with small |Re(s)| | is_dominant flag |

## L2: Core Concepts
- Dominant pole approximation
- Model order reduction (MOR)
- Balanced realization
- Singular perturbation theory
- Modal analysis and truncation
- Time-scale separation
- Pole-zero cancellation
- Moment matching
- Stability preservation in reduction
- H-infinity error bounds

## L3: Mathematical Structures
- Linear algebra: matrix multiply, transpose, inverse
- Eigenvalue decomposition via QR algorithm
- Lyapunov equation solver (Bartels-Stewart)
- Matrix exponential (scaling and squaring)
- LU decomposition with pivoting
- Cholesky decomposition
- Singular value decomposition
- Gramian computation
- Transfer function evaluation

## L4: Fundamental Laws
- Routh-Hurwitz stability criterion
- Lyapunov stability theorem
- Moore's balanced truncation theorem
- Glover's H-infinity error bound
- Pernebo-Silverman minimality condition
- Singular perturbation stability bounds
- Pade approximation theory
- Moment matching theory

## L5: Computational Methods
- Modal truncation (slowest, DC gain, HSV-based)
- Balanced truncation algorithm
- Singular perturbation reduction
- Davison method (DC gain correction)
- Marshall method (residue correction)
- Arnoldi/Krylov subspace methods
- Pade approximation
- Routh approximation method
- Proper Orthogonal Decomposition (POD)
- Frequency-weighted balanced truncation

## L6: Canonical Problems
- Boeing 747 lateral dynamics (4th order)
- DC motor with flexible shaft (4th order)
- Multi-room thermal system (8th order)
- Mass-spring-damper chain (variable order)
- Routh-Hurwitz stability analysis
- Dominant pole identification

## L7: Applications
- Boeing 747 flight control model reduction
- Servo motor control with flexible coupling
- Building thermal dynamics reduction
- Mechanical system model reduction

## L8: Advanced Topics
- Time-weighted balanced truncation
- Frequency-weighted balanced truncation
- Positive real balanced truncation
- Bounded real balanced truncation
- Second-order structure-preserving MOR
- Rational Krylov methods
- Multi-point Pade approximation

## L9: Research Frontiers
- Structure-preserving model reduction
- Data-driven MOR (POD)
- Error bounds for nonlinear MOR
- Machine learning for model reduction
