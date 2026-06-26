# Gap Report - mini-higher-order-reduction

## Missing Knowledge Points

### L7: Applications (Priority: Medium)
- More real-world industrial case studies
- Process control applications
- Power system model reduction

### L8: Advanced Topics (Priority: Low)
- Full frequency-weighted balanced truncation with weighting filters
- Nonlinear model reduction via empirical Gramians
- Parameter-varying model reduction

### L9: Research Frontiers (Priority: Low)
- Machine learning-based model reduction
- Structure-preserving MOR for port-Hamiltonian systems
- Data-driven MOR beyond POD

## Known Limitations
- Eigenvalue decomposition uses simplified QR algorithm
- Lyapunov solver falls back to Kronecker product for non-diagonalizable systems
- Singular perturbation requires balanced realization first
- Large-scale systems (>512 states) not supported by direct methods
