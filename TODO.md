# RoadMap to adaptive integrals
- [ ] Apply barycentric coords to set of points with matmul.
- [ ] Collect baryncentric coords in 1 array (allocate or use buffer as argument).
- [ ] Apply integrals using points collection: rewrite existing integrals and add interface of precalced points as aguments.
- [ ] Adaptive integration interface with fixed stop critetion.
## 2d integration
- [ ] When function $F(x,y)$ evals over matrix of point pairs, collect integral with $y^T A x$.
- [ ] Adaptive integration interface over 2 domains.

# Additional improvements
- [ ] Custom stop criterion: ideas about interface.
- [ ] Interface which accepts Funcs which can eval vector(matrix) of arguments more efficentrly.
