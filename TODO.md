# RoadMap to adaptive integrals

- [x] Apply barycentric coords to set of points with matmul.
- [x] Collect baryncentric coords in 1 array (allocate or use buffer as argument).
- [ ] If domain has func which transfer set of points then use it instead of naive cycle (first point)
- [ ] Apply integrals using points collection: rewrite existing integrals and add interface of precalced points as aguments.
- [ ] Adaptive integration interface with fixed stop critetion.

## 2d integration

- [ ] When function $F(x,y)$ evals over matrix of point pairs, collect integral with $y^T A x$.
- [ ] Adaptive integration interface over 2 domains.

# Additional improvements

- [ ] Custom stop criterion: ideas about interface.
- [ ] Interface which accepts Funcs which can eval vector(matrix) of arguments more efficentrly.
