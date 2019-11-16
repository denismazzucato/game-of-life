#report

## gol-seq
- *O(N ⋅ M ⋅ T)*

## gol-par ~80%
- row data partitioning
- MPI buffered mode
  to avoid unexpected deadlock
- communication over edge rows with neighbors
- *O((N/P) ⋅ M ⋅ T)*

## gol-par-nbl ~100%
- row data partitioning
- MPI non-blocking send and receive
  to overlap communication computation
- communication over edge rows with neighbors
- communication *O(T ⋅ 2M)*
- *O((N/P) ⋅ M ⋅ T)*

## gol-par-opt ~100%
this case has been done before checking the speedup of gol-par-nbl
- reduce communication overhead exploiting matrix structure
- row data partitioning
- MPI non-blocking send and receive
  to overlap communication computation
- communication over edge rows with neighbors
- communication *O(T ⋅ 2(M/31))*
- computational *O((N/P) ⋅ M ⋅ T ⋅ 2M)*