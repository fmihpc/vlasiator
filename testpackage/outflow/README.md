# Outflow population regression

Build with `make ARCH=<CPU architecture> vlasiator outflow-test` and run
`bash testpackage/outflow/run.sh`. Use the same precision/build flags for all
production and test objects; clean first when switching flags.

Six configurations run with one and two MPI ranks, using the production parser
and Outflow implementation. Each rank has a small DCCRG grid on MPI_COMM_SELF;
the parser and failure reduction use MPI_COMM_WORLD. In the two-rank schemes
case an additional grid is partitioned on MPI_COMM_WORLD. It transfers block
lists and VDFs through the production communication path, checks received
values and requires genuinely remote source cells on each rank.

The behavioral checks call vlasovBoundaryCondition on faces, edges and corners,
using a seeded closest-interior-cell cache. They compare every value in a
velocity block and population density before and after Copy, None and omitted
faces, and verify that the other population and source cell are unchanged.
Opposite schemes on all six faces exercise both populations independently.
Both L1 and L2 and R/V density moments are checked. L2 must copy moments without
copying the VDF. In a corner, Copy on one configured face is not vetoed by None
on another. The production restart helper is also checked for L1 VDF and L2
R/V moment updates; this is not a restart-file write/read round trip.

The distributed grid tests inter-rank exchange, not necessarily inter-node
transport. Full field/acceleration physics and scientific convergence remain
outside this regression's scope.

The solver fixtures run a short translation-only case and require step 2 without
bailout. The incomplete fixture must fail during configuration validation with
population and face named. Periodic faces count as covered. Face intersections
use the same stable precedence order as classification: higher precedence wins,
and later registration wins ties. Only the effective boundary supplies coverage;
Maxwellian applies to all populations, while Outflow requires an explicitly
configured face on that intersection for each population.
An explicitly configured None scheme counts as an explicit boundary choice.
This is a new admission policy requiring maintainer review, not just a storage fix.

The mixed solver tests cover both precedence orders, both registration orders at
equal precedence, explicit None, and overlapping x/z faces with corners. Valid
cases must reach step 2 without bailout; shadowed coverage must be rejected with
the missing population, face and effective boundary identified. Maxwellian input
data are population-specific and constant in time. These are short runtime
regressions, not physical-convergence tests. Inner-boundary geometry remains
outside the outer-face admission contract.

Other INI fixtures are partial parser/unit-test inputs. Different/empty face
fixtures deliberately test isolated Outflow parsing; the solver's stricter
cross-boundary coverage validation rejects such incomplete overall configs.
