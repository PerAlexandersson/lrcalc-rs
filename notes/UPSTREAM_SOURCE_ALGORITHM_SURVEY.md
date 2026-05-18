# Upstream Source And Algorithm Survey

Scope:

- upstream reference: `/tmp/lrcalc-upstream` or
  `/workspace/references/lrcalc-upstream`
- upstream commit: `8705a16e1575351684ed692f88552f97da3724f4`
- Rust state checked: current `src/` implementation and benchmark scripts
- documentation updated after the Rust survey

## Feature Surface

The upstream project exposes three user-facing layers.

Command-line programs:

- `lrcalc coef`: one Littlewood-Richardson coefficient.
- `lrcalc skew`: Schur expansion of a skew Schur function.
- `lrcalc mult`: Schur product, with row/column bounds.
- `lrcalc mult -q`: quantum product for a Grassmannian.
- `lrcalc mult -f`: fusion-ring notation for the same quantum data.
- `lrcalc coprod`: Schur coproduct.
- `lrcalc tab`: iterate LR tableaux; current upstream also accepts a weight.
- `schubmult`: Schubert polynomial product, including string notation.

Library/API surface:

- C ABI in installed headers under `include/lrcalc`.
- Python Cython binding in `python/lrcalc.pyx` and `python/liblrcalc.pxd`.
- Maple wrapper generated from `lrcalc.maple.src`.

External consumers:

- Sage documents `lrcalc` as a C library for LR coefficients, Schur products,
  coproducts, skew Schur functions, fusion products, and Schubert products.
- Buch's web page describes the implementation as C code for efficiency and
  notes a major speedup for single LR coefficients due to Pierre-Emmanuel
  Chaput.

## Upstream Architecture

The C code is small and allocation-conscious.

- `ivector` is the central flexible-array integer vector:
  `uint32_t length; int32_t array[1];`.
- `ivlincomb` is a custom hash table from `ivector *` to `int32_t`.
- `ilist` and `ivlist` are simple growable arrays from template headers.
- `part_iter` enumerates partitions in a box, under an outer partition, above an
  inner partition, or with a fixed size.
- Most algorithms mutate temporary `ivector.length` fields for chopped
  partitions; ABI `iv_free` must not depend on the current length.
- Debug-memory code can force allocation failures and check leaks, but release
  builds use libc allocation macros.

Important source paths:

- `src/schur.c`: top-level Schur, fusion, coproduct, and coefficient functions.
- `src/lriter.c`: LR tableau iterator and expansion counter.
- `src/lrcoef.c`: optimized single-coefficient counter.
- `src/optshape.c`: shape reduction and preprocessing.
- `src/schublib.c`: Schubert transition and multiplication.
- `src/perm.c`: permutation/string conversions and Bruhat helpers.
- `src/part.h`: partition predicates and iterator.

## Main Algorithms

### Schur Product

Path:

- `schur_mult`
- `optim_mult`
- `lrit_expand`
- `lrit_new`
- `lrit_count`

The product `s_sh1 * s_sh2` is turned into an LR tableau enumeration problem.
`optim_mult` chooses which factor should become the skew shape and which should
become the content. It removes full rows and columns imposed by row/column
bounds, then puts the smaller remaining shape on the left side of the LR count.

The expansion is then computed by iterating LR tableaux and accumulating the
content vector in an `ivlincomb`.

### Skew Schur Expansion

Path:

- `schur_skew`
- `optim_skew`
- `lrit_expand`

`optim_skew` validates containment, removes full-height columns, splits
independent components, and chooses the largest partition-shaped or
anti-partition-shaped component as the tableau content. Remaining components
are folded into a smaller skew shape.

This preprocessing is one of the main reasons the general expansion routines
are fast despite using tableau enumeration.

### Single LR Coefficient

Path:

- `schur_lrcoef`
- `optim_coef`
- `lrcoef_count`

`optim_coef` first handles cheap zeros and ones:

- invalid size or containment gives zero;
- empty residual shape gives one;
- horizontal and vertical compactification shrink the triple;
- removable rows and columns are stripped from alternating pairs of shapes.

Only if a nontrivial residual remains does `schur_lrcoef` call
`lrcoef_count`.

`lrcoef_count` builds a compact skew-tableau array. Each box stores:

- its current value;
- a maximum allowed value;
- pointers to the north and east constraint boxes;
- `se_supply`, `se_sz`, and `west_sz` pruning data.

The search fills boxes in a linear order and backtracks in place. The key
pruning test is the Chaput-style southeast supply bound:

- if the remaining larger letters cannot fill the required southeast region,
  the branch is skipped.

The upstream changelog says version 2.1 reimplemented single LR coefficients
and became more than twice as fast as version 2.0, mostly due to this pruning
idea.

### LR Tableau Iterator

Path:

- `lrit_new`
- `lrit_next`
- `lrit_count`

`lrit_new` creates a tableau array with sentinel boxes. For each real box it
precomputes:

- `above`: the box constraining column strictness;
- `right`: the box constraining row weak increase;
- `max`: the largest possible label after row/column and shape bounds.

It initializes the minimal LR tableau, stores the current content vector, and
`lrit_next` advances in place by increasing the first box that can be raised,
then refilling earlier boxes minimally.

This iterator is used for full expansions; `lrcoef_count` is separate because
one coefficient benefits from stronger pruning.

### Fusion And Quantum Products

Path:

- `schur_mult_fusion`
- `fusion_reduce`
- `optim_fusion`
- `lrit_expand`
- `fusion_reduce_lc`
- `part_qdegree`, `part_qentry`

The code reduces inputs by an affine sorting operation on shifted entries. A
zero sign kills terms with repeated reduced entries; otherwise the sign is
tracked. `optim_fusion` applies a Seidel-shift style rotation to choose a small
left-hand partition before ordinary LR enumeration.

The quantum command prints the same fusion output in quantum cohomology
notation using `part_qdegree` and `part_qentry`.

### Coproduct

Path:

- `schur_coprod`
- `schur_mult` with a rectangle
- `_schur_coprod_count`

The coproduct is reduced to multiplying by a rectangular partition. The result
is encoded as one `ivector` containing both partitions, then filtered so only
one representative of each symmetric pair is printed unless `all` is requested.

### Schubert Products

Path:

- `trans`
- `_trans`
- `monk`
- `_monk_add`
- `mult_poly_schubert`
- `_mult_ps`
- `mult_schubert`
- `mult_schubert_str`

`trans` recursively computes a Schubert polynomial in monomial/exponent-vector
form using a transition formula. `mult_poly_schubert` multiplies that
polynomial by a Schubert class by recursively stripping variables and applying
Monk's rule.

`mult_schubert` chooses the shorter permutation length first and uses a Bruhat
zero test when a finite rank is supplied. String notation is converted through
dimension vectors and permutations.

## Upstream Test Strategy

The upstream tests are worth copying conceptually.

- `test_partiter`: verifies partition enumeration counts and constraints.
- `test_lrmult`: compares Schur multiplication against Schubert string
  multiplication for Grassmannian permutations.
- `test_lrcoef`: compares each single coefficient against a full Schur product.
- `test_lrskew`: compares skew expansion against single coefficients.
- `test_fusion`: compares fusion multiplication against Schur multiplication
  followed by fusion reduction.
- `test_schubmult`: checks commutativity and rank truncation of Schubert
  products.
- `testsuite`: checks exact CLI output for representative examples.

For Rust, these are better than isolated sample values because they test
internal consistency among independent algorithms.

## Algorithm Context

Exact LR counting is not expected to have a general polynomial-time algorithm:
Narayanan proved that computing LR coefficients is #P-complete.

For deciding nonzero coefficients, hive and honeycomb models are more powerful.
Knutson-Tao hives are equivalent to the LR rule; Buch's exposition of the
saturation theorem records this equivalence. Burgisser-Ikenmeyer use the hive
model to decide positivity in polynomial time via flow methods.

For small exact coefficients, Ikenmeyer gives output-sensitive algorithms with
runtime polynomial in the input and quadratic in the coefficient. This suggests
that a Rust implementation could choose different methods depending on whether
it is proving zero/nonzero, finding a small coefficient, or expanding a large
product.

Narayanan's approximation work is relevant for large "deep" coefficients, but
it is not a drop-in replacement for exact `liblrcalc` semantics.

References:

- Buch lrcalc page: https://sites.math.rutgers.edu/~asbuch/lrcalc/
- Sage lrcalc interface:
  https://doc.sagemath.org/html/en/reference/libs/sage/libs/lrcalc/lrcalc.html
- Narayanan, #P-completeness: https://arxiv.org/abs/math/0501176
- Buch on Knutson-Tao saturation and hives:
  https://arxiv.org/abs/math/9810180
- Burgisser-Ikenmeyer positivity:
  https://arxiv.org/abs/1204.2484
- Ikenmeyer small coefficients: https://arxiv.org/abs/1209.1521
- Narayanan approximation: https://arxiv.org/abs/1306.4060
- Molev-Sagan factorial Schur LR rule:
  https://arxiv.org/abs/q-alg/9707028
- Molev Littlewood-Richardson polynomials:
  https://arxiv.org/abs/0704.0065

## LR Computation Techniques

No single method should be expected to dominate. Exact LR counting is
#P-complete, while nonvanishing is much easier. The Rust implementation should
make the coefficient engine swappable after a correct baseline exists.

### Direct LR Tableau Search

This is upstream `lrcalc`'s core approach.

- Counts skew semistandard tableaux with Yamanouchi/lattice-word constraints.
- Good fit for C ABI parity because it matches current semantics directly.
- Upstream adds strong shape reductions and Chaput-style southeast supply
  pruning for single coefficients.
- First Rust target: implement a clear version, then port the upstream pruning.

### LR Tableau Iterator For Expansions

Upstream keeps full expansion separate from one-coefficient counting.

- Iterates all LR tableaux and accumulates their contents in a linear
  combination.
- Good for `schur_mult` and `schur_skew`, where many coefficients are needed.
- Less ideal for a single sparse coefficient because the target content is
  known in advance.

### Horizontal-Strip / GT Dynamic Programming

This is the closest local Rust reference in `/workspace/rust/kostka/src/lr.rs`.

- Builds the skew shape by horizontal strips of sizes given by the content.
- Integrates Yamanouchi constraints into the DP state.
- Good first independent Rust oracle because it is easier to audit than the C
  backtracker.
- Potential risk: state explosion on shapes where the C pruning is very strong.

### GT-Polytope Single-Coefficient Route

The local `/workspace/rust/kostka` crate also has GT-polytope infrastructure:

- `/workspace/rust/kostka/src/gt_dim.rs` computes feasibility, propagated
  bounds, and dimension for `GT(lambda/mu, w)`.
- `/workspace/rust/kostka/GT_DIM_ALGORITHM.md` documents the chain model and
  the interval-propagation algorithm.
- `/workspace/rust/kostka/src/lr.rs` already adds Yamanouchi constraints for LR
  coefficients by augmenting the horizontal-strip DP state.

This suggests a serious single-coefficient strategy:

1. Model `c^lambda_{mu,nu}` as integer points in a GT chain polytope with
   additional Yamanouchi inequalities.
2. Run fast bound propagation before enumeration.
3. Return zero immediately when the constrained polytope is empty.
4. Return one when propagation proves a unique integral chain.
5. Otherwise enumerate with GT/horizontal-strip DP, using propagated bounds and
   Yamanouchi constraints.

This is a better fit for polytope tricks than Buch's flat tableau array. It can
use dimension, forced entries, Ehrhart data, and possibly chamber/fixed-rank
methods. The main caution is that the existing `gt_dim` dimension is for skew
Kostka polytopes before the full LR/Yamanouchi constraints; it needs either an
LR-constrained extension or a hive-equivalent feasibility/dimension layer before
it can certify `c = 0` or `c = 1`.

### Stretching And Ehrhart Interpolation

Stretching gives another route for large structured inputs.

For fixed GT/Kostka data, the local `kostka` crate already uses that

```text
n |-> K(n*lambda / n*mu, n*w)
```

is an Ehrhart polynomial. It computes the dimension first, samples enough
values, and interpolates over `Q`. It also uses Ehrhart-Macdonald reciprocity:
negative values of the polynomial are obtained from strict/interior GT counts,
which often vanish for small positive dilation and reduce the sampling cost.

The analogous LR plan would be:

1. Represent `c^{n*lambda}_{n*mu,n*nu}` or the affine family
   `c^{n*lambda+alpha}_{n*mu+beta,n*nu+gamma}` as lattice points in a
   constrained GT or hive polytope.
2. Determine the dimension/degree after the active chamber stabilizes.
3. Evaluate the exact coefficient at `d+1` values of `n` where polynomiality is
   known to hold.
4. Interpolate the polynomial and evaluate at the target large `n`.
5. Verify by checking one or two extra sample values outside the interpolation
   set.

This is attractive when many coefficients lie in the same stretched family, or
when the target dilation is much larger than the interpolation degree. It is
less attractive for one-off small ABI calls, where direct counting is cheaper.

The main technical questions:

- For pure dilation, LR stretching is an Ehrhart polynomial of the LR/hive
  polytope.
- For affine stretches `n*lambda+alpha`, eventual polynomiality depends on
  staying in one chamber of the relevant vector-partition/hive fan.
- A better framework for offsets is a two-parameter family
  `c^{n*lambda+m*alpha}_{n*mu+m*beta,n*nu+m*gamma}`. Rassart's chamber result
  says LR coefficients are polynomial functions of the boundary data on each
  cone of a vector-partition chamber complex. Restricting those chamber
  polynomials to the `(m,n)` plane gives a piecewise bivariate polynomial.
- The one-parameter affine stretch `n*lambda+alpha` is the slice `m = 1`.
  Eventual polynomiality is the statement that this ray eventually remains in
  one chamber.
- The whole quadrant in `(m,n)` has a single polynomial only when its image is
  contained in one chamber; otherwise it is subdivided by finitely many chamber
  rays or affine walls.
- The polynomiality threshold `n0` must be certified, not guessed, if this is
  used for production results.
- Ehrhart-Macdonald reciprocity should apply cleanly for pure dilations once
  the correct relative interior of the LR-constrained polytope is implemented.
- For the two-parameter/chamber-wise setting, the reciprocity question should
  be phrased as vector-dilated or parametric Ehrhart reciprocity. The exact
  strict/interior model has to be derived chamber by chamber before relying on
  negative interpolation points.

For Rust, this belongs in a native experimental path rather than the first C ABI
implementation. It can reuse `/workspace/rust/kostka/src/ehrhart.rs`,
`/workspace/rust/kostka/src/gt_dim.rs`, and strict GT counting ideas after the
LR/Yamanouchi constraints are incorporated into the polytope model.

### Near-Stretched Coefficients

A more practical version is to exploit triples close to a pure stretch:

```text
c^{N*lambda + alpha}_{N*mu + beta, N*nu + gamma}
```

with fixed small offsets and large `N`.

What can be exploited:

- If `N*(lambda,mu,nu)` lies in the interior of one Rassart chamber and the
  offset is bounded, then all sufficiently large `N` stay in the same chamber.
  The coefficient is then an exact one-variable polynomial in `N`.
- The pure stretched polynomial gives the `alpha=beta=gamma=0` slice of the
  chamber polynomial. It gives leading/asymptotic information, but it does not
  determine nearby offset slices by itself.
- If many offsets around the same ray are queried, sample enough nearby
  triples to reconstruct the local chamber polynomial in the offset variables
  and `N`; later queries become polynomial evaluation.
- If only one nearby coefficient is queried, use the stretched data as a
  prepass: estimate degree/size, predict likely active constraints, warm-start
  GT/hive bound propagation, and test whether the offset has crossed a chamber
  wall.

Safe exact workflow:

1. Compute or cache the pure stretched polynomial for the ray.
2. Use chamber inequalities, hive active sets, or finite-difference checks to
   test whether the offset family stays in the same chamber for the target
   range.
3. If stable, interpolate the offset slice
   `N -> c^{N*lambda+alpha}_{N*mu+beta,N*nu+gamma}` from exact small-`N`
   evaluations and validate with extra points.
4. If unstable or near a wall, fall back to direct GT/hive counting.

This is most useful for large `N`, repeated calls near one ray, or benchmark
families. It is not a replacement for the basic ABI counter because the chamber
certificate is the hard part.

### Kostka Matrix Inversion

Also present in `/workspace/rust/kostka/src/lr.rs`.

- Uses `K(lambda/mu, alpha) = sum_nu c^lambda_{mu,nu} K(nu, alpha)`.
- Solves by back-substitution in dominance order.
- Useful for testing and small examples.
- Not a likely production engine for `liblrcalc` because it computes many
  auxiliary Kostka numbers.

### Signed Kostka Expansion For LR

Shrivastava gives a signed expression of LR coefficients in terms of Kostka
numbers. In one convenient form, for fixed rank with
`rho = (ell-1, ell-2, ..., 0)`, terms are indexed by permutations and have
weights of the form

```text
weight_tau = nu + rho - tau(lambda + rho).
```

After discarding invalid weights, the coefficient is a signed sum of Kostka
numbers `K_{mu, weight_tau}`. This is different from the SymCat
`kostkaFromLR` anchor, which records the reverse fact that skew Kostka numbers
are special LR coefficients.

This route may be practical when the second lower partition `mu` is fixed and
the rank `ell` is small:

- Normalize first: use LR symmetry to swap the two lower partitions and use
  conjugation when helpful. Choose the Kostka shape side to be small in
  dominance order, so the dominance test rejects more candidate weights.
- Precompute all needed `K_{mu, alpha}` once.
- Kostka numbers are invariant under permutation of the weight, so canonicalize
  each `weight_tau` by sorting it before lookup.
- Accumulate the signed multiplicity of each canonical weight first; many
  permutation terms can cancel before any Kostka computation.
- If many LR calls share the same `mu`, the Kostka cache is reused across all
  calls.
- A recursive permutation generator can prune early using nonnegativity, total
  weight, and dominance `sort(weight_tau) <=_dom mu`.

Expected tradeoff:

- Good when `ell!` is small, `mu` is reused, and Kostka lookup is very fast.
- Bad when rank is large or cancellation requires many large Kostka values.
- Useful as a benchmark/oracle against Buch-style counting and GT/Yamanouchi
  DP.

### LR Symmetry Normalization

Before choosing a counting engine, Rust should generate a small set of
equivalent triples and pick the cheapest one for that engine.

Safe general symmetries:

- Swap the lower partitions:
  `c^lambda_{mu,nu} = c^lambda_{nu,mu}`.
- Conjugate all partitions:
  `c^lambda_{mu,nu} = c^{lambda'}_{mu',nu'}`.
- Determinant/rectangle translation in fixed rank `ell`: for integer shifts
  `a,b`, tensoring by determinant gives
  `c^lambda_{mu,nu} =
   c^{lambda-(a+b)^ell}_{mu-a^ell,nu-b^ell}`
  whenever all shifted weights remain partitions.
- Grassmann/Poincare complement symmetry: if all shapes fit in a rectangle `R`
  and `theta^vee` denotes complement in `R`, then
  `c^lambda_{mu,nu} = c^{nu^vee}_{mu,lambda^vee}
                   = c^{mu^vee}_{nu,lambda^vee}`.

Algorithmic use:

- For Buch-style counting, choose the variant with the smallest residual skew
  shape after row/column stripping.
- For GT/Yamanouchi DP, choose the variant with fewer free GT-chain entries and
  lower estimated peak states.
- For the signed Kostka expansion, choose the variant where the fixed Kostka
  shape is smallest in dominance order and where the permutation filter leaves
  the fewest canonical weights.
- Cache the chosen normalization with the original triple so ABI output remains
  unchanged.

There are more specialized or hidden symmetries, especially in small rank and
hive formulations. Treat those as research candidates until each is expressed
as an explicit, tested transformation on partition triples.

### Hive / Horn / Flow Methods

Hive inequalities are the right tool for nonzero checks.

- Positivity of LR coefficients can be decided in polynomial time.
- A hive feasibility prepass can cheaply reject many impossible triples.
- Counting hives is another exact coefficient model, but integer-point counting
  is not automatically faster than tableaux.
- Best Rust role: optional zero filter and benchmark corpus generator.

### BZ Patterns And Vector Partition Functions

Berenstein-Zelevinsky triangles and Steinberg/Kostant formulas are alternative
integer-point models.

- Good for theoretical cross-checking and fixed-rank experiments.
- Steinberg's formula expresses `c^nu_{lambda,mu}` as a signed double sum over
  `S_k x S_k` of Kostant partition-function values.
- The `k!^2` outer sum grows very quickly with rank.
- Signed Kostant sums can suffer heavy cancellation.
- It may be useful when `k` is tiny, when many summands vanish cheaply, or when
  a specialized/vectorized type-A Kostant partition function is already
  available.
- Less attractive as the first ABI engine unless a local implementation already
  exists.

### Puzzle / Equivariant / Factorial Models

Knutson-Tao puzzles and factorial Schur rules are natural for Grassmannian and
equivariant variants.

- They connect well to quantum/fusion and factorial-Schur generalizations.
- They are not the narrowest route to ordinary `schur_lrcoef`.
- Worth tracking for later compatibility with Schubert and equivariant code.

### Molev-Sagan / Shifted-Schur Recursion

The shifted-Schur recursion is worth studying as a separate engine.

- It computes a larger triangular family of coefficients.
- Its top homogeneous slice recovers ordinary LR coefficients.
- It gives a compact recursive oracle if shifted Schur evaluations are
  implemented exactly.
- It should be benchmarked against direct LR and horizontal-strip DP before it
  becomes part of the ABI path.

### Reductions, Symmetries, And Special Rules

These are low-risk improvements around any core counter.

- Swap the two lower partitions.
- Conjugate all partitions when height/width improves.
- Strip forced rows and columns, as upstream already does.
- Use Pieri and empty-shape cases before invoking a general engine.
- Add known reduction formulas only when tests show they reduce real workloads.

### Recommended Rust Evaluation Order

1. Correct direct LR or horizontal-strip DP for `schur_lrcoef`.
2. Upstream shape reductions and single-coefficient pruning.
3. Independent oracle: local Kostka inverse and/or shifted-Schur recursion.
4. Hive positivity prepass.
5. Expansion iterator and product/skew accumulation.
6. Later experimental engines: BZ, puzzles, factorial/shifted families.

## Shifted-Schur Comparison

The local `symmetricfunctions.com` source has the Molev-Sagan style recursion
for shifted Schur structure constants. With

```text
s*_lambda s*_mu = sum_nu c^nu_{lambda,mu} s*_nu,
```

and `mu, nu <= lambda`, it records

```text
c^lambda_{mu,nu}
  = (sum_{nu -> nu+} c^lambda_{mu,nu+}
     - sum_{lambda- -> lambda} c^{lambda-}_{mu,nu})
    / (|lambda| - |nu|).
```

The boundary value is

```text
c^lambda_{mu,lambda} = s*_mu(lambda),
```

computed from the shifted Jacobi-Trudi identity or an equivalent tableau sum.
These shifted coefficients are a strict extension of classical LR data: the
top homogeneous case recovers ordinary LR coefficients when
`|lambda| = |mu| + |nu|`.

Implementation implications:

- This gives a compact recursive specification and a useful independent oracle
  for small and medium LR coefficients.
- The recursion is triangular over pairs `(lambda', nu')`; memoization is
  essential.
- The base case requires exact evaluation of `s*_mu(lambda)`, so determinant or
  tableau evaluation needs integer arithmetic and overflow policy.
- It should not be the first ABI engine without benchmarks: ordinary LR has
  stronger direct tableau, horizontal-strip, and upstream `lrcoef_count`
  pruning methods.
- It is a good comparison target for the Rust roadmap because it computes a
  larger family whose top-degree slice is the classical coefficient required by
  `lrcalc`.

Local evidence:

- `/home/paxinum/Dropbox/webpages/symmetricfunctions.com/tex-source/schurShifted.tex`
  contains the shifted Schur definition, shifted Jacobi-Trudi identities, and
  this recursion.
- `/home/paxinum/Dropbox/webpages/symmetricfunctions.com/tex-source/littlewoodRichardson.tex`
  points to shifted Schur functions as a source of a simple LR recursion and to
  the Fomin-Greene/Molev-Sagan tableau model.
- `/workspace/rust/kostka/src/lr.rs` does not use this recursion. It implements
  an augmented GT/horizontal-strip DP with Yamanouchi constraints plus a Kostka
  matrix inverse check.

One source note needs cleanup before using it verbatim. The size relation
`|nu| = |lambda| + |mu|` is the right top-degree relation for the earlier
product coefficient `c^nu_{lambda,mu}`. Once the theorem reindexes the
coefficient as `c^lambda_{mu,nu}`, the ordinary LR extraction should read
`|lambda| = |mu| + |nu|`. The current webpage sentence is therefore reversed
under the notation used immediately above it.

## Current Rust State

The Rust crate now has several substantial native engines.

- `src/abi.rs`: C-facing ABI surface for vectors, lists, linear combinations,
  partition helpers, LR-tableau iterators, Schur/fusion operations, Schubert
  operations, Maple printing helpers, and low-level coefficient functions.
- `src/lrcoef.rs`: primary Buch-style single LR coefficient implementation,
  including upstream-style compactification, branch-pruned tableau counting,
  relative-interior counts, dimension, and stretch-cache helpers.
- `src/lr_gt.rs`: independent GT-chain/Yamanouchi DP for LR coefficients,
  stats, dimension, and relative interiors.
- `src/kostka_fast.rs`: packed `u128` ordinary and skew Kostka DP, plus
  interior counts.
- `src/lr_signed.rs`: signed Kostka expansion for LR coefficients.
- `src/lr_ehrhart.rs`: pure stretched LR h-vector interpolation using full and
  relative-interior sample points.
- `src/main.rs`: CLI entry points for coefficient, Kostka, GT, signed-Kostka,
  Buch-interior, and stretch experiments.
- `src/bin/` and `scripts/`: benchmark harnesses for in-process and upstream-C
  comparisons.

What remains for release hardening:

- broader C smoke tests for iterator, ownership-transfer, and printing edge
  cases;
- a full Sage rebuild recipe, beyond the current `LD_PRELOAD` smoke test;
- a final decision on whether CLI term order must exactly follow upstream
  hash-table iteration.

The concrete ABI issue from the survey has been fixed in the staged headers:
`ivlc_iter.index` and `ivlc_iter.i` now follow the upstream C header and use
`size_t`.  Python/Sage smoke tests currently pass through the staged Rust
prefix.

## Room For Rust Improvements

### Benchmark Engine Selection

The project now has enough independent engines to benchmark before adding new
algorithms.  The first comparison matrix should include:

- Buch-port `lrcoef`;
- GT-chain `lrcoef_gt_u128`;
- signed-Kostka `lrcoef_signed_kostka`;
- stretched/Ehrhart paths when the input is a pure dilation;
- upstream C `lrcalc coef`.

Use the same triples for correctness and timing.  Track value, time, peak state
counts, memo states, cache hits, and which compactification path was used.

### Keep The ABI Boundary Thin

Use C-compatible structs only at the ABI edge. Convert immediately to native
Rust data for computation, then convert back to `ivector` and `ivlincomb` for
return values.

Do not let native Rust data structures leak into installed C headers unless the
project intentionally drops full struct-layout compatibility for C callers.

### Use Wider Internal Counters

The current implementation already uses checked `u128` for nonnegative native
counts and downcasts at ABI boundaries. Keep this policy:

- `u128` for nonnegative LR counts;
- `i128` for signed linear-combination accumulation;
- explicit downcasts at ABI boundaries;
- compatibility mode returning upstream-style failure when the result does not
  fit the required C type.

Optional `BigInt` support belongs in native tools first, not in the ABI path.

### Add Symmetry And Zero Prepasses

Good low-risk improvements around the existing engines:

- choose between the two lower partitions before Buch or GT counting;
- conjugate triples when height/width improves the selected engine;
- add determinant/rectangle translations only with explicit tests;
- add a hive/Horn feasibility prepass for hard zero cases;
- cache the selected normalization for repeated benchmark families.

Each rewrite must be tested against upstream examples before being trusted.

### Reuse The Expansion Split

Upstream separates full expansions from one coefficient.  Rust should keep that
shape:

- one fast single-coefficient engine for `schur_lrcoef`;
- one `ivlincomb`-producing iterator/accumulator path for `schur_skew`,
  `schur_mult`, and `schur_coprod`;
- coefficient-by-coefficient fallbacks only when benchmarks show sparse support
  makes them faster.

### Reuse Ideas, Not Dependencies, From Local Rust Crates

Local candidates:

- `/workspace/rust/kostka/src/lr.rs`: LR DP reference algorithm.
- `/workspace/rust/kostka/src/kostka_dp.rs`: horizontal-strip DP.
- `/workspace/rust/combinatoric-core/src/partition.rs`: rich partition helpers.
- `/workspace/rust/sym-poly/core/src/tableau.rs`: tableau types and iterators.
- `/workspace/rust/sym-poly/sym/src/symmetric_function.rs`: skew Schur
  Jacobi-Trudi reference implementation.

Avoid direct path dependencies for now because `lrcalc-rs` should remain a
standalone ABI replacement. Copy small ideas with attribution or reimplement
locally.

Important caveat: the current `/workspace/rust/kostka` code is optimized around
stretched Kostka/Ehrhart experiments and can use substantial memory for dynamic
programming maps. Treat it as:

- a correctness oracle;
- a source of GT-polytope dimension/bounds ideas;
- a source of interpolation/reciprocity code;
- not automatically a production engine for one-off LR coefficients.

The signed Kostka expansion needs a separate memory benchmark: even if many
weights canonicalize or cancel, each surviving `K_{mu,alpha}` call may still be
expensive with the current DP.

However, DP may still be the right backend if the signed expansion is evaluated
as a batch for one fixed Kostka shape. Instead of calling
`K_{mu,alpha}` independently for many `alpha`, compute all requested weights
together:

- canonicalize and collect the requested weights first;
- share horizontal-strip transition data for the fixed shape `mu`;
- reuse intermediate DP states across weights with common prefixes or sorted
  weight multisets;
- stop once all requested weights are resolved;
- track peak state count and memory per requested weight, not just wall time.

This turns the signed expansion from many independent Kostka calls into a
multi-target fixed-shape DP problem. It is worth benchmarking separately.

### Preserve Upstream's Internal Consistency Tests

Port the relations, not just examples:

- `schur_lrcoef` equals lookup in `schur_mult`.
- `schur_skew` coefficients equal `schur_lrcoef`.
- fusion equals ordinary Schur multiplication followed by fusion reduction.
- Grassmannian Schur products match Schubert string products.

These tests will catch mistakes that fixed golden examples miss.

### Use Rust To Improve Safety Around Known C Hazards

Rust can improve:

- checked conversion between `int32_t`, `usize`, and `long long`;
- overflow reporting for coefficients above `i64::MAX`;
- explicit ownership for hash-table keys;
- iterator invalidation rules;
- no mutation of allocation-critical length metadata;
- deterministic cleanup on allocation failure.

The ABI still needs C semantics: null on allocation failure for pointer-returning
functions and `-1` for single-coefficient allocation/overflow failure.

## Suggested Milestones

1. Benchmark current coefficient engines against upstream C on a documented
   small/mixed/large corpus.
2. Implement `ivlincomb` and iterator ownership, then add C smoke tests for
   layout and lifetime behavior.
3. Implement `schur_skew` using an expansion accumulator, because it is the
   narrowest route to real `ivlincomb` output.
4. Implement `schur_mult` and `schur_coprod`, keeping output-order tests
   separate from coefficient-value tests.
5. Add `part_qdegree`, `part_qentry`, `fusion_reduce`, and
   `schur_mult_fusion`.
6. Add LR tableau iterator ABI.
7. Add Schubert transitions, Monk multiplication, string notation, and
   `schubmult` compatibility.
8. Rebuild upstream Python bindings against Rust `liblrcalc` and add Sage-style
   smoke tests.

The immediate performance target is to choose when each already-implemented
coefficient engine wins.  The immediate compatibility target is `ivlincomb`,
because all Schur expansion ABI functions depend on it.
