//! Certified exact affine hulls of rational polyhedra.
//!
//! For a polyhedron
//!
//! ```text
//! P = { x in Q^n : A x = b, C x <= d },
//! ```
//!
//! an inequality `C_i x <= d_i` is an *implicit equality* if it holds with
//! equality at every point of `P`.  The affine hull of a nonempty `P` is
//! `{ A x = b, C_I x = d_I }`, where `I` is the set of implicit equalities, and
//! its dimension is `n - rank [A; C_I]`.  A lattice point of `N P` lies in the
//! relative interior exactly when every inequality outside `I` is strict.
//!
//! # Method
//!
//! We use the single homogenized linear program of Freund, Roundy and Todd
//! (1985):
//!
//! ```text
//! maximize   sum_i t_i
//! subject to A x = b tau,   C x + t <= d tau,   0 <= t_i <= 1,   tau >= 1.
//! ```
//!
//! If `P` is empty, then so is this program, since `x / tau` would lie in `P`.
//! Otherwise choose `z` in the relative interior of `P`, let `s > 0` be its
//! least slack on a non-implicit inequality, and put `tau = max(1, 1/s)`,
//! `x = tau z`.  Then every non-implicit slack of `(x, tau)` is at least one,
//! so the optimum equals the number of non-implicit inequalities.  Each
//! implicit `t_i` is zero on the feasible region, hence every optimum has
//! `t_i = 1` exactly on the non-implicit inequalities, and `x / tau` is a
//! relative-interior point.  No sampling or denominator bound is involved.
//!
//! The program is solved by an exact two-phase primal simplex method over
//! `BigRational`, using Dantzig pricing and switching to Bland's rule during
//! long degenerate runs; Bland's rule cannot cycle.
//!
//! # Certificates
//!
//! The simplex result is not trusted directly.  Every answer is checked
//! against the original data, and a failed check panics as an internal bug:
//!
//! * Nonempty: the primal point `p = x / tau` satisfies every constraint
//!   exactly.  The inequalities with positive slack at `p` are therefore not
//!   implicit.  The optimal duals give `y >= 0` and `w` with
//!   `C^T y + A^T w = 0`, `y . d + w . b <= 0`, and `y_i > 0` on every
//!   inequality that is tight at `p`.  On `P`,
//!   `sum_i y_i (d_i - C_i x) = y . d + w . b <= 0` is a nonnegative
//!   combination of slacks, so every inequality with `y_i > 0` is implicit.
//! * Empty: the phase-one duals give a Farkas certificate `y >= 0`, `w` with
//!   `C^T y + A^T w = 0` and `y . d + w . b < 0`.
//!
//! # Cost
//!
//! The tableau has `p + 2m` rows and `2n + 3m + 1` structural columns for `n`
//! variables, `p` equalities and `m` inequalities.  Rows are sparse, and pivot
//! updates touch only rows that are nonzero in the pivot column.  This is much
//! more expensive than interval propagation, so callers should first contract
//! variables that cheap sound propagation has already identified, and should
//! not recompute an affine hull merely to prune a counting recursion.
//!
//! Related Rust: this is a standalone copy of Ehrcalc's
//! `crates/ehrcalc-kostka-engine/src/affine_hull.rs` (repair branch
//! `fix/correctness-audit-2026-10-08`, 2026-10-08), because lrcalc-rs and
//! Ehrcalc have no shared dependency.  Keep the two copies synchronized.

use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::{One, Signed, ToPrimitive, Zero};

/// A sparse linear form `sum coefficient * x_variable` with its right side.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct LinearConstraint {
    pub coefficients: Vec<(usize, BigRational)>,
    pub rhs: BigRational,
}

/// The polyhedron `{ x : A x = b, C x <= d }` with free rational variables.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct RationalPolyhedron {
    variables: usize,
    equalities: Vec<LinearConstraint>,
    inequalities: Vec<LinearConstraint>,
}

/// Exact affine-hull data for a nonempty polyhedron.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AffineHullData {
    /// Dimension of the affine hull.
    pub dimension: usize,
    /// `implicit_equalities[i]` is true when inequality `i` is tight on the
    /// whole polyhedron.
    pub implicit_equalities: Vec<bool>,
    /// A rational point that is strict on every non-implicit inequality.
    pub relative_interior_point: Vec<BigRational>,
}

/// Result of [`RationalPolyhedron::affine_hull`].
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum AffineHull {
    Empty,
    Nonempty(AffineHullData),
}

impl RationalPolyhedron {
    /// Start a system in `variables` free rational variables.
    pub fn new(variables: usize) -> Self {
        Self {
            variables,
            equalities: Vec::new(),
            inequalities: Vec::new(),
        }
    }

    pub fn variables(&self) -> usize {
        self.variables
    }

    pub fn inequality_count(&self) -> usize {
        self.inequalities.len()
    }

    /// Add `sum coefficient * x_variable = rhs` and return its index.
    pub fn add_equality<I>(&mut self, coefficients: I, rhs: BigRational) -> usize
    where
        I: IntoIterator<Item = (usize, BigRational)>,
    {
        let constraint = self.constraint(coefficients, rhs);
        self.equalities.push(constraint);
        self.equalities.len() - 1
    }

    /// Add `sum coefficient * x_variable <= rhs` and return its index.
    pub fn add_inequality<I>(&mut self, coefficients: I, rhs: BigRational) -> usize
    where
        I: IntoIterator<Item = (usize, BigRational)>,
    {
        let constraint = self.constraint(coefficients, rhs);
        self.inequalities.push(constraint);
        self.inequalities.len() - 1
    }

    fn constraint<I>(&self, coefficients: I, rhs: BigRational) -> LinearConstraint
    where
        I: IntoIterator<Item = (usize, BigRational)>,
    {
        let mut dense = std::collections::BTreeMap::<usize, BigRational>::new();
        for (variable, coefficient) in coefficients {
            assert!(
                variable < self.variables,
                "variable {variable} is outside 0..{}",
                self.variables
            );
            *dense.entry(variable).or_insert_with(BigRational::zero) += coefficient;
        }
        LinearConstraint {
            coefficients: dense
                .into_iter()
                .filter(|(_, coefficient)| !coefficient.is_zero())
                .collect(),
            rhs,
        }
    }

    /// Compute the exact affine hull, or report that the polyhedron is empty.
    ///
    /// The result is verified by exact primal and dual certificates; see the
    /// module documentation.
    pub fn affine_hull(&self) -> AffineHull {
        let mut program = HomogenizedProgram::new(self);
        match program.solve() {
            Solution::Infeasible { duals } => {
                self.verify_farkas(&duals);
                AffineHull::Empty
            }
            Solution::Optimal { point, duals } => {
                let implicit = self.verify_hull(&point, &duals);
                let mut equations = self
                    .equalities
                    .iter()
                    .map(|constraint| self.dense(constraint))
                    .collect::<Vec<_>>();
                equations.extend(
                    self.inequalities
                        .iter()
                        .zip(&implicit)
                        .filter(|(_, &is_implicit)| is_implicit)
                        .map(|(constraint, _)| self.dense(constraint)),
                );
                let dimension = self.variables - exact_rank(equations);
                AffineHull::Nonempty(AffineHullData {
                    dimension,
                    implicit_equalities: implicit,
                    relative_interior_point: point,
                })
            }
        }
    }

    fn dense(&self, constraint: &LinearConstraint) -> Vec<BigRational> {
        let mut row = vec![BigRational::zero(); self.variables];
        for (variable, coefficient) in &constraint.coefficients {
            row[*variable] = coefficient.clone();
        }
        row
    }

    fn evaluate(constraint: &LinearConstraint, point: &[BigRational]) -> BigRational {
        constraint
            .coefficients
            .iter()
            .map(|(variable, coefficient)| coefficient * &point[*variable])
            .fold(BigRational::zero(), |sum, term| sum + term)
    }

    /// Check `C^T y + A^T w = 0` and return `y . d + w . b`.
    fn certificate_value(&self, duals: &Duals) -> BigRational {
        assert!(
            duals.inequalities.iter().all(|value| !value.is_negative()),
            "affine-hull certificate has a negative inequality multiplier"
        );
        let mut combination = vec![BigRational::zero(); self.variables];
        let mut value = BigRational::zero();
        for (constraint, multiplier) in self
            .inequalities
            .iter()
            .zip(&duals.inequalities)
            .chain(self.equalities.iter().zip(&duals.equalities))
        {
            if multiplier.is_zero() {
                continue;
            }
            for (variable, coefficient) in &constraint.coefficients {
                combination[*variable] += multiplier * coefficient;
            }
            value += multiplier * &constraint.rhs;
        }
        assert!(
            combination.iter().all(Zero::is_zero),
            "affine-hull certificate does not cancel the variables"
        );
        value
    }

    fn verify_farkas(&self, duals: &Duals) {
        assert!(
            self.certificate_value(duals).is_negative(),
            "emptiness certificate does not prove infeasibility"
        );
    }

    fn verify_hull(&self, point: &[BigRational], duals: &Duals) -> Vec<bool> {
        for constraint in &self.equalities {
            assert!(
                Self::evaluate(constraint, point) == constraint.rhs,
                "relative-interior witness violates an equality"
            );
        }
        let implicit = self
            .inequalities
            .iter()
            .map(|constraint| {
                let value = Self::evaluate(constraint, point);
                assert!(
                    value <= constraint.rhs,
                    "relative-interior witness violates an inequality"
                );
                value == constraint.rhs
            })
            .collect::<Vec<_>>();
        assert!(
            !self.certificate_value(duals).is_positive(),
            "affine-hull certificate has a positive constant"
        );
        for (is_implicit, multiplier) in implicit.iter().zip(&duals.inequalities) {
            assert!(
                !*is_implicit || multiplier.is_positive(),
                "affine-hull certificate omits a tight inequality"
            );
        }
        implicit
    }
}

/// Exact rank of a dense rational matrix.
pub fn exact_rank(mut matrix: Vec<Vec<BigRational>>) -> usize {
    let columns = matrix.first().map_or(0, Vec::len);
    let mut rank = 0;
    for column in 0..columns {
        let Some(pivot) = (rank..matrix.len()).find(|&row| !matrix[row][column].is_zero()) else {
            continue;
        };
        matrix.swap(rank, pivot);
        let pivot_row = matrix[rank].clone();
        for row in matrix.iter_mut().skip(rank + 1) {
            if row[column].is_zero() {
                continue;
            }
            let factor = &row[column] / &pivot_row[column];
            for (entry, pivot_entry) in row[column..].iter_mut().zip(&pivot_row[column..]) {
                *entry -= &factor * pivot_entry;
            }
        }
        rank += 1;
        if rank == matrix.len() {
            break;
        }
    }
    rank
}

struct Duals {
    equalities: Vec<BigRational>,
    inequalities: Vec<BigRational>,
}

enum Solution {
    Infeasible {
        duals: Duals,
    },
    Optimal {
        point: Vec<BigRational>,
        duals: Duals,
    },
}

/// Exact rational used inside the tableau.
///
/// Most tableau entries of these combinatorial programs are small.  `Small`
/// stores a reduced fraction with `i64` parts and computes in `i128`, which
/// cannot overflow for one product or sum of two such fractions; results that
/// do not fit fall back to `BigRational`.  This is only a representation
/// change: every operation is exact.
#[derive(Clone, Debug)]
enum Q {
    Small(i64, i64),
    Big(BigRational),
}

impl Q {
    fn zero() -> Self {
        Q::Small(0, 1)
    }

    fn one() -> Self {
        Q::Small(1, 1)
    }

    fn from_big(value: &BigRational) -> Self {
        match (value.numer().to_i64(), value.denom().to_i64()) {
            (Some(numerator), Some(denominator)) => Q::Small(numerator, denominator),
            _ => Q::Big(value.clone()),
        }
    }

    fn to_big(&self) -> BigRational {
        match self {
            Q::Small(numerator, denominator) => {
                BigRational::new_raw(BigInt::from(*numerator), BigInt::from(*denominator))
            }
            Q::Big(value) => value.clone(),
        }
    }

    fn from_i128(mut numerator: i128, mut denominator: i128) -> Self {
        debug_assert!(denominator != 0);
        if denominator < 0 {
            numerator = -numerator;
            denominator = -denominator;
        }
        let divisor = gcd_u128(numerator.unsigned_abs(), denominator as u128);
        if divisor > 1 {
            numerator /= divisor as i128;
            denominator /= divisor as i128;
        }
        match (i64::try_from(numerator), i64::try_from(denominator)) {
            (Ok(numerator), Ok(denominator)) => Q::Small(numerator, denominator),
            _ => Q::Big(BigRational::new(
                BigInt::from(numerator),
                BigInt::from(denominator),
            )),
        }
    }

    fn normalize_big(value: BigRational) -> Self {
        Self::from_big(&value)
    }

    fn is_zero(&self) -> bool {
        match self {
            Q::Small(numerator, _) => *numerator == 0,
            Q::Big(value) => value.is_zero(),
        }
    }

    fn is_positive(&self) -> bool {
        match self {
            Q::Small(numerator, _) => *numerator > 0,
            Q::Big(value) => value.is_positive(),
        }
    }

    fn is_negative(&self) -> bool {
        match self {
            Q::Small(numerator, _) => *numerator < 0,
            Q::Big(value) => value.is_negative(),
        }
    }

    fn add(&self, other: &Q) -> Q {
        match (self, other) {
            (Q::Small(a, b), Q::Small(c, d)) => {
                let (a, b, c, d) = (*a as i128, *b as i128, *c as i128, *d as i128);
                if b == d {
                    Q::from_i128(a + c, b)
                } else {
                    Q::from_i128(a * d + c * b, b * d)
                }
            }
            _ => Q::normalize_big(self.to_big() + other.to_big()),
        }
    }

    fn sub(&self, other: &Q) -> Q {
        self.add(&other.neg())
    }

    fn mul(&self, other: &Q) -> Q {
        match (self, other) {
            (Q::Small(a, b), Q::Small(c, d)) => {
                Q::from_i128(*a as i128 * *c as i128, *b as i128 * *d as i128)
            }
            _ => Q::normalize_big(self.to_big() * other.to_big()),
        }
    }

    fn div(&self, other: &Q) -> Q {
        match (self, other) {
            (Q::Small(a, b), Q::Small(c, d)) => {
                Q::from_i128(*a as i128 * *d as i128, *b as i128 * *c as i128)
            }
            _ => Q::normalize_big(self.to_big() / other.to_big()),
        }
    }

    fn neg(&self) -> Q {
        match self {
            Q::Small(numerator, denominator) => match numerator.checked_neg() {
                Some(numerator) => Q::Small(numerator, *denominator),
                None => Q::Big(-self.to_big()),
            },
            Q::Big(value) => Q::normalize_big(-value.clone()),
        }
    }
}

impl PartialEq for Q {
    fn eq(&self, other: &Self) -> bool {
        self.cmp(other) == std::cmp::Ordering::Equal
    }
}

impl Eq for Q {}

impl PartialOrd for Q {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Q {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        match (self, other) {
            (Q::Small(a, b), Q::Small(c, d)) => {
                (*a as i128 * *d as i128).cmp(&(*c as i128 * *b as i128))
            }
            _ => self.to_big().cmp(&other.to_big()),
        }
    }
}

fn gcd_u128(mut left: u128, mut right: u128) -> u128 {
    while right != 0 {
        (left, right) = (right, left % right);
    }
    left
}

type SparseRow = Vec<(usize, Q)>;

/// Standard-form tableau `max c z`, `M z = r`, `z >= 0`.
///
/// Columns are `x+`, `x-`, `sigma = tau - 1`, `t`, `s`, `u`, then
/// artificials.  Rows are the equalities, the inequalities
/// `C x - d sigma + t + s = d`, and the bounds `t + u = 1`.  A row may be
/// negated to make its right side nonnegative; `row_sign` records this.
struct HomogenizedProgram {
    variables: usize,
    equalities: usize,
    inequalities: usize,
    rows: Vec<SparseRow>,
    rhs: Vec<Q>,
    basis: Vec<usize>,
    row_sign: Vec<bool>,
    identity_column: Vec<usize>,
    columns: usize,
    artificial_start: usize,
    objective: Vec<Q>,
    costs: Vec<Q>,
}

impl HomogenizedProgram {
    fn new(system: &RationalPolyhedron) -> Self {
        let n = system.variables;
        let p = system.equalities.len();
        let m = system.inequalities.len();
        let sigma = 2 * n;
        let t = |i: usize| 2 * n + 1 + i;
        let s = |i: usize| 2 * n + 1 + m + i;
        let u = |i: usize| 2 * n + 1 + 2 * m + i;
        let artificial_start = 2 * n + 1 + 3 * m;

        let mut rows = Vec::with_capacity(p + 2 * m);
        let mut rhs = Vec::with_capacity(p + 2 * m);
        let mut basis = Vec::with_capacity(p + 2 * m);
        let mut row_sign = Vec::with_capacity(p + 2 * m);
        let mut identity_column = Vec::with_capacity(p + 2 * m);
        let mut artificial = artificial_start;

        let structural = |constraint: &LinearConstraint| {
            let mut row = Vec::with_capacity(2 * constraint.coefficients.len() + 3);
            for (variable, coefficient) in &constraint.coefficients {
                let coefficient = Q::from_big(coefficient);
                row.push((n + variable, coefficient.neg()));
                row.push((*variable, coefficient));
            }
            if !constraint.rhs.is_zero() {
                row.push((sigma, Q::from_big(&constraint.rhs).neg()));
            }
            row
        };

        for constraint in &system.equalities {
            let mut row = structural(constraint);
            let mut value = Q::from_big(&constraint.rhs);
            let negate = value.is_negative();
            if negate {
                negate_row(&mut row);
                value = value.neg();
            }
            row.push((artificial, Q::one()));
            rows.push(row);
            rhs.push(value);
            basis.push(artificial);
            row_sign.push(negate);
            identity_column.push(artificial);
            artificial += 1;
        }
        for (index, constraint) in system.inequalities.iter().enumerate() {
            let mut row = structural(constraint);
            row.push((t(index), Q::one()));
            row.push((s(index), Q::one()));
            let mut value = Q::from_big(&constraint.rhs);
            if value.is_negative() {
                negate_row(&mut row);
                value = value.neg();
                row.push((artificial, Q::one()));
                basis.push(artificial);
                identity_column.push(artificial);
                row_sign.push(true);
                artificial += 1;
            } else {
                basis.push(s(index));
                identity_column.push(s(index));
                row_sign.push(false);
            }
            rows.push(row);
            rhs.push(value);
        }
        for index in 0..m {
            rows.push(vec![(t(index), Q::one()), (u(index), Q::one())]);
            rhs.push(Q::one());
            basis.push(u(index));
            row_sign.push(false);
            identity_column.push(u(index));
        }
        for row in &mut rows {
            row.sort_unstable_by_key(|(column, _)| *column);
        }

        let columns = artificial;
        Self {
            variables: n,
            equalities: p,
            inequalities: m,
            rows,
            rhs,
            basis,
            row_sign,
            identity_column,
            columns,
            artificial_start,
            objective: vec![Q::zero(); columns],
            costs: vec![Q::zero(); columns],
        }
    }

    fn solve(&mut self) -> Solution {
        if self.columns > self.artificial_start {
            self.costs = vec![Q::zero(); self.columns];
            for column in self.artificial_start..self.columns {
                self.costs[column] = Q::one().neg();
            }
            self.reset_objective();
            self.optimize(true);
            let infeasible = self
                .basis
                .iter()
                .zip(&self.rhs)
                .any(|(&column, value)| column >= self.artificial_start && value.is_positive());
            if infeasible {
                return Solution::Infeasible {
                    duals: self.duals(),
                };
            }
            self.drive_out_artificials();
        }

        self.costs = vec![Q::zero(); self.columns];
        for index in 0..self.inequalities {
            self.costs[2 * self.variables + 1 + index] = Q::one();
        }
        self.reset_objective();
        self.optimize(false);

        let mut values = vec![Q::zero(); self.columns];
        for (row, &column) in self.basis.iter().enumerate() {
            values[column] = self.rhs[row].clone();
        }
        let tau = BigRational::one() + values[2 * self.variables].to_big();
        let point = (0..self.variables)
            .map(|variable| {
                (values[variable].to_big() - values[self.variables + variable].to_big()) / &tau
            })
            .collect();
        Solution::Optimal {
            point,
            duals: self.duals(),
        }
    }

    /// Simplex multipliers of the original, unnegated constraint rows.
    fn duals(&self) -> Duals {
        let multiplier = |row: usize| {
            let column = self.identity_column[row];
            let value = self.costs[column].sub(&self.objective[column]).to_big();
            if self.row_sign[row] {
                -value
            } else {
                value
            }
        };
        Duals {
            equalities: (0..self.equalities).map(multiplier).collect(),
            inequalities: (self.equalities..self.equalities + self.inequalities)
                .map(multiplier)
                .collect(),
        }
    }

    fn reset_objective(&mut self) {
        self.objective = self.costs.clone();
        for (row, &column) in self.basis.iter().enumerate() {
            let cost = self.costs[column].clone();
            if cost.is_zero() {
                continue;
            }
            for (entry_column, value) in &self.rows[row] {
                self.objective[*entry_column] = self.objective[*entry_column].sub(&cost.mul(value));
            }
        }
    }

    fn entry(&self, row: usize, column: usize) -> Option<&Q> {
        let entries = &self.rows[row];
        entries
            .binary_search_by_key(&column, |(entry_column, _)| *entry_column)
            .ok()
            .map(|index| &entries[index].1)
    }

    fn optimize(&mut self, allow_artificial: bool) {
        const DEGENERATE_LIMIT: usize = 16;
        let column_limit = if allow_artificial {
            self.columns
        } else {
            self.artificial_start
        };
        let mut degenerate_run = 0;
        loop {
            let bland = degenerate_run >= DEGENERATE_LIMIT;
            let mut entering: Option<usize> = None;
            for column in 0..column_limit {
                if !self.objective[column].is_positive() {
                    continue;
                }
                if bland {
                    entering = Some(column);
                    break;
                }
                if entering.is_none_or(|best| self.objective[column] > self.objective[best]) {
                    entering = Some(column);
                }
            }
            let Some(entering) = entering else {
                return;
            };

            let mut leaving: Option<(usize, Q)> = None;
            for row in 0..self.rows.len() {
                let Some(value) = self.entry(row, entering) else {
                    continue;
                };
                if !value.is_positive() {
                    continue;
                }
                let ratio = self.rhs[row].div(value);
                let better = match &leaving {
                    None => true,
                    Some((best_row, best_ratio)) => {
                        ratio < *best_ratio
                            || (ratio == *best_ratio && self.basis[row] < self.basis[*best_row])
                    }
                };
                if better {
                    leaving = Some((row, ratio));
                }
            }
            let (leaving, ratio) = leaving.expect("the homogenized affine-hull program is bounded");
            if ratio.is_zero() {
                degenerate_run += 1;
            } else {
                degenerate_run = 0;
            }
            self.pivot(leaving, entering);
        }
    }

    /// Remove zero-valued artificial variables from the basis when possible.
    ///
    /// A row whose structural entries all vanish is redundant; its artificial
    /// remains basic at zero and is unaffected by later pivots.
    fn drive_out_artificials(&mut self) {
        for row in 0..self.rows.len() {
            if self.basis[row] < self.artificial_start {
                continue;
            }
            let replacement = self.rows[row]
                .iter()
                .find(|(column, value)| *column < self.artificial_start && !value.is_zero())
                .map(|(column, _)| *column);
            if let Some(column) = replacement {
                self.pivot(row, column);
            }
        }
    }

    fn pivot(&mut self, pivot_row: usize, pivot_column: usize) {
        let pivot = self
            .entry(pivot_row, pivot_column)
            .expect("pivot entry is nonzero")
            .clone();
        let mut normalized = std::mem::take(&mut self.rows[pivot_row]);
        for (_, value) in &mut normalized {
            *value = value.div(&pivot);
        }
        self.rhs[pivot_row] = self.rhs[pivot_row].div(&pivot);
        let normalized_rhs = self.rhs[pivot_row].clone();

        for row in 0..self.rows.len() {
            if row == pivot_row {
                continue;
            }
            let Some(factor) = self.entry(row, pivot_column).cloned() else {
                continue;
            };
            let current = std::mem::take(&mut self.rows[row]);
            self.rows[row] = subtract_multiple(current, &factor, &normalized);
            self.rhs[row] = self.rhs[row].sub(&factor.mul(&normalized_rhs));
        }
        let factor = self.objective[pivot_column].clone();
        if !factor.is_zero() {
            for (column, value) in &normalized {
                self.objective[*column] = self.objective[*column].sub(&factor.mul(value));
            }
        }
        self.rows[pivot_row] = normalized;
        self.basis[pivot_row] = pivot_column;
    }
}

fn negate_row(row: &mut SparseRow) {
    for (_, value) in row {
        *value = value.neg();
    }
}

/// Return `row - factor * pivot` for rows sorted by column.
fn subtract_multiple(row: SparseRow, factor: &Q, pivot: &SparseRow) -> SparseRow {
    let mut result = Vec::with_capacity(row.len() + pivot.len());
    let mut left = row.into_iter().peekable();
    let mut right = pivot.iter().peekable();
    loop {
        match (left.peek(), right.peek()) {
            (None, None) => break,
            (Some(_), None) => result.push(left.next().expect("peeked")),
            (None, Some(_)) => {
                let (column, value) = right.next().expect("peeked");
                result.push((*column, factor.mul(value).neg()));
            }
            (Some((left_column, _)), Some((right_column, _))) => {
                if left_column < right_column {
                    result.push(left.next().expect("peeked"));
                } else if right_column < left_column {
                    let (column, value) = right.next().expect("peeked");
                    result.push((*column, factor.mul(value).neg()));
                } else {
                    let (column, value) = left.next().expect("peeked");
                    let (_, pivot_value) = right.next().expect("peeked");
                    let value = value.sub(&factor.mul(pivot_value));
                    if !value.is_zero() {
                        result.push((column, value));
                    }
                }
            }
        }
    }
    result
}

/// Convert an integer into an exact rational coefficient.
pub fn rational(value: impl Into<BigInt>) -> BigRational {
    BigRational::from_integer(value.into())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn r(value: i64) -> BigRational {
        rational(value)
    }

    fn data(hull: AffineHull) -> AffineHullData {
        match hull {
            AffineHull::Nonempty(data) => data,
            AffineHull::Empty => panic!("expected a nonempty polyhedron"),
        }
    }

    #[test]
    fn interval_and_point() {
        let mut interval = RationalPolyhedron::new(1);
        interval.add_inequality([(0, r(-1))], r(0));
        interval.add_inequality([(0, r(1))], r(1));
        let hull = data(interval.affine_hull());
        assert_eq!(hull.dimension, 1);
        assert_eq!(hull.implicit_equalities, vec![false, false]);

        let mut point = RationalPolyhedron::new(1);
        point.add_inequality([(0, r(-1))], r(-2));
        point.add_inequality([(0, r(1))], r(2));
        let hull = data(point.affine_hull());
        assert_eq!(hull.dimension, 0);
        assert_eq!(hull.implicit_equalities, vec![true, true]);
        assert_eq!(hull.relative_interior_point, vec![r(2)]);
    }

    #[test]
    fn jointly_forced_upper_bounds() {
        // b <= 3, c <= a, d <= 2 - a, b + c + d = 5, 1 <= a <= 2:
        // the three upper bounds sum to the equality, so all are tight.
        let (a, b, c, d) = (0, 1, 2, 3);
        let mut system = RationalPolyhedron::new(4);
        system.add_equality([(b, r(1)), (c, r(1)), (d, r(1))], r(5));
        let upper_b = system.add_inequality([(b, r(1))], r(3));
        let upper_c = system.add_inequality([(c, r(1)), (a, r(-1))], r(0));
        let upper_d = system.add_inequality([(d, r(1)), (a, r(1))], r(2));
        let lower_a = system.add_inequality([(a, r(-1))], r(-1));
        let upper_a = system.add_inequality([(a, r(1))], r(2));
        let hull = data(system.affine_hull());
        assert_eq!(hull.dimension, 1);
        for index in [upper_b, upper_c, upper_d] {
            assert!(hull.implicit_equalities[index]);
        }
        for index in [lower_a, upper_a] {
            assert!(!hull.implicit_equalities[index]);
        }
    }

    #[test]
    fn rational_only_directions_are_detected() {
        // x + y = 1, 0 <= x, 0 <= y, and 2x <= 1: its only integer point is
        // (0, 1), but the rational polytope is a segment.
        let mut system = RationalPolyhedron::new(2);
        system.add_equality([(0, r(1)), (1, r(1))], r(1));
        system.add_inequality([(0, r(-1))], r(0));
        system.add_inequality([(1, r(-1))], r(0));
        system.add_inequality([(0, r(2))], r(1));
        let hull = data(system.affine_hull());
        assert_eq!(hull.dimension, 1);
        assert_eq!(hull.implicit_equalities, vec![false, false, false]);
    }

    #[test]
    fn infeasible_systems_are_empty() {
        let mut system = RationalPolyhedron::new(2);
        system.add_equality([(0, r(1)), (1, r(1))], r(4));
        system.add_inequality([(0, r(1))], r(1));
        system.add_inequality([(1, r(1))], r(2));
        assert_eq!(system.affine_hull(), AffineHull::Empty);

        let mut contradictory = RationalPolyhedron::new(0);
        contradictory.add_inequality([], r(-1));
        assert_eq!(contradictory.affine_hull(), AffineHull::Empty);
    }

    #[test]
    fn redundant_equalities_and_unbounded_directions() {
        let mut system = RationalPolyhedron::new(3);
        system.add_equality([(0, r(1)), (1, r(1))], r(2));
        system.add_equality([(0, r(2)), (1, r(2))], r(4));
        system.add_inequality([(0, r(-1))], r(0));
        let hull = data(system.affine_hull());
        // x2 is free and x0 + x1 = 2 with x0 >= 0.
        assert_eq!(hull.dimension, 2);
        assert_eq!(hull.implicit_equalities, vec![false]);
    }

    #[test]
    fn large_coefficients_fall_back_to_big_rationals() {
        // A segment translated by 2^130 exceeds the small-rational fast path.
        let offset = BigInt::from(1) << 130usize;
        let mut segment = RationalPolyhedron::new(1);
        segment.add_inequality([(0, r(-1))], -rational(offset.clone()));
        segment.add_inequality([(0, r(1))], rational(offset.clone() + 1));
        let hull = data(segment.affine_hull());
        assert_eq!(hull.dimension, 1);
        assert_eq!(hull.implicit_equalities, vec![false, false]);
        let mut point = RationalPolyhedron::new(2);
        point.add_equality([(0, r(1)), (1, r(1))], rational(offset.clone() * 2));
        point.add_inequality([(0, r(-1)), (1, r(1))], r(0));
        point.add_inequality([(0, r(1)), (1, r(-1))], r(0));
        let hull = data(point.affine_hull());
        assert_eq!(hull.dimension, 0);
        assert_eq!(
            hull.relative_interior_point,
            vec![rational(offset.clone()), rational(offset)]
        );
    }

    #[test]
    fn empty_system_is_full_dimensional() {
        let hull = data(RationalPolyhedron::new(3).affine_hull());
        assert_eq!(hull.dimension, 3);
    }
}
