import Mathlib

/- L1: a type-checked proposition, not a theorem or a formal proof. -/
namespace ElasticBall
noncomputable section

def a (u : ℝ) : ℝ := -(9/2 : ℝ) - 65/u^4 + 70/u^2
def b (u : ℝ) : ℝ := (191/2 : ℝ) - 130*u^2 + 35*u^4
def ap (u : ℝ) : ℝ := 260/u^5 - 140/u^3
def bp (u : ℝ) : ℝ := -260*u + 140*u^3
def rho (d u : ℝ) : ℝ := a u*d^2 + b u
def pr (d u : ℝ) : ℝ := a u*d^2 - b u
def q (d u : ℝ) : ℝ := u*(ap u*d^2 + bp u)/2
def pt (d u : ℝ) : ℝ := pr d u + q d u
def D (s z : ℝ) : ℝ := 1-2*s*z
def G (s d u z : ℝ) : ℝ := d/(u^2*Real.sqrt (D s z))-u
def DT (d u : ℝ) : ℝ := -(9/2 : ℝ)*d^2 - 65*u^2 + 35*(u^4+d^2/u^2)
def DL (d u : ℝ) : ℝ := a u*d^2
def K (u : ℝ) : ℝ := -65*u^2+35*u^4
def causal (w : ℝ) : Prop := 0 ≤ w ∧ w ≤ 1

def Conditions (s d u z : ℝ) : Prop :=
  0 < d ∧ 0 < u ∧ 0 < D s z ∧ 0 < 2*a u*d ∧ 0 < DT d u ∧ 0 < DL d u ∧
  0 ≤ pr d u ∧ 0 < pt d u ∧
  0 ≤ rho d u - pr d u ∧ 0 ≤ rho d u + pr d u ∧
  0 ≤ rho d u - pt d u ∧ 0 ≤ rho d u + pt d u ∧
  causal 1 ∧ causal 1 ∧ causal (K u / DT d u) ∧
  causal ((d/u^3)^2*K u / DL d u) ∧ causal ((-65*u^2+35*d^2/u^2)/DT d u)

def StaticODE (d u z : ℝ → ℝ) (s : ℝ) : Prop :=
  HasDerivAt d ((q (d s) (u s) - (ap (u s)*(d s)^2-bp (u s))*G s (d s) (u s) (z s)/2
    - s*(rho (d s) (u s)+pr (d s) (u s))*(z s+pr (d s) (u s))/(2*D s (z s)))
    / (2*a (u s)*d s) / s) s ∧
  HasDerivAt u (G s (d s) (u s) (z s)/(2*s)) s ∧
  HasDerivAt z ((rho (d s) (u s)-3*z s)/(2*s)) s

/-- Existence and quantitative core of the article. Schwarzschild matching is
    supplied by the article's analytic argument, not formalized here. -/
def MainClaim : Prop :=
  ∃ (S : ℝ) (d u z : ℝ → ℝ), 0 < S ∧
    ContDiffOn ℝ 1 d (Set.Icc 0 S) ∧ ContDiffOn ℝ 1 u (Set.Icc 0 S) ∧
    ContDiffOn ℝ 1 z (Set.Icc 0 S) ∧
    d 0 = 8 ∧ u 0 = 2 ∧ z 0 = (1415/6 : ℝ) ∧
    (∀ s ∈ Set.Ioc 0 S, StaticODE d u z s) ∧
    (∀ s ∈ Set.Icc 0 S, Conditions s (d s) (u s) (z s)) ∧
    (∀ s ∈ Set.Ico 0 S, 0 < pr (d s) (u s)) ∧ pr (d S) (u S) = 0 ∧
    (∃ v : ℝ, v < 0 ∧ HasDerivAt (fun s => pr (d s) (u s)) v S) ∧
    (1258970/100000000 : ℝ) < Real.sqrt (S/(4*Real.pi)) ∧
    Real.sqrt (S/(4*Real.pi)) < (1259407/100000000 : ℝ) ∧
    (91089951/100000000 : ℝ) < 2*S*z S ∧
    2*S*z S < (91180552/100000000 : ℝ) ∧
    (2201062/100000000 : ℝ) < 2*S*z S - (8/9 : ℝ)

#check MainClaim
end
end ElasticBall
