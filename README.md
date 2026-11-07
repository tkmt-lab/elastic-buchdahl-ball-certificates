# An interval-certified relativistic elastic ball

Version **1.0.0 — 2026-11-09**. DOI: —. [Article](paper.pdf).

For one Alho–Natário–Pani–Raposo material and central value δc=8, interval certificates establish a regular finite static ball with nonnegative pressures, the dominant energy condition and five principal-direction squared wave speeds in [0,1].
Its surface compactness satisfies **0.91089951 < 2M/R < 0.91180552**, exceeding 8/9.
The certificates give a global rigorous check of one anisotropic elastic example in a family whose numerical exceedance was established by those authors; Buchdahl's perfect-fluid assumptions differ from this model.

Python 3.10+, macOS/Linux; from the extracted packet root:
```sh
python3 -I -B verify.py
```
Expected final JSON: `"status": "PASS"`, `"full_replay": true`.
Default limits: A 3300 seconds, total 3600 seconds; see [REPRODUCE.md](REPRODUCE.md).

## Verification and mathematical content

The command checks both manifests, copies the certificates to a temporary directory, recomputes all 1592 B tubes and 465 A tubes (including all 4185 C1 entries), and compares the exact results with the article bounds. It prints a final PASS only after all checks pass. Both certificates use a shared integer interval engine with different continuation schemes. The check runs offline using Python's standard library. Allow about 8 minutes on the previously measured machine; slower machines may need more. [REPRODUCE.md](REPRODUCE.md) describes the time limits and failure conditions.

Run `python3 -I -B negative_controls.py` for nine rejection tests; full verification requires the command above as well.

The radius in the fixed normalization satisfies 0.01258970 < R < 0.01259407. B checks that the center lemma's endpoint enclosure lies in its first initial box and that the starting time equals x₀. A uses Picard self-inclusion for the entire initial box to connect its tubes to the same central solution. Section 3.2 gives all six exact dyadic endpoints of B's first initial box and the rational center constant B, so the inclusion in equation (7) can be checked from the paper alone. Section 3.3 defines the nonautonomous flow coefficient c₁₃(τ,v) and the variational flow V(t) up to the current time.

## Scope

The contribution is a global interval certificate for one fixed material and central value. The reference shear modulus is μ=−60. Physical material viability, arbitrary-direction causality, dynamical stability and general constitutive families are outside the claim. No priority claim is made for the phenomenon or for the first rigorous treatment of this family. The Lean file is a type-checked L1 proposition, with no formal proof. This article has not been externally peer-reviewed. The related computations by Merkushev cited below were not independently checked here. Approximate margin displays below are not exact decimal lower bounds; the integer endpoints and replay determine the proof. Review history and the scope of static AI review are recorded in [VERIFICATION-RECORD.md](VERIFICATION-RECORD.md).

## Obtain the matching version

The release version is [v1.0.0](https://github.com/tkmt-lab/elastic-buchdahl-ball-certificates/releases/tag/v1.0.0).

The article identifies the numerical and checker snapshot by `CERTIFICATE-SHA256SUMS`, SHA-256 **`431b5dc0ad9a67b3f17388b0560fbf6e1556431801ddd0ecc0e51d0e82fd15fa`**. That index fixes both certificates, their integer engine and checkers, the model specification, claims and entry scripts. `MANIFEST.sha256` identifies all files in this packet, including the article and documentation. From the extracted packet root, `shasum -a 256 CERTIFICATE-SHA256SUMS` (macOS) or `sha256sum CERTIFICATE-SHA256SUMS` (Linux) must match the digest printed in the article. The verifier checks every indexed file as well as the closed full manifest.

## Displayed global margins

The saved B-check record includes these approximate displays of whole-ball interval minima, including the center (see Scope):

| Quantity | Approximate display |
|---|---:|
| D | 0.08819448972521368 |
| ρ−pt | 55.23286197508267 |
| Third transverse squared speed | 0.10558269515105372 |
| 1 minus first transverse squared speed | 0.1596350304276371 |

## Contents

- `certificates/stage1-tubes.jsonl`: independently generated 1592 tubes (certificate B).
- `certificates/stage2-inputs/pilot-steps.jsonl`: 465 Taylor-flow tubes (certificate A).
- `certificates/independent.py`, `check_own_tubes.py`, `replay_producer.py`: integer interval engine and checkers.
- `spec.json`, `FORMAT.md`: fixed constitutive model, interval and certificate formats.
- `MainStatement.lean`, `lean-check.json`: existence-claim statement and L1 type-check record.
- `paper.tex`, `paper.pdf`: mathematical argument and final quantitative bounds, with the center-to-B inclusion and A tube/remainder definitions explicit.
- `CERTIFICATE-SHA256SUMS`: article-identified frozen numerical and checker snapshot; `MANIFEST.sha256`: complete packet inventory.
- `REVIEW-RESPONSE.md`: dispositions of all three blind-review requests.
- `VERIFICATION-RECORD.md`, `provenance.json`: evidence, input lineage and limits of independence.

## Related work and acknowledgment

Related certificate-supported equilibria for a different constitutive law are reported by [Merkushev (2026)](https://doi.org/10.5281/zenodo.20836296), at central value 1/2 and M/R approximately 0.11749. That model differs from the specified ACS material and the present high-compactness example.

We are grateful to the original authors for the material, numerical examples and detailed elasticity framework used here; see [sources and precise parameter correspondence](SOURCES.md).

## License and contact

Code is Apache-2.0; original article, documentation and certificate data are CC BY 4.0. See [LICENSE](LICENSE). External articles are available through the links in [SOURCES.md](SOURCES.md). Citation metadata are in [CITATION.cff](CITATION.cff).

Tatsuhiko Takimoto — research contact: [research@tatsuhikotakimoto.com](mailto:research@tatsuhikotakimoto.com).
