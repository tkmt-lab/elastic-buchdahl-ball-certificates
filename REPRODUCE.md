# Reproduction

## Version and retrieval

Version 1.0.0 (2026-11-09) is `v1.0.0` at <https://github.com/tkmt-lab/elastic-buchdahl-ball-certificates/releases/tag/v1.0.0>. The scientific/checker snapshot is pinned by `CERTIFICATE-SHA256SUMS`, whose SHA-256 is `431b5dc0ad9a67b3f17388b0560fbf6e1556431801ddd0ecc0e51d0e82fd15fa` and is also printed in the article. No public upload is performed by these instructions.

Compare the index fingerprint after extraction:

```sh
# macOS
shasum -a 256 CERTIFICATE-SHA256SUMS
# Linux
sha256sum CERTIFICATE-SHA256SUMS
```

`verify.py` then checks the thirteen indexed numerical/checker files against that index, in addition to the full closed `MANIFEST.sha256`. The latter covers revised paper/documentation too. A refreshed full manifest cannot silently replace a scientific input while retaining the article's fixed identity. Neither inventory is a signature or an assertion of publication.

## Required offline check

From the packet root on macOS or Linux with Python >=3.10:

```sh
python3 -I -B verify.py
python3 -I -B negative_controls.py
```

Only the standard library is required. Python `resource` is POSIX-specific. The verifier creates and removes a temporary directory and does not change the packet. The expected inventory is exactly the paths in `MANIFEST.sha256`, plus the manifest itself; Git's `.git` directory is excluded. A stray `.DS_Store`, editor backup or build output causes `manifest inventory mismatch`. Keep such files outside the packet and build TeX elsewhere. Run without `-O` or `-OO`: the numerical checkers use assertions, and the entry point explicitly refuses optimized Python. `-I -B` isolates Python configuration and avoids bytecode files.

The complete pass consists of both manifest/identity checks, 1592 center-to-surface tube checks and 465 full Taylor-flow checks (all 4185 Jacobian entries). Previously saved PASS JSON is never used instead of recomputation. The recomputed exact endpoints must match the shipped record, and the outward decimals in the article must enclose them. Incomplete replay, false booleans, wrong row ordering, timeout, or subprocess failure are failures.

The final JSON must report `status: PASS`, `full_replay: true`, 1592 B tubes, 465 A tubes, 4185 C1 entries, and the matching `certificate_identity_sha256`. The negative-control command must report `status: PASS` and `negative_controls: 9`. Intermediate progress or a saved PASS file is not the final success condition.

The outer default budget is 3600 seconds (`--budget-seconds` accepts values in (0,3600]). The A producer replay internal cutoff is 3300 seconds. Both are bounded and fail closed: all 465 A and 1592 B tubes remain mandatory. No partial run is accepted. On the earlier 480-second limit, the 2026-11-05 Linux audit completed only 346 A tubes; that historical failure motivated this revision.

## Optional article build

With an already provisioned TeX/Tectonic environment, build into a directory outside this packet:

```sh
tectonic --only-cached --keep-logs --outdir /path/to/build paper.tex
```

The build uses article, T1 fontenc, lmodern, microtype, amsmath, amssymb, amsthm, booktabs, geometry, xurl and hyperref. Cached fonts/packages are needed for an offline build. The shipped PDF is sufficient for the required verification.

## Optional Lean statement check

With Lean 4.35.0 and a matching existing mathlib build available, run `lake env lean MainStatement.lean` from that environment, passing the file path as needed. No Lean dependency is needed by the integer checker. `lean-check.json` records the 2026-10-27 type check, including the mathlib manifest hash. No theorem, `sorry` or added axiom occurs in the statement file. Schwarzschild matching and the analytic justification of the interval algorithms remain prose arguments.

## Clean clone

To check a clean clone, export only this packet into a local Git repository and clone it with `git clone --no-local`. Run the same two Python commands in that clone; rebuild TeX outside it. The checked packet manifest identifies the exact payload. This local preflight is distinct from remote CI, publication and DOI registration. The exported snapshot commit and results are retained with the dated verification report, outside the packet to avoid self-referential hashes.

## Center connections and A remainder

In `certificates/check_own_tubes.py`, `center_certificate()` recomputes the degree-seven center and its endpoint enclosure. On B's first row, `rx == x` fixes the start at x₀=1415/2097152, and `all(a.contains(b) for a,b in zip(yi,y))` contains that central endpoint in the first Y0 box. The same assertion on later rows chains the previous endpoint. No coordinate change occurs in B.

In `certificates/replay_producer.py`, `source_center()` checks the central polynomial and enclosure, and `main()` contains `local(own)` in the first v_start box. `verify.py` separately checks its first x_exact. Each call to `replay()` checks `start+[0,h]*vf(X,tube)` inside the stored tube, uniformly for every initial point in start. It then computes the field Jacobian, verifies the variational tube, and checks degree-13 remainder and all nine derivative entries. Article §3.3 defines c₁₃(τ,v) at elapsed time zero for a local flow starting at absolute time τ, and V(t) as the derivative of that flow from the initial point to time t. These are the quantities used by `coefficients(X,tube,13,True)` and the stored variational_tube. B arithmetic is unchanged; A changes only its runtime cutoff.
