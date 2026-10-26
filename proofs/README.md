# Proofs and their verification

`build.py` builds Online Resources 1–3 from the Markdown files below. Their PDFs are in `../online-resources/`.

| Online Resource | Source | Contents |
|---|---|---|
| 1 | `theorem1/symmetry_breaking_proof.md` | Complete proof of Theorem 1 (maximum likelihood, exchangeable design), the zero-uniqueness boundary (Section 4), the unrestricted ULS proof (Section 11) and an independent algebraic classification of the ML stationary points |
| 2 | `theorem1/weighted_least_squares_proof.md` | Unrestricted global proof of the weighted least-squares result for all positive pair-class weights, including the threshold and the equality cases |
| 3 | `unequal_blocks/unequal_blocks_proof.md` | Block reduction (Lemma 1), the exact two-parameter problem and its quartic, the cusp and fold computations (Theorem 3), mirrored heterogeneous designs, and the verification run |

Each proof file ends with its remaining gaps and exact scope. The numerical checks support the algebra and the implementation. The global claims and remainder bounds rest on the written arguments.

## theorem1/

| File | Contents |
|---|---|
| `section4.md` | Section 4 of Online Resource 1 (Heywood boundary points) as a separate file |
| `lemma71_add.md` | Addendum: the symbolic verification of equation (7.3) of Online Resource 1 |
| `verify_symmetry_breaking.py` | Checks the closed forms against independent optimization and verifies the symbolic identities (4,224 checks) |
| `verify_results.json`, `verify_results.log` | Saved output of its last run (0 failures) |
| `check_saddle_7_3_sympy.py`, `check73.py`, `numcheck.py`, `rotation.py` | SymPy check of the saddle expansion, equation (7.3), and the rotation argument |

Run from this folder: `python3 theorem1/verify_symmetry_breaking.py --output results.json`. It needs NumPy, SciPy, SymPy and mpmath.

## unequal_blocks/

| File | Contents |
|---|---|
| `verify_unequal_blocks.py` | Independent verification of Lemma 1, the exact quartic and Theorem 3 (133 checks); it does not import the reference computations |
| `verify_results.json` | Saved output of its last run (0 failures) |
| `reference_section31.py`, `reference_section31.json` | Reference computations that the verifier is compared against |

Run from this folder: `python3 unequal_blocks/verify_unequal_blocks.py --output results.json`.
Neither command overwrites the saved output.

`build.py --checks`, run from the repository root, reruns the exact algebra, this verifier and the interval certificate of `../checks/unequal_example_certificate.py`. It writes the results to `../build/checks/` and leaves the saved output here unchanged.
