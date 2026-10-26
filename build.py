#!/usr/bin/env python3
"""Single build entry point. Optional focused checks never rerun the simulation."""
from pathlib import Path
import argparse
import os
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
BUILD = ROOT / 'build'
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--checks', action='store_true', help='Also rerun focused Python/R checks')
args = parser.parse_args()
env = os.environ.copy()
env.update(OPENBLAS_NUM_THREADS='1', OMP_NUM_THREADS='1', MKL_NUM_THREADS='1',
           PYTHONDONTWRITEBYTECODE='1', NOT_CRAN='true', TZ='UTC')
for sub in ('manuscript', 'proofs', 'checks'):
    (BUILD / sub).mkdir(parents=True, exist_ok=True)

def run(command, log, cwd=ROOT):
    print('Running:', ' '.join(map(str, command)), flush=True)
    with (BUILD / log).open('w') as output:
        result = subprocess.run(list(map(str, command)), cwd=cwd, env=env,
                                stdout=output, stderr=subprocess.STDOUT)
    if result.returncode:
        raise SystemExit(f'Failed ({result.returncode}); see {BUILD / log}')

run(['latexmk', '-pdf', '-interaction=nonstopmode', '-halt-on-error',
     '-outdir=../build/manuscript', 'main.tex'], 'manuscript/build-console.log', ROOT/'manuscript')
run(['latexmk', '-pdf', '-interaction=nonstopmode', '-halt-on-error',
     '-outdir=../build/manuscript', 'cover_page.tex'], 'manuscript/cover-console.log', ROOT/'manuscript')
# Named preprint: the same sources with author details and declarations.
run(['latexmk', '-pdf', '-interaction=nonstopmode', '-halt-on-error', '-jobname=preprint',
     r'-usepretex=\def\preprint{}', '-outdir=../build/manuscript', 'main.tex'],
    'manuscript/preprint-console.log', ROOT/'manuscript')
placeholders = [f'{name}: {line.strip()}'
                for name in ('author.tex', 'declarations.tex', 'preprint_front.tex')
                for line in (ROOT/'manuscript'/name).read_text().splitlines()
                if ('AUTHOR MUST CONFIRM' in line
                    or (not line.lstrip().startswith('%') and '[' in line and 'to be' in line))]
print(f'Preprint placeholders to fill before posting: {len(placeholders)}')
for p in placeholders:
    print('  ' + p)
proofs = [('S1-exchangeable-ML-ULS', 'theorem1/symmetry_breaking_proof.md'),
          ('S2-weighted-least-squares', 'theorem1/weighted_least_squares_proof.md'),
          ('S3-unequal-blocks', 'unequal_blocks/unequal_blocks_proof.md')]
for name, source in proofs:
    run(['pandoc', ROOT/'proofs'/source, '-f', 'markdown+tex_math_single_backslash',
         '-t', 'latex', '-s', '-V', 'geometry:margin=1in', '-V', 'fontsize=10pt',
         '--pdf-engine=xelatex', '-o', BUILD/'proofs'/f'{name}.pdf'], f'proofs/{name}.log')
run([sys.executable, 'checks/verify_portable_evidence.py'], 'checks/portable-evidence.log')
run([sys.executable, 'checks/source_consistency.py'], 'checks/source-consistency.log')
run([sys.executable, 'checks/table2_ratios.py'], 'checks/table2-ratios.log')
if args.checks:
    run([sys.executable, 'checks/exact_algebra.py'], 'checks/exact-algebra.log')
    run([sys.executable, 'checks/theorem1_exact.py'], 'checks/theorem1-exact.log')
    run([sys.executable, 'checks/unequal_example_certificate.py'], 'checks/unequal-example-certificate.log')
    run([sys.executable, 'proofs/unequal_blocks/verify_unequal_blocks.py',
         '--output', BUILD/'checks/unequal-rerun.json'], 'checks/unequal-rerun.log')
    run(['Rscript', '-e', 'pkgload::load_all("undeRfactor",quiet=TRUE); '
         'testthat::test_dir("undeRfactor/tests/testthat",reporter="summary",stop_on_failure=TRUE)'],
        'checks/package-tests.log')
    run(['Rscript', 'checks/order_ctrl.R'], 'checks/order-console.log')
    run(['Rscript', 'checks/score_share.R'], 'checks/score-share.log')
    run(['Rscript', 'checks/phq_table5.R'], 'checks/phq-table5.log')
print(f'Built manuscript and S1–S3 in {BUILD}. See README.md for what the checks cover.')
