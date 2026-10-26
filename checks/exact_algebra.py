"""Focused exact identities; this is not a complete proof audit."""
from pathlib import Path
import runpy

root = Path(__file__).resolve().parents[1]
ns = runpy.run_path(str(root / 'proofs/unequal_blocks/verify_unequal_blocks.py'))
count = 0
for name in ('scalar_symbolic_checks', 'cusp_symbolic_checks', 'extension_symbolic_checks'):
    result = ns[name]()
    assert all(value is True for value in result.values()), (name, result)
    count += len(result)
    print(name, result)
print(f'{count} unequal-block exact checks passed')
text = (root / 'proofs/theorem1/weighted_least_squares_proof.md').read_text()
code = text.split('```python\n', 1)[1].split('```', 1)[0]
exec(compile(code, 'weighted_least_squares_proof.md::checks', 'exec'), {})
