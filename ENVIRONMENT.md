# Environment Compatibility

## Validation environment

All published disease-day results were produced on:

- **Hardware:** Raspberry Pi 5 (8GB RAM, ARM64)
- **OS:** Raspberry Pi OS (Debian bookworm), kernel 6.12.x
- **Python:** 3.13.x
- **Package manager:** pip with piwheels (ARM binary wheels)
- **Exact pins:** See `requirements-lock.txt`

## NumPy version note

`pyproject.toml` pins `numpy>=1.24,<2`. The validation environment uses
`numpy==1.26.4` via piwheels ARM wheels.

NumPy 1.26.4 officially supports Python 3.9-3.12. Python 3.13 support was
added in NumPy 2.1.0. The validation environment works because piwheels
provides a compatible ARM binary build of numpy 1.26.4 for Python 3.13.

**This may not reproduce on x86_64 with Python 3.13.** If numpy 1.26.4
fails to install on Python 3.13 in your environment, use Python 3.12.

## Recommended reproduction environment

For public reproduction, we recommend:

- **Python 3.11 or 3.12** (officially supported by all pinned dependencies)
- `pip install -e ".[clustering]"` for full pipeline functionality
- Or `pip install -r requirements-lock.txt` for exact validation match

## hdbscan ARM compatibility

The `numpy<2` pin exists because the piwheels ARM binary for `hdbscan`
is compiled against numpy 1.x. Importing hdbscan with numpy 2.x causes
an ABI crash on ARM. This constraint does not apply to x86_64 environments
where hdbscan wheels are compiled for numpy 2.x.

## CI test matrix

GitHub Actions tests on Python 3.11, 3.12, and 3.13 (x86_64).
The numpy<2 pin resolves successfully on all three via PyPI wheels.
