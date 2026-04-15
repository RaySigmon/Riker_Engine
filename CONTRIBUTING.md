# Contributing to Riker Engine

Thank you for your interest in contributing to the Riker Engine.

## How to Contribute

1. **Report bugs** — Open an issue on GitHub with steps to reproduce.
2. **Suggest features** — Open an issue describing the use case and proposed solution.
3. **Submit code** — Fork the repo, create a branch, make your changes, and open a pull request.

## Development Setup

```bash
git clone https://github.com/RaySigmon/Riker_Engine.git
cd Riker_Engine
pip install -e ".[clustering]"
python -m pytest tests/ -q  # 300 tests must pass
```

## Running Tests

```bash
# Full suite
python -m pytest tests/ -q

# Specific layers
python -m pytest tests/test_stats.py -v       # Statistical primitives
python -m pytest tests/test_ingestion.py -v    # Data loading
python -m pytest tests/test_phases.py -v       # Pipeline phases
python -m pytest tests/test_qc.py -v           # QC gates

# With coverage
python -m pytest tests/ --cov=riker --cov-report=term-missing
```

All tests must pass before submitting a pull request.

## Code Style

- Python 3.11+ with type hints
- No unnecessary dependencies
- Tests for all new functionality
- Keep functions focused and well-named — comments only where logic isn't self-evident

## License

By contributing, you agree that your contributions will be licensed under the [AGPL-3.0 License](LICENSE).

A Contributor License Agreement (CLA) may be required for substantial contributions. This ensures the project maintainer can manage the licensing terms if needed (e.g., for dual-licensing to support commercial integrations while keeping the open-source version AGPL).

## Questions?

Open an issue or reach out via GitHub.
