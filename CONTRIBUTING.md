# Contributing to GenomeViz

Thank you for your interest in contributing to GenomeViz! This document provides guidelines for contributing to the project.

## Code of Conduct

Be respectful and constructive in all interactions. We want GenomeViz to be welcoming to contributors of all backgrounds and skill levels.

## How Can I Contribute?

### Reporting Bugs

Before submitting a bug report:
1. Check the [Issues](https://github.com/Aaron-Thiel/GenomeViz/issues) page
2. Make sure you're using the latest version
3. Collect information about the bug:
   - Your operating system
   - Python version (`python --version`)
   - GenomeViz version
   - Full error message
   - Input file formats and sizes

Submit bugs via GitHub Issues with the "bug" label.

### Suggesting Enhancements

Enhancement suggestions are welcome! Please:
1. Check existing suggestions first
2. Explain your use case clearly
3. Describe the proposed solution
4. Consider implementation complexity

Submit via GitHub Issues with the "enhancement" label.

### Pull Requests

1. Fork the repository
2. Create a new branch (`git checkout -b feature/amazing-feature`)
3. Make your changes
4. Test thoroughly
5. Commit your changes (`git commit -m 'Add amazing feature'`)
6. Push to your branch (`git push origin feature/amazing-feature`)
7. Open a Pull Request

#### Pull Request Guidelines

- Follow existing code style
- Add docstrings to new functions
- Update documentation if needed
- Test with example data
- Keep PRs focused on a single feature/fix

## Development Setup

```bash
# Clone your fork
git clone https://github.com/Aaron-Thiel/GenomeViz.git
cd GenomeViz

# Install dependencies
pip install -r requirements.txt

# Run tests (if available)
python -m pytest tests/

# Try the example
python genomeViz.py \
  --reference examples/input/reference.fna \
  --assembly examples/input/sample.fna \
  --gff examples/input/reference.gff3 \
  --output test_output/
```

## Code Style

- Follow PEP 8 guidelines
- Use meaningful variable names
- Add comments for complex logic
- Keep functions focused and short
- Use type hints where helpful

Example:
```python
def calculate_quality_score(coverage_pct: float, identity: float) -> float:
    """
    Calculate gene quality score.
    
    Args:
        coverage_pct: Coverage percentage (0-100)
        identity: Sequence identity (0-100)
    
    Returns:
        Quality score (0-100)
    """
    return (coverage_pct * 0.3) + (identity * 0.7)
```

## Testing

- Test with various genome sizes (small plasmids to large chromosomes)
- Test with different assemblers (SPAdes, Flye, etc.)
- Test with different annotators (Bakta, Prokka, PGAP)
- Check edge cases (no alignments, all gaps, etc.)

## Documentation

- Update README.md for new features
- Update docs/USAGE.md for usage changes
- Add examples when appropriate
- Keep CHANGELOG.md current

## Priority Areas

We especially welcome contributions in these areas:

1. **Performance optimization**
   - Faster processing for large genomes
   - Memory optimization

2. **Visualization improvements**
   - Additional plot types
   - Customizable colors/themes
   - Export to more formats

3. **Feature additions**
   - Batch comparison mode
   - Additional quality metrics
   - Integration with databases

4. **Testing**
   - Unit tests
   - Integration tests
   - CI/CD setup

5. **Documentation**
   - Tutorial videos
   - More examples
   - Troubleshooting guides

## Questions?

- Open a Discussion on GitHub
- Email [aaron.chris.thiel@gmail.com]
- Check existing Issues

## License

By contributing, you agree that your contributions will be licensed under the MIT License.

---

Thank you for helping make GenomeViz better! 🧬
