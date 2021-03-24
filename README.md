# PurduePHYS324-LSST
## SETTING UP THE ENVIRONMENT

### Module Dependency

```bash
pip install -r requirement.txt
```

### IDE

VSCode with Python Extension is optimal, because it allows separating blocks using `#%%`, making the code easy to manipulate and test.

Jupyter Notebook is great for preserving outputs. The Interactive Shell in VSCode can be exported to a jupyter notebook.

### Testing

No test cases for now. Maybe you can try running `example.py` as a test.

## Development

### Code

Each task (for example, the GPR model) is realized in a single `.py` file in the root directory.
By doing this, the package `Supernova` can be directly imported without having to install it.

If something is replicable, it can be refactored into the package `Supernova`.

### Assets

All assets are organized in `src`.

## Documentation

Documentation for package `Supernova` is at `.doc`.
