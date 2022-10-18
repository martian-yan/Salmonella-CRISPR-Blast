# Salmonella-WGCT-pipeline

## Dependencies

- blast+
- bedtools
- Python >= 3.9
- Python modules: (will be install by `pip` automatically)
    - biopython
    - pandas

## Installation

### Install from source

Prepare a proper environment with the necessary dependencies, clone this project and install:

```
git clone https://github.com/martian-yan/Salmonella-WGCT-pipeline.git
cd Salmonella-WGCT-pipeline
python -m build
pip install dist/SWGCT-1.0.0-py3-none-any.whl
```

If anything goes wrong, check if you have the latest version of Pyhton `setuptools` and `build`:

```
pip install --upgrade setuptools
pip install --upgrade build
```

## Usage

TBC