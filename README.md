# Transforming VCFs for Data Science


List of participants and affiliations:
- Nezar Abdennur, UMass Chan Medical School (Team Leader)
- David Adeleke
- Clark Cucinell
- Se-Ran Jun
- Mehmet Kuscuoglu
- Lei Ma
- Garrett Ng, UMass Chan Medical School
- Thomas Reimonn, UMass Chan Medical School

## Project Goals
Facilitate data analysis by mapping VCF data into a flat schema where the information nested within columns is parsed into new columns

## Approach
We will use [oxbow](https://github.com/abdenlab/oxbow) to read VCF files into Apache Arrow IPC and then parse out nested data into new columnar representations. The resulting data can be manipulated interactively in a Jupyter Notebook environment or dumped into other data formats such as Parquet.

[Daily Logs](https://hackmd.io/hG1YtJvcSiiHr7HplDAJQA)

## Results
- Proof-of-concept accomplished
- Implemented extraction for VCF INFO fields
- Implemented extraction for VCF sample genotype fields (FORMAT + samples)
- Support for extracting into Pandas and Polars DataFrames in Python
- Refactored notebook code into a Python module

## Future Work
- Support for multiple backends (oxbow, pysam, cyvcf2)
- Convenient integrations for out-of-core and distributed computing
  - Dask
  - Polars LazyFrame
- Release as an open-source Python package on PyPI
