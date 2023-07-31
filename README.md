# Team Project Name

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

## Results

## Future Work
