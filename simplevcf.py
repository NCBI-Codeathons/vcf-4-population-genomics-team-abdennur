from __future__ import annotations
from typing import Any

import numpy as np
import oxbow as ox
import pandas as pd
import polars as pl
import pyarrow as pa
import pysam
import bioframe


# Maps the VCF-derived input types to numpy/pandas dtypes
TYPE_MAP = {
    "Integer": "int64",
    "Float": "float64",
    "String": "object",
    "Flag": "bool",
}


def _read_info_schema(f: pysam.VariantFile) -> pd.DataFrame:
    return pd.DataFrame(
        [
            (obj.name, obj.number, obj.type, obj.description)
            for obj in f.header.info.values()
        ],
        columns=["name", "number", "type", "description"],
    )


def _read_sample_schema(f: pysam.VariantFile) -> pd.DataFrame:
    return pd.DataFrame(
        [
            (obj.name, obj.number, obj.type, obj.description)
            for obj in f.header.formats.values()
        ],
        columns=["name", "number", "type", "description"],
    )


def _read_vcf_as_records(
    f: pysam.VariantFile,
    query: str | None,
    info_fields: set[str],
    sample_fields: set[str],
    samples: list[str],
    include_unspecified: bool,
) -> list[dict[str, Any]]:
    if query is not None:
        record_iter = f.fetch(*bioframe.parse_region(query))
    else:
        record_iter = f

    records = []
    for record in record_iter:
        # Main fields
        dct = {
            "chrom": record.chrom,
            "pos": record.pos,
            "id": record.id,
            "ref": record.ref,
            "alts": list(record.alts),
            "qual": record.qual,
            "filters": list(record.filter.keys()),
        }

        # Info
        for key, value in record.info.items():
            if key in info_fields or include_unspecified:
                if isinstance(value, tuple):
                    dct[key] = list(value)
                else:
                    dct[key] = value

        # Samples
        for sample, genotype in record.samples.items():
            if sample in samples or include_unspecified:
                for key, value in genotype.items():
                    if key in sample_fields or include_unspecified:
                        if key == "GT":
                            dct[f"{sample}.GT"] = list(genotype[key])
                            dct[f"{sample}.phased"] = genotype.phased
                        elif isinstance(value, tuple):
                            dct[f"{sample}.{key}"] = list(value)
                        else:
                            dct[f"{sample}.{key}"] = value
        records.append(dct)

    return records


def read_info_schema(path: str):
    """
    Read the schema of the INFO column of a VCF.

    This currently uses pysam.

    Parameters
    ----------
    path : str
        Path to VCF file.

    Returns
    -------
    pd.DataFrame
        Dataframe with [name, number, type] columns

    Notes
    -----
    Possible values for `type` are: "Integer", "Float", "String", "Flag".

    Possible values for `number` are:
    - An integer (e.g. 0, 1, 2, 3, 4, etc.) - for fields where the number of
      values per VCF record is fixed. 0 means the field is a "Flag".
    - A string ("A", "G", "R") - for fields where the number of values per VCF
      record is determined by the number of alts, the total number of alleles,
      or the number of genotypes, respectively.
    - A dot (".") - for fields where the number of values per VCF record
      varies, is unknown, or is unbounded.
    """
    with pysam.VariantFile(path) as f:
        return _read_info_schema(f)


def read_sample_schema(path: str):
    """
    Read the schema of the genotype sample columns of a VCF.

    This currently uses pysam.

    Parameters
    ----------
    f : pysam.VariantFile

    Returns
    -------
    pd.DataFrame
        Dataframe with [name, number, type] columns

    Notes
    -----
    Possible values for `type` are: "Integer", "Float", and "String".
    """
    with pysam.VariantFile(path) as f:
        return _read_sample_schema(f)


def read_vcf_as_pandas(
    path: str,
    query: str | None = None,
    info_fields: list[str] | None = None,
    sample_fields: list[str] | None = None,
    samples: list[str] | None = None,
    include_unspecified: bool = False,
) -> pd.DataFrame:
    """
    Read a VCF into a pandas dataframe, extracting INFO and sample genotype fields.

    Parameters
    ----------
    path : str
        Path to VCF file.
    query : str, optional
        Genomic range query string. If None, all records will be read.
    info_fields : list[str], optional
        List of fields to extract from the INFO column. If None, all fields
        will be extracted.
    sample_fields: list[str], optional
        List of fields to extract from the sample genotype columns. If None,
        all fields will be extracted.
    samples : list[str], optional
        List of samples to extract. If None, all samples will be extracted.

    Returns
    -------
    pd.DataFrame
        Pandas DataFrame with columns corresponding to the requested fields.
    """
    with pysam.VariantFile(path) as f:
        info_schema = _read_info_schema(f)
        sample_schema = _read_sample_schema(f)

        if info_fields is None:
            info_fields = list(info_schema["name"])
        if sample_fields is None:
            sample_fields = list(sample_schema["name"])
        if samples is None:
            samples = list(f.header.samples)

        records = _read_vcf_as_records(
            f, query, set(info_fields), set(sample_fields), samples, include_unspecified
        )

    df = pd.DataFrame.from_records(records)

    columns = (
        ["chrom", "pos", "id", "ref", "alts", "qual", "filters"]
        + info_fields
        + [f"{sample}.{field}" for sample in samples
           for field in [f"phased", *sample_fields]]
    )
    columns = columns + list(set(df.columns) - set(columns))

    for col in columns:
        if col not in df.columns:
            df[col] = pd.NA

    return df[columns].fillna(np.nan)


def read_vcf_as_polars(
    path: str,
    query: str | None = None,
    info_fields: list[str] | None = None,
    sample_fields: list[str] | None = None,
    samples: list[str] | None = None,
    include_unspecified: bool = False,
) -> pd.DataFrame:
    """
    Read a VCF into a polars dataframe, extracting INFO and sample genotype fields.

    Parameters
    ----------
    path : str
        Path to VCF file.
    query : str, optional
        Genomic range query string. If None, all records will be read.
    info_fields : list[str], optional
        List of fields to extract from the INFO column. If None, all fields
        will be extracted.
    sample_fields: list[str], optional
        List of fields to extract from the sample genotype columns. If None,
        all fields will be extracted.
    samples : list[str], optional
        List of samples to extract. If None, all samples will be extracted.

    Returns
    -------
    pl.DataFrame
        Polars DataFrame with columns corresponding to the requested fields.
    """
    with pysam.VariantFile(path) as f:
        info_schema = _read_info_schema(f)
        sample_schema = _read_sample_schema(f)

        if info_fields is None:
            info_fields = list(info_schema["name"])
        if sample_fields is None:
            sample_fields = list(sample_schema["name"])
        if samples is None:
            samples = list(f.header.samples)

        records = _read_vcf_as_records(
            f, query, set(info_fields), set(sample_fields), samples, include_unspecified
        )

    df = pl.from_records(records)

    columns = (
        ["chrom", "pos", "id", "ref", "alts", "qual", "filters"]
        + info_fields
        + [f"{sample}.{field}" for sample in samples
           for field in [f"phased", *sample_fields]]
    )
    columns = columns + list(set(df.columns) - set(columns))

    for col in columns:
        if col not in df.columns:
            df = df.with_columns(pl.lit(None).alias(col))

    return df.select(columns)
    # return df[columns].fillna(np.nan)
