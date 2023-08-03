from __future__ import annotations
from io import BytesIO

import numpy as np
import oxbow as ox
import pandas as pd
import polars as pl
import pyarrow as pa
import pysam
import bioframe


# This dictionary maps the VCF-derived input types to numpy/pandas dtypes
TYPE_MAP = {
    "Integer": "int64",
    "Float": "float64",
    "String": "object",
    "Flag": "bool",
}


def read_info_schema(f: pysam.VariantFile) -> pd.DataFrame:
    """
    Read the schema of the INFO column of a VCF.

    This currently uses pysam.

    Parameters
    ----------
    f: pysam.VariantFile

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
    return pd.DataFrame(
        [(obj.name, obj.number, obj.type) for obj in f.header.info.values()],
        columns=["name", "number", "type"],
    )


def read_sample_schema(f: pysam.VariantFile) -> pd.DataFrame:
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
    """
    return pd.DataFrame(
        [(obj.name, obj.number, obj.type) for obj in f.header.formats.values()],
        columns=["name", "number", "type"],
    )


def _read_vcf_as_records(f, query, info_fields, sample_fields, samples, include_unspecified):
    if query is not None:
        query = f.fetch(*bioframe.parse_region(query))
    else:
        query = f

    records = []
    for record in query:
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

        # Items
        for key, value in record.info.items():
            if key in info_fields:
                if isinstance(value, tuple):
                    dct[key] = list(value)
                else:
                    dct[key] = value
            elif include_unspecified:
                dct[key] = value

        # Samples
        for sample, genotype in record.samples.items():
            if sample in samples:
                for key, value in genotype.items():
                    if key in sample_fields:
                        if key == "GT":
                            dct[sample] = list(genotype[key])
                            dct[f"{sample}.phased"] = genotype.phased
                        elif isinstance(value, tuple):
                            dct[f"{sample}.{key}"] = list(value)
                        else:
                            dct[f"{sample}.{key}"] = value
                    elif include_unspecified:
                        dct[f"{sample}.{key}"] = value
        records.append(dct)

    return records


def read_vcf_as_pandas(
    path: str,
    query: str | None = None,
    info_fields: list[str] | None = None,
    sample_fields: list[str] | None = None,
    samples: list[str] | None = None
) -> pd.DataFrame:

    with pysam.VariantFile(path) as f:

        info_schema = read_info_schema(f)
        sample_schema = read_sample_schema(f)

        if info_fields is None:
            info_fields = set(info_schema["name"])
        if sample_fields is None:
            sample_fields = set(sample_schema["name"])
        if samples is None:
            samples = list(f.header.samples)

        records = _read_vcf_as_records(f, query, info_fields, sample_fields, samples)

    return pd.DataFrame.from_records(records)


def read_vcf_as_polars(
    vcf_path: str,
    query: str | None = None,
    info_fields: list[str] | None = None,
    sample_fields: list[str] | None = None,
    samples: list[str] | None = None
) -> pd.DataFrame:

    with pysam.VariantFile(vcf_path) as f:

        info_schema = read_info_schema(f)
        sample_schema = read_sample_schema(f)

        if info_fields is None:
            info_fields = set(info_schema["name"])
        if sample_fields is None:
            sample_fields = set(sample_schema["name"])
        if samples is None:
            samples = list(f.header.samples)

        records = _read_vcf_as_records(f, query, info_fields, sample_fields, samples)

    return pl.from_records(records)
