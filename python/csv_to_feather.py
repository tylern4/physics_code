#!/usr/bin/env python

import pandas as pd
import pyarrow as pa
from pyarrow import csv
from pyarrow import feather

import sys
import time

if len(sys.argv) <= 1:
    exit(1)


file_name = str(sys.argv[1])

names = ["electron_sector", "w", "q2", "theta", "phi", "mm2", "helicty", "type", "hash"]
dtype = {
    "electron_sector": "int8",
    "helicty": "int8",
    "w": "float32",
    "q2": "float32",
    "theta": "float32",
    "phi": "float32",
    "mm2": "float32",
}

start = time.time()
# df = pd.read_csv(file_name, names=names, index_col=False, dtype=dtype)
pyTable = csv.read_csv(
    file_name,
    read_options=csv.ReadOptions(use_threads=True, column_names=names),
    convert_options=csv.ConvertOptions(column_types=dtype),
)
stop = time.time()
print(f"read time: {stop - start}")

output_name = file_name[:-3] + "feather"
print(output_name)

start = time.time()
# df.to_feather(output_name)
feather.write_feather(pyTable.to_pandas(), output_name)
stop = time.time()
print(f"write time: {stop - start}")
