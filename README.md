# crys-gsu

## src/split.py - script for splitting grid into several *.cry files

```
[Overview]:
split.py script reads *.dat file with surface grid and split it into several *.cry files

[Usage]:
split.py <grid-file> <output-dir> <split-strategy> <ais-zone-name-1> ... <ais-zone-name-n>
    <grid-file> - full name of grid file
    <output-dir> - directory name for output files
        if there is no such directory, it will be created
    <split-strategy>:
        h<n> - hierarchical split into 2^n zones in addition to AIS zones
    <ais-zone-name-i> - names of AIS zones

[Examples]:
bunny.dat -> bunny_<mpi_i>_000000000000.cry
bunny_000000000100.dat -> bunny_<mpi_i>_000000000100.cry
```

## src/merge.py - script for merging *.dat file with *.txt data files

```
[Overview]:
merge.py script merge initial *.dat grid file with *.txt data files
         and produces *.dat files with timestamps

[Usage]:
merge.py <grid-file> <txt-files-dir> <r-files-dir>
    <grid-file> - initial file name
    <txt-files-dir> - name of directory with data files
    <r-files-dir> - name of directory with output *_r_* files

[Examples]:
            [bunny_00000_000000000100.txt bunny_00001_000000000100.txt] -> [bunny_r_000000000100.dat]
bunny.dat + [bunny_00000_000000000200.txt bunny_00001_000000000200.txt] -> [bunny_r_000000000200.dat]
            [bunny_00000_000000000300.txt bunny_00001_000000000300.txt] -> [bunny_r_000000000300.dat]
```

## Building documentation

Documentation is foramted as [Sphinx reStructuredText](https://www.sphinx-doc.org/en/master/usage/restructuredtext/index.html).

To build documentation in HTML you should run:
```
cd docs
make html
```

## Testing

Run `python -m unittest` in project's top directory to execute tests.
