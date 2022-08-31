Little script to take MSMC2 file format and convert an unfolded SFS.
Assumes alleles from at least one outgroup are available at each input locus.  
Unit tested, but not field-tested. Use at your own risk. 

See `python msmc2_to_sfs.py --help` for usage details

For demo, run:

```bash
python msmc2_to_sfs.py \
  --msmc2_file example/example_msmc2.txt \
  --outgroup_index 24 27 \
  --allele_total 28
```

To confirm all unit tests pass, run `pytest msmc2_to_sfs.py`
