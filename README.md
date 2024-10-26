<a href="https://jinghuazhao.github.io/Caprion/"><img src="https://jinghuazhao.github.io/Caprion/qrcode.png" height=200 width=200 align="right"></img></a>
# Caprion data analysis

### Welcome!

This repository/site is dedicated to protein/peptide quantitative trait analysis using the Caprion platform, which is organised chonologically/logistically into the following sections.

## Pilot studies

- [autoencoder](pilot/autoencoder)
- [gwas2](pilot/gwas2)
- [Pilot studies](pilot/)

## Analysis

- [Protein analysis](progs/)
- [Peptide analysis](peptide_progs)
- [Miscellaneous analysis](misc/)

## Additional information

- [Caprion panel](https://jinghuazhao.github.io/pQTLdata/reference/caprion.html)
- [Notes](https://jinghuazhao.github.io/Caprion/Notes/)

## Mirror of SRCF

There are available in two ways, both from the following directory:

`/rds/project/rds-zuZwCZMsS0w/Caprion_proteomics/analysis/srcf`

### File exploration

This allows for examination of each file

```bash
module load ceuadmin/chrome
chrome index.html &
```

### Web-style navigation

This is more involved and via a port number, e.g., 8000,

```bash
export pn=8000
if lsof -i :${pn}; then
    echo "Port ${pn} is already in use."
else
    module load ceuadmin/chrome
    module load python/3.8.11/gcc/pqdmnzmw
    python -m http.server ${pn} &
    server_pid=$!
    chrome http://localhost:${pn} &
fi
```

where the port number can be released with `kill $server_pid`.
