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

## Local file/web browsing

A web-style navigation is furnised via a port number, e.g., 8000,

```bash
module load ceuadmin/edge

cd /rds/project/rds-zuZwCZMsS0w/Caprion_proteomics/analysis

export pn=8000
if lsof -i :${pn}; then
    echo "Port ${pn} is already in use' try another one."
else
    python3 -m http.server ${pn} &
    server_pid=$!
    edge http://localhost:${pn} &
fi
```

where the port number can be released with `kill $server_pid`. In case it does now show, use

`edge --user-data-dir=${TMPDIR} http://localhost:${pn} &`

where TMPDIR is a directory name.

One could browse files as well as mirrors of two web sites.

1. SRCF. The mirror is within the following subdirectory: `/srcf`.
2. Web site. This is from `/site` as above.
