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
cd /rds/project/rds-zuZwCZMsS0w/Caprion_proteomics/analysis

module load ceuadmin/edge

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
3. Isotopes associated with >1 proteins, srcf/peptides/dup/json/dup.htm

## Non-CSD3 browser(s)

This approach seems less problematic with `user-data-dir` mentioned above. We can again set up tunneling from CSD3 with

```bash
python3 -m http.server 8000 &
hostname
```

Once succeeded, we establish the connection elsewhere.

```bash
ssh -4 -L 8080:127.0.0.1:8000 -fN jhz22@${hostname}.hpc.cam.ac.uk
```

where hostname from CSD3 and ${hostname} have to be the same. We can then browse `http://127.0.0.1:8080`.
