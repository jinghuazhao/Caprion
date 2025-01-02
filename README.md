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

A web-style navigation is furnised as follows,

```bash
cd /rds/project/rds-zuZwCZMsS0w/Caprion_proteomics/analysis
module load ceuadmin/edge
python3 -m http.server &
edge http://127.0.0.1:8000 &
```
where 8000 is a port number (pn).

In case the browser does not show, use

`edge --user-data-dir=/tmp/edge http://127.0.0.1:8000 &`

'/tmp/edge' directory in replace of `~/.config/microsoft-edge`.

One could browse mirrors of two web sites as well as files,

1. Caprion site. This is from `/site` as above.
2. SRCF site. The mirror is within the following subdirectory: `/srcf`.
3. Colocalisation. See /json/coloc.html. Chromosomal positions are in hg19.
4. Multiprotein mapping isotopes, /dup/json/dup.htm
5. Supplementary tables.

To facilitate navigation, an `index.html` is created in place, so `python3 -m http.server 8000 &` starts a home page for 1-5 above.

In case pn is already in use, a different one can be chosen as follows,

```bash
export pn=8000
if lsof -i :${pn}; then
    echo "Port ${pn} is already in use' try another one."
else
    python3 -m http.server ${pn} &
    server_pid=$!
    edge http://127.0.0.1:${pn} &
fi
```

and pn can be released with `kill $server_pid` (can be checked with `ps`).

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
