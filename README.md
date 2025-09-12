<a href="https://jinghuazhao.github.io/Caprion/">
  <img src="https://jinghuazhao.github.io/Caprion/qrcode.png" height="200" width="200" align="right" alt="Caprion QR Code">
</a>

# Caprion Data Analysis

### Welcome!

This repository/site is dedicated to protein/peptide quantitative trait analysis using the **Caprion platform**, organized chronologically and logistically into the following sections.

## 📊 Pilot Studies

- [Autoencoder](pilot/autoencoder)
- [GWAS Round 2](pilot/gwas2)
- [General Pilot Overview](pilot)

## 🔬 Analysis

- [Protein-level Analysis](progs)
- [Peptide-level Analysis](peptide_progs)
- [Miscellaneous Analysis](misc)

## 📎 Additional Information

- [Caprion Panel Description](https://jinghuazhao.github.io/pQTLdata/reference/caprion.html)
- [Project Notes](https://jinghuazhao.github.io/Caprion/Notes/)

## 🌐 Local File/Web Browsing

A local web-style navigation can be set up as follows:

```bash
cd /rds/project/rds-zuZwCZMsS0w/Caprion_proteomics/analysis
# As of 07/09/2025
module load ceuadmin/firefox/nightly
python3 -m http.server &
firefox http://127.0.0.1:8000 &
```

📌 **Note:** 8000 is a port number[^port]. One can also use other browsers such as Microsoft Edge or Google Chrome[^browsers].

## 🗂️ Web Page Structure

One can browse local mirrors and resources via the home page (`index.html`) served at the chosen port:

- **SRCF mirror** — under `/srcf`
- **Colocalisation view** — `/json/coloc.html` (hg19 positions)
- **Multiprotein isotope mappings** — `/dup/json/dup.htm`
- **Supplementary tables**
- **Caprion site** — from `/site`

## 🔐 Remote (Non-CSD3) Access

Set up SSH tunneling to access the local web server from another machine:

### 1. On CSD3

```bash
python3 -m http.server 8000 &
hostname
```

### 2. On **local machine**

```bash
ssh -4 -L 8080:127.0.0.1:8000 -fN CRSid@${hostname}.hpc.cam.ac.uk
```

Then visit: <http://127.0.0.1:8080>

Ensure `${hostname}` matches the result from CSD3 `hostname`.

[^port]: **Use of specific port number**

    One can replace port number `8000` with any free port whose availability can be handled as follows,

        export pn=8000
        if lsof -i :${pn}; then
            echo "Port ${pn} is already in use. Try another one."
        else
            python3 -m http.server ${pn} &
            server_pid=$!
            firefox http://127.0.0.1:${pn} &
        fi

    Release the port using:

        kill $server_pid

    Check active processes with `ps`.

[^browsers]: **Other browsers**

    One can launch Edge using a module or a temporary user data directory if needed:

        module load ceuadmin/edge
        # ~/.config/microsoft-edge
        edge http://127.0.0.1:${pn} &

    If the Edge config is not available, use a temporary directory:

        edge --user-data-dir=/tmp/edge http://127.0.0.1:${pn} &

    Google Chrome follows the Microsoft Edge syntax.
