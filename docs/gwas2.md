## gwas2

This is a promising alternative showing through RCN3/FCGRN with [gwas2.sh](utils/gwas2.sh).
```mermaid
graph TB;
utils/gwas2.sh --> utils/gwas.do
utils/gwas2.sh --> utils/gwas2.do
```
where `gwas.do` (`caprion.dat` also contains _invn data) and `gwas2.do` (`gwas2_invn.do` for _invn data) are for the pilot and batch 2 data, respectively.

See [gwas2](https://jinghuazhao.github.io/gwas2/) repository for additional information.
