# 4-3-2020 JHZ

# recycle tensorQTL venv
# ln -s $HOME/tensorqtl/venv
# pip install tensorqtl

module load python/3.6
virtualenv --system-site-package venv
source venv/bin/activate
python -m pip install rpy2
hostname
jupyter notebook --ip=127.0.0.1 --no-browser --port 8087
# ssh -4 -L 8087:127.0.0.1:8087 -fN login-e-10.hpc.cam.ac.uk
firefox http://127.0.0.1:8087/?token=d991ea12ce42b216d3aacd3c573e95280b6cd30d4b4aeeed &

