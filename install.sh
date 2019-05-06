set -e

python setup.py install
python --version
python -c "\
try:
    import tgirt_smRNA_correction
except ImportError:
    pass
"
python -c "import multiprocessing as mp; print('%d CPUs' % mp.cpu_count())"
