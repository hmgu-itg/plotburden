from pathlib import Path


__version__ = '0.0.1'


DEPENDENCIES = {
    'tabix',
    'plink'
}

_BASE_DIR = Path(__file__).parent

GET_LD_SH = str(_BASE_DIR / 'getld.sh')