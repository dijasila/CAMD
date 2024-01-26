import subprocess
import sys
from pathlib import Path
from time import time


def main(folder: Path) -> None:  # pragma: no cover
    while True:
        t = time()
        subprocess.run(f'python3 -m camdweb.web {folder}/*/', shell=True)
        t = time() - t
        if t < 20:
            return
        subprocess.run('git pull'.split())


if __name__ == '__main__':  # pragma: no cover
    folder = Path(sys.argv[1])
    main(folder)
