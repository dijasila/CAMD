import subprocess


def coverage():  # pragma: no cover
    subprocess.run(
        'coverage run --branch -m pytest -v && '
        'coverage report && '
        'coverage html', shell=True)
