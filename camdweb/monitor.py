from urllib.request import urlopen


def check_cmr():
    CMR = 'https://cmrdb.fysik.dtu.dk'
    with urlopen(f'{CMR}/agau309/material/2/download/xyz') as fd:
        xyz = fd.read().decode()
    assert xyz.startswith('309\n')


if __name__ == '__main__':
    check_cmr()
