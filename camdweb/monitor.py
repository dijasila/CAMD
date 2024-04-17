from urllib.request import urlopen
import json

CMR = 'https://cmrdb.fysik.dtu.dk'
C2DB = 'https://c2db.fysik.dtu.dk'


def check_cmr():
    with urlopen(f'{CMR}/agau309/material/2/download/xyz') as fd:
        xyz = fd.read().decode()
    assert xyz.startswith('309\n')


def check_optimade():
    with urlopen(f'{C2DB}/optimade/structures?filter=nelements=4') as fd:
        dct = json.loads(fd.read().decode())
    assert dct['meta']['data_returned'] == 20


if __name__ == '__main__':
    ok = True
    for check in [check for name, check in globals().items()
                  if name.startswith('check_')]:
        print(check)
        try:
            check()
        except Exception as ex:
            print(ex)
            ok = False
    assert ok
