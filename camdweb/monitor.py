from urllib.request import urlopen
import json

CMR = 'https://cmrdb.fysik.dtu.dk'
C2DB = 'https://c2db.fysik.dtu.dk'
BiDB = 'https://bidb.fysik.dtu.dk'
CRYSP = 'https://crysp.fysik.dtu.dk/'
OQMD12345 = 'https://oqmd12345.fysik.dtu.dk/'


def check_cmr():
    with urlopen(f'{CMR}/agau309/material/1') as fd:
        html = fd.read().decode()
    assert 'Unique ID' in html


def check_c2db():
    with urlopen(f'{C2DB}/material/1MoS2-1') as fd:
        html = fd.read().decode()
    assert 'Structure origin' in html


def check_bidb():
    with urlopen(f'{BiDB}/material/1MoS2-1') as fd:
        html = fd.read().decode()
    assert 'Space group number' in html


# def check_crysp():
#     with urlopen(f'{CRYSP}/material/{UID}') as fd:
#         html = fd.read().decode()
#     assert 'Atoms: {ATOMS}}' in html


def check_OQMD12345():
    with urlopen(f'{OQMD12345}/material/12B-620379') as fd:
        html = fd.read().decode()
    assert 'OQMD Ref. ID' in html


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
