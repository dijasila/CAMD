from pathlib import Path

fro = Path('/home/niflheim/joaso/paper_calcs/htp-spiral_paper/n1-magnets/tree')
to = Path('/home/niflheim2/cmr/C2DB-ASR/tree')


def copy_spinspiral_data_to_tree():
    for p1 in fro.glob('A*/*/*/'):
        p2 = to.joinpath(p1.parts[-3:])
        assert p2.is_dir()
        for name in ['collect_spiral',
                     'dmi',
                     'spinorbit',
                     'spinorbit@calculate']:
            f1 = p1 / f'results-asr.{name}.json'
            if f1.is_file():
                f2 = p2 / f1.name
                assert not f2.is_file()
                if 0:
                    f2.write_bytes(f1.read_bytes())
                else:
                    print(f1, f2)
