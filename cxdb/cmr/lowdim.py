from collections import namedtuple
from itertools import combinations
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
import numpy as np

from cxdb.utils import FormPart, table
from cxdb.panels.panel import Panel
from cxdb.material import Material, Materials


Desc = namedtuple('Desc', ['short', 'long', 'unit'])


def get_combinations(dims: list[int]) -> Iterable[tuple[int, ...]]:
    return (x for y in (combinations(dims, i)
                        for i in range(1, len(dims) + 1)) for x in y)


def sab_key_descriptions() -> dict[str, Desc]:
    combs = list(get_combinations([0, 1, 2, 3]))
    keys = []
    for x in ['s', 'a', 'b']:
        for comb in combs:
            combstr = ''.join([str(d) for d in comb])
            keys.append('{}_{}'.format(x, combstr))

    key_descriptions = {}
    # add description of s, a, and b keys
    labels = {'s': ('', 'score'),
              'a': ('Start of', 'k-interval'),
              'b': ('End of', 'k-interval')}
    for k in keys:
        x, numbers = k.split('_')
        beg, end = labels.get(x, (None, None))
        if beg is None and end is None:
            continue
        else:
            assert beg is not None
            assert end is not None
            beg = beg + ' '
            end = ' ' + end
        desc = beg + '+'.join([n + 'D' for n in numbers]) + end
        key_descriptions[k] = Desc(short=desc, long='', unit='')

    return key_descriptions


def score(material, path: Path) -> None:
    vs = []
    for k in keysforfigure:
        v = material.get(k, 0)
        vs.append(v)
    x = np.arange(len(keysforfigure))
    fig, ax = plt.subplots()
    ax.bar(x, vs)
    ax.set_xticks(x)
    prettykeys = ['${}_{{{}}}$'.format(*k.split('_')) for k in keysforfigure]
    ax.set_xticklabels(prettykeys, rotation=90)
    ax.set_ylabel('Score')
    plt.tight_layout()
    plt.savefig(path)
    plt.close()


def all_key_descriptions():
    key_descriptions = sab_key_descriptions()
    # add description of number of dimensional components
    for d in range(4):
        desc = 'Number of {}D components'.format(d)
        key = 'numc_{}'.format(d)
        key_descriptions[key] = Desc(short=desc, long='', unit='')
    # add description of external db ids
    key_descriptions['source'] = Desc('Source', 'Database source', '')
    key_descriptions['dbid'] = Desc('ID #', 'Database ID Number', '')
    key_descriptions['doi'] = Desc('DOI', '', '')
    key_descriptions['publication'] = Desc('Publication', '', '')
    key_descriptions['spacegroup_number'] = Desc('Space group #',
                                                 'Space group Number', '')

    key_descriptions['dimtype'] = Desc(
        'Dimensionality',
        'Dimensionality with highest score', '')
    key_descriptions['h'] = Desc(
        'Component count',
        '# components with each dimensionality type', '')
    key_descriptions['warning'] = Desc('Warning',
                                       'Potential problems with structure \
    (comma separated)', '')
    return key_descriptions


all_keydescs = all_key_descriptions()
keysfortable0 = ['source', 'spacegroup_number', 'publication']
keysfortable2 = [k for k in all_keydescs if k.startswith('numc')]
keysforfigure = [k for k in all_keydescs if k.startswith('s_')]


HTML = """
<div class="row">
  <div class="col-6">
    {}
  </div>
  <div class="col-6">
    {}
  </div>
</div>
"""


class LowDimPanel(Panel):
    title = 'Dimensionality analysis'

    def get_html(self,
                 material: Material,
                 materials: Materials) -> tuple[str, str]:
        uid = material.uid
        path = material.folder / f'lowdim/{uid}.png'
        path.parent.mkdir(exist_ok=True)
        if not path.is_file():
            score(material, path)

        col1 = table(['Item', ''],
                     materials.table(material, keysfortable2))
        col2 = f'<img alt="Dim. analysis for {uid}" src="/abs3/png/{uid}" />'
        return HTML.format(col1, col2), ''


class LowDimRange(FormPart):
    def __init__(self):
        super().__init__('Score range', 's')

    def render(self, query: dict) -> str:
        fro = query.get('from_s', '')
        to = query.get('to_s', '')
        parts = [
            f'<label class="form-label">{self.text}</label>',
            '<input',
            '  class="form-control"',
            '  type="text"',
            '  name="from_s"',
            f'  value="{fro}" />',
            '<input',
            '  class="form-control"',
            '  type="text"',
            '  name="to_s"',
            f'  value="{to}" />']
        selection = query.get('s')
        parts += [f'<label class="form-label">{self.text}</label>\n'
                  f'<select name="s" class="form-select">']
        for txt, val in [('0D', 's_0'),
                         ('1D', 's_1'),
                         ('2D', 's_2'),
                         ('3D', 's_3'),
                         ('0D+1D', 's_01'),
                         ('0D+2D', 's_02'),
                         ('0D+3D', 's_03'),
                         ('1D+2D', 's_12'),
                         ('1D+3D', 's_13'),
                         ('2D+3D', 's_23'),
                         ('0D+1D+2D', 's_012'),
                         ('0D+1D+3D', 's_013'),
                         ('1D+2D+3D', 's_123'),
                         ('0D+1D+2D+3D', 's_0123')]:
            selected = ' selected' if selection == val else ''
            parts.append(f'  <option value="{val}"{selected}>{txt}</option>')
        parts.append('</select>')
        return '\n'.join(parts)

    def get_filter_strings(self, query: dict) -> list[str]:
        filters = []
        fro = query.get('from_s', '')
        if fro:
            filters.append(f's>={fro}')
        to = query.get('to_s', '')
        if to:
            filters.append(f's<={to}')
        return filters
