from __future__ import annotations

import abc
from pathlib import Path
from typing import Iterable, Sequence
from ase.formula import Formula


def table(header: list[str] | None, rows: Sequence[Iterable],
          responsive: bool = True) -> str:
    """Create HTML table.

    Example:

    === =====
    A   B
    === =====
    1.2 hello
    2.0 hi!
    === =====

    >>> print(table(['A', 'B'], [[1.2, 'hello'], [2.0, 'hi!']]))
    <div class="table-responsive">
     <table class="table table-striped">
      <thead>
       <tr>
        <th>A</th>
        <th>B</th>
       </tr>
      </thead>
      <tbody>
       <tr>
        <td>1.2</td>
        <td>hello</td>
       </tr>
       <tr>
        <td>2.0</td>
        <td>hi!</td>
       </tr>
      </tbody>
     </table>
    </div>
    """
    if header is None:
        head = ''
    else:
        head = ('  <thead>\n   <tr>\n    <th>' +
                '</th>\n    <th>'.join(header) +
                '</th>\n   </tr>\n  </thead>\n')

    html_table = (
        f' <table class="table table-striped">\n{head}' +
        '  <tbody>\n   <tr>\n    ' +
        '\n   </tr>\n   <tr>\n    '.join(
            '\n    '.join(f'<td>{x}</td>' for x in row)
            for row in rows) +
        '\n   </tr>\n  </tbody>\n </table>')

    return f'<div class="table-responsive">\n{html_table}\n</div>' \
        if responsive else html_table


def image(path: Path | str, alt=None) -> str:
    """Create <img> tag.

    >>> image('abc/def.png', alt='Short description')
    '<img alt="Short description" src="/png/abc/def.png" class="img-fluid">'
    """
    return f'<img alt="{alt or path}" src="/png/{path}" class="img-fluid">'


class FormPart(abc.ABC):
    def __init__(self, text: str, name: str):
        self.text = text
        self.name = name

    @abc.abstractmethod
    def render(self) -> str:
        raise NotImplementedError

    def get_filter_strings(self, query: dict) -> list[str]:
        val = query.get(self.name, '')
        if val:
            return [f'{self.name}={val}']
        return []


class Select(FormPart):
    def __init__(self, text, name, options, names=None, default=None):
        super().__init__(text, name)
        self.options = options
        self.names = names
        self.default = default or options[0]

    def render(self) -> str:
        """Render select block.

        >>> s = Select('Bla-bla', 'xyz', ['A', 'B', 'C'])
        >>> html = s.render()
        >>> s.get_filter_strings({'xyz': 'C'})
        ['xyz=C']
        """
        parts = [
            '<form>',
            '<div class="form-group row pb-1">',
            f'<label for="form_{self.name}" class="col-4 col-form-label-sm">'
            f'  {self.text}</label>'
            '<div class="col d-flex align-items-center">',
            f'<select class="form-select" name="{self.name}" id="{self.name}">'
        ]
        names = self.names or self.options
        for val, txt in zip(self.options, names):
            selected = ' selected' if val == self.default else ''
            parts.append(f'  <option value="{val}"{selected}>{txt}</option>')
        parts.append('</select></div></div></form>')
        return '\n'.join(parts)


class Input(FormPart):
    def __init__(self, text, name, placeholder='...'):
        super().__init__(text, name)
        self.placeholder = placeholder

    def render(self) -> str:
        """Render input block.

        >>> s = Input('Bla-bla', 'xyz')
        >>> html = s.render()
        """
        parts = [
            '<form>',
            '<div class="form-group row pb-1">',
            '<label class="col-4 col-form-label-sm">',
            f'  {self.text}',
            '</label>',
            '<div class="col d-flex align-items-center">',
            '<input',
            '  class="form-control"',
            '  type="text"',
            f'  name="{self.name}"',
            '  value=""',
            f'  placeholder="{self.placeholder}" /></div></div></form>']
        return '\n'.join(parts)


class StoichiometryInput(Input):
    def __init__(self):
        super().__init__('Stoichiometry:', 'stoichiometry', 'A, AB2, ABC, ...')

    def get_filter_strings(self, query: dict) -> list[str]:
        """Make sure A2B and AB2 both work.

        >>> s = StoichiometryInput()
        >>> s.get_filter_strings({'stoichiometry': 'A2B'})
        ['stoichiometry=AB2']
        >>> s.get_filter_strings({})
        []
        >>> s.get_filter_strings({'stoichiometry': 'garbage'})
        ['stoichiometry=garbage']
        """
        val = query.get(self.name, '')
        if not val:
            return []
        # Reduce A2B2 to AB and so on.
        try:
            f = Formula(val)
        except ValueError:
            pass
        else:
            val = f.reduce()[0].stoichiometry()[0].format('ab2')
        return [f'{self.name}={val}']


class Range(FormPart):
    def __init__(self,
                 text: str,
                 name: str,
                 nonnegative: bool = False,
                 default: tuple[str, str] = ('', '')):
        super().__init__(text, name)
        self.nonnegative = nonnegative
        self.default = default

    def render(self) -> str:
        """Render range block.

        >>> s = Range('Band gap', 'gap')
        >>> html = s.render()
        """
        v1, v2 = self.default
        parts = [
            '<form>',
            '<div class="form-group row pb-1">',
            f'<label for="form_{self.name}"',
            '  class="col-4 col-form-label-sm">',
            f'  {self.text}',
            '</label>',
            '<div class="col d-flex align-items-center">',
            '<input class="form-control"',
            '  type="text"',
            f'  name="from_{self.name}"',
            f'  id="form_{self.name}"',
            f'  value="{v1}" />\n</div>',
            f'<label for="to_{self.name}" class="col-1 align-items-center ',
            'd-flex justify-content-center">-</label>',
            '<div class="col d-flex align-items-center">',
            '<input class="form-control"',
            '  type="text"',
            f'  id="to_{self.name}" ',
            f'  name="to_{self.name}"',
            f'  value="{v2}" /></div></div></form>']
        return '\n'.join(parts)

    def get_filter_strings(self, query: dict) -> list[str]:
        filters = []
        fro = query.get(f'from_{self.name}', '')
        if fro:
            limit = float(fro)
            if not self.nonnegative or limit > 0.0:
                filters.append(f'{self.name}>={fro}')
        to = query.get(f'to_{self.name}', '')
        if to:
            filters.append(f'{self.name}<={to}')
        return filters


class RangeX(Range):
    def __init__(self, text, name, options, names=None):
        super().__init__(text, name)
        self.options = options
        self.names = names

    def render(self) -> str:
        parts = [
            '<form>',
            '<div class="form-group row pb-1">',
            f'<label for="form_{self.name}"',
            '  class="col-4 col-form-label-sm">',
            f'  {self.text}',
            '</label>',
            '<div class="col d-flex align-items-center">',
            '<input',
            '  class="form-control"',
            '  type="text"',
            f'  name="from_{self.name}"',
            '  value="" />',
            '</div>',
            f'<label for="to_{self.name}" class="col-1 align-items-center ',
            'd-flex justify-content-center">-</label>',
            '<div class="col d-flex align-items-center">',
            '<input',
            '  class="form-control"',
            '  type="text"',
            f'  name="to_{self.name}"',
            '  value="" />',
            '</div>'
            '<div class="col d-flex align-items-center">',
            f'<select name="{self.name}" class="form-select">']
        names = self.names or self.options
        for val, txt in zip(self.options, names):
            selected = ' selected' if val == '' else ''
            parts.append(f'  <option value="{val}"{selected}>{txt}</option>')
        parts.append('</select></div></div></form>')
        return '\n'.join(parts)

    def get_filter_strings(self, query: dict) -> list[str]:
        filters = []
        x = query.get(self.name)
        fro = query.get(f'from_{self.name}', '')
        if fro:
            filters.append(f'{x}>={fro}')
        to = query.get(f'to_{self.name}', '')
        if to:
            filters.append(f'{x}<={to}')
        return filters


class RangeS(Range):
    def __init__(self, text, name, options, names=None):
        super().__init__(text, name)
        self.options = options
        self.names = names

    def render(self) -> str:
        names = self.names or self.options

        parts = [
            '<form>',
            '<div class="form-group row pb-1">',
            f'<select name="from_{self.name}" class="form-select">']
        for val, txt in zip(self.options, names):
            selected = ' selected' if val == '' else ''
            parts.append(f'  <option value="{val}"{selected}>{txt}</option>')
        parts.append('</select></div>')

        parts += [
            '<div class="col">',
            f'<select name="to_{self.name}" class="form-select">']
        for val, txt in zip(self.options, names):
            selected = ' selected' if val == '' else ''
            parts.append(f'  <option value="{val}"{selected}>{txt}</option>')
        parts.append('</select></div></div></form>')

        return '\n'.join(parts)

    def get_filter_strings(self, query: dict) -> list[str]:
        filters = []
        fro = query.get(f'from_{self.name}', '')
        if fro:
            filters.append(f'{self.name}>={fro}')
        to = query.get(f'to_{self.name}', '')
        if to:
            filters.append(f'{self.name}<={to}')
        return filters
