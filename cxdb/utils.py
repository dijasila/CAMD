import abc
from __future__ import annotations
from typing import Iterable, Sequence


def table(header: list[str] | None, rows: Sequence[Iterable]) -> str:
    """Create HTML table.

    Example:

    === =====
    A   B
    === =====
    1.2 hello
    2.0 hi!
    === =====

    >>> print(table(['A', 'B'], [[1.2, 'hello'], [2.0, 'hi!']]))
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
    """
    if header is None:
        head = ''
    else:
        head = (' <thead>\n  <tr>\n   <th>' +
                '</th>\n   <th>'.join(header) +
                '</th>\n  </tr>\n </thead>\n')
    return (
        f'<table class="table table-striped">\n{head} <tbody>\n  <tr>\n   ' +
        '\n  </tr>\n  <tr>\n   '.join(
            '\n   '.join(f'<td>{x}</td>' for x in row)
            for row in rows) +
        '\n  </tr>\n </tbody>\n</table>')


class FormPart(abc.ABC):
    def __init__(self, text: str, name: str):
        self.text = text
        self.name = name

    @abc.abstractmethod
    def render(self, query: dict) -> str:
        raise NotImplementedError

    def get_filter_strings(self, query: dict) -> list[str]:
        val = query.get(self.name, '')
        if val:
            return [f'{self.name}={val}']
        return []


class Select(FormPart):
    def __init__(self, text, name, options):
        super().__init__(text, name)
        self.options = options

    def render(self, query: dict) -> str:
        """Render select block.

        >>> s = Select('Bla-bla', 'xyz', ['A', 'B', 'C'])
        >>> print(s.render({'xyz': 'B'}))
        <label class="form-label">Bla-bla</label>
        <select name="xyz" class="form-select">
          <option value="A">A</option>
          <option value="B" selected>B</option>
          <option value="C">C</option>
        </select>
        """
        selection = query.get(self.name)
        parts = [f'<label class="form-label">{self.text}</label>\n'
                 f'<select name="{self.name}" class="form-select">']
        for val in self.options:
            selected = ' selected' if selection == val else ''
            parts.append(f'  <option value="{val}"{selected}>{val}</option>')
        parts.append('</select>')
        return '\n'.join(parts)


class Input(FormPart):
    def __init__(self, text, name, placeholder='...'):
        super().__init__(text, name)
        self.placeholder = placeholder

    def render(self, query: dict) -> str:
        """Render input block.

        >>> s = Input('Bla-bla', 'xyz')
        >>> print(s.render({'xyz': 'abc'}))
        <label class="form-label">Bla-bla</label>
        <input
          class="form-control"
          type="text"
          name="xyz"
          value="abc"
          placeholder="..." />
        """
        value = query.get(self.name, '')
        parts = [
            f'<label class="form-label">{self.text}</label>',
            '<input',
            '  class="form-control"',
            '  type="text"',
            f'  name="{self.name}"',
            f'  value="{value}"',
            f'  placeholder="{self.placeholder}" />']
        return '\n'.join(parts)
