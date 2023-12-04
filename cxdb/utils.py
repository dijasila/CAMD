from typing import Iterable


def table(header: list[str] | None, rows: list[Iterable]) -> str:
    """Create HTML table.

    >>> print(table(['A', 'B'], [[1.2, 'hello'], [2.0, 'hi!']]))
    <table>
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
        f'<table>\n{head} <tbody>\n  <tr>\n   ' +
        '\n  </tr>\n  <tr>\n   '.join(
            '\n   '.join(f'<td>{x}</td>' for x in row)
            for row in rows) +
        '\n  </tr>\n </tbody>\n</table>')
