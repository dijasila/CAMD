from typing import Iterable


def table(header: list[str] | None, rows: list[Iterable]) -> str:
    if header is None:
        head = ''
    else:
        head = (' <thead>\n  <tr>\n   <th>' +
                '</th>\n   '.join(header) +
                '</th>\n  </tr> </thead>\n')
    return (
        f'<table>\n{head} <tbody>\n  <tr>\n' +
        '</tr>\n  <tr>'.join(
            '\n   '.join(f'<td>{x}</td>' for x in row)
            for row in rows) +
        '\n  </tr>\n </tbody>\n</table>')
