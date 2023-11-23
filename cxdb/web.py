from __future__ import annotations

import pandas as pd
import json
import plotly
import plotly.express as px
from bottle import run, route, request, template, Bottle

INDEX = """
<html>
  <head>
  <link rel="stylesheet" href="https://cdn.simplecss.org/simple.min.css">
  <title>C2DB</title>
  </head>
  <body>
    <form>
      <input
       type="text"
       name="query"
       value="{{ query }}"
       placeholder="chemical symbol"
       size="60">
      <button type="submit">Submit</button>
    </form>
    <br/>
    <table>
      <tr>
      % for column in header:
        <th>{{ column }}</th>
      % end
      </tr>
      % for id, row in rows:
      <tr>
        % for column in row:
        <td>{{ column }}</td>
        % end
        <td><a href=/material/{{ id }}>{{ id }}</a></td>
      </tr>
      % end
    </table>
  </body>
</html>
"""

MATERIAL = """
<html>
<head>
<script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
<script src="https://ajax.googleapis.com/ajax/libs/jquery/3.5.1/jquery.min.js">
</script>
<script>
  function cb(selection, what) {
      $.getJSON(
          {
              url: "/callback",
              data: { 'data': selection, 'what': what },
              success: function (result) {
                  Plotly.newPlot(what, result, {staticPlot: true});;
                  }
              });
      }
</script>
<title>C2DB-{{ name }}</title>
</head>
<body>
<a href="/">Home</a>
% for title, html in sections:
<h1>{{ title }}</h1>
{{ !html }}
% end
</html>
"""


class C2DB:
    def __init__(self, materials):
        self.materials = materials
        self.app = Bottle()
        self.app.route('/')(self.index)
        self.app.route('/material/<id:int>')(self.material)
        self.app.route('/callback')(self.callback)

    def index(self) -> str:
        query = request.query.get('query', '')
        return template(INDEX,
                        query=query,
                        rows=[(12345, ['a', 'b']),
                              (42, ['H2O', '234'])],
                        header=['xxx', 'yyy'])

    def material(self, id: int) -> str:
        return template(MATERIAL,
                        name=f'I{id}',
                        atoms_json=plot_atoms())

    def callback(self):
        request.query.data
        name = request.query.name
        id = request.query.id
        return self.callbacks[name](id, request)




C2DB(42).app.run(host='localhost', port=8080, debug=True)
