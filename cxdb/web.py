from __future__ import annotations

import pandas as pd
import json
import plotly
import plotly.express as px
from bottle import run, route, request, template, Bottle

INDEX = """
<html>
  <head>
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
        n = int(request.query.data)
        print(request.query.what)
        return plot_atoms(min(n, 100))


def plot_atoms(n=5):
    df = pd.DataFrame(
        {'Fruit': ['Apples', 'Oranges', 'Bananas', 'Apples', 'Oranges',
                   'Bananas'],
         'Amount': [4, 1, 2, 2, 4, n],
         'City': ['SF', 'SF', 'SF', 'Montreal', 'Montreal', 'Montreal']})
    fig = px.bar(df, x='Fruit', y='Amount', color='City',
                 barmode='group')
    return json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)


C2DB(42).app.run(host='localhost', port=8080, debug=True)
