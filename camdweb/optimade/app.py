from __future__ import annotations

import json
import os
from datetime import datetime
from functools import lru_cache
from pathlib import Path
from typing import Any, Callable, Iterable

from ase.db import connect
from ase.db.row import AtomsRow
from bottle import Bottle, request, template  # type: ignore
from bottle import response as bottle_response  # type: ignore
from optimade.server.exceptions import BadRequest, NotImplementedResponse  # type: ignore # noqa

from oase import version
from oase.asedata import extract_data_from_db_row, row2dict
from oase.indexdb import IndexDB
from oase.parser import create_parse_function
from oase.tests import make_index_db


def get_config() -> dict:
    if 'OASE_CONFIG_FILE' in os.environ:
        cfg = json.loads(Path(os.environ['OASE_CONFIG_FILE']).read_text())
    else:
        cfg = {'url': 'http://localhost:8080/v1',
               'provider': {
                   'name': 'OASE',
                   'description': 'OPTIMADE for ASE',
                   'prefix': 'cmr',
                   'homepage': 'https://gitlab.com/jensj/ase-optimade'}}
    return cfg


API_VERSION = (1, 1, 0)
API_VERSION_STR = ".".join(str(d) for d in API_VERSION)


INFO = {
    'id': '/',
    'type': 'info',
    'attributes': {
        'api_version': API_VERSION_STR,
        'available_api_versions': [
            {'version': API_VERSION_STR},
        ],
        'formats': ['json'],
        'available_endpoints': ['info',
                                'links',
                                'structures'],
        'entry_types_by_format': {'json': ['structures']},
        'is_index': False}}

META = {
    'api_version': API_VERSION_STR,
    'more_data_available': False,
    'implementation': {
        'name': 'oase',
        'version': version,
        'source_url': 'https://gitlab.com/jensj/ase-optimade',
        'maintainer': {'email': 'jjmo@dtu.dk'}}}

PROPERTIES = [
    ('id', 'string'),
    ('type', 'string'),
    ('elements', 'list'),
    ('nelements', 'integer'),
    ('chemical_formula_reduced', 'string'),
    ('chemical_formula_hill', 'string'),
    ('chemical_formula_anonymous', 'string'),
    ('dimension_types', 'list'),
    ('nperiodic_dimensions', 'integer'),
    ('lattice_vectors,Å', 'list'),
    ('cartesian_site_positions,Å', 'list'),
    ('nsites', 'integer'),
    ('species', 'list'),
    ('species_at_sites', 'list'),
    ('structure_features', 'list')]


INDEX = """
<html>
  <head>
    <title>OASE</title>
  </head>
  <body>
    <form>
      <input
       type="text"
       name="query"
       value="{{ query }}"
       placeholder="OPTIMADE expression"
       size="60">
      <button type="submit">Submit</button>
    </form>
    <br/>
    % for page, text in pages:
    %     if page == -1:
    ...
    %     else:
    <a href="?query={{ query }}&page={{ page }}">{{ text }}</a>
    %     end
    % end
    <br/>
    <table>
      <tr>
      % for column in header:
        <th>{{ column }}</th>
      % end
      </tr>
      % for row in rows:
      <tr>
        % for column in row:
        <td>{{ column }}</td>
        % end
      </tr>
      % end
    </table>
  </body>
</html>
"""


def add_meta(method: Callable) -> Callable:
    """Decorator for adding "meta" stuff to JSON response."""
    def new_method(self: App, *args: Any, **kwargs: Any) -> Any:
        dct = method(self, *args, **kwargs)

        data_returned = dct.get("meta", {}).get("data_returned", 0)

        endpoint = request.urlparts.path.split('/')[-1]
        data_available = getattr(
            self, "meta", {}).get("entries_available", {}).get(endpoint, 1)

        dct['meta'] = META | {
            'time_stamp': f'{datetime.now().isoformat()}',
            'provider': self.cfg['provider'],
            'query': {
                'representation': "?".join((request.urlparts.path,
                                            request.urlparts.query or ""))},
            'data_available': data_available,
            'data_returned': data_returned}

        return dct

    return new_method


def add_error_handling(method: Callable) -> Callable:
    """Decorator for adding error handling to the API,
    based on the optimade errors.

    """
    def new_method(self: App, *args: Any, **kwargs: Any) -> Any:
        try:
            return method(self, *args, **kwargs)
        except BadRequest as exc:
            bottle_response.status = 400
            return {'errors': [{'status': '400',
                                'title': 'Bad Request',
                                'detail': str(exc)}]}

        except (NotImplementedError, NotImplementedResponse) as exc:
            bottle_response.status = 501
            return {'errors': [{'status': '501',
                                'title': 'Not Implemented',
                                'detail': str(exc)}]}
        except Exception:
            bottle_response.status = 500
            # omit any unhandled error details from the response
            return {'errors': [{'status': '500',
                                'title': 'Internal Server Error',
                                'detail': ""}]}

    return new_method


def make_test_app(dbfile: str = 'test.db') -> Bottle:
    """Create simple test app."""
    indexdb, asedb = make_index_db(Path(dbfile))
    return make_app(indexdb, asedb.get)


def make_app_from_file(path: Path) -> Bottle:
    """Create app from ase.db database file."""
    db = connect(path)
    indexdb = IndexDB()
    indexdb.insert([extract_data_from_db_row(row)
                    for row in db.select()])
    return make_app(indexdb, db.get)


def get_pages(page: int,
              nrows: int,
              limit: int,
              extra: int = 2) -> list[tuple[int, str]]:
    npages = nrows // limit + 1
    pages = set(range(extra))
    pages.update(range(page - extra, page + extra + 1))
    pages.update(range(npages - extra, npages))
    buttons = [(max(page - 1, 0), 'previous'),
               (min(page + 1, npages - 1), 'next')]
    prev = -1
    for page in sorted(pages):
        if not 0 <= page < npages:
            continue
        if page - prev > 1:
            buttons.append((-1, '...'))
        r1 = page * limit + 1
        r2 = min((page + 1) * limit, nrows)
        buttons.append((page, f'{r1}-{r2}'))
        prev = page
    return buttons


class App:
    def __init__(self,
                 indexdb: IndexDB,
                 get_row: Callable[[int], AtomsRow]):
        self.indexdb = indexdb
        self.get_row = get_row
        self.parse = create_parse_function()
        self._query = None
        self.cfg = get_config()
        self.get_row_ids = lru_cache(maxsize=64)(self._get_row_ids)

    @property
    def query(self) -> dict[str, str]:
        if self._query is not None:
            return self._query
        return request.query

    def index(self) -> str:
        query = self.query.get('query', '')
        page = int(self.query.get('page', '0'))
        limit = min(int(self.query.get('limit', '10')), 100)
        rowids = self.get_row_ids(query)
        pages = get_pages(page, len(rowids), limit=limit)
        rowids = rowids[page * limit:(page + 1) * limit]
        rows = []
        for id in rowids:
            row = self.get_row(id)
            rows.append([row.id, row.formula])
        return template(INDEX,
                        query=query,
                        rows=rows,
                        pages=pages,
                        header=['id', 'formula'])

    def versions(self) -> bottle_response:
        bottle_response.content_type = (
            "text/csv; header=present; charset=utf-8")
        return 'version\n1\n'

    @add_meta
    def info(self) -> dict:
        dct = INFO.copy()
        aav = dct['attributes']['available_api_versions']  # type: ignore
        aav[0]['url'] = self.cfg['url']
        return {'data': dct}

    @add_meta
    def links(self) -> dict:
        return {'data': []}

    @add_meta
    def info_structures(self) -> dict:
        properties = {}
        for name, type in PROPERTIES:
            name, _, unit = name.partition(',')
            dct = dict(type=type,
                       sortable=type != 'list',
                       description='...')
            if unit:
                dct['unit'] = unit
            properties[name] = dct,
        return {
            'data': {
                'formats': ['json'],
                'description': '...',
                'properties': properties,
                'output_fields_by_format': {
                    'json': list(properties)}}}

    def _get_row_ids(self, query: str) -> list[int]:
        if query:
            tree = self.parse(query)
            selection = self.indexdb.select(tree)
            rowids = self.indexdb.execute(selection)
        else:
            rowids = list(range(1, self.indexdb.nstructures + 1))
        return rowids

    @add_meta
    @add_error_handling
    def structures(self) -> dict:
        rowids = self.get_row_ids(self.query.get('filter', ''))
        offset = int(self.query.get('page_offset', '0'))
        limit = min(int(self.query.get('page_limit', '20')), 100)
        limited_rowids = rowids[offset:offset + limit]

        response = self.response(limited_rowids)
        response["meta"]["data_returned"] = len(rowids)
        return response

    def structures_id(self, id: str) -> dict:
        dct = self.response([int(id)])
        dct['data'] = dct['data'][0]
        return dct

    @add_meta
    def response(self,
                 rowids: Iterable[int]) -> dict:
        response_fields = self.query.get('response_fields')
        if response_fields:
            fields = set(response_fields.split(','))
        else:
            fields = None

        data = []
        for id in rowids:
            row = self.get_row(id)
            dct = row2dict(row)
            optimade_id = str(dct.pop("id"))
            if fields:
                dct = {key: value
                       for key, value in dct.items()
                       if key in fields}
            data.append({'id': optimade_id,
                         'attributes': dct,
                         'type': 'structures'})

        return {'data': data}


def make_app(indexdb: IndexDB,
             get_row: Callable[[int], AtomsRow]) -> Bottle:
    app = Bottle()
    a = App(indexdb, get_row)
    app.route('/')(a.index)
    app.route('/versions')(a.versions)

    a.meta = {}  # type: ignore[attr-defined]
    a.meta["entries_available"] = {  # type: ignore[attr-defined]
        "structures": indexdb.nstructures, "references": 0, "links": 0
    }

    for prefix in [
        "",
        f"/v{API_VERSION[0]}",
        f"/v{API_VERSION[0]}.{API_VERSION[1]}",
        f"/v{API_VERSION[0]}.{API_VERSION[1]}.{API_VERSION[2]}"
    ]:
        app.route(f'{prefix}/info')(a.info)
        app.route(f'{prefix}/links')(a.links)
        app.route(f'{prefix}/info/structures')(a.info_structures)
        app.route(f'{prefix}/structures')(a.structures)
        app.route(f'{prefix}/structures/<id>')(a.structures_id)

    return app
