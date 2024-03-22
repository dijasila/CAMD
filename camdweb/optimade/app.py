from __future__ import annotations

import importlib
from datetime import datetime
from typing import Any, Callable, Iterable

import numpy as np
from ase.formula import Formula
from bottle import request, response

from camdweb.optimade.filter import create_parse_function, select

CFG = {'url': 'http://localhost:8080/v1',
       'provider': {
           'name': 'CAMD',
           'description': 'OPTIMADE for CAMD',
           'prefix': 'cmr',
           'homepage': 'https://gitlab.com/camd/camd-web'}}

API_VERSION = (1, 1, 0)
API_VERSION_STR = ".".join(str(d) for d in API_VERSION)


INFO = {
    'id': '/',
    'type': 'info',
    'attributes': {
        'api_version': API_VERSION_STR,
        'available_api_versions': [
            {'version': API_VERSION_STR,
             'url': 'http://localhost:8080/v1'}],
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
        'version': importlib.metadata.version('camdweb'),
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


def add_meta(method: Callable) -> Callable:
    """Decorator for adding "meta" stuff to JSON response."""
    def new_method(self, *args: Any, **kwargs: Any) -> Any:
        dct = method(self, *args, **kwargs)

        data_returned = dct.get("meta", {}).get("data_returned", 0)

        endpoint = request.urlparts.path.split('/')[-1]
        data_available = getattr(
            self, "meta", {}).get("entries_available", {}).get(endpoint, 1)

        dct['meta'] = META | {
            'time_stamp': f'{datetime.now().isoformat()}',
            'provider': CFG['provider'],
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
    def new_method(self, *args: Any, **kwargs: Any) -> Any:
        try:
            return method(self, *args, **kwargs)
        except Exception:
            response.status = 500
            # omit any unhandled error details from the response
            return {'errors': [{'status': '500',
                                'title': 'Internal Server Error',
                                'detail': ""}]}

    return new_method


class Optimade:
    def __init__(self, materials):
        self.materials = materials
        self.parse = create_parse_function()
        self.meta = {
            "entries_available": {
                "structures": len(materials),
                "references": 0,
                "links": 0}}

    def versions(self):
        response.content_type = (
            "text/csv; header=present; charset=utf-8")
        return 'version\n1\n'

    @add_meta
    def info(self) -> dict:
        dct = INFO.copy()
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

    @add_meta
    @add_error_handling
    def structures(self) -> dict:
        rowids = select(self.parse(request.query.get('filter', '')),
                        self.materials.index)
        offset = int(request.query.get('page_offset', '0'))
        limit = min(int(request.query.get('page_limit', '20')), 100)
        limited_rowids = sorted(rowids)[offset:offset + limit]

        response = self.response(limited_rowids)
        response["meta"]["data_returned"] = len(limited_rowids)
        return response

    def structures_id(self, id: str) -> dict:
        dct = self.response([int(id)])
        dct['data'] = dct['data'][0]
        return dct

    @add_meta
    def response(self,
                 rowids: Iterable[int]) -> dict:
        response_fields = request.query.get('response_fields')
        if response_fields:
            fields = set(response_fields.split(','))
        else:
            fields = None

        data = []
        for id in rowids:
            uid = self.materials.i2uid[id]
            material = self.materials[uid]
            dct = material2dict(material)
            if fields:
                dct = {key: value
                       for key, value in dct.items()
                       if key in fields}
            data.append({'id': int(id),
                         'attributes': dct,
                         'type': 'structures'})

        return {'data': data}


def add_optimade(webapp) -> None:
    om = Optimade(webapp.materials)
    webapp.optimade = om
    app = webapp.app
    app.route('/optimade/versions')(om.versions)

    for prefix in ["",
                   f"/v{API_VERSION[0]}",
                   f"/v{API_VERSION[0]}.{API_VERSION[1]}",
                   f"/v{API_VERSION[0]}.{API_VERSION[1]}.{API_VERSION[2]}"]:
        prefix = f'/optimade{prefix}'
        app.route(f'{prefix}/info')(om.info)
        app.route(f'{prefix}/links')(om.links)
        app.route(f'{prefix}/info/structures')(om.info_structures)
        app.route(f'{prefix}/structures')(om.structures)
        app.route(f'{prefix}/structures/<id>')(om.structures_id)

    return app


def get_optimade_things(formula: Formula, pbc: np.ndarray) -> dict:
    """Collect some OPTIMADE stuff."""
    _, reduced, num = formula.stoichiometry()
    count = reduced.count()

    # Alphapetically sorted:
    reduced = Formula.from_dict({symbol: count[symbol]
                                 for symbol in sorted(count)})

    # Elements with highest proportion should appear first:
    c = ord('A')
    dct = {}
    for n in sorted(count.values(), reverse=True):
        dct[chr(c)] = n
        c += 1
    anonymous = Formula.from_dict(dct)

    return {
        'chemical_formula_descriptive': None,
        'chemical_formula_reduced': f'{reduced}',
        'chemical_formula_anonymous': f'{anonymous}',
        'chemical_formula_hill': f'{formula:hill}',
        'nsites': num * sum(count.values()),
        'nelements': len(count),
        'nperiodic_dimensions': int(sum(pbc))}


def material2dict(material) -> dict[str, Any]:
    """Create OPTIMADE dictionary."""
    count = material.count
    formula = Formula.from_dict(count)
    last_modified = '%Y-%m-%dT%H:%M:%SZ'  # ???
    atoms = material.atoms
    dct = {
        'cartesian_site_positions': atoms.positions.tolist(),
        'species_at_sites': atoms.get_chemical_symbols(),
        'species': [{'name': symbol,
                     'chemical_symbols': [symbol],
                     'concentration': [1.0]}
                    for symbol in count],
        'lattice_vectors': atoms.cell.tolist(),
        'dimension_types': atoms.pbc.astype(int).tolist(),
        'last_modified': last_modified,
        'elements': list(count),
        'elements_ratios': [n / len(atoms) for n in count.values()],
        'structure_features': []}
    dct |= get_optimade_things(formula, atoms.pbc)
    return dct
