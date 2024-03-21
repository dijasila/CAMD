from camdweb.session import Sessions, Session


def test_sessions():
    s = Sessions(['a', 'b'], max_sessions=5)
    print(s)
    for _ in range(20):
        s.get(-1)
    assert len(s.sessions) <= 5


def test_session():
    s = Session(1, ['a', 'b'])
    s.update({'toggle': 'a'})
    assert s.columns == ['b']
    s.update({'toggle': 'a'})
    assert s.columns == ['b', 'a']
    s.update({'sort': 'a'})
    s.update({'sort': 'a'})
    assert s.direction == -1
