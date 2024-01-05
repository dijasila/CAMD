from cxdb.session import Sessions, Session


def test_sessions():
    s = Sessions(['a', 'b'], 5)
    for _ in range(20):
        s.get(-1)
    assert len(s.sessions) <= 5


def test_session():
    s = Session(1, ['a', 'b'])
    s.update({'togle': 'a'})
    assert s.columns == ['a', 'b']
    s.update({'togle': 'a'})
    assert s.columns == ['a', 'b']
    s.update({'sort': 'a'})
    s.update({'sort': 'a'})
    assert s.direction == -1
