from cxdb.session import Sessions


def test_session():
    s = Sessions(['a', 'b'], 5)
    for _ in range(20):
        s.get(-1)
    assert len(s.sessions) <= 5
