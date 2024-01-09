class Sessions:
    def __init__(self,
                 columns: list[str],
                 max_sessions=200):
        self.columns = columns
        self.max_sessions = max_sessions
        self.sid = 0
        self.sessions: dict[int, Session] = {}

    def get(self, sid):
        if sid not in self.sessions:
            sid = self.sid
            self.sid += 1
            self.sessions[sid] = Session(sid, self.columns)
            if len(self.sessions) > self.max_sessions:
                self.sessions = {
                    sid: session
                    for sid, session
                    in list(self.sessions.items())[self.max_sessions // 2:]}
        return self.sessions[sid]


class Session:
    def __init__(self,
                 sid: int,
                 columns: list[str]):
        self.sid = sid
        self.columns = list(columns)
        self.filter = ''
        self.page = 0
        self.sort = ''
        self.direction = 1
        self.rows_per_page = 25

    def update(self,
               filter: str,
               query: dict) -> None:
        if filter != self.filter:
            self.filter = filter
            self.page = 0
            return

        column = query.get('toggle')
        if column:
            if column in self.columns:
                self.columns.remove(column)
            else:
                self.columns.append(column)

        column = query.get('sort')
        if column:
            if column == self.sort:
                self.direction *= -1
            else:
                self.sort = column
                self.direction = 1

        page = query.get('page')
        if page:
            self.page = int(page)
