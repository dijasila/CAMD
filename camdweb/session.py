"""Session handling.

Each client gets a unique session id (1, 2, ...) and a Session object.
This object keeps track of:

* name of columns to show
* filter string
* page number
* name of column used for sorting
* sort direction
* number of rows per page
"""
from __future__ import annotations


class Sessions:
    def __init__(self,
                 columns: list[str],
                 max_sessions: int = 200):
        """Session factory.

        Parameters
        ==========
        columns:
            Initial list of column names for new clients.
        max_sessions:
            Maximum number of session to keep in memory.
        """
        self.columns = columns
        self.max_sessions = max_sessions
        self.sid = 0
        self.sessions: dict[int, Session] = {}

    def get(self, sid: int) -> Session:
        """Return existing session object or create a new one."""
        if sid not in self.sessions:
            sid = self.sid
            self.sid += 1
            self.sessions[sid] = Session(sid, self.columns)
            # Remove old session objects:
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
        """Update session object.

        toggle:
            add/remove a column

        sort:
            change direction of sorting

        page:
            go to another page
        """
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
