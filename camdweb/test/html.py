from html.parser import HTMLParser


class HTMLCheckParser(HTMLParser):
    def __init__(self):
        super().__init__()
        self.tags = []

    def handle_starttag(self, tag, attrs):
        if tag in {'hr', 'input'}:
            return
        self.tags.append(tag)

    def handle_endtag(self, tag):
        assert tag == self.tags.pop()


def check_html(html):
    p = HTMLCheckParser()
    p.feed(html)
    p.close()
    assert not p.tags
