from html.parser import HTMLParser


class HTMLCheckParser(HTMLParser):
    def __init__(self):
        super().__init__()
        self.tags = []

    def handle_startendtag(self, tag, attrs):
        pass

    def handle_starttag(self, tag, attrs):
        if tag in {'hr', 'input', 'img'}:
            # We allow skipping </img>
            return
        self.tags.append(tag)

    def handle_endtag(self, tag):
        starttag = self.tags.pop()
        assert tag == starttag, f'<{starttag}> ... </{tag}>'


def check_html(html):
    p = HTMLCheckParser()
    p.feed(html)
    p.close()
    assert not p.tags
