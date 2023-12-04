def get_pages(page: int,
              nrows: int,
              limit: int = 25,
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
