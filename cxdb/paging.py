def get_pages(page: int,
              nrows: int,
              limit: int = 25,
              pad: int = 2) -> list[tuple[int, str]]:
    """Paging button helper.

    >>> for p in get_pages(5, 100, 9, 1):
    ...     print(p)
    (0, 'previous')
    (2, 'next')
    (0, '1-9')
    (1, '10-18')
    (2, '19-27')
    (3, '28-36')
    (-1, '...')
    (10, '91-99')
    (11, '100-100')
    """
    npages = nrows // limit + 1
    pages = set(range(pad))
    pages.update(range(page - pad, page + pad + 1))
    pages.update(range(npages - pad, npages))
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
