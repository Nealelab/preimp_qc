import abc


class Element:
    @abc.abstractmethod
    def to_markdown(self):
        pass


class Text(Element):
    def __init__(self, text):
        self.text = text

    def to_markdown(self):
        return self.text


class Table(Element):
    def __init__(self):
        self.header_row = None
        self.rows = None

    def add_header_row(self, row):
        if self.rows:
            assert len(row) == len(self.rows[0])
        self.header_row = row

    def add_row(self, row):
        if self.header_row:
            assert len(row) == len(self.header_row)
        self.rows.append(row)

    def __aenter__(self):
        return self

    def __aexit__(self, exc_type, exc_val, exc_tb):
        pass

    def to_markdown(self):
        md = []

        def row_md(row):
            return f'| {" | ".join(row)} |'

        if self.header_row:
            md.append(row_md(self.header_row))
            divider = " | ".join(['---' for item in self.header_row])
            md.append(f'| {divider} |')

        for row in self.rows:
            md.append(row_md(row))

        return '\n'.join(md)


class OrderedList(Element):
    def __init__(self):
        self.items = []

    def add_item(self, item):
        self.items.append(item)

    def __aenter__(self):
        return self

    def __aexit__(self, exc_type, exc_val, exc_tb):
        pass

    def to_markdown(self):
        md = []
        for i, item in enumerate(self.items):
            md.append(f'{i}. {item}')
        return '\n'.join(md)


class UnorderedList(Element):
    def __init__(self):
        self.items = []

    def add_item(self, item):
        self.items.append(item)

    def __aenter__(self):
        return self

    def __aexit__(self, exc_type, exc_val, exc_tb):
        pass

    def to_markdown(self):
        md = []
        for item in self.items:
            md.append(f'- {item}')
        return '\n'.join(md)


class Code(Element):
    def __init__(self, language=None, code=None):
        self.language = language if language else ''
        self.code = code if code else ''

    def to_markdown(self):
        return f'''
```{self.language}
{self.code}
```
'''


class Image(Element):
    def __init__(self, path, alt_text=None):
        self.path = path
        self.alt_text = alt_text

    def to_markdown(self):
        alt_text = self.alt_text if self.alt_text else ''
        return f'![{alt_text}]({self.path})'


class HLine(Element):
    def to_markdown(self):
        return '---'


class Header(Element):
    def __init__(self, level, header):
        assert 1 <= level <= 6
        self.level = level
        self.header = header

    def to_markdown(self):
        return f'{"#" * self.level} {self.header}'


class Quote(Element):
    def __init__(self):
        self.quotes = []

    def add_quote(self, quote):
        self.quotes.append(quote)

    def __aenter__(self):
        return self

    def __aexit__(self, exc_type, exc_val, exc_tb):
        pass

    def to_markdown(self):
        return '\n'.join(f'> {quote}' for quote in self.quotes)


class Inline(Element):
    def __init__(self):
        self.elements = []

    def add_item(self, element):
        self.elements.append(element)

    def __aenter__(self):
        return self

    def __aexit__(self, exc_type, exc_val, exc_tb):
        pass

    def to_markdown(self):
        md = [element.to_markdown() for element in self.elements]
        return ' '.join(md)


class Section:
    def __init__(self, level):
        assert 1 <= level <= 6, f'Too many nested sections. Max of 6.'
        self.header = None
        self.level = level
        self.contents = []

    def add_header(self, header):
        assert self.header is None
        self.header = header
        h = Header(self.level, header)
        self.contents.append(h)

    def add_table(self) -> Table:
        t = Table()
        self.contents.append(t)
        return t

    def add_ordered_list(self) -> OrderedList:
        ordered_list = OrderedList()
        self.contents.append(ordered_list)
        return ordered_list

    def add_unordered_list(self) -> UnorderedList:
        unordered_list = UnorderedList()
        self.contents.append(unordered_list)
        return unordered_list

    def add_text(self, text):
        t = Text(text)
        self.contents.append(t)

    def add_hline(self):
        hline = HLine()
        self.contents.append(hline)

    def add_inline(self) -> Inline:
        inline = Inline()
        self.contents.append(inline)
        return inline

    def add_image(self, path, alt_text=None):
        image = Image(path, alt_text=alt_text)
        self.contents.append(image)

    def add_quote(self) -> Quote:
        quote = Quote()
        self.contents.append(quote)
        return quote

    def add_section(self):
        s = Section(self.level + 1)
        self.contents.append(s)
        return s

    def to_markdown(self):
        md = []
        for content in self.contents:
            md.append(content.to_markdown())
        return '\n'.join(md)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        pass


class Report:
    def __init__(self, dest):
        self.dest = dest
        self.sections = []

    def add_section(self):
        s = Section(1)
        self.sections.append(s)
        return s

    def to_markdown(self):
        md = []
        for section in self.sections:
            md.append(section.to_markdown())
        return '\n\n'.join(md)

    def write(self):
        md = self.to_markdown()
        with open(self.dest, 'w') as f:
            f.write(md)
