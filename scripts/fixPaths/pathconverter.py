"""
Path Converter.

pymdownx.pathconverter
An extension for Python Markdown.

An extension to covert tag paths to relative or absolute:

Given an absolute base and a target relative path, this extension searches for file
references that are relative and converts them to a path relative
to the base path.

-or-

Given an absolute base path, this extension searches for file
references that are relative and converts them to absolute paths.

MIT license.

Copyright (c) 2014 - 2017 Isaac Muse <isaacmuse@gmail.com>

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated
documentation files (the "Software"), to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software,
and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions
of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED
TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.
"""
from __future__ import unicode_literals
import os
import re
import sys
import logging
from markdown import Extension
from markdown.postprocessors import Postprocessor
from . import util

RE_TAG_HTML = r'''(?xus)
    (?:
        (?P<comments>(\r?\n?\s*)<!--[\s\S]*?-->(\s*)(?=\r?\n)|<!--[\s\S]*?-->)|
        (?P<open><(?P<tag>(?:%s)))
        (?P<attr>(?:\s+[\w\-:]+(?:\s*=\s*(?:"[^"]*"|'[^']*'))?)*)
        (?P<close>\s*(?:\/?)>)
    )
    '''

RE_TAG_LINK_ATTR = re.compile(
    r'''(?xus)
    (?P<attr>
        (?:
            (?P<name>\s+(?:href|src)\s*=\s*)
            (?P<path>"[^"]*"|'[^']*')
        )
    )
    '''
)

log = logging.getLogger(__name__)


def pprint(*args):
    print(*args)
    sys.stdout.flush()


def repl_absolute(m, key, val):
    """Replace path with absolute path."""

    link = m.group(0)
    try:
        scheme, netloc, path, params, query, fragment, is_url, is_absolute = util.parse_url(m.group('path')[1:-1])
        new_path = m.group('path')[1:-1].replace('../{{ %s }}' % key, val)

        if (not is_absolute and not is_url):
            link = '%s"%s"' % (m.group('name'), new_path)

    except Exception:  # pragma: no cover
        # Parsing crashed and burned; no need to continue.
        pass

    return link


def repl(m, key, val):
    """Replace."""

    if m.group('comments'):
        tag = m.group('comments')
    else:
        tag = m.group('open')
        tag += RE_TAG_LINK_ATTR.sub(lambda m2: repl_absolute(m2, key, val), m.group('attr'))
        tag += m.group('close')
    return tag


class PathConverterPostprocessor(Postprocessor):
    """Post process to find tag links to convert."""

    def run(self, text):
        """Find and convert paths."""

        variables = self.config['variables']
        # relativepath = self.config['relative_path']
        # absolute = bool(self.config['absolute'])
        tags = re.compile(RE_TAG_HTML % '|'.join(self.config['tags'].split()))
        # pprint(absolute, basepath, relativepath)
        # if not absolute and basepath and relativepath:
        #     text = tags.sub(lambda m: repl(m, basepath, relativepath), text)
        for key, val in variables.items():
            text = tags.sub(lambda m, k=key, v=val: repl(m, k, v), text)
        return text


class PathConverterExtension(Extension):
    """PathConverter extension."""

    def __init__(self, *args, **kwargs):
        """Initialize."""

        self.config = {
            'variables': [{}, "Dict of variables to replace"],
            'tags': ["a link", "tags to convert src and/or href in - Default: 'img scripts a link'"]
        }

        super(PathConverterExtension, self).__init__(*args, **kwargs)

    def extendMarkdown(self, md, md_globals):
        """Add post processor to Markdown instance."""

        rel_path = PathConverterPostprocessor(md)
        rel_path.config = self.getConfigs()
        md.postprocessors.add("path-converter", rel_path, "_end")
        md.registerExtension(self)


def makeExtension(*args, **kwargs):
    """Return extension."""

    return PathConverterExtension(*args, **kwargs)
