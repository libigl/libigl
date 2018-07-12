"""
General utilities.

MIT license.

Copyright (c) 2017 Isaac Muse <isaacmuse@gmail.com>
"""
from __future__ import unicode_literals
import sys
import copy
import re

PY3 = sys.version_info >= (3, 0)
PY34 = sys.version_info >= (3, 4)

if PY3:
    uchr = chr  # noqa
    from urllib.request import pathname2url, url2pathname  # noqa
    from urllib.parse import urlparse, urlunparse, quote  # noqa
    from html.parser import HTMLParser  # noqa
    if PY34:
        import html  # noqa
        html_unescape = html.unescape  # noqa
    else:  # pragma: no cover
        html_unescape = HTMLParser().unescape  # noqa
else:
    uchr = unichr  # noqa
    from urllib import pathname2url, url2pathname, quote  # noqa
    from urlparse import urlparse, urlunparse  # noqa
    from HTMLParser import HTMLParser  # noqa
    html_unescape = HTMLParser().unescape  # noqa

RE_WIN_DRIVE_LETTER = re.compile(r"^[A-Za-z]$")
RE_WIN_DRIVE_PATH = re.compile(r"^[A-Za-z]:(?:\\.*)?$")
RE_URL = re.compile('(http|ftp)s?|data|mailto|tel|news')
IS_NARROW = sys.maxunicode == 0xFFFF

if IS_NARROW:
    def get_code_points(s):
        """Get the Unicode code points."""

        pt = []

        def is_full_point(p, point):
            """
            Check if we have a full code point.

            Surrogates are stored in point.
            """
            v = ord(p)
            if 0xD800 <= v <= 0xDBFF:
                del point[:]
                point.append(p)
                return False
            if point and 0xDC00 <= v <= 0xDFFF:
                point.append(p)
                return True
            del point[:]
            return True

        return [(''.join(pt) if pt else c) for c in s if is_full_point(c, pt)]

    def get_ord(c):
        """Get Unicode ord."""

        if len(c) == 2:
            high, low = [ord(p) for p in c]
            ordinal = (high - 0xD800) * 0x400 + low - 0xDC00 + 0x10000
        else:
            ordinal = ord(c)

        return ordinal

    def get_char(value):
        """Get the Unicode char."""
        if value > 0xFFFF:
            c = ''.join(
                [
                    uchr(int((value - 0x10000) / (0x400)) + 0xD800),
                    uchr((value - 0x10000) % 0x400 + 0xDC00)
                ]
            )
        else:
            c = uchr(value)
        return c

else:
    def get_code_points(s):
        """Get the Unicode code points."""

        return [c for c in s]

    def get_ord(c):
        """Get Unicode ord."""

        return ord(c)

    def get_char(value):
        """Get the Unicode char."""

        return uchr(value)


def escape_chars(md, echrs):
    """
    Add chars to the escape list.

    Don't just append as it modifies the global list permanently.
    Make a copy and extend **that** copy so that only this Markdown
    instance gets modified.
    """

    escaped = copy.copy(md.ESCAPED_CHARS)
    for ec in echrs:
        if ec not in escaped:
            escaped.append(ec)
    md.ESCAPED_CHARS = escaped


def parse_url(url):
    """
    Parse the URL.

    Try to determine if the following is a file path or
    (as we will call anything else) a URL.

    We return it slightly modified and combine the path parts.

    We also assume if we see something like c:/ it is a Windows path.
    We don't bother checking if this **is** a Windows system, but
    'nix users really shouldn't be creating weird names like c: for their folder.
    """

    is_url = False
    is_absolute = False
    scheme, netloc, path, params, query, fragment = urlparse(html_unescape(url))

    if RE_URL.match(scheme):
        # Clearly a url
        is_url = True
    elif scheme == '' and netloc == '' and path == '':
        # Maybe just a url fragment
        is_url = True
    elif scheme == 'file' and (RE_WIN_DRIVE_PATH.match(netloc)):
        # file://c:/path or file://c:\path
        path = '/' + (netloc + path).replace('\\', '/')
        netloc = ''
        is_absolute = True
    elif scheme == 'file' and netloc.startswith('\\'):
        # file://\c:\path or file://\\path
        path = (netloc + path).replace('\\', '/')
        netloc = ''
        is_absolute = True
    elif scheme == 'file':
        # file:///path
        is_absolute = True
    elif RE_WIN_DRIVE_LETTER.match(scheme):
        # c:/path
        path = '/%s:%s' % (scheme, path.replace('\\', '/'))
        scheme = 'file'
        netloc = ''
        is_absolute = True
    elif scheme == '' and netloc != '' and url.startswith('//'):
        # //file/path
        path = '//' + netloc + path
        scheme = 'file'
        netloc = ''
        is_absolute = True
    elif scheme != '' and netloc != '':
        # A non-filepath or strange url
        is_url = True
    elif path.startswith(('/', '\\')):
        # /root path
        is_absolute = True

    return (scheme, netloc, path, params, query, fragment, is_url, is_absolute)


class PymdownxDeprecationWarning(UserWarning):  # pragma: no cover
    """Deprecation warning for Pymdownx that is not hidden."""
