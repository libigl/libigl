#!/bin/bash
# Usage: cd $LIBIGL; scripts/update_gh-pages.sh
set -o xtrace

HEADER="title: libigl
author: Alec Jacobson and Daniele Panozzo and others
css: tutorial/style.css
html header:   <script type='text/javascript' src='http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML'></script>
<link rel='stylesheet' href='http://yandex.st/highlightjs/7.3/styles/default.min.css'>
<script src='http://yandex.st/highlightjs/7.3/highlight.min.js'></script>
<script>hljs.initHighlightingOnLoad();</script>

"

echo "$HEADER" \
  | cat - README.md | multimarkdown -o index.html

echo "$HEADER" \
  | cat - style-guidelines.md | multimarkdown -o style-guidelines.html 

HEADER="title: libigl
author: Alec Jacobson and Daniele Panozzo and others
css: ../tutorial/style.css
html header:   <script type='text/javascript' src='http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML'></script>
<link rel='stylesheet' href='http://yandex.st/highlightjs/7.3/styles/default.min.css'>
<script src='http://yandex.st/highlightjs/7.3/highlight.min.js'></script>
<script>hljs.initHighlightingOnLoad();</script>

"

echo "$HEADER" \
  | cat - optional/README.md | multimarkdown -o optional/index.html

multimarkdown tutorial/tutorial.md -o tutorial/tutorial.html
