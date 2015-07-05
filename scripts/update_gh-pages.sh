#!/bin/bash
# Usage: cd $LIBIGL; scripts/update_gh-pages.sh
set -o xtrace
git pull && git checkout gh-pages && git rebase master && git pull

HEADER="title: libigl
author: Alec Jacobson and Daniele Panozzo and others
css: tutorial/style.css
html header:   <script type='text/javascript' src='http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML'></script>
<link rel='stylesheet' href='http://yandex.st/highlightjs/7.3/styles/default.min.css'>
<script src='http://yandex.st/highlightjs/7.3/highlight.min.js'></script>
<script>hljs.initHighlightingOnLoad();</script>

"

echo "$HEADER" \
  | cat - README.md | multimarkdown -o index.html && \
  git commit -m "update index.html to match README.md" README.md index.html

HEADER="title: libigl
author: Alec Jacobson and Daniele Panozzo and others
css: ../tutorial/style.css
html header:   <script type='text/javascript' src='http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML'></script>
<link rel='stylesheet' href='http://yandex.st/highlightjs/7.3/styles/default.min.css'>
<script src='http://yandex.st/highlightjs/7.3/highlight.min.js'></script>
<script>hljs.initHighlightingOnLoad();</script>

"

echo "$HEADER" \
  | cat - google-soc/google-soc.md | multimarkdown -o google-soc/index.html && \
  git commit -m "update google-soc/index.html to match google-soc/google-soc.md" google-soc/google-soc.md google-soc/index.html 

echo "$HEADER" \
  | cat - optional/README.md | multimarkdown -o optional/index.html && \
  git commit -m "update index.html to match README.md" optional/README.md \
  optional/index.html

git push origin gh-pages && git checkout master && git merge gh-pages && git push
