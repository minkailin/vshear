#!/bin/bash
latex poster.tex
dvips -o poster.ps poster.dvi
###ps2pdf -sPAPERSIZE=a0 poster.ps
ps2pdf poster.ps
echo done
