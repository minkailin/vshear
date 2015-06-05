#!/bin/bash
latex proceed.tex
dvips -o proceed.ps proceed.dvi
ps2pdf -sPAPERSIZE=a4 proceed.ps
echo 'done'
