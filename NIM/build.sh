#!/bin/sh -f
cd NIM_TOF12
rm nim_tof12.dvi nim_tof12.dvi.ps  nim_tof12.dvi.pdf nim_tof12.spl nim_tof12.log nim_tof12.aux nim_tof12.pdf

#read -p "Do you wish to spell check? (y/n) " -n 1
#echo
#if [[ $REPLY =~ ^[Yy]$ ]]; then
#    aspell -t -c nick/nick_pmt_test.tex;  
#fi


pdflatex nim_tof12
bibtex nim_tof12
pdflatex nim_tof12
pdflatex nim_tof12
rm ralf/*.aux
rm evan/*.aux
rm gleb/*.aux
rm ye/*.aux
rm arjun/*.aux
rm gary/*.aux
rm nick/*.aux
rm *.aux *.log
#dvips nim_tof12.dvi
#ps2pdf nim_tof12.dvi.ps
#/Applications/Adobe\ Reader.app/Contents/MacOS/AdobeReader nim_tof12.pdf &
#acroread nim_tof12.dvi.pdf &
#acroread nim_tof12.pdf &
