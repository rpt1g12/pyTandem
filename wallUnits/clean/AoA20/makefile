# ============================================================================
# Name        : Makefile
# Author      : 
# Version     :
# Copyright   : Your copyright notice
# Description : Makefile for Fortran
# ============================================================================
.PHONY: all

# Change this line if you are using a different Fortran compiler
COMPILER = pdflatex
PDFCOMPILER = pdflatex
NAME = tCf

all: 
	$(COMPILER) -draftmode $(NAME).tex
	$(COMPILER) -draftmode $(NAME).tex >/dev/null
	$(PDFCOMPILER) $(NAME).tex >/dev/null
