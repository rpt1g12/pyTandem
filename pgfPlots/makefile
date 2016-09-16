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
NAME = pgfPlots

all: 
	$(COMPILER) -draftmode -enable-write18 $(NAME).tex
	$(COMPILER) -draftmode -enable-write18 $(NAME).tex >/dev/null
	$(PDFCOMPILER) -enable-write18 $(NAME).tex >/dev/null
