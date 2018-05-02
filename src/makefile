# Makefile for the SEING package for molecular fingerprinting
#
# 
# Author: Mardochee Reveil
# Date created: 9/29/17
#

CC = g++
DEBUG = -g
CFLAGS = -Wall -c $(DEBUG)
LFLAGS = -Wall $(DEBUG)


SRC = atom.cpp atomicsystem.cpp neighborlist.cpp fingerprintgenerator.cpp \
      genericlocalcalculator.cpp gaussiancalculator.cpp inputs.cpp \
      utilities.cpp bispectrumcalculator.cpp zernikecalculator.cpp \
      main.cpp periodictable.cpp agnicalculator.cpp

SEING_EXE = ../bin/seing

all : seing

seing: $(SRC)
	$(CC) -std=c++11 $(SRC) -o $(SEING_EXE) -fopenmp

clean :
	\rm *.o *~ p1
