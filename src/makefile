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

INCLUDES = -I /fs/home/mr937/Boost/boost_1_65_1

SRC = atom.cpp atomicsystem.cpp neighborlist.cpp fingerprintgenerator.cpp \
      genericlocalcalculator.cpp gaussiancalculator.cpp inputs.cpp \
      utilities.cpp bispectrumcalculator.cpp zernikecalculator.cpp \
      diamondcalculator.cpp main.cpp periodictable.cpp agnicalculator.cpp

SEING_EXE = ../bin/seing

all : seing

seing: $(SRC)
	$(CC) $(SRC) -o $(SEING_EXE) -fopenmp

clean :
	\rm *.o *~ p1
