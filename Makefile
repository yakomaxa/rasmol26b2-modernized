# Makefile
# RasMol2 Molecular Graphics
# Roger Sayle, March 1995
# Version 2.6
#
# Modernised build: targets a current gcc/clang + SDL2 instead of
# 1990s Unix compilers and raw X11/Xlib. See sdlwin.c for the display
# backend that replaced x11win.c.

CC = gcc-13
CFLAGS = -g -O2 -finline-functions $(shell sdl2-config --cflags)

# For Debugging use LFLAGS =
LFLAGS = -s

LIBS = -lm $(shell sdl2-config --libs)


SRCS = rasmol.c molecule.c infile.c transfor.c command.c abstree.c \
       render.c repres.c sdlwin.c pixutils.c outfile.c script.c

OBJS = rasmol.o molecule.o infile.o transfor.o command.o abstree.o \
       render.o repres.o sdlwin.o pixutils.o outfile.o script.o


rasmol:		$(OBJS)
		$(CC) -o rasmol $(LFLAGS) $(OBJS) $(LIBS)
		chmod 755 rasmol

rasmol.o:	rasmol.c rasmol.h molecule.h transfor.h command.h \
		abstree.h render.h graphics.h pixutils.h outfile.h
		$(CC) -c $(CFLAGS) rasmol.c

molecule.o:	molecule.c molecule.h rasmol.h command.h abstree.h \
		transfor.h render.h
		$(CC) -c $(CFLAGS) molecule.c

infile.o:	infile.c infile.h
		$(CC) -c $(CFLAGS) infile.c

transfor.o:	transfor.c transfor.h rasmol.h molecule.h command.h \
		abstree.h render.h graphics.h
		$(CC) -c $(CFLAGS) transfor.c

command.o:	command.c command.h rasmol.h tokens.h abstree.h \
		molecule.h transfor.h render.h graphics.h pixutils.h \
                outfile.h
		$(CC) -c $(CFLAGS) command.c

abstree.o:	abstree.c abstree.h rasmol.h molecule.h
		$(CC) -c $(CFLAGS) abstree.c

render.o:	render.c render.h rasmol.h molecule.h transfor.h \
		command.h abstree.h graphics.h pixutils.h
		$(CC) -c $(CFLAGS) render.c

repres.o:	repres.c repres.h rasmol.h
		$(CC) -c $(CFLAGS) repres.c

sdlwin.o:	sdlwin.c graphics.h rasmol.h command.h
		$(CC) -c $(CFLAGS) sdlwin.c

pixutils.o:	pixutils.c pixutils.h rasmol.h font.h molecule.h \
		transfor.h render.h graphics.h
		$(CC) -c $(CFLAGS) pixutils.c

outfile.o:	outfile.c outfile.h rasmol.h molecule.h command.h \
		abstree.h transfor.h render.h graphics.h pixutils.h \
		script.h
		$(CC) -c $(CFLAGS) outfile.c

script.o:	script.c script.h rasmol.h molecule.h command.h \
		abstree.h transfor.h render.h graphics.h pixutils.h
		$(CC) -c $(CFLAGS) script.c


cflow:
		cflow -I/usr/local/include $(SRCS)

clean:
		rm -f rasmol $(OBJS)
