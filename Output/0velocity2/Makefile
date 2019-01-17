all: MiniSimmons

MiniSimmons: MiniSimmons.c
	qcc -g -O2 MiniSimmons.c -o MiniSimmons -I${BASILISK} -L${BASILISK}/gl -lglutils -lfb_osmesa -lOSMesa -lGLU -lm

clean:
	rm -f MiniSimmons
