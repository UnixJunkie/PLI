CC = gcc
CFLAGS = -O2
LIBS = -lz -lm

OBJS =	genio.o		pdbio.o		protein.o	molecule.o	math.o \
	map.o		contacts.o	voronoi.o	settings.o 	atom_type.o \
	rings.o		atomgeom.o	mdlio.o		hbonds.o	molio.o \
	score.o		ff.o		histograms.o	mapio.o		field.o	\
	lebedev.o	memory.o	minimise.o	system.o	bfactor.o \
	resolve.o	pliff_score.o	sybylio.o	site.o		water.o	\
	genff.o		constraints.o	torsion.o	molint.o	sysmols.o \
	dofs.o 		params.o	modes.o		sfuncs.o	pli.o

obj/%.o: src/%.c src/pli.h | obj
	$(CC) $(CFLAGS) -c $< -o $@

pli: bin/pli

bin/pli: $(addprefix obj/,$(OBJS)) | bin
	$(CC) $(CFLAGS) -o bin/pli $(addprefix obj/,$(OBJS)) $(LIBS)

obj:
	mkdir -p obj

bin:
	mkdir -p bin

clean:
	/bin/rm -f bin/pli obj/*.o
