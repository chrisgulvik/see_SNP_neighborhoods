SHELL=	/bin/bash
CC=		gcc
CFLAGS=	-Wall -O2
PREFIX=	$(PWD)
export SAMDIR=$(PREFIX)/samtools-0.1.18


all:samtools-depend $(PREFIX)/bin/ttview $(PREFIX)/bin/bamsorted $(PREFIX)/bin/bam2wig


samtools-depend:
	cd $(PREFIX)/samtools-0.1.18/ && make && cd $(PREFIX)
	ln -sf $(PREFIX)/samtools-0.1.18/samtools-0.1.18 $(PREFIX)/bin/samtools-0.1.18


checksamenv:
	echo "Compiling with SAMDIR=$$SAMDIR"


$(PREFIX)/bin/ttview:$(PREFIX)/src/ttview.c $(PREFIX)/bin checksamenv
	$(CC) -o $@ -DSTANDALONE_VERSION  -I${SAMDIR} -L${SAMDIR} -L${SAMDIR}/bcftools  $<  ${SAMDIR}/bam2bcf.o  ${SAMDIR}/errmod.o  ${SAMDIR}/bam_color.o ${SAMDIR}/libbam.a -lbcf  -lm -lz
$(PREFIX)/bin/bamsorted:$(PREFIX)/src/bamsorted.c $(PREFIX)/bin checksamenv
	$(CC) ${CFLAGS} -o $@ -I ${SAMDIR} -L ${SAMDIR} $< -lbam -lz
$(PREFIX)/bin/bam2wig:$(PREFIX)/src/bam2wig.c $(PREFIX)/bin checksamenv
	$(CC) ${CFLAGS} -o $@ -I ${SAMDIR} -L ${SAMDIR} $< -lbam -lz


test:test-samtools $(PREFIX)/bin/bam2wig
	$(PREFIX)/bin/bam2wig ${SAMDIR}/examples/toy.bam
	$(PREFIX)/bin/bam2wig ${SAMDIR}/examples/toy.bam "ref2:10-20"

test-samtools:
	${SAMDIR}/samtools view -b ${SAMDIR}/examples/toy.sam -t ${SAMDIR}/examples/toy.fa -o ${SAMDIR}/examples/toy.bam
	${SAMDIR}/samtools index  ${SAMDIR}/examples/toy.bam


clean:
	rm -f $(PREFIX)/bin/samtools-0.1.18 $(PREFIX)/bin/ttview $(PREFIX)/bin/bamsorted $(PREFIX)/bin/bam2wig
	cd $(PREFIX)/samtools-0.1.18 && make clean
