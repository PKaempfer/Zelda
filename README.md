# The Zelda Genome Assembler

Motivation: Currently, the reconstruction of genomic information using short-read data utilizes two distinct
graph-based approaches, namely the overlap-layout-consensus concept (OLC) and de Bruijn graphs
(dBG). The overlap-layout-consensus concept is considered superior to de Bruijn graphs in which the unit
of assembly is a read as opposed to a small k-mer, causing the graph and its path structure to be simpler
and easier to disambiguate, together resulting in higher contig lengths. However, popular NGS assemblers
rely solely on de Bruijn graphs due to their superior runtime efficiency compared to the quadratic nature
of the OLC approach.

Results: Here we present Zelda, a new hybrid graph approach. Zelda links the fast, linear-time construction
of a dBG together with the higher contig resolution of the OLC approach using string graphs. This is
accomplished by touring the dBG and collecting read information while simultaneously constructing an
overlap graph directly from an expanded dBG which contains information on the original sequencing
reads prior their decomposition into k-mers. By further deriving a string graph from the transitively reduced
overlap graph, Zelda is able to reconstruct large unique contigs of the genome. The Zelda assembler can
be used with various sequencing technologies with low error rates. Tests on several short read data sets
have shown this hybrid approach is comparable, in time and space, with state of the art dBG assemblers,
but provides higher contig resolution as expected from OLC based assemblers, while avoiding the long
run-time and storage consumption of the dBG approach.

Availability: The Zelda assembler is implemented in C for x86-64bit processor architecture and UNIX
based operating systems. The source code and binaries are freely available under GNU Public License.
