/* infile.h
 * RasMol2 Molecular Graphics
 * Roger Sayle, August 1995
 * Version 2.6
 */

#ifdef INFILE

#else

#ifdef FUNCPROTO
int LoadPDBMolecule( FILE*, int );
int LoadCIFMolecule( FILE* );

int SavePDBMolecule( char* );

#else /* non-ANSI C compiler */
int LoadPDBMolecule();
int LoadCIFMolecule();

int SavePDBMolecule();

#endif
#endif

