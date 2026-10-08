/* outfile.h
 * RasMol2 Molecular Graphics
 * Roger Sayle, August 1995
 * Version 2.6
 */

#ifdef OUTFILE
int UseTransparent;
int UseOutLine;

#else
extern int UseTransparent;
extern int UseOutLine;

#ifdef FUNCPROTO
int WriteVectPSFile( char* );
int WritePPMFile( char*, int );
void InitialiseOutFile();

#else /* non-ANSI C compiler */
int WriteVectPSFile();
int WritePPMFile();
void InitialiseOutFile();

#endif
#endif

