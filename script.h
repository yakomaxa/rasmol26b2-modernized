/* script.h
 * RasMol2 Molecular Graphics
 * Roger Sayle, August 1995
 * Version 2.6
 */

#ifdef SCRIPT
int KinemageFlag;

#else
extern int KinemageFlag;

#ifdef FUNCPROTO
int WriteMolScriptFile( char* );
int WriteKinemageFile( char* );
int WriteScriptFile( char* );
int WritePOVRayFile( char* );

#else /* non-ANSI C compiler */
int WriteMolScriptFile();
int WriteKinemageFile();
int WriteScriptFile();
int WritePOVRayFile();

#endif
#endif

