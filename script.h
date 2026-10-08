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
int WriteScriptFile( char* );

#else /* non-ANSI C compiler */
int WriteScriptFile();

#endif
#endif

