/* command.h
 * RasMol2 Molecular Graphics
 * Roger Sayle, August 1995
 * Version 2.6
 */

#define MAXBUFFLEN   256
#define MAXLINELEN   256

/* Format values are related to Tokens */
#define Tok2Format(x) ((x)-349)
#define Format2Tok(x) ((x)+349)

#define FormatPDB        1
#define FormatMacroMod   2
#define FormatGaussian   3
#define FormatAlchemy    4
#define FormatNMRPDB     5
#define FormatCharmm     6
#define FormatBiosym     7
#define FormatMOPAC      8
#define FormatSHELX      9
#define FormatMol2      10
#define FormatFDAT      11
#define FormatMMDB      12
#define FormatMDL       13
#define FormatXYZ       14
#define FormatCIF       15
#define FormatCEX       16

#define IPC_Ok      0
#define IPC_Error   1
#define IPC_Exit    2
#define IPC_Quit    3

#ifdef COMMAND
int DataFileFormat;
char DataFileName[256];
char CurLine[MAXBUFFLEN];
int CurState,StateOption;
int CommandActive;
Long SelectCount;
int Interactive;
int FileDepth;
int IsPaused;

#else
extern int DataFileFormat;
extern char DataFileName[256];
extern char CurLine[MAXBUFFLEN];
extern int CurState,StateOption;
extern int CommandActive;
extern Long SelectCount;
extern int Interactive;
extern int FileDepth;
extern int IsPaused;

#ifdef FUNCPROTO
int ProcessCharacter( int );
int FetchFile( int, int, char* );
int GuessMoleculeFormat( char* );
void LoadScriptFile( FILE*, char* );
void ResetCommandLine( int );
void InitialiseCommand();
int ExecuteIPCCommand( char __huge* );
int ExecuteCommand();
void ZapDatabase();

#else /* non-ANSI C compiler */
int ProcessCharacter();
int FetchFile();
int GuessMoleculeFormat();
void LoadScriptFile();
void ResetCommandLine();
void InitialiseCommand();
int ExecuteIPCCommand();
int ExecuteCommand();
void ZapDatabase();

#endif
#endif

