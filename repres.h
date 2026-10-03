/* repres.h
 * RasMol2 Molecular Graphics
 * Roger Sayle, August 1995
 * Version 2.6
 */

#define DotMax    100
typedef struct _DotStruct {
        struct _DotStruct __far *next;
        short col[DotMax];
        Long xpos[DotMax];
        Long ypos[DotMax];
        Long zpos[DotMax];
        int count;
    } DotStruct;


typedef struct _Monitor {
        struct _Monitor *next;
        Atom __far *src;
        Atom __far *dst;
        unsigned short dist;
        short col;
    } Monitor;


typedef struct _Label {
        struct _Label *next;
        Long  refcount;
        char *label;
    } Label;



#ifdef REPRES
DotStruct __far *DotPtr;
Monitor *MonitList;
Label *LabelList;

int CartoonHeight;
int SolventDots;
int ProbeRadius;

int DrawDots,DrawLabels;
int DrawMonitDistance;
int DrawBetaArrows;

#else
extern DotStruct __far *DotPtr;
extern Monitor *MonitList;
extern Label *LabelList;

extern int CartoonHeight;
extern int ProbeRadius;
extern int SolventDots;

extern int DrawDots,DrawLabels;
extern int DrawMonitDistance;
extern int DrawBetaArrows;


#ifdef FUNCPROTO
int DeleteLabels();
void DeleteLabel( Label* );
Label *CreateLabel( char*, int );
void DefineLabels( char* );
void DefaultLabels( int );
void DisplayLabels();

void DeleteMonitors();
void AddMonitors( Atom __far*, Atom __far* );
void CreateMonitor( Long, Long );
void DisplayMonitors();

void DeleteSurface();
void CalculateSurface( int );
void DisplaySurface();

/* Ribbons & Cartoons */
void DisplayRibbon( Chain __far* );

void ResetRepres();
void InitialiseRepres();

#else /* non-ANSI C compiler */
int DeleteLabels();
void DeleteLabel();
Label *CreateLabel();
void DefineLabels();
void DefaultLabels();
void DisplayLabels();

void DeleteMonitors();
void AddMonitors();
void CreateMonitor();
void DisplayMonitors();

void DeleteSurface();
void CalculateSurface();
void DisplaySurface();

/* Ribbons & Cartoons */
void DisplayRibbon();

void ResetRepres();
void InitialiseRepres();

#endif
#endif

