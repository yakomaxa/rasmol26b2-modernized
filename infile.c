/* infile.c
 * RasMol2 Molecular Graphics
 * Roger Sayle, August 1995
 * Version 2.6
 */
#include "rasmol.h"

#ifdef IBMPC
#include <windows.h>
#include <malloc.h>
#endif
#ifdef APPLEMAC
#include <Types.h>
#endif
#ifndef sun386
#include <stdlib.h>
#endif

#include <string.h>
#include <ctype.h>
#include <stdio.h>
#include <math.h>

#define INFILE
#include "infile.h"
#include "molecule.h"
#include "abstree.h"
#include "command.h"
#include "transfor.h"

#ifndef APPLEMAC
#ifndef IBMPC
#ifndef VMS
#include <sys/types.h>
#endif
#include <sys/time.h>
#endif
#include <time.h>
#endif
 

#ifdef MMIO
#include "mmio.h"
#endif

#define GroupPool    8

#define FeatHelix    1
#define FeatSheet    2
#define FeatTurn     3


typedef struct {
	int init, term;
	char chain;
	char type;
	} FeatEntry;

#define FeatSize    32
typedef struct _Feature {
	struct _Feature __far *fnext;
        FeatEntry data[FeatSize];
        int count;
    } Feature;


typedef struct {
	  char src[4];
	  char dst[4];
	  } ConvTable;

#define MAXALCATOM   5
static ConvTable AlcAtomTable[MAXALCATOM] = {
    { { 'S', 'O', '2', ' ' }, { ' ', 'S', '2', ' ' } },  /*  1 */
    { { 'C', 'A', 'R', ' ' }, { ' ', 'C', ' ', ' ' } },  /*  2 */
    { { 'N', 'A', 'R', ' ' }, { ' ', 'N', ' ', ' ' } },  /*  3 */
    { { 'N', 'A', 'M', ' ' }, { ' ', 'N', ' ', ' ' } },  /*  4 */
    { { 'N', 'P', 'L', '3' }, { ' ', 'N', '3', ' ' } },  /*  5 */
				 };

static char PDBInsert;
static Feature __far *FeatList;
static char Record[202];
static FILE *DataFile;

/* Macros for commonly used loops */
#define ForEachAtom  for(chain=Database->clist;chain;chain=chain->cnext) \
		     for(group=chain->glist;group;group=group->gnext)    \
		     for(aptr=group->alist;aptr;aptr=aptr->anext)
#define ForEachBond  for(bptr=Database->blist;bptr;bptr=bptr->bnext)


/* Forward Reference */
void DestroyDatabase();


#ifdef APPLEMAC
/* External RasMac Function Declaration! */
void SetFileInfo( char*, OSType, OSType, short );
#endif


static void FatalInFileError(ptr)
    char *ptr;
{
    char buffer[80];

    sprintf(buffer,"InFile Error: %s!",ptr);
    RasMolFatalExit(buffer);
}


/*================================*/
/* File/String Handling Functions */
/*================================*/

static int FetchRecord()
{
    register char *ptr;
    register int ch;

    if( feof(DataFile) )
    {   *Record = '\0';
        return( False );
    }

    ptr = Record;
    do {
        ch = getc(DataFile);
        if( ch == '\n' )
        {   *ptr = 0;
            return( True );
        } else if( ch == '\r' )
        {   ch = getc(DataFile);
            if( ch != '\n' )
                ungetc(ch,DataFile);
            *ptr = 0;
            return( True );
        } else if( ch == EOF )
        {   *ptr = 0;
            return( ptr != Record+1 );
        } else *ptr++ = ch;
    } while( ptr < Record+200 );

    /* skip to the end of the line! */
    do { ch = getc(DataFile);
    } while( (ch!='\n') && (ch!='\r') && (ch!=EOF) );

    if( ch == '\r' )
    {   ch = getc(DataFile);
        if( ch != '\n' )
            ungetc(ch,DataFile);
    }
    *ptr = 0;
    return( True );
}


static void ExtractString( len, src, dst )
    int len;  char *src, *dst;
{
    register char *ptr;
    register char ch;
    register int i;

    ptr = dst;
    for( i=0; i<len; i++ )
    {   if( *src )
	{   ch = *src++;
            *dst++ = ch;
            if( ch != ' ' ) 
		ptr = dst;
	} else break;
    }
    *ptr = 0;
}


static Long ReadValue( pos, len )
    int pos, len;
{
    register Long result;
    register char *ptr;
    register char ch;
    register int neg;

    result = 0;
    neg = False;
    ptr = Record+pos;
    while( len-- )
    {   ch = *ptr++;
	if( (ch>='0') && (ch<='9') )
	{   result = (10*result)+(ch-'0');
	} else if( ch=='-' )
	    neg = True;
    }
    return( neg? -result : result );
}


/*===================================*/
/* File Format Independent Functions */
/*===================================*/

#ifdef FUNCPROTO
static void UpdateFeature( FeatEntry __far*, int );
#endif

static FeatEntry __far *AllocFeature()
{
    register Feature __far *ptr;
 
    if( !FeatList || (FeatList->count==FeatSize) )
    {   ptr = (Feature __far*)_fmalloc(sizeof(Feature));
        if( !ptr ) FatalInFileError("Memory allocation failed");
        /* Features are always deallocated! */
 
        ptr->fnext = FeatList;
        ptr->count = 0;
        FeatList = ptr;
    } else ptr = FeatList;
 
    return( &(ptr->data[ptr->count++]) );
}


static void UpdateFeature( ptr, mask )
    FeatEntry __far *ptr;  int mask;
{
    register Chain __far *chain;
    register Group __far *group;
 
    for( chain=Database->clist; chain; chain=chain->cnext )
        if( chain->ident == ptr->chain )
        {   group=chain->glist;
            while( group && (group->serno<ptr->init) )
                group = group->gnext;
 
            while( group && (group->serno<=ptr->term) )
            {   group->struc |= mask;
                group = group->gnext;
            }
 
            if( NMRModel )
            {  continue;
            } else return;
        }
}
 
 
static void ProcessFeatures()
{
    register Feature __far *next;
    register Feature __far *ptr;
    register int i;
 
    InfoTurnCount = 0;
    InfoHelixCount = 0;
    InfoLadderCount = 0;
    InfoStrucSrc = SourcePDB;
 
    for( ptr=FeatList; ptr; ptr=next )
    {    if( Database )
         {   for( i=0; i<ptr->count; i++ )
                 if( ptr->data[i].type==FeatHelix )
                 {   UpdateFeature( &ptr->data[i], HelixFlag );
                     InfoHelixCount++;
                 } else if( ptr->data[i].type==FeatSheet )
                 {   UpdateFeature( &ptr->data[i], SheetFlag );
                     InfoLadderCount++;
                 } else /* FeatTurn */
                 {   UpdateFeature( &ptr->data[i], TurnFlag );
                     InfoTurnCount++;
                 }
         }

         /* Deallocate Memory */
         next = ptr->fnext;
         _ffree( ptr );
    }
}
 
 


/*==============================*/
/* Molecule File Format Parsing */
/*==============================*/

static Long ReadPDBCoord( offset )
    int offset;
{
    register int len,neg;
    register Long result;
    register char *ptr;
    register char ch;
 
    result = 0;
    neg = False;
    len = 8;
 
    ptr = Record+offset;
    while( len-- )
    {   ch = *ptr++;
        if( (ch>='0') && (ch<='9') )
        {   result = (10*result)+(ch-'0');
        } else if( ch=='-' )
            neg = True;
    }
 
    /* Handle Chem3D PDB Files! */
    if( Record[offset+3]=='.' )
        result /= 10;
    return( neg? -result : result );
}


static void ProcessPDBGroup( heta, serno )
    int heta, serno;
{
    PDBInsert = Record[26];
    if( !CurChain || (CurChain->ident!=Record[21]) )
        CreateChain( Record[21] );
    CreateGroup( GroupPool );
 
    CurGroup->refno = FindResNo( Record+17 );
    CurGroup->serno = serno;
    ProcessGroup( heta );
}
 

static void ProcessPDBAtom( heta )
    int heta;
{
    register Bond __far *bptr;
    register Atom __far *ptr;
    register Long dx,dy,dz;
    register int temp,serno;
 
    dx = ReadPDBCoord(30);
    dy = ReadPDBCoord(38);
    dz = ReadPDBCoord(46);
 
    /* Process Pseudo Atoms Limits!! */
    if( (Record[13]=='Q') && (Record[12]==' ') )
    {   temp = (int)ReadValue(60,6);
        if( MMinMaxFlag )
        {   if( temp < MinMainTemp )
            {   MinMainTemp = temp;
            } else if( temp > MaxMainTemp )
                MaxMainTemp = temp;
        }
 
        /* Dummy co-ordinates! */
        if( (dx==dy) && (dx==dz) )
        {   if( !dx || (dx == 9999000L) )
                return;
        }
 
        if( HMinMaxFlag || MMinMaxFlag )
        {   if( dx < MinX )
            {   MinX = dx;
            } else if( dx > MaxX )
                MaxX = dx;
 
            if( dy < MinY )
            {   MinY = dy;
            } else if( dy > MaxY )
                MaxY = dy;
 
            if( dz < MinZ )
            {   MinZ = dz;
            } else if( dz > MaxZ )
                MaxZ = dz;
        }
        return;
    }
 
 
    /* Ignore XPLOR Pseudo Atoms!! */
    if( (dx==9999000L) && (dy==9999000L) && (dz==9999000L) )
        return;
 
    serno = (int)ReadValue(22,4);
    if( !CurGroup || (CurGroup->serno!=serno)
        || (CurChain->ident!=Record[21])
        || (PDBInsert!=Record[26]) )
        ProcessPDBGroup( heta, serno );
 
 
    ptr = CreateAtom();
    ptr->refno = ComplexAtomType(Record+12);
    ptr->serno = (int)ReadValue(6,5);
    ptr->temp = (int)ReadValue(60,6);
    ptr->altl = Record[16];
 
    ptr->xorg =  dx/4;
    ptr->yorg =  dy/4;
    ptr->zorg = -dz/4;
 
    if( heta ) ptr->flag |= HeteroFlag;
    ProcessAtom( ptr );
 
    /* Create biopolymer Backbone */
    if( IsAlphaCarbon(ptr->refno) && IsProtein(CurGroup->refno) )
    {   if( ConnectAtom )
        {   dx = ConnectAtom->xorg - ptr->xorg;
            dy = ConnectAtom->yorg - ptr->yorg;
            dz = ConnectAtom->zorg - ptr->zorg;
 
            /* Break backbone if CA-CA > 7.00A */
            if( dx*dx+dy*dy+dz*dz < (Long)1750*1750 )
            {   bptr = ProcessBond(ptr,ConnectAtom,NormBondFlag);
                bptr->bnext = CurChain->blist;
                CurChain->blist = bptr;
            } else ptr->flag |= BreakFlag;
        }
        ConnectAtom = ptr;
    } else if( IsSugarPhosphate(ptr->refno) && IsNucleo(CurGroup->refno) )
    {   if( ConnectAtom )
        {   bptr = ProcessBond(ConnectAtom,ptr,NormBondFlag);
            bptr->bnext = CurChain->blist;
            CurChain->blist = bptr;
        }
        ConnectAtom = ptr;
    }
}
 

static void ProcessPDBColourMask()
{
    register MaskDesc *ptr;
    register char *mask;
    register int i;
 
    if( MaskCount==MAXMASK )
        FatalInFileError("Too many COLOR records in file");
    ptr = &UserMask[MaskCount];
    mask = ptr->mask;
 
 
    ptr->flags = 0;
    for( i=6; i<11; i++ )
        if( (*mask++ = Record[i]) != '#' )
            ptr->flags |= SerNoFlag;
 
    for( i=12; i<20; i++ )
        *mask++ = Record[i];
    *mask++ = Record[21];
 
    for( i=22; i<26; i++ )
        if( (*mask++ = Record[i]) != '#' )
            ptr->flags |= ResNoFlag;
    *mask++ = Record[26];
 
    ptr->r = (int)(ReadPDBCoord(30)>>2) + 5;
    ptr->g = (int)(ReadPDBCoord(38)>>2) + 5;
    ptr->b = (int)(ReadPDBCoord(46)>>2) + 5;
    ptr->radius = (short)(5*ReadValue(54,6))>>1;
    MaskCount++;
}
 

int LoadPDBMolecule( fp, flag )
    FILE *fp;  int flag;
{
    register FeatEntry __far *ptr;
    register int srcatm, dstatm;
    register char *src, *dst;
    register int i,ignore;
 
    ignore = False;
    FeatList = (void __far*)0;
    DataFile = fp;
    NMRModel = 0;
 
    while( FetchRecord() )
    {   if( *Record == 'A' )
        {   if( !ignore && !strncmp("ATOM",Record,4) )
                ProcessPDBAtom( False );

        } else switch(*Record)
        {   case('C'):    if( !strncmp("CONE",Record,4) )
                          {   if( ignore || flag ) continue;
 
                              srcatm = (int)ReadValue(6,5);
                              if( srcatm )
                                  for( i=11; i<=36 && Record[i]; i+=5 )
                                  {   dstatm = (int)ReadValue(i,5);
                                      if( dstatm && (dstatm>srcatm) )
                                          CreateBondOrder(srcatm,dstatm);
                                  }
                               
                          } else if( !strncmp("COMP",Record,4) )
                          {   /* First or MOLECULE: COMPND record */
                              if( (Record[9]==' ') && 
                                  strncmp(Record+10,"MOL_ID:",7) )  
                              {   ExtractString(60,Record+10,InfoMoleculeName);
                              } else if( !InfoMoleculeName[0] &&
                                         !strncmp(Record+11,"MOLECULE: ",10) )
                                  ExtractString(49,Record+21,InfoMoleculeName);
                          } else if( !strncmp("CRYS",Record,4) )
                          {   dst = InfoSpaceGroup;  src=Record+55;
                              while( *src && src<Record+66 )
                                  if( *src!=' ' ) 
                                  {   *dst++ = *src++;
                                  } else src++;
                              *dst = 0;
 
                              InfoCellA = ReadValue( 6,9)/1000.0;
                              InfoCellB = ReadValue(15,9)/1000.0;
                              InfoCellC = ReadValue(24,9)/1000.0;
 
                              InfoCellAlpha = Deg2Rad*(ReadValue(33,7)/100.0);
                              InfoCellBeta =  Deg2Rad*(ReadValue(40,7)/100.0);
                              InfoCellGamma = Deg2Rad*(ReadValue(47,7)/100.0);
 
                          } else if( !strncmp("COLO",Record,4) )
                              ProcessPDBColourMask();
                          break;

            case('E'):    if( !strncmp("ENDM",Record,4) )
                          {   /* break after single model??? */
                              if( flag )
                              {   ConnectAtom = (void __far*)0;
                                  CurGroup = (void __far*)0;
                                  CurChain = (void __far*)0;
                              } else ignore = True;
 
                          } else if( !strncmp("END",Record,3) )
                              if( !Record[4] || (Record[4]==' ') )
                              {   /* Treat END same as TER! */
                                  ConnectAtom = (void __far*)0;
                                  CurGroup = (void __far*)0;
                                  CurChain = (void __far*)0;
                              }
                          break;

            case('H'):    if( !strncmp("HETA",Record,4) )
                          {   if( !ignore ) ProcessPDBAtom(True);
                          } else if( !strncmp("HELI",Record,4) )
                          {   if( ignore ) continue;
 
                              /* Remaining HELIX record fields   */
                              /* 38-39 .... Helix Classification */
                              /* 31 ....... Same Chain as 19?    */
                              ptr = AllocFeature();
                              ptr->type = FeatHelix;
                              ptr->chain = Record[19];
                              ptr->init = (int)ReadValue(21,4);
                              ptr->term = (int)ReadValue(33,4);
                              
                          } else if( !strncmp("HEAD",Record,4) )
                          {   ExtractString(40,Record+10,InfoClassification);
                              ExtractString( 4,Record+62,InfoIdentCode);
                          }
                          break;

            case('M'):    if( !strncmp("MODE",Record,4) )
                              if( flag ) NMRModel++;
                          break;
 
            case('S'):    if( !strncmp("SHEE",Record,4) )
                          {   if( ignore ) break;
                              /* Remaining SHEET record fields   */
                              /* 38-39 .... Strand Parallelism   */
                              /* 32 ....... Same Chain as 21?    */
                              ptr = AllocFeature();
                              ptr->type = FeatSheet;
                              ptr->chain = Record[21];
                              ptr->init = (int)ReadValue(22,4);
                              ptr->term = (int)ReadValue(33,4);
                          }
                          break;

            case('T'):    if( !strncmp("TURN",Record,4) )
                          {   if( ignore ) continue;
 
                              ptr = AllocFeature();
                              ptr->type = FeatTurn;
                              ptr->chain = Record[19];
                              ptr->init = (int)ReadValue(20,4);
                              ptr->term = (int)ReadValue(31,4);
                          } else if( !strncmp("TER",Record,3) )
                          {   if( !Record[3] || (Record[3]==' ') )
                              {   ConnectAtom = (void __far*)0;
                                  CurGroup = (void __far*)0;
                                  CurChain = (void __far*)0;
                              }
                          }
                          break;
        }
    }
 
    if( Database )
        strcpy(InfoFileName,DataFileName);
    if( FeatList ) ProcessFeatures();
    return( True );
}


int LoadMDLMolecule( fp )
    FILE *fp;
{
    register Bond __far *bptr;
    register Atom __far *src;
    register Atom __far *dst;
    register Atom __far *ptr;
 
    register int i,type;
    register int atoms, bonds;
    register int srcatm,dstatm;
    register Long dx, dy, dz;
    register Card dist2;
    register Real scale;
    register char *cptr;
 
    DataFile = fp;
 
    FetchRecord(); /* Molecule Name */
    ExtractString(78,Record,InfoMoleculeName);
 
    FetchRecord(); /* Program Details */
    FetchRecord(); /* Comments */
 
    FetchRecord();
    atoms = (int)ReadValue(0,3);
    bonds = (int)ReadValue(3,3);
 
    if( !atoms )
        return( False );
 
    CreateMolGroup();
    for( i=1; i<=atoms; i++ )
    {   FetchRecord();
        ptr = CreateAtom();
 
        cptr = Record+31;
        while( *cptr == ' ' ) cptr++;
        ptr->refno = SimpleAtomType(cptr);
 
        switch( (int)ReadValue(36,3) )
        {   case(1):  ptr->temp =  300;  break;
            case(2):  ptr->temp =  200;  break;
            case(3):  ptr->temp =  100;  break;
            case(5):  ptr->temp = -100;  break;
            case(6):  ptr->temp = -200;  break;
            case(7):  ptr->temp = -300;  break;
            default:  ptr->temp = 0;
        }
        ptr->serno = i;
 
        ptr->xorg =  ReadValue( 0,10)/40;
        ptr->yorg =  ReadValue(10,10)/40;
        ptr->zorg = -ReadValue(20,10)/40;
        ProcessAtom( ptr );
    }
 
    for( i=0; i<bonds; i++ )
    {   FetchRecord();
        srcatm = (int)ReadValue(0,3);
        dstatm = (int)ReadValue(3,3);
        type =   (int)ReadValue(6,3);
 
        if( type==2 )                 /* DOUBLE */
        {   CreateBond(srcatm,dstatm,DoubBondFlag);
        } else if( type==3 )          /* TRIPLE */
        {   CreateBond(srcatm,dstatm,TripBondFlag);
        } else if( type==4 )          /* AROMATIC */
        {   CreateBond(srcatm,dstatm,AromBondFlag);
        } else                        /* SINGLE */
            CreateBond(srcatm,dstatm,NormBondFlag);
    }
 
    for( bptr=Database->blist; bptr; bptr=bptr->bnext )
        if( bptr->flag & NormBondFlag )
        {   src = bptr->srcatom;
            dst = bptr->dstatom;
            if( (src->refno==2) && (dst->refno==2) )
            {   dx = dst->xorg - src->xorg;
                dy = dst->yorg - src->yorg;
                dz = dst->zorg - src->zorg;
                if( dx || dy || dz )
                {   dist2 = dx*dx + dy*dy + dz*dz;
                    scale = 385.0/sqrt(dist2);
                    break;
                }
            }
        }
 
    if( bptr )
    {   for( ptr=CurGroup->alist; ptr; ptr=ptr->anext )
        {   ptr->xorg = (Long)(ptr->xorg*scale);
            ptr->yorg = (Long)(ptr->yorg*scale);
            ptr->zorg = (Long)(ptr->zorg*scale);
        }
        MinX = (Long)(MinX*scale);  MaxX = (Long)(MaxX*scale);
        MinY = (Long)(MinY*scale);  MaxY = (Long)(MaxY*scale);
        MinZ = (Long)(MinZ*scale);  MaxZ = (Long)(MaxZ*scale);
    }
    return( True );
}
 

int LoadXYZMolecule( fp )
    FILE *fp;
{
    auto char type[12];
    auto double xpos, ypos, zpos;
    auto double charge, u, v, w;
    auto int atoms;
 
    register Atom __far *ptr;
    register char *src,*dst;
    register int i,count;
 
 
    DataFile = fp;
    /* Number of Atoms */
    FetchRecord();
    sscanf(Record,"%d",&atoms);
 
    /* Molecule (step) Description */
    FetchRecord();
    src = Record;
    while( *src == ' ' )
        src++;
 
    dst = InfoMoleculeName;
    for( i=0; i<78; i++ )
        if( *src ) *dst++ = *src++;
    *dst = '\0';
 
    if( atoms )
    {   CreateMolGroup();
        for( i=0; i<atoms; i++ )
        {   FetchRecord();
            ptr = CreateAtom();
            ptr->serno = i;
 
            xpos = ypos = zpos = 0.0;
            count = sscanf(Record,"%s %lf %lf %lf %lf %lf %lf %lf",
                           type, &xpos, &ypos, &zpos, &charge, &u, &v, &w );
 
            ptr->refno = SimpleAtomType(type);
            ptr->xorg =  (Long)(250.0*xpos);
            ptr->yorg =  (Long)(250.0*ypos);
            ptr->zorg = -(Long)(250.0*zpos);
 
            if( (count==5) || (count==8) )
            {   ptr->temp = (short)(100.0*charge);
            } else ptr->temp = 0;
            ProcessAtom( ptr );
        }
    }
    return( True );
}
 

static int FindSybylRefNo( ptr )
    char *ptr;
{
    register char *src,*dst;
    auto char name[4];
 
    src = ptr;
    dst = name;
    if( ptr[1] && (ptr[1]!='.') )
    {   *dst++ = ToUpper(*src);  src++;
        *dst++ = ToUpper(*src);  src++;
    } else
    {   *dst++ = ' ';
        *dst++ = ToUpper(*src);
        src++;
    }
 
    if( *src )
    {   src++;
 
        if( *src == 'a' )
        {   *dst++ = ' ';
            *dst = ' ';
        } else if( *src == 'p' )
        {   *dst++ = '3';
            *dst = ' ';
        } else
        {   *dst++ = *src++;
            if( *src && (*src!='+') )
            {   *dst = *src;
            } else *dst = ' ';
        }
    } else
    {   *dst++ = ' ';
        *dst = ' ';
    }
    return( NewAtomType(name) );
}
 
 
int LoadMol2Molecule( fp )
    FILE *fp;
{
    double xpos, ypos, zpos;
    int features, sets, serno;
    int atoms, bonds, structs;
    int srcatm, dstatm;
 
    char name[8];
    char type[8];
 
    register Atom __far *ptr;
    register char *src, *dst;
    register int i;
 
 
    DataFile = fp;
    while( FetchRecord() )
    {   if( !*Record || *Record=='#' )
            continue;
 
        if( !strncmp("@<TRIPOS>MOLECULE",Record,17) ||
            !strncmp("@MOLECULE",Record,9) )
        {   FetchRecord();  /* Molecule Name */
            src = Record;
            while( *src==' ' )
                src++;
 
            dst = InfoMoleculeName;
            while( (*dst++ = *src++) );
 
            FetchRecord();
            atoms = bonds = structs = features = sets = 0;
            sscanf(Record,"%d %d %d %d %d", &atoms, &bonds, &structs,
                                            &features, &sets );
 
            FetchRecord();  /* Molecule Type  */
            FetchRecord();  /* Charge Type    */
 
        } else if( !strncmp("@<TRIPOS>ATOM",Record,13) ||
                   !strncmp("@ATOM",Record,5) )
        {   if( !atoms ) continue;
 
            CreateMolGroup();
            for( i=0; i<atoms; i++ )
            {    FetchRecord();
                 ptr = CreateAtom();
 
                 sscanf(Record,"%d %s %lf %lf %lf %s", &serno, name,
                                &xpos, &ypos, &zpos, type );
 
                 ptr->refno = FindSybylRefNo( type );
                 ptr->serno = serno;
                 /* ptr->serno = i; */
 
                 ptr->xorg =  (Long)(250.0*xpos);
                 ptr->yorg =  (Long)(250.0*ypos);
                 ptr->zorg = -(Long)(250.0*zpos);
                 ProcessAtom( ptr );
            }
 
        } else if( !strncmp("@<TRIPOS>BOND",Record,13) ||
                   !strncmp("@BOND",Record,5) )
            for( i=0; i<bonds; i++ )
            {   FetchRecord();
                sscanf(Record,"%d %d %d %s",&serno,&srcatm,&dstatm,type);
                if( !strncmp(type,"ar",2) )
                {   CreateBond(srcatm,dstatm,AromBondFlag);
                } else if( *type == '2' )
                {   CreateBond(srcatm,dstatm,DoubBondFlag);
                } else /* *type == '1' */
                    CreateBond(srcatm,dstatm,NormBondFlag);
            }
    }
    return( True );
}


static int FindAlchemyRefNo()
{
    register char *ptr;
    register int i;
    char name[4];
 
    ptr = Record+6;
    if( !isalpha(ptr[1]) )
    {   name[0] = ' ';
        for( i=0; i<3; i++ )
            name[i+1] = ToUpper(ptr[i]);
        ptr = name;
    } else
    {   for( i=0; i<4; i++ )
            ptr[i] = ToUpper(ptr[i]);
 
        for( i=0; i<MAXALCATOM; i++ )
            if( !strncmp(AlcAtomTable[i].src,ptr,4) )
            {   ptr = AlcAtomTable[i].dst;
                break;
            }
    }
    return( NewAtomType(ptr) );
}
 

int LoadAlchemyMolecule( fp )
    FILE *fp;
{
    auto int serno,srcatm,dstatm;
    register Atom __far *ptr;
    register int atoms, bonds;
    register char *chptr;
    register int i;
 
    DataFile = fp;
    FetchRecord();
    atoms = (int)ReadValue(0,5);
    bonds = (int)ReadValue(13,5);
    ExtractString(38,Record+41,InfoMoleculeName);
 
    if( !atoms )
        return( False );
 
    CreateMolGroup();
    for( i=0; i<atoms; i++ )
    {   FetchRecord();
        ptr = CreateAtom();
 
        ptr->refno = FindAlchemyRefNo();
        ptr->temp = (int)ReadValue(40,8);
        ptr->serno = (int)ReadValue(0,5);
        /* ptr->serno = i+1; */
 
        ptr->xorg =  ReadValue(12,7)/4;
        ptr->yorg =  ReadValue(21,7)/4;
        ptr->zorg = -ReadValue(30,7)/4;
        ProcessAtom( ptr );
    }
 
    for( i=0; i<bonds; i++ )
    {   FetchRecord();
        sscanf(Record,"%d %d %d",&serno,&srcatm,&dstatm);
 
        chptr = Record;
        while( *chptr && !isalpha(*chptr) )
            chptr++;
 
        if( *chptr =='A' )             /* AROMATIC */
        {   CreateBond(srcatm,dstatm,AromBondFlag);
        } else if( *chptr == 'D' )     /* DOUBLE */
        {   CreateBond(srcatm,dstatm,DoubBondFlag);
        } else if( *chptr == 'T' )     /* TRIPLE */
        {   CreateBond(srcatm,dstatm,TripBondFlag);
        } else if( *chptr == 'H' )     /* HYDROGEN */
        {   CreateBond(srcatm,dstatm,HydrBondFlag);
        } else /* (*chptr == 'S') */   /* SINGLE */
            CreateBond(srcatm,dstatm,NormBondFlag);
    }
    return( True );
}
 
 
int LoadCharmmMolecule( fp )
    FILE *fp;
{
    auto char buffer[4];
    register Atom __far *ptr;
    register int resno,serno;
    register int atoms,i;
    register int chain;

    DataFile = fp;
 
    do {
        FetchRecord();
    } while( *Record=='*' );
    atoms = (int)ReadValue(0,5);
    if( !atoms ) return False;
 
    MinHetaRes = MaxHetaRes = 0;
    strcpy(InfoFileName,DataFileName);
    MainGroupCount = 0;
 
    chain = 0;
    CurChain = NULL;
    for( serno=0; serno<atoms; serno++ )
    {   FetchRecord();
 
        if( !CurChain || strncmp(Record+51,buffer,4) )
        {   for( i=0; i<4; i++ )
                buffer[i] = Record[51+i];
            CreateChain(chain+49);
            chain++;
        }

        /* Non-standard Sequential Residue Numbering */
        /* resno = (int)ReadValue(5,5); */

        resno = (int)ReadValue(56,5);
        if( !CurGroup || (CurGroup->serno!=resno) )
        {   CreateGroup( GroupPool );
            CurGroup->refno = FindResNo(Record+11);
            CurGroup->serno = resno;
            ProcessGroup( False );
        }
 
        ptr = CreateAtom();
        ptr->refno = ComplexAtomType(Record+15);
        ptr->temp = (int)ReadValue(60,9);
        ptr->serno = (int)ReadValue(0,5);
        /* ptr->serno = serno+1; */
 
        ptr->xorg =  ReadValue(20,8)/4;
        ptr->yorg =  ReadValue(30,8)/4;
        ptr->zorg = -ReadValue(40,8)/4;
        ProcessAtom( ptr );
    }
    return( True );
}
 

static int MOPACAtomType( type )
    char *type;
{
    auto char name[4];
    register char ch1,ch2;
    register int i;
 
    if( *type == ' ' )
        type++;
 
    name[2] = name[3] = ' ';
    if( isdigit(type[0]) )
    {   i = *type++ - '0';
        while( isdigit(*type) )
            i = (10*i) + (*type++ - '0');
 
        /* Last Atom in File! */
        if( i == 0 )
        {   return( -1 );
        } else if( i >= 99 )
            return( 1 );
 
        /* Construct Name */
        ch1 = Element[i].symbol[0];
        ch2 = Element[i].symbol[1];
        if( ch2 == ' ' )
        {   name[1] = ch1;
            name[0] = ' ';
        } else
        {   name[1] = ToUpper(ch2);
            name[0] = ch1;
        }
 
    } else
    {   ch1 = ToUpper(type[0]);
        ch2 = ToUpper(type[1]);
        if( (ch1=='X') || (ch1=='+') || (ch1=='-') )
        {   return( 1 );
        } else if( (ch1=='T') && (ch2=='V') )
            return( 1 );
 
        if( ch2 && (ch2!=' ') && (ch2!='(') && !isdigit(ch2) )
        {   name[0] = ch1;
            name[1] = ch2;
        } else
        {   name[1] = ch1;
            name[0] = ' ';
        }
    }
    return( NewAtomType(name) );
}
 
 
static int ReadMOPACOutputFile()
{
    register Atom __far *atm;
    register int i,init;
    register char *ptr;
 
    init = False;
    while( FetchRecord() )
    {   ptr = Record;
        while( *ptr == ' ' )
            ptr++;
 
        if( !strncmp(ptr,"CARTESIAN COORDINATES",21) )
        {   for( i=0; i<3; i++ )
                FetchRecord();
 
            if( Database )
            {   atm = CurGroup->alist;
                MMinMaxFlag = False;
                HasHydrogen = False;
                MainAtomCount = 0;
            }
 
            while( FetchRecord() && *Record && isdigit(Record[5]) )
            {   if( !Database )
                {   atm = (Atom __far*)0;
                    CreateMolGroup();
                    init = True;
                }
 
                if( !atm )
                {   atm = CreateAtom();
                    atm->serno = (int)ReadValue(0,6);
                    atm->refno = MOPACAtomType(Record+14);
                    atm->temp = 0;
 
                    atm->xorg = ReadValue(20,10)/40;
                    atm->yorg = ReadValue(30,10)/40;
                    atm->zorg = ReadValue(40,10)/40;
                } else
                {   atm->xorg = ReadValue(30,10)/40;
                    atm->yorg = ReadValue(40,10)/40;
                    atm->zorg = ReadValue(50,10)/40;
                }
                ProcessAtom(atm);
                atm = atm->anext;
            }
             
        } else if( !strncmp(ptr,"NET ATOMIC CHARGES",18) )
        {   FetchRecord();
            FetchRecord();
 
            if( Database )
            {   atm = CurGroup->alist;
                MMinMaxFlag = False;
                HasHydrogen = False;
                MainAtomCount = 0;
            }
 
            while( FetchRecord() && strncmp(Record," DIPOLE",7) )
            {   if( !Database )
                {   atm = (Atom __far*)0;
                    CreateMolGroup();
                }
 
                if( !atm )
                {   atm = CreateAtom();
                    atm->serno = (int)ReadValue(0,12);
                    atm->refno = MOPACAtomType(Record+21);
                    atm->temp = (int)(ReadValue(27,13)/100);
                    atm->xorg = atm->yorg = atm->zorg = 0;
                } else
                    atm->temp = (int)(ReadValue(27,13)/100);
                ProcessAtom(atm);
                atm = atm->anext;
            }
        }
    }
 
    if( !init )
    {   if( Database )
            DestroyDatabase();
        return( False );
    } else return( True );
}
 
 
static int MoreMOPACKeywords()
{
    register char *ptr;
 
    ptr = Record;
    while( *ptr )
    {   if( *ptr == '+' )
            if( !ptr[1] || (ptr[1]==' ') )
                return(True);
 
        /* Skip Next Keyword */
        while( *ptr && *ptr!=' ' ) ptr++;
        while( *ptr == ' ' ) ptr++;
    }
    return( False );
}
 
 
int LoadMOPACMolecule( fp )
    FILE *fp;
{
    static int na,nb,nc,lopt;
    static double dist,angle,dihed;
    register IntCoord __far *coord;
    register Atom __far *aptr;
    register int count,refno;
    register int cartflag;
    register char *ptr;
 
    DataFile = fp;
    FetchRecord();
 
    /* Test for MOPAC output file */
    if( !strncmp(Record," ***",4) )
        return( ReadMOPACOutputFile() );
 
    /* MOPAC Keywords */
    while( MoreMOPACKeywords() )
        FetchRecord();
 
    FetchRecord(); /* Molecule Name */
    ExtractString(78,Record,InfoMoleculeName);
    FetchRecord(); /* Comments */
 
    count = 0;
    cartflag = False;
    InitInternalCoords();
    while( FetchRecord() )
    {   /* Process Record! */
        for( ptr=Record; *ptr; ptr++ )
            if( (*ptr==',') || (*ptr==0x09) )
                *ptr = ' ';
 
        ptr = Record;
        while( *ptr == ' ' )
           ptr++;
 
        if( !*ptr ) break;
        refno = MOPACAtomType(ptr);
        if( refno == -1 ) break;
 
        while( *ptr && (*ptr!=' ') )
            if( *ptr == '(' )
            {   /* Skip Atom Label */
                while( *ptr && (*ptr!=')') )
                    ptr++;
            } else ptr++;
 
        na = nb = nc = 0;
        dist = angle = dihed = 0.0;
        sscanf(ptr,"%lf %*d %lf %*d %lf %d %d %d %d",
               &dist, &angle, &dihed, &lopt, &na, &nb, &nc );
        count++;
 
        if( count == 3 )
        {   /* Handle missing dihedral */
            if( lopt == 2 )
            {   na = 1;  nb = 2;
                dihed = 0.0;
            } else if( lopt == 1)
            {   /* Safe FP comparison for Cartesian */
                if( (dihed>=1.9999) && (dihed<=2.0001) )
                {   na = 2;  nb = 1;
                    dihed = 0.0;
                }
            }
        } else if( count == 4 )
            cartflag = (na == 0);
 
        coord = AllocInternalCoord();
        coord->na = na; coord->nb = nb; coord->nc = nc;
        coord->refno = refno;
        coord->angle = angle;
        coord->dihed = dihed;
        coord->dist = dist;
    }
 
    if( !count )
        return( False );
 
    /* Co-ordinate conversion! */
    if( !cartflag )
        if( !ConvertInternal2Cartesian() )
        {   FreeInternalCoords();
            if( CommandActive )
                WriteChar('\n');
            WriteString("Error: Invalid MOPAC z-matrix file!\n\n");
            CommandActive=False;
            return( False );
        }
 
    count = 0;
    for( coord=IntList; coord; coord=coord->inext )
        if( coord->refno != 1 )
        {   if( !Database )
                CreateMolGroup();
 
            aptr = CreateAtom();
            aptr->refno = (Byte)coord->refno;
            aptr->serno = ++count;
            aptr->temp = 0;
 
            aptr->xorg = (Long)(250.0*coord->dist);
            aptr->yorg = (Long)(250.0*coord->angle);
            aptr->zorg = (Long)(250.0*coord->dihed);
            ProcessAtom(aptr);
        } else count++;
 
    FreeInternalCoords();
    return( True );
}
 

int LoadMacroModelMolecule( fp )
    FILE *fp;
{
#ifdef MMIO
    return( True );
#else
    register char *src,*dst;
    register int i;
    auto int atoms;
 
    DataFile = fp;
 
    /* Number of Atoms & Description */
    FetchRecord();
    sscanf(Record,"%d",&atoms);
 
    src = Record;
    dst = InfoMoleculeName;
    for( i=0; i<78; i++ )
        if( *src ) *dst++ = *src++;
    *dst = '\0';
 
    for( i=0; i<atoms; i++ )
    {    FetchRecord();
    }
    return( True );
#endif
}


int LoadBiosymMolecule( fp )
    FILE *fp;
{
    return( True );
}
 
 
int LoadSHELXMolecule( fp )
    FILE *fp;
{
    return( True );
}
 

int LoadFDATMolecule( fp )
    FILE *fp;
{
    return( True );
}


/*====================================*/
/* mmCIF (Macromolecular CIF) Parsing */
/*====================================*/

/* Minimal, dependency-free parser for the subset of mmCIF used by the
 * Protein Data Bank: "data_" blocks containing scalar "_category.item
 * value" pairs and "loop_" tables of "_category.item" tags followed by
 * whitespace/quote-delimited rows of values. Only the "_atom_site" loop
 * is interpreted structurally; everything else is tokenised and skipped,
 * except for a handful of scalar items used to populate the same summary
 * fields the PDB loader fills in from HEADER/COMPND/CRYST1.
 */

#define CIFTokMax   512
#define CIFMaxCols   48

static char CIFTokBuf[CIFTokMax];
static int  CIFTokPushed;
static int  CIFAtBOL;
static char CIFInsert;

enum { ACOL_GROUP, ACOL_ID, ACOL_TYPE, ACOL_ATOM, ACOL_ALT, ACOL_COMP,
       ACOL_ASYM, ACOL_SEQ, ACOL_INS, ACOL_X, ACOL_Y, ACOL_Z, ACOL_TEMP,
       ACOL_AUTHSEQ, ACOL_AUTHCOMP, ACOL_AUTHASYM, ACOL_AUTHATOM,
       ACOL_MODEL, ACOL_MAX };

static char *CIFAtomSiteItem[ACOL_MAX] = {
    "group_PDB", "id", "type_symbol", "label_atom_id", "label_alt_id",
    "label_comp_id", "label_asym_id", "label_seq_id", "pdbx_PDB_ins_code",
    "Cartn_x", "Cartn_y", "Cartn_z", "B_iso_or_equiv", "auth_seq_id",
    "auth_comp_id", "auth_asym_id", "auth_atom_id", "pdbx_PDB_model_num" };


/* Case insensitive comparison of two NUL terminated strings */
static int CIFEqual( a, b )
    char *a, *b;
{
    while( *a && *b )
    {   if( ToUpper(*a) != ToUpper(*b) )
            return False;
        a++;  b++;
    }
    return( *a == *b );
}


/* Case insensitive comparison of the first n characters */
static int CIFEqualN( a, b, n )
    char *a, *b;  int n;
{
    while( n-- )
    {   if( !*a || !*b )
            return( *a == *b );
        if( ToUpper(*a) != ToUpper(*b) )
            return False;
        a++;  b++;
    }
    return True;
}


/* Return the item part of a "_category.item" tag */
static char *CIFItemName( tag )
    char *tag;
{
    register char *ptr;

    ptr = strchr(tag,'.');
    return( ptr? ptr+1 : tag );
}


/* True if tok cannot be a loop/scalar data value, i.e. it starts a new
 * tag, loop, data block, save frame or other top level construct.
 */
static int CIFIsReserved( tok )
    char *tok;
{
    if( tok[0] == '_' )                 return True;
    if( CIFEqual(tok,"loop_") )         return True;
    if( CIFEqual(tok,"stop_") )         return True;
    if( CIFEqual(tok,"global_") )       return True;
    if( CIFEqualN(tok,"data_",5) )      return True;
    if( CIFEqualN(tok,"save_",5) )      return True;
    return False;
}


/* Fetch the next raw CIF token from DataFile into CIFTokBuf, handling
 * '#' comments, '...'/"..." quoted strings (where the closing quote
 * must be followed by whitespace, per the CIF grammar) and ';'
 * semicolon delimited multi-line text fields. Returns False on EOF.
 */
static int CIFRawToken()
{
    register int ch;
    register int len;
    register int q;
    register int atLineStart;

    for(;;)
    {   for(;;)
        {   ch = getc(DataFile);
            if( ch == EOF )
                return False;
            if( ch == '\n' )
            {   CIFAtBOL = True;
                continue;
            }
            if( (ch==' ') || (ch=='\t') || (ch=='\r') )
                continue;
            if( ch == '#' )
            {   do { ch = getc(DataFile);
                } while( (ch!=EOF) && (ch!='\n') );
                if( ch == '\n' ) CIFAtBOL = True;
                if( ch == EOF )  return False;
                continue;
            }
            break;
        }

        if( CIFAtBOL && (ch==';') )
        {   /* Semicolon delimited text field */
            CIFAtBOL = False;
            len = 0;
            atLineStart = False;
            for(;;)
            {   ch = getc(DataFile);
                if( ch == EOF )
                {   CIFAtBOL = True;
                    break;
                }
                if( atLineStart && (ch==';') )
                {   while( (ch=getc(DataFile))!=EOF && (ch!='\n') )
                        ;
                    CIFAtBOL = True;
                    break;
                }
                if( len < CIFTokMax-1 )
                    CIFTokBuf[len++] = (char)ch;
                atLineStart = (ch=='\n');
            }
            CIFTokBuf[len] = '\0';
            return True;
        }

        CIFAtBOL = False;

        if( (ch=='\'') || (ch=='"') )
        {   register int nc;

            q = ch;
            len = 0;
            for(;;)
            {   ch = getc(DataFile);
                if( (ch==EOF) || (ch=='\n') )
                    break;
                if( ch == q )
                {   nc = getc(DataFile);
                    if( (nc==EOF) || (nc==' ') || (nc=='\t') ||
                        (nc=='\n') || (nc=='\r') )
                    {   if( nc == '\n' ) CIFAtBOL = True;
                        break;
                    } else
                    {   if( len < CIFTokMax-1 )
                            CIFTokBuf[len++] = (char)ch;
                        ungetc(nc,DataFile);
                        continue;
                    }
                }
                if( len < CIFTokMax-1 )
                    CIFTokBuf[len++] = (char)ch;
            }
            CIFTokBuf[len] = '\0';
            return True;
        }

        /* Bare (unquoted) token */
        len = 0;
        if( len < CIFTokMax-1 )
            CIFTokBuf[len++] = (char)ch;
        for(;;)
        {   ch = getc(DataFile);
            if( (ch==EOF) || (ch==' ') || (ch=='\t') ||
                (ch=='\n') || (ch=='\r') )
            {   if( ch == '\n' ) CIFAtBOL = True;
                break;
            }
            if( len < CIFTokMax-1 )
                CIFTokBuf[len++] = (char)ch;
        }
        CIFTokBuf[len] = '\0';
        return True;
    }
}


static int CIFNextToken()
{
    if( CIFTokPushed )
    {   CIFTokPushed = False;
        return True;
    }
    return CIFRawToken();
}


static void CIFPushToken()
{
    CIFTokPushed = True;
}


/* Copy CIFTokBuf into a CIFTokMax sized row slot, NUL terminated */
static void CIFStoreTok( dst )
    char *dst;
{
    register int len;

    len = strlen(CIFTokBuf);
    if( len > CIFTokMax-1 )
        len = CIFTokMax-1;
    memcpy(dst,CIFTokBuf,len);
    dst[len] = '\0';
}


/* Build the 4 character fixed-column atom name RasMol's PDB parser
 * would have seen, from a free-form mmCIF atom name and its explicit
 * element symbol. The element occupies column 13 only when it takes
 * two characters (e.g. "CA" Calcium, "ZN", "FE"); single character
 * elements (C, N, O, S, P, H, ...) leave column 13 blank unless the
 * full atom name needs all four columns (e.g. hydrogen "HG21"). This
 * matches the convention ComplexAtomType()/GetElemNumber() assume.
 */
static void FormatCIFAtomName( name, elem, buf )
    char *name, *elem;  char buf[4];
{
    char n[5];
    register int i, nlen, elen;

    nlen = 0;
    for( i=0; name[i] && (nlen<4); i++ )
        n[nlen++] = ToUpper(name[i]);
    for( i=nlen; i<4; i++ )
        n[i] = ' ';

    elen = 0;
    for( i=0; elem[i]; i++ )
        if( elem[i] != ' ' )
            elen++;
    if( !elen || (elen>2) )
        elen = 1;

    if( (elen==2) || (nlen==4) )
    {   for( i=0; i<4; i++ )
            buf[i] = n[i];
    } else
    {   buf[0] = ' ';
        buf[1] = n[0];
        buf[2] = n[1];
        buf[3] = n[2];
    }
}


/* Map an mmCIF label_comp_id/auth_comp_id into the 3 character
 * residue field RasMol's residue table expects, including the legacy
 * right justified single letter nucleotide codes ("  A", "  U", ...)
 * and the modern two letter DNA codes (DA/DC/DG/DT).
 */
static void NormalizeResidueName( compId, out3 )
    char *compId;  char out3[4];
{
    char u[8];
    register int i, len;

    len = 0;
    for( i=0; compId[i] && (len<7); i++ )
        u[len++] = ToUpper(compId[i]);
    u[len] = '\0';

    if( len >= 3 )
    {   out3[0] = u[0];  out3[1] = u[1];  out3[2] = u[2];
        out3[3] = '\0';
        return;
    }

    if( (len==1) && strchr("ACGUI",u[0]) )
    {   out3[0] = ' ';  out3[1] = ' ';  out3[2] = u[0];
        out3[3] = '\0';
        return;
    }

    if( (len==2) && (u[0]=='D') && strchr("ACGT",u[1]) )
    {   out3[0] = ' ';  out3[1] = ' ';  out3[2] = u[1];
        out3[3] = '\0';
        return;
    }

    out3[0] = (len>0)? u[0] : ' ';
    out3[1] = (len>1)? u[1] : ' ';
    out3[2] = ' ';
    out3[3] = '\0';
}


/* Process one row of the _atom_site loop, given the token values for
 * this row and the column index of each recognised item (-1 if the
 * file did not provide that item). Mirrors ProcessPDBAtom() closely,
 * including shapely backbone tracing for cartoons/ribbons.
 */
static void ProcessCIFAtom( rowVal, col )
    char rowVal[][CIFTokMax];  int col[];
{
    register Bond __far *bptr;
    register Atom __far *ptr;
    register Long dx, dy, dz;
    register int heta, serno;
    char *compVal, *asymVal, *seqVal, *atomVal, *insVal, *altVal, *typeVal;
    char name3[4], name4[4];
    int chainId, insCode;

    heta = (col[ACOL_GROUP]>=0) &&
           !CIFEqualN(rowVal[col[ACOL_GROUP]],"ATOM",4);

    atomVal = (col[ACOL_ATOM]>=0)?     rowVal[col[ACOL_ATOM]] :
              (col[ACOL_AUTHATOM]>=0)? rowVal[col[ACOL_AUTHATOM]] : "";

    compVal = (col[ACOL_AUTHCOMP]>=0)? rowVal[col[ACOL_AUTHCOMP]] :
              (col[ACOL_COMP]>=0)?     rowVal[col[ACOL_COMP]] : "UNK";

    asymVal = (col[ACOL_AUTHASYM]>=0)? rowVal[col[ACOL_AUTHASYM]] :
              (col[ACOL_ASYM]>=0)?     rowVal[col[ACOL_ASYM]] : " ";

    seqVal  = (col[ACOL_AUTHSEQ]>=0)?  rowVal[col[ACOL_AUTHSEQ]] :
              (col[ACOL_SEQ]>=0)?      rowVal[col[ACOL_SEQ]] : "0";

    insVal  = (col[ACOL_INS]>=0)?  rowVal[col[ACOL_INS]]  : ".";
    altVal  = (col[ACOL_ALT]>=0)?  rowVal[col[ACOL_ALT]]  : ".";
    typeVal = (col[ACOL_TYPE]>=0)? rowVal[col[ACOL_TYPE]] : "";

    chainId = (*asymVal && (*asymVal!='.') && (*asymVal!='?'))?
              ToUpper(*asymVal) : ' ';
    insCode = (*insVal && (*insVal!='.') && (*insVal!='?'))?
              *insVal : ' ';
    serno = atoi(seqVal);

    if( !CurGroup || (CurGroup->serno!=serno) ||
        (CurChain->ident!=chainId) || (CIFInsert!=insCode) )
    {   CIFInsert = (char)insCode;
        if( !CurChain || (CurChain->ident!=chainId) )
            CreateChain( chainId );
        CreateGroup( GroupPool );

        NormalizeResidueName( compVal, name3 );
        CurGroup->refno = FindResNo( name3 );
        CurGroup->serno = serno;
        ProcessGroup( heta );
    }

    ptr = CreateAtom();

    FormatCIFAtomName( atomVal, typeVal, name4 );
    ptr->refno = ComplexAtomType( name4 );
    ptr->serno = (short)((col[ACOL_ID]>=0)? atoi(rowVal[col[ACOL_ID]]) : 0);
    ptr->temp  = (short)((col[ACOL_TEMP]>=0)?
                          (int)(100.0*atof(rowVal[col[ACOL_TEMP]])) : 0);
    ptr->altl  = (*altVal && (*altVal!='.') && (*altVal!='?'))?
                 ToUpper(*altVal) : ' ';

    ptr->xorg =  (Long)(250.0*atof(rowVal[col[ACOL_X]]));
    ptr->yorg =  (Long)(250.0*atof(rowVal[col[ACOL_Y]]));
    ptr->zorg = -(Long)(250.0*atof(rowVal[col[ACOL_Z]]));

    if( heta ) ptr->flag |= HeteroFlag;
    ProcessAtom( ptr );

    /* Create biopolymer Backbone */
    if( IsAlphaCarbon(ptr->refno) && IsProtein(CurGroup->refno) )
    {   if( ConnectAtom )
        {   dx = ConnectAtom->xorg - ptr->xorg;
            dy = ConnectAtom->yorg - ptr->yorg;
            dz = ConnectAtom->zorg - ptr->zorg;

            /* Break backbone if CA-CA > 7.00A */
            if( dx*dx+dy*dy+dz*dz < (Long)1750*1750 )
            {   bptr = ProcessBond(ptr,ConnectAtom,NormBondFlag);
                bptr->bnext = CurChain->blist;
                CurChain->blist = bptr;
            } else ptr->flag |= BreakFlag;
        }
        ConnectAtom = ptr;
    } else if( IsSugarPhosphate(ptr->refno) && IsNucleo(CurGroup->refno) )
    {   if( ConnectAtom )
        {   bptr = ProcessBond(ConnectAtom,ptr,NormBondFlag);
            bptr->bnext = CurChain->blist;
            CurChain->blist = bptr;
        }
        ConnectAtom = ptr;
    }
}


/* Read and discard, or (for "_atom_site") structurally parse, a single
 * "loop_" construct. Returns False only on a premature EOF.
 */
static int LoadCIFAtomSiteAndSkipLoop()
{
    static char rowVal[CIFMaxCols][CIFTokMax];

    int col[ACOL_MAX];
    int tagCount, atomSite, i;
    int firstModel, haveModel;

    for( i=0; i<ACOL_MAX; i++ )
        col[i] = -1;
    tagCount = 0;
    atomSite = False;

    /* Read the loop's column tags */
    while( CIFNextToken() )
    {   if( CIFTokBuf[0] != '_' )
        {   CIFPushToken();
            break;
        }

        if( CIFEqualN(CIFTokBuf,"_atom_site.",11) )
            atomSite = True;

        if( tagCount < CIFMaxCols )
        {   char *item = CIFItemName(CIFTokBuf);
            for( i=0; i<ACOL_MAX; i++ )
                if( CIFEqual(item,CIFAtomSiteItem[i]) )
                {   col[i] = tagCount;
                    break;
                }
        }
        tagCount++;
    }

    if( !tagCount )
        return True;

    if( !atomSite ||
        (col[ACOL_X]<0) || (col[ACOL_Y]<0) || (col[ACOL_Z]<0) ||
        ((col[ACOL_ATOM]<0) && (col[ACOL_AUTHATOM]<0)) ||
        ((col[ACOL_COMP]<0) && (col[ACOL_AUTHCOMP]<0)) )
    {   /* Not a loop we understand: skip over its data rows */
        for(;;)
        {   if( !CIFNextToken() ) return False;
            if( CIFIsReserved(CIFTokBuf) )
            {   CIFPushToken();
                return True;
            }
            for( i=1; i<tagCount; i++ )
                if( !CIFNextToken() ) return False;
        }
    }

    firstModel = 0;
    haveModel = False;
    for(;;)
    {   if( !CIFNextToken() )
            return True;
        if( CIFIsReserved(CIFTokBuf) )
        {   CIFPushToken();
            return True;
        }

        if( 0 < CIFMaxCols )
            CIFStoreTok(rowVal[0]);

        for( i=1; i<tagCount; i++ )
        {   if( !CIFNextToken() )
                return True;   /* Truncated final row! */
            if( i < CIFMaxCols )
                CIFStoreTok(rowVal[i]);
        }

        if( col[ACOL_MODEL] >= 0 )
        {   int model = atoi(rowVal[col[ACOL_MODEL]]);
            if( !haveModel )
            {   firstModel = model;
                haveModel = True;
            } else if( model != firstModel )
                continue;   /* Only the first NMR/multi-state model! */
        }

        ProcessCIFAtom( rowVal, col );
    }
}


/* Process a single scalar "_category.item value" pair, extracting the
 * handful of summary fields the PDB loader fills in from its HEADER,
 * COMPND and CRYST1 records.
 */
static void ProcessCIFScalarItem( tag )
    char *tag;
{
    char tagbuf[64];
    char *item;
    register int i;
    register char *src, *dst;

    for( i=0; tag[i] && (i<63); i++ )
        tagbuf[i] = tag[i];
    tagbuf[i] = '\0';

    if( !CIFNextToken() )
        return;
    item = CIFItemName(tagbuf);

    if( CIFEqualN(tagbuf,"_entry.",7) && CIFEqual(item,"id") )
    {   for( i=0; CIFTokBuf[i] && (i<5); i++ )
            InfoIdentCode[i] = CIFTokBuf[i];
        InfoIdentCode[i] = '\0';

    } else if( CIFEqualN(tagbuf,"_struct.",8) && CIFEqual(item,"title") )
    {   if( !InfoMoleculeName[0] )
        {   for( i=0; CIFTokBuf[i] && (i<79); i++ )
                InfoMoleculeName[i] = CIFTokBuf[i];
            InfoMoleculeName[i] = '\0';
        }

    } else if( CIFEqualN(tagbuf,"_cell.",6) )
    {   if( CIFEqual(item,"length_a") )
        {   InfoCellA = atof(CIFTokBuf);
        } else if( CIFEqual(item,"length_b") )
        {   InfoCellB = atof(CIFTokBuf);
        } else if( CIFEqual(item,"length_c") )
        {   InfoCellC = atof(CIFTokBuf);
        } else if( CIFEqual(item,"angle_alpha") )
        {   InfoCellAlpha = Deg2Rad*atof(CIFTokBuf);
        } else if( CIFEqual(item,"angle_beta") )
        {   InfoCellBeta = Deg2Rad*atof(CIFTokBuf);
        } else if( CIFEqual(item,"angle_gamma") )
            InfoCellGamma = Deg2Rad*atof(CIFTokBuf);

    } else if( CIFEqualN(tagbuf,"_symmetry.",10) &&
               CIFEqual(item,"space_group_name_H-M") )
    {   dst = InfoSpaceGroup;
        for( src=CIFTokBuf; *src && (dst<InfoSpaceGroup+10); src++ )
            if( *src != ' ' )
                *dst++ = *src;
        *dst = '\0';
    }
}


int LoadCIFMolecule( fp )
    FILE *fp;
{
    DataFile = fp;
    CIFAtBOL = True;
    CIFTokPushed = False;
    NMRModel = 0;

    while( CIFNextToken() )
    {   if( CIFEqual(CIFTokBuf,"loop_") )
        {   if( !LoadCIFAtomSiteAndSkipLoop() )
                break;
        } else if( CIFTokBuf[0] == '_' )
            ProcessCIFScalarItem(CIFTokBuf);
        /* "data_", "save_", "stop_", "global_" and anything else
         * outside a loop or scalar item is simply ignored.
         */
    }

    if( Database )
        strcpy(InfoFileName,DataFileName);
    return( Database? True : False );
}


/*=================================*/
/* Molecule File Format Generation */
/*=================================*/

int SavePDBMolecule( filename )
    char *filename;
{
    register double x, y, z;
    register Group __far *prev;
    register Chain __far *chain;
    register Group __far *group;
    register Atom __far *aptr;
    register char *ptr;
    register int count;
    register char ch;
    register int i;
 
    if( !Database )
        return( False );
 
    DataFile = fopen( filename, "w" );
    if( !DataFile )
    {   if( CommandActive )
            WriteChar('\n');
        WriteString("Error: Unable to create file!\n\n");
        CommandActive=False;
        return( False );
    }
 
    if( *InfoClassification || *InfoIdentCode )
    {   fputs("HEADER    ",DataFile);
 
        ptr = InfoClassification;
        for( i=11; i<=50; i++ )
            putc( (*ptr ? *ptr++ : ' '), DataFile );
        fprintf(DataFile,"13-JUL-93   %.4s\n",InfoIdentCode);
    }
 
    if( *InfoMoleculeName )
        fprintf(DataFile,"COMPND    %.60s\n",InfoMoleculeName);
 
    prev = (void __far*)0;
    ch = ' ';

    count = 1;
    ForEachAtom
        if( aptr->flag&SelectFlag )
        {   if( prev && (chain->ident!=ch) )
                fprintf( DataFile, "TER   %5d      %.3s %c%4d \n",
                         count++, Residue[prev->refno], ch, prev->serno);
 
            if( aptr->flag&HeteroFlag )
            {      fputs("HETATM",DataFile);
            } else fputs("ATOM  ",DataFile);
            fprintf( DataFile, "%5d %.4s %.3s %c%4d    ",
                     count++, ElemDesc[aptr->refno], Residue[group->refno],
                     chain->ident, group->serno );
 
            x = (double)aptr->xorg/250.0;
            y = (double)aptr->yorg/250.0;
            z = (double)aptr->zorg/250.0;
 
#ifdef INVERT
            fprintf(DataFile,"%8.3f%8.3f%8.3f",x,-y,-z);
#else
            fprintf(DataFile,"%8.3f%8.3f%8.3f",x,y,-z);
#endif
            fprintf(DataFile,"  1.00%6.2f\n",aptr->temp/100.0);
 
            ch = chain->ident;
            prev = group;
        }
 
    if( prev )
        fprintf( DataFile, "TER   %5d      %.3s %c%4d \n",
                 count, Residue[prev->refno], ch, prev->serno);
 
    fputs("END   \n",DataFile);
    fclose( DataFile );
#ifdef APPLEMAC
    SetFileInfo(filename,'RSML','TEXT',131);
#endif
    return( True );
}
 
 
int SaveMDLMolecule( filename )
    char *filename;
{
    register Chain __far *chain;
    register Group __far *group;
    register Atom __far *aptr;
    register Bond __far *bptr;
    register int atoms,bonds;
    register double x,y,z;
    register int ch,temp;
    register int atomno;

#ifndef APPLEMAC
    register struct tm *ptr;
    static time_t curtime;
#endif
 
    if( !Database )
        return( False );

    DataFile = fopen( filename, "w" );
    if( !DataFile )
    {   if( CommandActive )
            WriteChar('\n');
        WriteString("Error: Unable to create file!\n\n");
        CommandActive=False;
        return( False );
    }

    /* Write Mol file header */
    fprintf(DataFile,"%.80s\n  RasMol  ",InfoMoleculeName);

#ifndef APPLEMAC
    curtime = time((time_t*)0);
    ptr = localtime(&curtime);
    fprintf(DataFile,"%02d%02d%02d%02d%02d",ptr->tm_mon+1,ptr->tm_mday,
                             ptr->tm_year%100,ptr->tm_hour,ptr->tm_min);
#else
    fputs("0101951200",DataFile);
#endif

    atoms = 0;
    ForEachAtom
        if( aptr->flag & SelectFlag )
            atoms++;

    bonds = 0;
    ForEachBond
        if( bptr->srcatom->flag & bptr->dstatom->flag & SelectFlag )
             bonds++;

    fprintf(DataFile,"%cD\n\n", (MinZ||MaxZ)?'3':'2' );
    fprintf(DataFile,"%3d%3d",atoms,bonds);
    fputs("  0        0              1 V2000\n",DataFile);

    atomno = 1;
    ForEachAtom
        if( aptr->flag & SelectFlag )
        {   x = (double)aptr->xorg/250.0;
            y = (double)aptr->yorg/250.0;
            z = (double)aptr->zorg/250.0;
#ifdef INVERT
            fprintf(DataFile,"%10.4lf%10.4lf%10.4lf ",x,-y,-z);
#else
            fprintf(DataFile,"%10.4lf%10.4lf%10.4lf ",x,y,-z);
#endif
            fputc(Element[aptr->elemno].symbol[0],DataFile);
            fputc(Element[aptr->elemno].symbol[1],DataFile);
            fputs("  0  ",DataFile);

            ch = '0';
            /* Test for charges or b-factors */
            if( (MinMainTemp<0) || (MinHetaTemp<0) )
            {   temp = aptr->temp;
                if( temp > 50 )
                {   if( temp > 250 )
                    {   ch = '1';
                    } else if( temp > 150 )
                    {   ch = '2';
                    } else ch = '3';
                } else if( temp < -50 )
                {   if( temp < -250 )
                    {   ch = '7';
                    } else if( temp < -150 )
                    {   ch = '6';
                    } else ch = '5';
                } else ch = '0';
            } 
            fputc(ch,DataFile);

            fputs("  0  0  0  0  0  0  0  0  0  0\n",DataFile);
            aptr->mbox = atomno++;
        }

    ForEachBond
        if( bptr->srcatom->flag & bptr->dstatom->flag & SelectFlag )
        {   fprintf(DataFile,"%3d%3d  ",bptr->srcatom->mbox,
                                        bptr->dstatom->mbox);
            if( bptr->flag & AromBondFlag )
            {      fputc('4',DataFile);
            } else if( bptr->flag & TripBondFlag )
            {      fputc('3',DataFile);
            } else if( bptr->flag & DoubBondFlag )
            {      fputc('2',DataFile);
            } else fputc('1',DataFile);
            fputs("  0\n",DataFile);
        }

    fputs("M  END\n",DataFile);
    fclose( DataFile );
#ifdef APPLEMAC
    SetFileInfo(filename,'RSML','mMOL',131);
#endif
    return( True );
}


int SaveAlchemyMolecule( filename )
    char *filename;
{
    register Real x, y, z;
    register float xpos, ypos, zpos;
    register Chain __far *chain;
    register Group __far *group;
    register Atom __far *aptr;
    register Bond __far *bptr;
    register char *ptr;
    register int atomno;
    register int bondno;
    register int num;
 
    if( !Database )
        return( False );
 
    DataFile = fopen( filename, "w" );
    if( !DataFile )
    {   if( CommandActive )
            WriteChar('\n');
        WriteString("Error: Unable to create file!\n\n");
        CommandActive=False;
        return( False );
    }
 
    atomno = 0;
    ForEachAtom
        if( aptr->flag & SelectFlag )
            aptr->mbox = 0;
 
    ForEachBond
        if( ((bptr->srcatom->flag&bptr->dstatom->flag)&SelectFlag) &&
           !((bptr->srcatom->flag|bptr->dstatom->flag)&HydrogenFlag) )
        {   if( bptr->flag&AromBondFlag )
            {   bptr->srcatom->mbox = -1;
                bptr->dstatom->mbox = -1;
            } else
            {   num = (bptr->flag&DoubBondFlag)? 2 : 1;
                if( bptr->srcatom->mbox>0 )
                    bptr->srcatom->mbox += num;
                if( bptr->dstatom->mbox>0 )
                    bptr->dstatom->mbox += num;
            }
        }
 
    fprintf(DataFile,"%5d ATOMS, ",MainAtomCount+HetaAtomCount);
    fprintf(DataFile,"%5d BONDS, ",InfoBondCount);
    fprintf(DataFile,"    0 CHARGES, %s\n", InfoMoleculeName );
 
    atomno = 1;
    ForEachAtom
        if( aptr->flag & SelectFlag )
        {   aptr->mbox = atomno;
            fprintf(DataFile,"%5d ",atomno++);
 
            switch( aptr->elemno )
            {   case( 6 ):  if( aptr->mbox == -1 )
                            {   ptr = "CAR ";
                            } else if( aptr->mbox == 1 )
                            {   ptr = "C3  ";
                            } else ptr = "C2  ";
                            fputs( ptr, DataFile );
                            break;
 
                case( 7 ):  if( aptr->mbox == -1 )
                            {   ptr = "NAR ";
                            } else ptr = "N2  ";
                            fputs( ptr, DataFile );
                            break;
 
                case( 8 ):  if( aptr->mbox == 2 )
                            {   ptr = "O2  ";
                            } else ptr = "O3  ";
                            fputs( ptr, DataFile );
                            break;
 
                case( 1 ):  fputs( "H   ", DataFile );  break;
 
                default:    ptr = ElemDesc[aptr->refno];
                            if( *ptr==' ' )
                            {   fprintf(DataFile,"%.3s ",ptr+1);
                            } else fprintf(DataFile,"%.4s",ptr);
            }
 
            x = aptr->xorg/250.0;
            y = aptr->yorg/250.0;
            z = aptr->zorg/250.0;
 
            /* Apply Current Viewpoint Rotation Matrix */
            xpos = (float)(x*RotX[0] + y*RotX[1] + z*RotX[2]);
            ypos = (float)(x*RotY[0] + y*RotY[1] + z*RotY[2]);
            zpos = (float)(x*RotZ[0] + y*RotZ[1] + z*RotZ[2]);
 
#ifdef INVERT
            fprintf(DataFile,"  %8.4f %8.4f %8.4f",xpos,-ypos,-zpos);
#else
            fprintf(DataFile,"  %8.4f %8.4f %8.4f",xpos,ypos,-zpos);
#endif
            fprintf(DataFile,"    %7.4f\n",aptr->temp/1000.0);
        }
 
    bondno = 1;
    ForEachBond
        if( (bptr->srcatom->flag&bptr->dstatom->flag) & SelectFlag )
        {   fprintf(DataFile,"%5d %5d %5d  ", bondno++,
                        bptr->srcatom->mbox, bptr->dstatom->mbox );
            if( bptr->flag & AromBondFlag )
            {   ptr = "AROMATIC\n";
            } else if( bptr->flag & TripBondFlag )
            {   ptr = "TRIPLE\n";
            } else if( bptr->flag & DoubBondFlag )
            {   ptr = "DOUBLE\n";
            } else ptr = "SINGLE\n";
            fputs( ptr, DataFile );
        }
 
    ForEachAtom
            if( aptr->flag & SelectFlag )
                aptr->mbox = 0;
    fclose( DataFile );
#ifdef APPLEMAC
    SetFileInfo(filename,'RSML','TEXT',131);
#endif
    return( True );
}


int SaveXYZMolecule( filename )
    char *filename;
{
    if( !Database )
        return( False );
    return( True );
}
 
 
int SaveCIFMolecule( filename )
    char *filename;
{
    if( !Database )
        return( False );
    return( True );
}

