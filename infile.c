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


/* True if chain's identifier is exactly the single raw PDB column
 * character ch. Fixed-column PDB records (ATOM/HETATM/HELIX/SHEET/
 * TURN/TER) only ever create or reference single character chain
 * identifiers, so this is the only comparison those code paths need,
 * even though Chain.ident is a general string.
 */
static int ChainIdentIs( chain, ch )
    Chain __far *chain;  int ch;
{
    return( (chain->ident[0]==ch) && !chain->ident[1] );
}


static void UpdateFeature( ptr, mask )
    FeatEntry __far *ptr;  int mask;
{
    register Chain __far *chain;
    register Group __far *group;

    for( chain=Database->clist; chain; chain=chain->cnext )
        if( ChainIdentIs(chain,ptr->chain) )
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
    char buf[2];

    PDBInsert = Record[26];
    if( !CurChain || !ChainIdentIs(CurChain,Record[21]) )
    {   buf[0] = Record[21];  buf[1] = '\0';
        CreateChain( buf );
    }
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
        || !ChainIdentIs(CurChain,Record[21])
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
    char name3[4], name4[4], chainBuf[MAXCHAINID];
    int insCode, k;

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

    if( *asymVal && (*asymVal!='.') && (*asymVal!='?') )
    {   for( k=0; asymVal[k] && (k<MAXCHAINID-1); k++ )
            chainBuf[k] = ToUpper(asymVal[k]);
        chainBuf[k] = '\0';
    } else
    {   chainBuf[0] = ' ';  chainBuf[1] = '\0'; }

    insCode = (*insVal && (*insVal!='.') && (*insVal!='?'))?
              *insVal : ' ';
    serno = atoi(seqVal);

    if( !CurGroup || (CurGroup->serno!=serno) ||
        strcmp(CurChain->ident,chainBuf) || (CIFInsert!=insCode) )
    {   CIFInsert = (char)insCode;
        if( !CurChain || strcmp(CurChain->ident,chainBuf) )
            CreateChain( chainBuf );
        CreateGroup( GroupPool );

        NormalizeResidueName( compVal, name3 );
        CurGroup->refno = FindResNo( name3 );
        CurGroup->serno = serno;
        ProcessGroup( heta );
    }

    ptr = CreateAtom();

    FormatCIFAtomName( atomVal, typeVal, name4 );
    ptr->refno = ComplexAtomType( name4 );
    ptr->serno = (col[ACOL_ID]>=0)? atoi(rowVal[col[ACOL_ID]]) : 0;
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
    {   /* Classic Brookhaven codes are 4 characters ("4HHB"); newer
         * "extended" PDB identifiers are 12 ("pdb_00004hhb") now that
         * 4-letter codes are running out. Take whichever _entry.id
         * actually is, up to MAXIDENTCODE-1, rather than assuming 4.
         */
        for( i=0; CIFTokBuf[i] && (i<MAXIDENTCODE-1); i++ )
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
        {   /* PDB's chain column is a single character; this is the
             * best any PDB writer can do for a chain whose real
             * (mmCIF) identifier is longer than one character.
             */
            if( prev && (chain->ident[0]!=ch) )
                fprintf( DataFile, "TER   %5d      %.3s %c%4d \n",
                         count++, Residue[prev->refno], ch, prev->serno);

            if( aptr->flag&HeteroFlag )
            {      fputs("HETATM",DataFile);
            } else fputs("ATOM  ",DataFile);
            fprintf( DataFile, "%5d %.4s %.3s %c%4d    ",
                     count++, ElemDesc[aptr->refno], Residue[group->refno],
                     chain->ident[0], group->serno );
 
            x = (double)aptr->xorg/250.0;
            y = (double)aptr->yorg/250.0;
            z = (double)aptr->zorg/250.0;
 
#ifdef INVERT
            fprintf(DataFile,"%8.3f%8.3f%8.3f",x,-y,-z);
#else
            fprintf(DataFile,"%8.3f%8.3f%8.3f",x,y,-z);
#endif
            fprintf(DataFile,"  1.00%6.2f\n",aptr->temp/100.0);
 
            ch = chain->ident[0];
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
 
 
