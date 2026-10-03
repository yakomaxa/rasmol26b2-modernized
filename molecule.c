/* molecule.c
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

#define MOLECULE
#include "molecule.h"
#include "command.h"
#include "abstree.h"
#include "transfor.h"
#include "render.h"

#define HBondPool   32
#define BondPool    32
#define AtomPool    32

#define NoLadder     0x00
#define ParaLadder   0x01
#define AntiLadder   0x02

#define Cos70Deg     0.34202014332567

#define MaxHBondDist   ((Long)300*300)
#define MaxBondDist    ((Long)475*475)
#define MinBondDist    ((Long)100*100)
#define AbsMaxBondDist 600

#ifdef APPLEMAC
#define AllocSize   256
typedef struct _AllocRef {
	struct _AllocRef *next;
	void *data[AllocSize];
	int count;
	} AllocRef;
static AllocRef *AllocList;  
#endif


typedef struct {
          char name[4];
          int code;
      } SynonymTable;

#define RESSYNMAX 16
static SynonymTable ResSynonym[RESSYNMAX] = {
    { "ADE", 24 },  /*   A : Adenosine   */
    { "CPR", 11 },  /* PRO : Cis-proline */
    { "CSH", 17 },  /* CYS : Cystine     */
    { "CSM", 17 },  /* CYS : Cystine     */
    { "CYH", 17 },  /* CYS : Cystine     */
    { "CYT", 25 },  /*   C : Cytosine    */
    { "D2O", 47 },  /* DOD : Heavy Water */
    { "GUA", 26 },  /*   G : Guanosine   */
    { "H2O", 46 },  /* HOH : Solvent     */
    { "SOL", 46 },  /* HOH : Solvent     */
    { "SUL", 48 },  /* SO4 : Sulphate    */
    { "THY", 27 },  /*   T : Thymidine   */
    { "TIP", 46 },  /* HOH : Water       */
    { "TRY", 20 },  /* TRP : Tryptophan  */
    { "URI", 28 },  /*   U : Uridine     */
    { "WAT", 46 }   /* HOH : Water       */
        };


static Molecule __far *FreeMolecule;
static HBond __far *FreeHBond;
static Chain __far *FreeChain;
static Group __far *FreeGroup;
static Atom __far *FreeAtom;
static Bond __far *FreeBond;

static IntCoord __far *IntPrev;
static HBond __far * __far *CurHBond;
static int MemSize;


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


static void FatalDataError(ptr)
    char *ptr;
{
    char buffer[80];

    sprintf(buffer,"Database Error: %s!",ptr);
    RasMolFatalExit(buffer);
}



void DescribeMolecule()
{
    char buffer[40];

    if( CommandActive )
	WriteChar('\n');
    CommandActive=False;

    if( *InfoMoleculeName )
    {   WriteString("Molecule name ....... ");
	WriteString(InfoMoleculeName);
	WriteChar('\n');
    }

    if( *InfoClassification )
    {   WriteString("Classification ...... ");
	WriteString(InfoClassification);
	WriteChar('\n');
    }

    if( Database && (MainGroupCount>1) )
    {   WriteString("Secondary Structure . ");
	if( InfoStrucSrc==SourceNone )
	{   WriteString("No Assignment\n");
	} else if( InfoStrucSrc==SourcePDB )
	{   WriteString("PDB Data Records\n");
	} else WriteString("Calculated\n");
    }


    if( *InfoIdentCode )
    {   WriteString("Brookhaven Code ..... ");
	WriteString(InfoIdentCode);
	WriteChar('\n');
    }

    if( InfoChainCount>1 )
    {   sprintf(buffer,"Number of Chains .... %d\n",InfoChainCount);
	WriteString(buffer);
    }

    sprintf(buffer,"Number of Groups .... %d",MainGroupCount);
    WriteString(buffer);
    if( HetaAtomCount )
    {   sprintf(buffer," (%d)\n",HetaGroupCount);
	WriteString(buffer);
    } else WriteChar('\n');

    sprintf(buffer,"Number of Atoms ..... %ld",(long)MainAtomCount);
    WriteString(buffer);
    if( HetaAtomCount )
    {   sprintf(buffer," (%d)\n",HetaAtomCount);
	WriteString(buffer);
    } else WriteChar('\n');

    if( InfoBondCount )
    {   sprintf(buffer,"Number of Bonds ..... %ld\n",(long)InfoBondCount);
	WriteString(buffer);
    }

    if( InfoSSBondCount != -1 )
    {   WriteString("Number of Bridges ... ");
	sprintf(buffer,"%d\n\n",InfoSSBondCount);
	WriteString(buffer);
    }

    if( InfoHBondCount != -1 )
    {   WriteString("Number of H-Bonds ... ");
	sprintf(buffer,"%d\n",InfoHBondCount);
	WriteString(buffer);
    }

    if( InfoHelixCount != -1 )
    {   WriteString("Number of Helices ... ");
	sprintf(buffer,"%d\n",InfoHelixCount);
	WriteString(buffer);

	WriteString("Number of Strands ... ");
	sprintf(buffer,"%d\n",InfoLadderCount);
	WriteString(buffer);

	WriteString("Number of Turns ..... ");
	sprintf(buffer,"%d\n",InfoTurnCount);
	WriteString(buffer);
    }
}


#ifdef APPLEMAC
/* Avoid System Memory Leaks! */
static void RegisterAlloc( data )
    void *data;
{
    register AllocRef *ptr;
    
    if( !AllocList || (AllocList->count==AllocSize) )
    {   ptr = (AllocRef *)_fmalloc( sizeof(AllocRef) );
	if( !ptr ) FatalDataError("Memory allocation failed");
	
	ptr->next = AllocList;
	ptr->data[0] = data;
	ptr->count = 1;
	AllocList = ptr;
    } else AllocList->data[AllocList->count++] = data;
}
#else
#define RegisterAlloc(x)
#endif


/*==================================*/
/* Group & Chain Handling Functions */
/*==================================*/

void CreateChain( ident )
    char *ident;
{
    register Chain __far *prev;
    register int i;

    if( !CurMolecule )
    {   if( !(CurMolecule = FreeMolecule) )
	{   MemSize += sizeof(Molecule);
	    CurMolecule = (Molecule __far *)_fmalloc(sizeof(Molecule));
	    if( !CurMolecule ) FatalDataError("Memory allocation failed");
	    RegisterAlloc( CurMolecule );
	} else FreeMolecule = (void __far*)0;

	CurChain = (void __far*)0;
	CurMolecule->slist = (void __far*)0;
	CurMolecule->hlist = (void __far*)0;
	CurMolecule->blist = (void __far*)0;
	CurMolecule->clist = (void __far*)0;
	Database = CurMolecule;
    }

    /* Handle chain breaks! */
    if( !(prev=CurChain) )
	if( (prev=CurMolecule->clist) )
	    while( prev->cnext )
		prev = prev->cnext;

    if( !(CurChain = FreeChain) )
    {   MemSize += sizeof(Chain);
	CurChain = (Chain __far *)_fmalloc(sizeof(Chain));
	if( !CurChain ) FatalDataError("Memory allocation failed");
	RegisterAlloc( CurChain );
    } else FreeChain = FreeChain->cnext;

    if( prev )
    {   prev->cnext = CurChain;
    } else CurMolecule->clist = CurChain;
    CurChain->cnext = (void __far*)0;
     
    i = 0;
    while( ident[i] && (i<MAXCHAINID-1) )
    {   CurChain->ident[i] = ident[i];  i++; }
    CurChain->ident[i] = '\0';
    CurChain->identval = ChainIdentValue( CurChain->ident );

    CurChain->model = NMRModel;
    CurChain->glist = (void __far*)0;
    CurChain->blist = (void __far*)0;
    ConnectAtom = (void __far*)0;
    CurGroup = (void __far*)0;
    InfoChainCount++;
}


void CreateGroup( pool )
    int pool;
{
    register Group __far *ptr;
    register int i;

    if( !(ptr = FreeGroup) )
    {   MemSize += pool*sizeof(Group);
	ptr = (Group __far *)_fmalloc( pool*sizeof(Group) );
	if( !ptr ) FatalDataError("Memory allocation failed");
	RegisterAlloc( ptr );
	for( i=1; i<pool; i++ )
	{   ptr->gnext = FreeGroup;
	    FreeGroup = ptr++;
	} 
    } else FreeGroup = ptr->gnext;
    
    if( CurGroup )
    {   ptr->gnext = CurGroup->gnext;
        CurGroup->gnext = ptr;
    } else 
    {   ptr->gnext = CurChain->glist;
        CurChain->glist = ptr;
    }
    CurGroup = ptr;

    CurAtom = (void __far*)0;
    ptr->alist = (void __far*)0;
    ptr->insert = ' ';
    ptr->struc = 0;
    ptr->flag = 0;
    ptr->col1 = 0;
    ptr->col2 = 0;
}


int FindResNo( ptr )
    char *ptr;
{
    register int hi,lo;
    register int refno;
    register int flag;
    register int mid;

    for( refno=0; refno<ResNo; refno++ )
	if( !strncmp(Residue[refno],ptr,3) )
	    return( refno );

    lo = 0;
    hi = RESSYNMAX;
    while( lo < hi )
    {   mid = (hi+lo)>>1;
        flag = strncmp(ResSynonym[mid].name,ptr,3);
        if( !flag ) return( ResSynonym[mid].code );

        /* Binary Search */
        if( flag<0 )
        {   lo = mid+1;
        } else hi = mid;
    }

    if( ResNo++ == MAXRES )
	FatalDataError("Too many new residues");
    Residue[refno][0] = *ptr++;
    Residue[refno][1] = *ptr++;
    Residue[refno][2] = *ptr;
    return( refno );
}


void ProcessGroup( heta )
    int heta;
{
    register int serno;

    serno = CurGroup->serno;
    if( IsSolvent(CurGroup->refno) )
        heta = True;

    if( heta )
    {   HetaGroupCount++;
        if( HMinMaxFlag )
        {   if( serno > MaxHetaRes )
            {   MaxHetaRes = serno;
            } else if( serno < MinHetaRes )
                MinHetaRes = serno;
        } else MinHetaRes = MaxHetaRes = serno;
    } else
    {   MainGroupCount++;
        if( MMinMaxFlag )
        {   if( serno > MaxMainRes )
            {   MaxMainRes = serno;
            } else if( serno < MinMainRes )
                MinMainRes = serno;
        } else MinMainRes = MaxMainRes = serno;
    }
}


void CreateMolGroup()
{
    strcpy(InfoFileName,DataFileName);

    CreateChain( " " );
    CreateGroup( 1 );

    CurGroup->refno = FindResNo( "MOL" );
    CurGroup->serno = 1;
	
    MinMainRes = MaxMainRes = 1;
    MinHetaRes = MaxHetaRes = 0;
    MainGroupCount = 1;
}


/*=========================*/
/* Atom Handling Functions */
/*=========================*/


Atom __far *CreateAtom()
{
    register Atom __far *ptr;
    register int i;

    if( !(ptr = FreeAtom) )
    {   MemSize += AtomPool*sizeof(Atom);
	ptr = (Atom __far *)_fmalloc( AtomPool*sizeof(Atom) );
	if( !ptr ) FatalDataError("Memory allocation failed");
	RegisterAlloc( ptr );
	for( i=1; i<AtomPool; i++ )
	{   ptr->anext = FreeAtom;
	    FreeAtom = ptr++;
	} 
    } else FreeAtom = ptr->anext;

    if( CurAtom )
    {   ptr->anext = CurAtom->anext;
        CurAtom->anext = ptr;
    } else 
    {   ptr->anext = CurGroup->alist;
        CurGroup->alist = ptr;
    }
    CurAtom = ptr;

    SelectCount++;
    ptr->flag = SelectFlag | NonBondFlag;
    ptr->label = (void*)0;
    ptr->radius = 375;
    ptr->altl = ' ';
    ptr->mbox = 0;
    ptr->col = 0;

    return( ptr );
}


void ProcessAtom( ptr )
    Atom __far *ptr;
{
    ptr->elemno = GetElemNumber(CurGroup,ptr);
    if( ptr->elemno == 1 )
    {   ptr->flag |= HydrogenFlag;
	HasHydrogen = True;
    }

    if( !IsSolvent(CurGroup->refno) )
    {   if( !(ptr->flag&(HydrogenFlag|HeteroFlag)) )
	    ptr->flag |= NormAtomFlag;
    } else ptr->flag |= HeteroFlag;

#ifdef INVERT
    ptr->yorg = -ptr->yorg;
#endif

    if( HMinMaxFlag || MMinMaxFlag )
    {   if( ptr->xorg < MinX ) 
	{   MinX = ptr->xorg;
	} else if( ptr->xorg > MaxX ) 
	    MaxX = ptr->xorg;

	if( ptr->yorg < MinY ) 
	{   MinY = ptr->yorg;
	} else if( ptr->yorg > MaxY ) 
	    MaxY = ptr->yorg;

	if( ptr->zorg < MinZ ) 
	{   MinZ = ptr->zorg;
	} else if( ptr->zorg > MaxZ ) 
	    MaxZ = ptr->zorg;
    } else 
    {   MinX = MaxX = ptr->xorg;
	MinY = MaxY = ptr->yorg;
	MinZ = MaxZ = ptr->zorg;
    }
	    
    if( ptr->flag & HeteroFlag )
    {   if( HMinMaxFlag )
	{   if( ptr->temp < MinHetaTemp ) 
	    {   MinHetaTemp = ptr->temp;
	    } else if( ptr->temp > MaxHetaTemp ) 
		MaxHetaTemp = ptr->temp;
	} else MinHetaTemp = MaxHetaTemp = ptr->temp;
	HMinMaxFlag = True;
	HetaAtomCount++;
    } else
    {   if( MMinMaxFlag )
	{   if( ptr->temp < MinMainTemp ) 
	    {   MinMainTemp = ptr->temp;
	    } else if( ptr->temp > MaxMainTemp ) 
		MaxMainTemp = ptr->temp;
	} else MinMainTemp = MaxMainTemp = ptr->temp;
	MMinMaxFlag = True;
	MainAtomCount++;
    }
}


Atom __far *FindGroupAtom( group, n )
    Group __far *group;  Byte n;
{
    register Atom __far *ptr;

    for( ptr=group->alist; ptr; ptr=ptr->anext )
        if( ptr->refno == n ) return( ptr );
    return( (Atom __far*)0 );
}


int NewAtomType( ptr )
    char *ptr;
{
    register int refno;
    register int i;

    for( refno=0; refno<ElemNo; refno++ )
	if( !strncmp(ElemDesc[refno],ptr,4) )
	    return(refno);

    if( ElemNo++ == MAXELEM )
        FatalDataError("Too many new atom types");

    for( i=0; i<4; i++ )
        ElemDesc[refno][i] = ptr[i];
    return( refno );
}


int SimpleAtomType( type )
    char *type;
{
    char name[4];

    name[2] = name[3] = ' ';
    if( type[1] && (type[1]!=' ') )
    {   name[0] = ToUpper(type[0]);
	name[1] = ToUpper(type[1]);
    } else
    {   name[1] = ToUpper(type[0]);
	name[0] = ' ';
    }
    return( NewAtomType(name) );
}


int ComplexAtomType( ptr )
    char *ptr;
{
    static char name[4];
    register int i;

    if( isdigit(ptr[1]) )
    {   name[1] = ToUpper(*ptr);
	name[2] = name[3] = ' ';
	name[0] = ' ';
    } else for( i=0; i<4; i++ )
	name[i] = ToUpper(ptr[i]);


    /* Handle Unconventional Naming */
    if( IsProtein(CurGroup->refno) )
    {   if( name[0]=='H' )
        {   name[0]=' ';
            name[1]='H';
        }
    } else if( IsNucleo(CurGroup->refno) )
    {   if( name[3]=='\'' )
            name[3] = '*';

        if( (name[1]=='O') && (name[2]=='P') )
        {   if( !strncmp(name," OP1",4) ||
                !strncmp(name,"1OP ",4) )
                return( 8 );
            if( !strncmp(name," OP2",4) ||
                !strncmp(name,"2OP ",4) )
                return( 9 );
        }
    }
    return( NewAtomType(name) );
}





/*===============================*/
/* Z-Matrix Conversion Functions */
/*===============================*/

#ifdef FUNCPROTO
static IntCoord __far* GetInternalCoord( int );
#endif

void InitInternalCoords()
{
    IntList = (IntCoord __far*)0;
    IntPrev = (IntCoord __far*)0;
}


IntCoord __far* AllocInternalCoord()
{
    register IntCoord __far *ptr;

    ptr = (IntCoord __far*)_fmalloc(sizeof(IntCoord));
    if( !ptr ) FatalDataError("Memory allocation failed");
    ptr->inext = (IntCoord __far*)0;

    if( IntPrev )
    {   IntPrev->inext = ptr;
    } else IntList = ptr;
    IntPrev = ptr;
    return( ptr );
}


static IntCoord __far* GetInternalCoord( index )
    int index;
{
    register IntCoord __far *ptr;

    ptr = IntList;
    while( (index>1) && ptr->inext )
    {   ptr = ptr->inext;
        index--;
    }
    return( ptr );
}


void FreeInternalCoords()
{
    register IntCoord __far *ptr;

    while( (ptr = IntList) )
    {    IntList = ptr->inext;
         _ffree( ptr );
    }
}


int ConvertInternal2Cartesian()
{
    register IntCoord __far *ptr;
    register IntCoord __far *na;
    register IntCoord __far *nb;
    register IntCoord __far *nc;
    register double cosph,sinph,costh,sinth,coskh,sinkh;
    register double cosa,sina,cosd,sind;
    register double dist,angle,dihed;

    register double xpd,ypd,zpd,xqd,yqd,zqd;
    register double xa,ya,za,xb,yb,zb;
    register double rbc,xyb,yza,temp;
    register double xpa,ypa,zqa;
    register double xd,yd,zd;
    register int flag;


    /* Atom #1 */
    ptr = IntList;
    ptr->dist  = 0.0;
    ptr->angle = 0.0;
    ptr->dihed = 0.0;

    if( !(ptr=ptr->inext) )
        return( True );

    /* Atom #2 */
    ptr->angle = 0.0;
    ptr->dihed = 0.0;

    if( !(ptr=ptr->inext) )
        return( True );

    /* Atom #3 */
    dist = ptr->dist;
    angle = Deg2Rad*ptr->angle;
    cosa = cos(angle);
    sina = sin(angle);
    if( ptr->na == 1 )
    {   na = IntList;
        ptr->dist = na->dist + cosa*dist;
    } else /* ptr->na == 2 */
    {   na = IntList->inext;
        ptr->dist = na->dist - cosa*dist;
    }
    ptr->angle = sina*dist;
    ptr->dihed = 0.0;

    while( ptr=ptr->inext )
    {   dist = ptr->dist;
        angle = Deg2Rad*ptr->angle;
        dihed = Deg2Rad*ptr->dihed;

        /* Optimise this access?? */
        na = GetInternalCoord(ptr->na);
        nb = GetInternalCoord(ptr->nb);
        nc = GetInternalCoord(ptr->nc);

        xb = nb->dist  - na->dist;
        yb = nb->angle - na->angle;
        zb = nb->dihed - na->dihed;

        rbc = xb*xb + yb*yb + zb*zb;
        if( rbc < 0.0001 )
            return( False );
        rbc= 1.0/sqrt(rbc);

        cosa = cos(angle);
        sina = sin(angle);


        if( fabs(cosa) >= 0.999999 )
        {   /* Colinear */
            temp = dist*rbc*cosa;
            ptr->dist  = na->dist  + temp*xb;
            ptr->angle = na->angle + temp*yb;
            ptr->dihed = na->dihed + temp*zb;
        } else
        {   xa = nc->dist  - na->dist;
            ya = nc->angle - na->angle;
            za = nc->dihed - na->dihed;

            sind = -sin(dihed);
            cosd = cos(dihed);

            xd = dist*cosa;
            yd = dist*sina*cosd;
            zd = dist*sina*sind;

            xyb = sqrt(xb*xb + yb*yb);
            if( xyb < 0.1 )
            {   /* Rotate about y-axis! */
                temp = za; za = -xa; xa = temp;
                temp = zb; zb = -xb; xb = temp;
                xyb = sqrt(xb*xb + yb*yb);
                flag = True;
            } else flag = False;

            costh = xb/xyb;
            sinth = yb/xyb;
            xpa = costh*xa + sinth*ya;
            ypa = costh*ya - sinth*xa;

            sinph = zb*rbc;
            cosph = sqrt(1.0 - sinph*sinph);
            zqa = cosph*za  - sinph*xpa;

            yza = sqrt(ypa*ypa + zqa*zqa);

            if( yza > 1.0E-10 )
            {   coskh = ypa/yza;
                sinkh = zqa/yza;

                ypd = coskh*yd - sinkh*zd;
                zpd = coskh*zd + sinkh*yd;
            } else
            {   /* coskh = 1.0; */
                /* sinkh = 0.0; */
                ypd = yd;
                zpd = zd;
            }

            xpd = cosph*xd  - sinph*zpd;
            zqd = cosph*zpd + sinph*xd;
            xqd = costh*xpd - sinth*ypd;
            yqd = costh*ypd + sinth*xpd;

            if( flag )
            {   /* Rotate about y-axis! */
                ptr->dist  = na->dist  - zqd;
                ptr->angle = na->angle + yqd;
                ptr->dihed = na->dihed + xqd;
            } else
            {   ptr->dist  = na->dist  + xqd;
                ptr->angle = na->angle + yqd;
                ptr->dihed = na->dihed + zqd;
            }
        }
    }
    return( True );
}



/*=========================*/
/* Bond Handling Functions */
/*=========================*/


#ifdef FUNCPROTO
Bond __far *ProcessBond( Atom __far*, Atom __far*, int );
static void CreateHydrogenBond( Atom __far*, Atom __far*,
				Atom __far*, Atom __far*, int, int );
#endif

Bond __far *ProcessBond( src, dst, flag )
    Atom __far *src, __far *dst;
    int flag;
{
    register Bond __far *ptr;
    register int i;

    if( flag & (DoubBondFlag|TripBondFlag) )
        DrawDoubleBonds = True;

    if( !(ptr = FreeBond) )
    {   MemSize += BondPool*sizeof(Bond);
	ptr = (Bond __far *)_fmalloc( BondPool*sizeof(Bond) );
	if( !ptr ) FatalDataError("Memory allocation failed");
	RegisterAlloc( ptr );
	for( i=1; i<BondPool; i++ )
	{   ptr->bnext = FreeBond;
	    FreeBond = ptr++;
	} 
    } else FreeBond = ptr->bnext;

    ptr->flag = flag | SelectFlag;
    ptr->srcatom = src;
    ptr->dstatom = dst;
    ptr->radius = 0;
    ptr->col = 0;

    return( ptr );
}

static void CreateHydrogenBond( srcCA, dstCA, src, dst, energy, offset )
    Atom __far *srcCA, __far *dstCA;
    Atom __far *src, __far *dst;
    int energy, offset;
{
    register HBond __far *ptr;
    register int i,flag;

    if( !(ptr = FreeHBond) )
    {   MemSize += HBondPool*sizeof(HBond);
	ptr = (HBond __far *)_fmalloc( HBondPool*sizeof(HBond) );
	if( !ptr ) FatalDataError("Memory allocation failed");
	RegisterAlloc( ptr );
	for( i=1; i<HBondPool; i++ )
	{   ptr->hnext = FreeHBond;
	    FreeHBond = ptr++;
	} 
    } else FreeHBond = ptr->hnext;

    if( (offset>=-128) && (offset<127) )
    {   ptr->offset = (Char)offset;
    } else ptr->offset = 0;

    flag = ZoneBoth? src->flag&dst->flag : src->flag|dst->flag;
    ptr->flag = flag & SelectFlag;

    ptr->src = src;
    ptr->dst = dst;
    ptr->srcCA = srcCA;
    ptr->dstCA = dstCA;
    ptr->energy = energy;
    ptr->col = 0;

    *CurHBond = ptr;
    ptr->hnext = (void __far*)0;
    CurHBond = &ptr->hnext;
    InfoHBondCount++;
}


void CreateBond( src, dst, flag )
    int src, dst, flag;
{
    register Chain __far *chain;
    register Group __far *group;
    register Atom __far *aptr;
    register Atom __far *sptr;
    register Atom __far *dptr;
    register Bond __far *bptr;
    register int done;


    if( src == dst )
	return;

    sptr = (void __far*)0;
    dptr = (void __far*)0;

    done = False;
    for( chain=Database->clist; chain && !done; chain=chain->cnext )
	for( group=chain->glist; group && !done; group=group->gnext )
	    for( aptr=group->alist; aptr; aptr=aptr->anext )
	    {   if( aptr->serno == src )
		{   sptr = aptr;
		    if( dptr )
		    {   done = True;
			break;
		    }
		} else if( aptr->serno == dst )
		{   dptr = aptr;
		    if( sptr )
		    {   done = True;
			break;
		    }
		}
	    }


    /* Both found! */
    if( done ) 
    {   if( flag )
        {   /* Reset Non-bonded flags! */
	    sptr->flag &= ~NonBondFlag;
	    dptr->flag &= ~NonBondFlag;

	    bptr = ProcessBond( sptr, dptr, flag );
	    bptr->bnext = CurMolecule->blist;
	    CurMolecule->blist = bptr;
	    InfoBondCount++;

        } else /* Hydrogen Bond! */
        {   if( InfoHBondCount<0 ) 
            {   CurHBond = &CurMolecule->hlist;
                InfoHBondCount = 0;
            }
            CreateHydrogenBond( NULL, NULL, sptr, dptr, 0, 0 );
        }
    }
}


void CreateBondOrder( src, dst )
    int src, dst;
{
    register Bond __far *bptr;
    register int bs,bd;

    ForEachBond
    {   bs = bptr->srcatom->serno;
	bd = bptr->dstatom->serno;

	if( ((bs==src)&&(bd==dst)) || ((bs==dst)&&(bd==src)) )
        {   DrawDoubleBonds = True;
	    if( bptr->flag & NormBondFlag )
	    {  /* Convert Single to Double */
	       bptr->flag &= ~(NormBondFlag);
	       bptr->flag |= DoubBondFlag;
	    } else if( bptr->flag & DoubBondFlag )
	    {  /* Convert Double to Triple */
	       bptr->flag &= ~(DoubBondFlag);
	       bptr->flag |= TripBondFlag;
	    }
	    return;
	}
    }
    CreateBond( src, dst, NormBondFlag );
}


static void TestBonded( sptr, dptr, flag )
    Atom __far *sptr, __far *dptr; 
    int flag;
{
    register Bond __far *bptr;
    register Long dx, dy, dz;
    register Long max, dist;

    if( flag )
    {    /* Sum of covalent radii with 0.56A tolerance */
         dist = Element[sptr->elemno].covalrad + 
                Element[dptr->elemno].covalrad + 140;
         max = dist*dist;  
    } else 
    {    /* Fast Bio-Macromolecule Bonding Calculation */
         if( (sptr->flag|dptr->flag) & HydrogenFlag )
	 {      max = MaxHBondDist;
         } else max = MaxBondDist;
    }

    dx = sptr->xorg-dptr->xorg;   if( (dist=dx*dx)>max ) return;
    dy = sptr->yorg-dptr->yorg;   if( (dist+=dy*dy)>max ) return;
    dz = sptr->zorg-dptr->zorg;   if( (dist+=dz*dz)>max ) return;

    if( dist > MinBondDist )
    {   /* Reset Non-bonded flags! */
	sptr->flag &= ~NonBondFlag;
	dptr->flag &= ~NonBondFlag;

	bptr = ProcessBond(sptr,dptr,NormBondFlag);
	bptr->bnext = CurMolecule->blist;
	CurMolecule->blist = bptr;
	InfoBondCount++;
    }
}


static void ReclaimHBonds( ptr )
    HBond __far *ptr;
{
    register HBond __far *temp;

    if( (temp = ptr) )
    {   while( temp->hnext )
	    temp=temp->hnext;
	temp->hnext = FreeHBond;
	FreeHBond = ptr;
    }
}


static void ReclaimBonds( ptr )
    Bond __far *ptr;
{
    register Bond __far *temp;

    if( (temp = ptr) )
    {   while( temp->bnext )
	    temp=temp->bnext;
	temp->bnext = FreeBond;
	FreeBond = ptr;
    }
}


static Bond __far *ExtractBonds( ptr )
    Bond __far *ptr;
{
    register Bond __far *result;
    register Bond __far *temp;

    result = (Bond __far*)0;

    while( (temp = ptr) )
    {   ptr = temp->bnext;
        if( temp->flag & NormBondFlag )
        {   temp->bnext = FreeBond;
	    FreeBond = temp;
        } else /* Double or Triple! */
        {   temp->bnext = result;
            result = temp;
        }
    }
    return( result );
}


static void InsertBonds( list, orig )
    Bond __far **list;  Bond __far *orig;
{
    register Atom __far *src;
    register Atom __far *dst;
    register Bond __far *temp;
    register Bond __far *ptr;

    while( (ptr=orig) )
    {   orig = ptr->bnext;
        src = ptr->srcatom;
        dst = ptr->dstatom;
        for( temp=*list; temp; temp=temp->bnext )
            if( ((temp->srcatom==src)&&(temp->dstatom==dst)) ||
                ((temp->srcatom==dst)&&(temp->dstatom==src)) )
                break;

        if( temp )
        {   temp->flag = ptr->flag;
            ptr->bnext = FreeBond;
            FreeBond = ptr;
        } else
        {   ptr->bnext = *list;
            *list = ptr;
        }
    }
}


void CreateMoleculeBonds( info, flag )
    int info, flag;
{
    register int i, x, y, z;
    register Long tx, ty, tz;
    register Long mx, my, mz; 
    register Long dx, dy, dz;
    register int lx, ly, lz, hx, hy, hz;
    register Atom __far *aptr, __far *dptr;
    register Chain __far *chain;
    register Group __far *group;
    register Bond __far *list;
    char buffer[40];


    if( !Database ) 
	return;

    dx = (MaxX-MinX)+1;
    dy = (MaxY-MinY)+1;
    dz = (MaxZ-MinZ)+1;

    /* Save Explicit Double and Triple Bonds! */
    list = ExtractBonds( CurMolecule->blist );
    CurMolecule->blist = (Bond __far*)0;
    InfoBondCount = 0;

    ResetVoxelData();

    for( chain=Database->clist; chain; chain=chain->cnext )
    {   ResetVoxelData();
	for( group=chain->glist; group; group=group->gnext )
	    for( aptr=group->alist; aptr; aptr=aptr->anext )
	    {   /* Initially non-bonded! */
		aptr->flag |= NonBondFlag;

		mx = aptr->xorg-MinX;
		my = aptr->yorg-MinY;
		mz = aptr->zorg-MinZ;

		tx = mx-AbsMaxBondDist;  
		ty = my-AbsMaxBondDist;  
		tz = mz-AbsMaxBondDist;  

		lx = (tx>0)? (int)((VOXORDER*tx)/dx) : 0;
		ly = (ty>0)? (int)((VOXORDER*ty)/dy) : 0;
		lz = (tz>0)? (int)((VOXORDER*tz)/dz) : 0;

		tx = mx+AbsMaxBondDist;  
		ty = my+AbsMaxBondDist;  
		tz = mz+AbsMaxBondDist;

		hx = (tx<dx)? (int)((VOXORDER*tx)/dx) : VOXORDER-1;
		hy = (ty<dy)? (int)((VOXORDER*ty)/dy) : VOXORDER-1;
		hz = (tz<dz)? (int)((VOXORDER*tz)/dz) : VOXORDER-1;
	
		for( x=lx; x<=hx; x++ )
		{   i = VOXORDER2*x + VOXORDER*ly;
		    for( y=ly; y<=hy; y++ )
		    {   for( z=lz; z<=hz; z++ )
			    if( (dptr = (Atom __far*)HashTable[i+z]) )
				do { TestBonded(aptr,dptr,flag);
				} while( (dptr = dptr->next) );
			i += VOXORDER;
		    }
		}
		
		x = (int)((VOXORDER*mx)/dx);
		y = (int)((VOXORDER*my)/dy);
		z = (int)((VOXORDER*mz)/dz);

		i = VOXORDER2*x + VOXORDER*y + z;
		aptr->next = (Atom __far*)HashTable[i];
		HashTable[i] = (void __far*)aptr;
	    }
	VoxelsClean = False;
    }

    /* Replace Double & Triple Bonds! */
    InsertBonds(&CurMolecule->blist,list);

    if( info )
    {   if( CommandActive )
	    WriteChar('\n');
	CommandActive=False;
	sprintf(buffer,"Number of Bonds ..... %ld\n\n",(long)InfoBondCount);
	WriteString(buffer);
    }
}


/*===============================*/
/* Disulphide bridging functions */
/*===============================*/

#ifdef FUNCPROTO
static Atom __far *FindCysSulphur(  Group __far* );
#endif

static Atom __far *FindCysSulphur( group )
    Group __far *group;
{
    register Atom __far *ptr;
    register char *elem;

    for( ptr=group->alist; ptr; ptr=ptr->anext )
    {   elem = ElemDesc[ptr->refno];
        if( (elem[1]=='S') && (elem[0]==' ')  )
            return( ptr );
    }
    return( (Atom __far*)0 );
}


static void TestDisulphideBridge( group1, group2, cys1 )
    Group __far *group1, __far *group2;  
    Atom __far *cys1;
{
    register HBond __far *ptr;
    register Atom __far *cys2;
    register int dx, dy, dz;
    register Long max,dist;

    if( !(cys2=FindCysSulphur(group2)) )
	return;

    max = (Long)750*750;
    dx = (int)(cys1->xorg-cys2->xorg);   if( (dist=(Long)dx*dx)>max ) return;
    dy = (int)(cys1->yorg-cys2->yorg);   if( (dist+=(Long)dy*dy)>max ) return;
    dz = (int)(cys1->zorg-cys2->zorg);   if( (dist+=(Long)dz*dz)>max ) return;

    if( !(ptr = FreeHBond) )
    {   MemSize += sizeof(HBond);
	ptr = (HBond __far*)_fmalloc(sizeof(HBond));
	if( !ptr ) FatalDataError("Memory allocation failed");
	RegisterAlloc( ptr );
    } else FreeHBond = ptr->hnext;

    ptr->dst = cys1;
    if( !(ptr->dstCA=FindGroupAtom(group1,1)) )
	ptr->dstCA = cys1;

    ptr->src = cys2;
    if( !(ptr->srcCA=FindGroupAtom(group2,1)) )
	ptr->srcCA = cys2;

    ptr->offset = 0;
    ptr->energy = 0;
    ptr->flag = 0;
    ptr->col = 0;

    ptr->hnext = CurMolecule->slist;
    CurMolecule->slist = ptr;

    group1->flag |= CystineFlag;
    group2->flag |= CystineFlag;
    InfoSSBondCount++;
}


void FindDisulphideBridges()
{
    register Chain __far *chn1;
    register Chain __far *chn2;
    register Group __far *group1;
    register Group __far *group2;
    register Atom __far *cys;
    char buffer[40];

    if( !Database ) return;
    ReclaimHBonds( CurMolecule->slist );
    InfoSSBondCount = 0;

    for(chn1=Database->clist;chn1;chn1=chn1->cnext)
	for(group1=chn1->glist;group1;group1=group1->gnext)
	    if( IsCysteine(group1->refno) && (cys=FindCysSulphur(group1)) )
	    {   for(group2=group1->gnext;group2;group2=group2->gnext)
		    if( IsCysteine(group2->refno) )
			TestDisulphideBridge(group1,group2,cys);

		for(chn2=chn1->cnext;chn2;chn2=chn2->cnext)
		    for(group2=chn2->glist;group2;group2=group2->gnext)
			if( IsCysteine(group2->refno) )
			    TestDisulphideBridge(group1,group2,cys);
	    }

    if( FileDepth == -1 )
    {   if( CommandActive )
            WriteChar('\n');
        CommandActive=False;
    
        sprintf(buffer,"Number of Bridges ... %d\n\n",InfoSSBondCount);
        WriteString(buffer);
    }
}


/*=========================================*/
/* Kabsch & Sander Structure Determination */
/*=========================================*/

#ifdef FUNCPROTO
static int CalculateBondEnergy( Group __far* );
static void CalcProteinHBonds( Chain __far* );
static void CalcNucleicHBonds( Chain __far* );
static int IsHBonded( Atom __far*, Atom __far*, HBond __far* );
static void TestLadder( Chain __far* );
#endif


/* Coupling constant for Electrostatic Energy   */
/* QConst = -332 * 0.42 * 0.2 * 1000.0 [*250.0] */
#define QConst (-6972000.0)
#define MaxHDist ((Long)2250*2250)
#define MinHDist ((Long)125*125)


/* Protein Donor Atom Coordinates */
static int hxorg,hyorg,hzorg;
static int nxorg,nyorg,nzorg;
static Atom __far *best1CA;
static Atom __far *best2CA;
static Atom __far *best1;
static Atom __far *best2;
static Atom __far *optr;
static int res1,res2;
static int off1,off2;


static int CalculateBondEnergy( group )
    Group __far *group;
{
    register double dho,dhc;
    register double dnc,dno;

    register Atom __far *cptr;
    register Long dx,dy,dz;
    register Long dist;
    register int result;

    if( !(cptr=FindGroupAtom(group,2)) )  return(0);
    if( !(optr=FindGroupAtom(group,3)) )  return(0);

    dx = hxorg-optr->xorg;  
    dy = hyorg-optr->yorg;  
    dz = hzorg-optr->zorg;
    dist = dx*dx+dy*dy+dz*dz;
    if( dist < MinHDist ) 
	return( -9900 );
    dho = sqrt((double)dist);

    dx = hxorg-cptr->xorg;  
    dy = hyorg-cptr->yorg;  
    dz = hzorg-cptr->zorg;
    dist = dx*dx+dy*dy+dz*dz;
    if( dist < MinHDist ) 
	return( -9900 );
    dhc = sqrt((double)dist);

    dx = nxorg-cptr->xorg;  
    dy = nyorg-cptr->yorg;  
    dz = nzorg-cptr->zorg;
    dist = dx*dx+dy*dy+dz*dz;
    if( dist < MinHDist ) 
	return( -9900 );
    dnc = sqrt((double)dist);

    dx = nxorg-optr->xorg;  
    dy = nyorg-optr->yorg;  
    dz = nzorg-optr->zorg;
    dist = dx*dx+dy*dy+dz*dz;
    if( dist < MinHDist ) 
	return( -9900 );
    dno = sqrt((double)dist);

    result = (int)(QConst/dho - QConst/dhc + QConst/dnc - QConst/dno);

    if( result<-9900 ) 
    {   return( -9900 );
    } else if( result>-500 ) 
	return( 0 );

    return( result );
}


static void CalcProteinHBonds( chn1 )
    Chain __far *chn1;
{
    register int energy, offset;
    register Chain __far *chn2;
    register Group __far *group1;
    register Group __far *group2;
    register Atom __far *ca1;
    register Atom __far *ca2;
    register Atom __far *pc1;
    register Atom __far *po1;
    register Atom __far *n1;
    register int pos1,pos2;
    register int dx,dy,dz;
    register double dco;
    register Long dist;

    pos1 = 0;
    pc1 = po1 = (void __far*)0;
    for(group1=chn1->glist;group1;group1=group1->gnext)
    {   pos1++;
        if( pc1 && po1 )
        {   dx = (int)(pc1->xorg - po1->xorg);
            dy = (int)(pc1->yorg - po1->yorg);
            dz = (int)(pc1->zorg - po1->zorg);
        } else
        {   pc1 = FindGroupAtom(group1,2);
            po1 = FindGroupAtom(group1,3);
	    continue;
	}

        pc1 = FindGroupAtom(group1,2);
        po1 = FindGroupAtom(group1,3);

        if( !IsAmino(group1->refno) || IsProline(group1->refno) )
            continue;

        if( !(ca1=FindGroupAtom(group1,1)) ) continue;
        if( !(n1=FindGroupAtom(group1,0)) )  continue;

	dist = (Long)dx*dx + (Long)dy*dy + (Long)dz*dz;
	dco = sqrt( (double)dist )/250.0;

	nxorg = (int)n1->xorg;   hxorg = nxorg + (int)(dx/dco);
	nyorg = (int)n1->yorg;   hyorg = nyorg + (int)(dy/dco);
	nzorg = (int)n1->zorg;   hzorg = nzorg + (int)(dz/dco);
	res1 = res2 = 0;

	/* Only Hydrogen Bond within a single chain!       */
	/* for(chn2=Database->clist;chn2;chn2=chn2->cnext) */

	chn2 = chn1;
	{   /* Only consider non-empty peptide chains! */
	    /* if( !chn2->glist || !IsProtein(chn2->glist->refno) ) */
	    /*     continue;                                        */

	    pos2 = 0;
	    for(group2=chn2->glist;group2;group2=group2->gnext)
	    {   pos2++;
		if( (group2==group1) || (group2->gnext==group1) )
		    continue;

		if( !IsAmino(group2->refno) ) 
		    continue;
		if( !(ca2=FindGroupAtom(group2,1)) ) 
		    continue;

		dx = (int)(ca1->xorg-ca2->xorg);
		if( (dist=(Long)dx*dx) > MaxHDist )
		    continue;

		dy = (int)(ca1->yorg-ca2->yorg);
		if( (dist+=(Long)dy*dy) > MaxHDist )
		    continue;

		dz = (int)(ca1->zorg-ca2->zorg);
		if( (dist+=(Long)dz*dz) > MaxHDist )
		    continue;

		if( (energy = CalculateBondEnergy(group2)) )
		{   if( chn1 == chn2 )
		    {   offset = pos1 - pos2;
		    } else offset = 0;

		    if( energy<res1 )
		    {   best2CA = best1CA;  best1CA = ca2;
			best2 = best1;      best1 = optr;
			res2 = res1;        res1 = energy;
			off2 = off1;        off1 = offset;
		    } else if( energy<res2 )
		    {   best2CA = ca2;
			best2 = optr;
			res2 = energy;
			off2 = offset;
		    }
		}
	    }  /* group2 */
	}      /* chn2 */

	if( res1 ) 
	{   if( res2 ) 
		CreateHydrogenBond(ca1,best2CA,n1,best2,res2,off2);
	    CreateHydrogenBond(ca1,best1CA,n1,best1,res1,off1);
	}
    }
}


static void CalcNucleicHBonds( chn1 )
    Chain __far *chn1;
{
    register Chain __far *chn2;
    register Group __far *group1;
    register Group __far *group2;
    register Group __far *best;
    register Atom __far *ca1;
    register Atom __far *ca2;
    register Atom __far *n1;
    register Long max,dist;
    register int dx,dy,dz;
    register int refno;


    for(group1=chn1->glist;group1;group1=group1->gnext)
    {   if( !IsPurine(group1->refno) ) continue;
	/* Find N1 of Purine Group */
	if( !(n1=FindGroupAtom(group1,21)) )
	    continue;

	/* Maximum N1-N3 distance 5A */
	refno = NucleicCompl(group1->refno);
	max = (Long)1250*1250;
	best = (void __far*)0;

	for(chn2=Database->clist;chn2;chn2=chn2->cnext)
	{   /* Only consider non-empty nucleic acid chains! */
	    if( (chn1==chn2) || !chn2->glist || 
		!IsDNA(chn2->glist->refno) )
		continue;

	    for(group2=chn2->glist;group2;group2=group2->gnext)
		if( group2->refno == (Byte)refno )
		{   /* Find N3 of Pyramidine Group */
		    if( !(ca1=FindGroupAtom(group2,23)) )
			continue;

		    dx = (int)(ca1->xorg - n1->xorg);
		    if( (dist=(Long)dx*dx) >= max ) 
			continue;

		    dy = (int)(ca1->yorg - n1->yorg);
		    if( (dist+=(Long)dy*dy) >= max ) 
			continue;

		    dz = (int)(ca1->zorg - n1->zorg);
		    if( (dist+=(Long)dz*dz) >= max )
			continue;

		    best1 = ca1;
		    best = group2;
		    max = dist;
		}
	}

	if( best )
	{   /* Find the sugar phosphorous atoms */
	    ca1 = FindGroupAtom( group1, 7 );
	    ca2 = FindGroupAtom( best, 7 );

	    CreateHydrogenBond( ca1, ca2, n1, best1, 0, 0 );
	    if( IsGuanine(group1->refno) )
	    {   /* Guanine-Cytosine */
		if( (ca1=FindGroupAtom(group1,22)) &&  /* G.N2 */
		    (ca2=FindGroupAtom(best,26)) )     /* C.O2 */
		    CreateHydrogenBond( (void __far*)0, (void __far*)0,
					ca1, ca2, 0, 0 );

		if( (ca1=FindGroupAtom(group1,28)) &&  /* G.O6 */
		    (ca2=FindGroupAtom(best,24)) )     /* C.N4 */
		    CreateHydrogenBond( (void __far*)0, (void __far*)0,
					ca1, ca2, 0, 0 );

	    } else /* Adenine-Thymine */
		if( (ca1=FindGroupAtom(group1,25)) &&  /* A.N6 */
		    (ca2=FindGroupAtom(best,27)) )     /* T.O4 */
		    CreateHydrogenBond( (void __far*)0, (void __far*)0,
					ca1, ca2, 0, 0 );
	}
    }
}


void CalcHydrogenBonds()
{
    register Chain __far *chn1;
    char buffer[40];

    if( !Database ) return;
    ReclaimHBonds( CurMolecule->hlist );
    CurMolecule->hlist = (void __far*)0;
    CurHBond = &CurMolecule->hlist;
    InfoHBondCount = 0;

    if( MainAtomCount > 10000 )
    {   if( CommandActive )
	    WriteChar('\n');
	WriteString("Please wait... ");
	CommandActive=True;
    }

    for(chn1=Database->clist; chn1; chn1=chn1->cnext)
        if( chn1->glist )
        {   if( IsProtein(chn1->glist->refno) )
	    {   CalcProteinHBonds(chn1);
	    } else if( IsDNA(chn1->glist->refno) )
                CalcNucleicHBonds(chn1);
        }

    if( FileDepth == -1 )
    {   if( CommandActive )
	    WriteChar('\n');
        CommandActive=False;
    
        sprintf(buffer,"Number of H-Bonds ... %d\n",InfoHBondCount);
        WriteString(buffer);
    }
}


static int IsHBonded( src, dst, ptr )
    Atom __far *src, __far *dst;
    HBond __far *ptr;
{
    while( ptr && (ptr->srcCA==src) )
	if( ptr->dstCA == dst )
	{   return( True );
	} else ptr=ptr->hnext;
    return( False );
}


static void FindAlphaHelix( pitch, flag )
    int pitch, flag;
{
    register HBond __far *hbond;
    register Chain __far *chain;
    register Group __far *group;
    register Group __far *first;
    register Group __far *ptr;
    register Atom __far *srcCA;
    register Atom __far *dstCA;
    register int res,dist,prev;

    /* Protein chains only! */
    hbond = Database->hlist;
    for( chain=Database->clist; chain; chain=chain->cnext )
    if( (first=chain->glist) && IsProtein(first->refno) )
    {   prev = False; dist = 0;
	for( group=chain->glist; hbond && group; group=group->gnext )
	{   if( IsAmino(group->refno) && (srcCA=FindGroupAtom(group,1)) )
	    {   if( dist==pitch )
		{   res = False;
		    dstCA=FindGroupAtom(first,1);

		    while( hbond && hbond->srcCA == srcCA )
		    {   if( hbond->dstCA==dstCA ) res=True;
			hbond = hbond->hnext;
		    }

		    if( res )
		    {   if( prev )
			{   if( !(first->struc&HelixFlag) ) 
				InfoHelixCount++;

			    ptr = first;
			    do {
				ptr->struc |= flag;
				ptr = ptr->gnext;
			    } while( ptr != group );
			} else prev = True;
		    } else prev = False;
		} else while( hbond && hbond->srcCA==srcCA )
		    hbond = hbond->hnext;
	    } else prev = False;

	    if( group->struc&HelixFlag )
	    {   first = group; prev = False; dist = 1;
	    } else if( dist==pitch )
	    {   first = first->gnext;
	    } else dist++;
	}
    } else if( first && IsNucleo(first->refno) )
	while( hbond && !IsAminoBackbone(hbond->src->refno) )
	    hbond = hbond->hnext;
}


static Atom __far *cprevi, __far *ccurri, __far *cnexti;
static HBond __far *hcurri, __far *hnexti;
static Group __far *curri, __far *nexti;



static void TestLadder( chain )
    Chain __far *chain;
{
    register Atom __far *cprevj, __far *ccurrj, __far *cnextj;
    register HBond __far *hcurrj, __far *hnextj;
    register Group __far *currj, __far *nextj;
    register int count, result, found;

    /* Already part of atleast one ladder */
    found = curri->flag & SheetFlag;
    nextj = nexti->gnext;

    hnextj = hnexti;
    while( hnextj && hnextj->srcCA==cnexti )
	hnextj = hnextj->hnext;

    while( True )
    {   if( nextj )
	    if( IsProtein(chain->glist->refno) )
	    {   count = 1;
		do {
		    cnextj = FindGroupAtom(nextj,1);
		    if( count == 3 )
		    {   if( IsHBonded(cnexti,ccurrj,hnexti) &&
			    IsHBonded(ccurrj,cprevi,hcurrj) )
			{   result = ParaLadder;
			} else if( IsHBonded(cnextj,ccurri,hnextj) &&
				   IsHBonded(ccurri,cprevj,hcurri) )
			{   result = ParaLadder;
			} else if( IsHBonded(cnexti,cprevj,hnexti) &&
				   IsHBonded(cnextj,cprevi,hnextj) )
			{   result = AntiLadder;
			} else if( IsHBonded(ccurrj,ccurri,hcurrj) &&
				   IsHBonded(ccurri,ccurrj,hcurri) )
			{   result = AntiLadder;
			} else result = NoLadder;

			if( result )
			{   curri->struc |= SheetFlag;
			    currj->struc |= SheetFlag;
			    if( found ) return;
			    found = True;
			}
		    } else count++;

		    cprevj = ccurrj; ccurrj = cnextj; 
		    currj = nextj;   hcurrj = hnextj;

		    while( hnextj && hnextj->srcCA==cnextj )
			hnextj = hnextj->hnext;
		} while( (nextj = nextj->gnext) );

	    } else if( IsNucleo(chain->glist->refno) )
		while( hnextj && !IsAminoBackbone(hnextj->src->refno) )
		    hnextj = hnextj->hnext;

	if( (chain = chain->cnext) ) 
	{   nextj = chain->glist;
	} else return;
    }
}


static void FindBetaSheets()
{
    register Chain __far *chain;
    register int ladder;
    register int count;

    hnexti = Database->hlist;
    for( chain=Database->clist; chain; chain=chain->cnext )
	if( (nexti = chain->glist) )
	    if( IsProtein(nexti->refno) )
	    {   count = 1;
		ladder = False;
		do {
		    cnexti = FindGroupAtom(nexti,1);

		    if( count == 3 )
		    {   TestLadder( chain );
			if( curri->struc & SheetFlag )
			{   if( !ladder )
			    {   InfoLadderCount++;
				ladder = True;
			    }
			} else ladder = False;
		    } else count++;

		    cprevi = ccurri; ccurri = cnexti; 
		    curri = nexti;   hcurri = hnexti;
		    while( hnexti && hnexti->srcCA==cnexti )
			hnexti = hnexti->hnext;
		} while( (nexti = nexti->gnext) );

	    } else if( IsNucleo(nexti->refno) )
		while( hnexti && !IsAminoBackbone(hnexti->src->refno) )
		    hnexti = hnexti->hnext;
}


static void FindTurnStructure()
{
    static Atom __far *aptr[5];
    register Chain __far *chain;
    register Group __far *group;
    register Group __far *prev;
    register Atom __far *ptr;
    register Long ux,uy,uz,mu;
    register Long vx,vy,vz,mv;
    register int i,found,len;
    register Real CosKappa;

    for( chain=Database->clist; chain; chain=chain->cnext )
	if( chain->glist && IsProtein(chain->glist->refno) )
	{   len = 0;  found = False;
	    for( group=chain->glist; group; group=group->gnext )
	    {    ptr = FindGroupAtom(group,1);
		 if( ptr && (ptr->flag&BreakFlag) )
		 {   found = False;
		     len = 0;
		 } else if( len==5 )
		 {   for( i=0; i<4; i++ )
			 aptr[i] = aptr[i+1];
		     len = 4;
		 } else if( len==2 )
		     prev = group;

		 aptr[len++] = ptr;
		 if( len==5 ) 
		 {   if( !(prev->struc&(HelixFlag|SheetFlag)) &&
			 aptr[0] && aptr[2] && aptr[4] )
		     {   ux = aptr[2]->xorg - aptr[0]->xorg;
			 uy = aptr[2]->yorg - aptr[0]->yorg;
			 uz = aptr[2]->zorg - aptr[0]->zorg;

			 vx = aptr[4]->xorg - aptr[2]->xorg;
			 vy = aptr[4]->yorg - aptr[2]->yorg;
			 vz = aptr[4]->zorg - aptr[2]->zorg;

			 mu = ux*ux + uy*uy + uz*uz;
			 mv = vx*vx + vz*vz + vy*vy;
			 if( mu && mv )
			 {   CosKappa = (Real)(ux*vx + uy*vy + uz*vz);
			     CosKappa /= sqrt( (Real)mu*mv );
			     if( CosKappa<Cos70Deg )
			     {   if( !found )
				     InfoTurnCount++;
				 prev->struc |= TurnFlag;
			     }
			 }
		     }
		     found = prev->struc&TurnFlag;
		     prev = prev->gnext;
		 } /* len==5 */
	    }
	}
}


static void FindBetaTurns()
{
    static Atom __far *aptr[4];
    register Chain __far *chain;
    register Group __far *group;
    register Group __far *prev;
    register Group __far *next;
    register Atom __far *ptr;
    register Long dx,dy,dz;
    register int found,len;
    register int flag;

    for( chain=Database->clist; chain; chain=chain->cnext )
        if( chain->glist && IsProtein(chain->glist->refno) )
        {   prev = chain->glist;  
            len = 0;  found = False;
            for( next=chain->glist; next; next=next->gnext )
            {   ptr = FindGroupAtom(next,1);
                if( ptr && (ptr->flag&BreakFlag) )
                {   found = False;
                    prev = next;
                    len = 0;
                } else if( len==4 )
                {   aptr[0] = aptr[1];
                    aptr[1] = aptr[2];
                    aptr[2] = aptr[3];
                    aptr[3] = ptr;

                } else aptr[len++] = ptr;
                if( len==4 ) 
                {   flag = False;
                    if( aptr[0] && aptr[3] )
                    {   dx = aptr[3]->xorg - aptr[0]->xorg;
                        dy = aptr[3]->yorg - aptr[0]->yorg;
                        dz = aptr[3]->zorg - aptr[0]->zorg;
                        if( dx*dx + dy*dy + dz*dz < (Long)1750*1750 )
                        {   group = prev;
                            while( group!=next->gnext )
                            {   if( !(group->struc&(HelixFlag|SheetFlag)) )
                                {   group->struc |= TurnFlag;
                                    flag = True;
                                }
                                group = group->gnext;
                            }
                            if( !found && flag ) 
                                InfoTurnCount++;
                        }
                    }
                    prev = prev->gnext;   
                    found = flag;
                } /* len==4 */
            }
        }
}


void DetermineStructure( flag )
    int flag;
{
    register Chain __far *chain;
    register Group __far *group;
    char buffer[40];

    if( !Database )
	return;

    if( InfoHBondCount<0 )
	CalcHydrogenBonds();

    if( InfoHelixCount>=0 )
	for( chain=Database->clist; chain; chain=chain->cnext )
	    for( group=chain->glist; group; group=group->gnext )
		group->struc = 0;

    InfoStrucSrc = SourceCalc;
    InfoLadderCount = 0;
    InfoHelixCount = 0;
    InfoTurnCount = 0;

    if( InfoHBondCount )
    {   FindAlphaHelix(4,Helix4Flag);
	FindBetaSheets();
	FindAlphaHelix(3,Helix3Flag);
	FindAlphaHelix(5,Helix5Flag);

        if( !flag )
	{   FindTurnStructure();
        } else FindBetaTurns();
    }

    if( FileDepth == -1 )
    {   if( CommandActive )
	    WriteChar('\n');
        CommandActive=False;

        sprintf(buffer,"Number of Helices ... %d\n",InfoHelixCount);
        WriteString(buffer);
        sprintf(buffer,"Number of Strands ... %d\n",InfoLadderCount);
        WriteString(buffer);
        sprintf(buffer,"Number of Turns ..... %d\n",InfoTurnCount);
        WriteString(buffer);
    }
}


void RenumberMolecule( start )
    int start;
{
    register Chain __far *chain;
    register Group __far *group;
    register int hinit, minit;
    register int resno;

    if( !Database )
	return;

    hinit = minit = False;
    for( chain=Database->clist; chain; chain=chain->cnext )
    {   resno = start;
	for( group=chain->glist; group; group=group->gnext )
	{   if( group->alist->flag & HeteroFlag )
	    {   if( hinit )
		{   if( resno > MaxHetaRes )
		    {   MaxHetaRes = resno;
		    } else if( resno < MinHetaRes )
			MinHetaRes = resno;
		} else MinHetaRes = MaxHetaRes = resno;
		hinit = True;
	    } else
	    {   if( minit )
		{   if( resno > MaxMainRes )
		    {   MaxMainRes = resno;
		    } else if( resno < MinMainRes )
			MinMainRes = resno;
		} else MinMainRes = MaxMainRes = resno;
		minit = True;
	    }
	    group->serno = resno++;
	}
    }
}


/*===============================*/
/* Molecule Database Maintenance */
/*===============================*/

static void ReclaimAtoms( ptr )
    Atom __far *ptr;
{
    register Atom __far *temp;

    if( (temp = ptr) )
    {   while( temp->anext )
	    temp=temp->anext;
	temp->anext = FreeAtom;
	FreeAtom = ptr;
    }
}

static void ResetDatabase()
{
    Database = CurMolecule = (void __far*)0;
    MainGroupCount = HetaGroupCount = 0;
    InfoChainCount = HetaAtomCount = 0;
    MainAtomCount = InfoBondCount = 0;  
    SelectCount = 0;

    InfoStrucSrc = SourceNone;
    InfoSSBondCount = InfoHBondCount = -1;
    InfoHelixCount = InfoLadderCount = -1;
    InfoTurnCount = -1;

    CurGroup = (void __far*)0;
    CurChain = (void __far*)0;
    CurAtom = (void __far*)0;

    MinX = MinY = MinZ = 0;
    MaxX = MaxY = MaxZ = 0;

    MinMainTemp = MaxMainTemp = 0;
    MinHetaTemp = MaxHetaTemp = 0;
    MinMainRes = MaxMainRes = 0;
    MinHetaRes = MaxHetaRes = 0;

    *InfoMoleculeName = 0;
    *InfoClassification = 0;
    *InfoIdentCode = 0;
    *InfoSpaceGroup = 0;
    *InfoFileName = 0;

    VoxelsClean = False;
    HMinMaxFlag = False;
    MMinMaxFlag = False;
    HasHydrogen = False;
    ElemNo = MINELEM;
    ResNo = MINRES;
    MaskCount = 0;
    NMRModel = 0;
}


void DestroyDatabase()
{
    register void __far *temp;
    register Group __far *gptr;

    if( Database )
    {   ReclaimHBonds( Database->slist );
	ReclaimHBonds( Database->hlist );
	ReclaimBonds( Database->blist );

	while( Database->clist )
	{   ReclaimBonds(Database->clist->blist);
	    if( (gptr = Database->clist->glist) )
	    {   ReclaimAtoms(gptr->alist);
		while( gptr->gnext )
		{   gptr = gptr->gnext;
		    ReclaimAtoms(gptr->alist);
		}
		gptr->gnext = FreeGroup;
		FreeGroup = Database->clist->glist;
	    }
	    temp = Database->clist->cnext;
	    Database->clist->cnext = FreeChain;
	    FreeChain = Database->clist;
	    Database->clist = temp;
	}

	FreeMolecule = Database;
	Database = (void __far*)0;
    }
    ResetDatabase();
}


void PurgeDatabase()
{
#ifdef APPLEMAC
    register AllocRef *ptr;
    register AllocRef *tmp;
    register int i;
    
    /* Avoid Memory Leaks */
    for( ptr=AllocList; ptr; ptr=tmp )
    {   for( i=0; i<ptr->count; i++ )
	    _ffree( ptr->data[i] );
	tmp = ptr->next;
	_ffree( ptr );
    }
#endif
}


void InitialiseDatabase()
{
    FreeMolecule = (void __far*)0;
    FreeHBond = (void __far*)0;
    FreeChain = (void __far*)0;
    FreeGroup = (void __far*)0;
    FreeAtom = (void __far*)0;
    FreeBond = (void __far*)0;

#ifdef APPLEMAC
    AllocList = (void*)0;
#endif

    ResetDatabase();
}

