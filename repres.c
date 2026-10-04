/* repres.c
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
#include <Errors.h>
#ifdef __CONDITIONALMACROS__
#include <Printing.h>
#else
#include <PrintTraps.h>
#endif
#endif
#ifndef sun386
#include <stdlib.h>
#endif

#include <string.h>
#include <ctype.h>
#include <stdio.h>
#include <math.h>

#define REPRES
#include "molecule.h"
#include "graphics.h"
#include "repres.h"
#include "render.h"
#include "command.h"
#include "abstree.h"
#include "transfor.h"
#include "pixutils.h"


typedef struct { int dx, dy, dz; } DotVector;
typedef struct {
        DotVector __far *probe;
        DotVector __far *dots;
        int count;
    } ElemDotStruct;
 

static ElemDotStruct __far *ElemDots;
static Atom __far *Exclude;
static Monitor *FreeMonit;
static Label *FreeLabel;



#define ForEachAtom  for(chain=Database->clist;chain;chain=chain->cnext) \
                     for(group=chain->glist;group;group=group->gnext)    \
                     for(aptr=group->alist;aptr;aptr=aptr->anext)
#define ForEachBond  for(bptr=Database->blist;bptr;bptr=bptr->bnext)
#define ForEachBack  for(chain=Database->clist;chain;chain=chain->cnext) \
                     for(bptr=chain->blist;bptr;bptr=bptr->bnext)
 

static void FatalRepresError(ptr)
    char *ptr;
{
    char buffer[80];
 
    sprintf(buffer,"Renderer Error: Unable to allocate %s!",ptr);
    RasMolFatalExit(buffer);
}
 

/*============================*/
/*  Label Handling Functions  */
/*============================*/


static void ResetLabels()
{
    register Label *ptr;
 
    while( LabelList )
    {   ptr = LabelList;
        LabelList = ptr->next;
        ptr->next = FreeLabel;
        free(ptr->label);
        FreeLabel = ptr;
    }
}
 

void DeleteLabel( label )
    Label *label;
{
    register Label **ptr;
 
    if( label->refcount == 1 )
    {   ptr = &LabelList;
        while( *ptr != label )
        ptr = &(*ptr)->next;
 
        *ptr = label->next;
        label->next = FreeLabel;
        free(label->label);
        FreeLabel = label;
    } else label->refcount--;
}
 

int DeleteLabels()
{
    register Chain __far *chain;
    register Group __far *group;
    register Atom __far *aptr;
    register int result;
 
    if( !Database )
        return( True );
 
    result = True;
 
    ForEachAtom
        if( aptr->flag & SelectFlag )
        {   if( aptr->label )
            {   DeleteLabel( (Label*)aptr->label );
                aptr->label = (void*)0;
            }
            result = False;
        }
    DrawLabels = LabelList? True : False;
    return( result );
}


Label *CreateLabel( text, len )
    char *text;  int len;
{
    register Label *ptr;
 
    /* Test for existing label */
    for( ptr=LabelList; ptr; ptr=ptr->next )
        if( !strcmp(ptr->label,text) )
            return( ptr );
 
    if( FreeLabel )
    {   ptr = FreeLabel;  FreeLabel = ptr->next;
    } else if( !(ptr=(Label*)malloc(sizeof(Label))) )
        FatalRepresError("label");
 
    ptr->label = (char*)malloc(len+1);
    if( !ptr->label ) FatalRepresError("label");
    strcpy(ptr->label,text);
 
    ptr->next = LabelList;
    ptr->refcount = 0;
    LabelList = ptr;
    return( ptr );
}
 
 
void DefineLabels( label )
    char *label;
{
    register Chain __far *chain;
    register Group __far *group;
    register Atom __far *aptr;
    register Label *ptr;
    register char *cptr;
    register int len;
 
 
    if( !Database ) return;
    if( DeleteLabels() )
        return;
 
    len = 0;
    for( cptr=label; *cptr; cptr++ )
        len++;
 
    /* Strip trailing spaces */
    while( len && cptr[-1]==' ' )
    {   cptr--;  len--;
        *cptr = '\0';
    }
 
    if( !len )
        return;
 
    ptr = CreateLabel(label,len);
    DrawLabels = True;
 
    ForEachAtom
        if( aptr->flag & SelectFlag )
        {   aptr->label = ptr;
            ptr->refcount++;
        }
}
 
 
void DefaultLabels( enable )
    int enable;
{
    register Chain __far *chain;
    register Group __far *group;
    register Atom __far *aptr;
    register Label *label1;
    register Label *label2;
 
    if( !Database )
        return;
 
    label1 = (Label*)0;
    label2 = (Label*)0;
 
    if( MainGroupCount > 1 )
    {   ForEachAtom
            if( IsAlphaCarbon(aptr->refno) || IsSugarPhosphate(aptr->refno) )
            {   if( aptr->flag & SelectFlag )
                {   if( enable )
                    {   if( InfoChainCount > 1 )
                        {   if( chain->ident[1] || isdigit(chain->ident[0]) )
                            {   if( !label1 )
                                    label1 = CreateLabel("%n%r:%c",7);
                                aptr->label = label1;
                                label1->refcount++;
                            } else
                            {   if( !label2 )
                                    label2 = CreateLabel("%n%r%c",6);
                                aptr->label = label2;
                                label2->refcount++;
                            }
                        } else
                        {   if( !label1 )
                                label1 = CreateLabel("%n%r",4);
                            aptr->label = label1;
                            label1->refcount++;
                        }
                    } else if( aptr->label )
                    {   DeleteLabel( (Label*)aptr->label );
                        aptr->label = (Label*)0;
                    }
                    ReDrawFlag |= RFRefresh;
                }
                break;
            }
 
    } else /* Small Molecule! */
        ForEachAtom
            if( (aptr->flag&SelectFlag) && (aptr->elemno!=6)
                                        && (aptr->elemno!=1) )
            {   if( enable )
                {   if( !label1 )
                        label1 = CreateLabel("%e",2);
                    aptr->label = label1;
                    label1->refcount++;
                } else if( aptr->label )
                {   DeleteLabel( (Label*)aptr->label );
                    aptr->label = (Label*)0;
                }
                ReDrawFlag |= RFRefresh;
            }
 
    DrawLabels = LabelList? True : False;
}


void DisplayLabels()
{
    register Chain __far *chain;
    register Group __far *group;
    register Atom __far *aptr;
    register Label *label;
    register int col,z;
 
    auto char buffer[256];
 
 
    if( !Database )
        return;
 
    if( !UseSlabPlane )
    {   z = ImageRadius + ZOffset;
    } else z = SlabValue - 1;
 
    ForEachAtom
        if( aptr->label )
        {   /* Peform Label Slabbing! */
            if( !ZValid(aptr->z) )
                continue;
 
            label = (Label*)aptr->label;
            FormatLabel(chain,group,aptr,label->label,buffer);
 
            if( !UseLabelCol )
            {   /* Depth-cue atom labels */
                /* col = aptr->col + (ColorDepth*                  */
                /*       (aptr->z+ImageRadius-ZOffset))/ImageSize; */
                col = aptr->col + (ColourMask>>1);
            } else col = LabelCol;
 
            /* (aptr->z+2) + ((aptr->flag & SphereFlag)?aptr->irad:0); */
            DisplayString(aptr->x+4,aptr->y,z,buffer,col);
        }
}


/*==============================*/
/*  Monitor Handling Functions  */
/*==============================*/

#ifdef FUNCPROTO
/* Function Prototype */
void AddMonitors( Atom __far*, Atom __far* );
#endif
 
void DeleteMonitors()
{
    register Monitor *ptr;
 
    while( MonitList )
    {   ptr = MonitList;
        if( ptr->col )
            Shade[Colour2Shade(ptr->col)].refcount--;
 
        MonitList = ptr->next;
        ptr->next = FreeMonit;
        FreeMonit = ptr;
    }
}
 

void AddMonitors( src, dst )
    Atom __far *src, __far *dst;
{
    register Monitor **prev;
    register Monitor *ptr;
    register Long dx,dy,dz;
    register Long dist;
 
    /* Delete an already existing monitor! */
    for( prev=&MonitList; (ptr=*prev); prev=&ptr->next )
         if( ((ptr->src==src) && (ptr->dst==dst)) ||
             ((ptr->src==dst) && (ptr->dst==src)) )
         {   if( ptr->col )
                 Shade[Colour2Shade(ptr->col)].refcount--;
 
             *prev = ptr->next;
             ptr->next = FreeMonit;
             FreeMonit = ptr;
             return;
         }
 
 
    /* Create a new monitor! */
    if( FreeMonit )
    {   ptr = FreeMonit;  FreeMonit = ptr->next;
    } else if( !(ptr=(Monitor*)malloc(sizeof(Monitor))) )
        FatalRepresError("monitor");
 
    dx = src->xorg - dst->xorg;
    dy = src->yorg - dst->yorg;
    dz = src->zorg - dst->zorg;
 
    /* ptr->dist = 100.0*CalcDistance(src,dst) */
    dist = isqrt( dx*dx + dy*dy + dz*dz );
    ptr->dist = (unsigned short)((dist<<1)/5);
 
    ptr->src = src;
    ptr->dst = dst;
    ptr->col = 0;
 
    ptr->next = MonitList;
    MonitList = ptr;
}


void CreateMonitor( src, dst )
    Long src, dst;
{
    register Chain __far *chain;
    register Group __far *group;
    register Atom __far *aptr;
    register Atom __far *sptr;
    register Atom __far *dptr;
    register int done;
    char buffer[20];
 
    if( src == dst )
    {   if( !CommandActive )
            WriteChar('\n');
        WriteString("Error: Duplicate atom serial numbers!\n");
        CommandActive = False;
        return;
    }
 
    done = False;
    sptr = (Atom __far*)0;
    dptr = (Atom __far*)0;
 
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
 
    if( !done )
    {   if( !CommandActive )
            WriteChar('\n');
        WriteString("Error: Atom serial number");
        if( sptr )
        {   sprintf(buffer," %d",dst);
        } else if( dptr )
        {   sprintf(buffer," %d",src);
        } else sprintf(buffer,"s %d and %d",src,dst);
        WriteString(buffer); WriteString(" not found!\n");
        CommandActive = False;
 
    } else AddMonitors( sptr, dptr );
}
 
 
void DisplayMonitors()
{
    register Atom __far *s;
    register Atom __far *d;
    register Monitor *ptr;
    register int x,y,z;
    register int sc,dc;
    register int col;
 
    register char *cptr;
    register int dist;
    char buffer[10];
 
    if( !Database )
        return;
 
    if( !UseSlabPlane )
    {   z = ImageRadius + ZOffset;
    } else z = SlabValue-1;
    buffer[9] = '\0';
    buffer[6] = '.';
 
    for( ptr=MonitList; ptr; ptr=ptr->next )
    {   s = ptr->src;
        d = ptr->dst;
 
        if( !ptr->col )
        {   sc = s->col;
            dc = d->col;
        } else sc = dc = ptr->col;
 
        ClipDashVector(s->x,s->y,s->z,d->x,d->y,d->z,sc,dc);
 
        if( DrawMonitDistance )
            if( ZValid( (s->z+d->z)/2 ) )
            {   x = (s->x+d->x)/2;
                y = (s->y+d->y)/2;
 
                if( !UseLabelCol )
                {   /* Use Source atom colour! */
                    col = sc + (ColourMask>>1);
                } else col = LabelCol;
 
                dist = ptr->dist;
                buffer[8] = (dist%10)+'0';  dist /= 10;
                buffer[7] = (dist%10)+'0';
                cptr = &buffer[5];
 
                if( dist > 9 )
                {   do {
                       dist /= 10;
                       *cptr-- = (dist%10)+'0';
                    } while( dist > 9 );
                    cptr++;
                } else *cptr = '0';
 
                DisplayString(x+4,y,z,cptr,col);
            }
    }
}
 
 
/*=========================*/
/*  Dot Surface Functions  */
/*=========================*/

#ifdef FUNCPROTO
/* Function Prototype */
static void AddDot( Long, Long, Long, int );
static void CheckVDWDot( Long, Long, Long, int );
static int TestSolventDot( Long, Long, Long );
#endif
 
 
void DeleteSurface()
{
    register DotStruct __far *ptr;
    register int shade;
    register int i;
 
    while( DotPtr )
    {   for( i=0; i<DotPtr->count; i++ )
        {   shade = Colour2Shade(DotPtr->col[i]);
            Shade[shade].refcount--;
        }
 
        ptr = DotPtr->next;
        _ffree( DotPtr );
        DotPtr = ptr;
    }
    DrawDots = False;
}
 
 
static void AddDot( x, y, z, col )
    Long x, y, z; int col;
{
    register DotStruct __far *ptr;
    register int i, shade;
 
    if( !DotPtr || (DotPtr->count==DotMax) )
    {   ptr = (DotStruct __far*)_fmalloc(sizeof(DotStruct));
        if( !ptr ) FatalRepresError("dot surface");
 
        ptr->next = DotPtr;
        ptr->count = 0;
        DotPtr = ptr;
    } else ptr = DotPtr;
 
    shade = Colour2Shade(col);
    Shade[shade].refcount++;
 
    i = ptr->count++;
    ptr->col[i] = col;
    ptr->xpos[i] = x;
    ptr->ypos[i] = y;
    ptr->zpos[i] = z;
    DrawDots = True;
}
 
 
static void CheckVDWDot( x, y, z, col )
    Long x, y, z; int col;
{
    register Item __far *item;
    register Atom __far *aptr;
    register int ix,iy,iz;
    register int dx,dy,dz;
    register Long dist;
    register Long rad;
    register int i;
 
 
    ix = (int)((x+Offset)*IVoxRatio);
    iy = (int)((y+Offset)*IVoxRatio);
    iz = (int)((z+Offset)*IVoxRatio);
 
    i = VOXORDER2*ix + VOXORDER*iy + iz;
    for( item=HashTable[i]; item; item=item->list )
        if( item->data != Exclude )
        {   aptr = item->data;
            if( !ProbeRadius )
            {   rad = ElemVDWRadius(aptr->elemno);
            } else rad = ProbeRadius;
            rad = rad*rad;
 
            /* Optimized Test! */
            dx = (int)(aptr->xorg - x);
            if( (dist=(Long)dx*dx) < rad )
            {   dy = (int)(aptr->yorg - y);
                if( (dist+=(Long)dy*dy) < rad )
                {   dz = (int)(aptr->zorg - z);
                    if( (dist+=(Long)dz*dz) < rad )
                        return;
                }
            }
        }
    AddDot( x, y, z, col );
}
 
 
static int TestSolventDot( x, y, z )
    Long x, y, z;
{
    register Item __far *item;
    register Atom __far *aptr;
    register int lx,ly,lz;
    register int ux,uy,uz;
    register int dx,dy,dz;
    register int ix,iy,iz;
    register Long dist;
    register Long rad;
    register int i;
 
    dist = Offset-ProbeRadius;
    lx = (int)((x+dist)*IVoxRatio);
    if( lx >= VOXORDER ) return( True );
    ly = (int)((y+dist)*IVoxRatio);
    if( ly >= VOXORDER ) return( True );
    lz = (int)((z+dist)*IVoxRatio);
    if( lz >= VOXORDER ) return( True );
 
    dist = Offset+ProbeRadius;
    ux = (int)((x+dist)*IVoxRatio);
    if( ux < 0 ) return( True );
    uy = (int)((y+dist)*IVoxRatio);
    if( uy < 0 ) return( True );
    uz = (int)((z+dist)*IVoxRatio);
    if( uz < 0 ) return( True );
 
    if( lx < 0 ) lx = 0;  if( ux >= VOXORDER ) ux = VOXORDER-1;
    if( ly < 0 ) ly = 0;  if( uy >= VOXORDER ) uy = VOXORDER-1;
    if( lz < 0 ) lz = 0;  if( uz >= VOXORDER ) uz = VOXORDER-1;
 
    for( ix=lx; ix<=ux; ix++ )
       for( iy=ly; iy<=uy; iy++ )
          for( iz=lz; iz<=uz; iz++ )
          {   i = VOXORDER2*ix + VOXORDER*iy + iz;
              for( item=HashTable[i]; item; item=item->list )
                  if( item->data != Exclude )
                  {   aptr = item->data;
                      rad = ElemVDWRadius(aptr->elemno);
                      rad = (rad+ProbeRadius)*(rad+ProbeRadius);
 
                      /* Optimized Test! */
                      dx = (int)(aptr->xorg - x);
                      if( (dist=(Long)dx*dx) < rad )
                      {   dy = (int)(aptr->yorg - y);
                          if( (dist+=(Long)dy*dy) < rad )
                          {   dz = (int)(aptr->zorg - z);
                              if( (dist+=(Long)dz*dz) < rad )
                                  return( False );
                          }
                      }
                  }
          }
    return( True );
}
 
 
static void InitElemDots()
{
    register int i,size;
 
    size = MAXELEMNO*sizeof(ElemDotStruct);
    ElemDots = (ElemDotStruct __far*)_fmalloc(size);
    if( !ElemDots ) FatalRepresError("dot vector table");
 
    for( i=0; i<MAXELEMNO; i++ )
        ElemDots[i].count = 0;
}
 
 
static void AddElemDots( elem, density )
    int elem, density;
{
    register DotVector __far *ptr;
    register DotVector __far *probe;
    register Real x, y, z, p, q, xy;
    register int equat,vert,horz;
    register int rad,count;
    register int i,j,k;
    register int temp;
 
 
    if( SolventDots || !ProbeRadius )
    {   rad = ElemVDWRadius(elem);
    } else rad = ProbeRadius;
 
    count = (int)(((Long)density*rad*rad)/((Long)250*250));
    ptr = (DotVector __far*)_fmalloc(count*sizeof(DotVector));
    if( !ptr ) FatalRepresError("dot vectors");
 
    if( SolventDots )
    {   probe = (DotVector __far*)_fmalloc(count*sizeof(DotVector));
        if( !probe ) FatalRepresError("probe vectors");
        temp = rad + ProbeRadius;
    } else probe = NULL;
 
    equat = (int)sqrt(PI*count);
    if( !(vert=equat>>1) )
        vert = 1;
 
    i = 0;
    for( j=0; (i<count) && (j<vert); j++ )
    {   p = (PI*j)/(Real)vert;
        z = cos(p);  xy = sin(p);
        horz = (int)(equat*xy);
        if( !horz ) horz = 1;
 
        for( k=0; (i<count) && (k<horz); k++ )
        {   q = (2.0*PI*k)/(Real)horz;
            x = xy*sin(q);
            y = xy*cos(q);
 
            ptr[i].dx = (int)(rad*x);
            ptr[i].dy = (int)(rad*y);
            ptr[i].dz = (int)(rad*z);
            if( probe )
            {   probe[i].dx = (int)(temp*x);
                probe[i].dy = (int)(temp*y);
                probe[i].dz = (int)(temp*z);
            }
            i++;
        }
    }
    ElemDots[elem].probe = probe;
    ElemDots[elem].dots = ptr;
    ElemDots[elem].count = i;
}
 
 
static void FreeElemDots()
{
    register int i;
 
    for( i=0; i<MAXELEMNO; i++ )
        if( ElemDots[i].count )
            _ffree( ElemDots[i].dots );
    _ffree( ElemDots );
}
 
 
void CalculateSurface( density )
    int density;
{
    register DotVector __far *probe;
    register DotVector __far *ptr;
    register Chain __far *chain;
    register Group __far *group;
    register Atom __far *aptr;
    register int i,count;
    register int elem;
 
    if( !Database )
        return;
 
    DeleteSurface();
    ResetVoxelData();
 
    InitElemDots();
    CreateVoxelData( AllAtomFlag );
    VoxelsClean = False;
 
    ForEachAtom
        if( aptr->flag & SelectFlag )
        {   elem = aptr->elemno;
            if( !ElemDots[elem].count )
                AddElemDots(elem,density);
 
            Exclude = aptr;
            ptr = ElemDots[elem].dots;
            probe = ElemDots[elem].probe;
            count = ElemDots[elem].count;
            if( SolventDots )
            {   for( i=0; i<count; i++ )
                    if( TestSolventDot( aptr->xorg + probe[i].dx,
                                        aptr->yorg + probe[i].dy,
                                        aptr->zorg + probe[i].dz ) )
                        AddDot( aptr->xorg + ptr[i].dx,
                                aptr->yorg + ptr[i].dy,
                                aptr->zorg + ptr[i].dz,
                                aptr->col );
            } else
                for( i=0; i<count; i++ )
                    CheckVDWDot( aptr->xorg + ptr[i].dx,
                                 aptr->yorg + ptr[i].dy,
                                 aptr->zorg + ptr[i].dz,
                                 aptr->col);
        }
 
    FreeElemDots();
}
 
 
void DisplaySurface()
{
    register DotStruct __far *ptr;
    register int xi,yi,zi;
    register Real x,y,z;
    register int i;
 
    for( ptr=DotPtr; ptr; ptr=ptr->next )
        for( i=0; i<ptr->count; i++ )
        {   x = ptr->xpos[i];
            y = ptr->ypos[i];
            z = ptr->zpos[i];
 
            xi = (int)(x*MatX[0]+y*MatX[1]+z*MatX[2])+XOffset;
            if( XValid(xi) )
            {   yi = (int)(x*MatY[0]+y*MatY[1]+z*MatY[2])+YOffset;
                if( YValid(yi) )
                {   zi = (int)(x*MatZ[0]+y*MatZ[1]+z*MatZ[2])+ZOffset;
                    if( ZValid(zi) )
                        PlotDeepPoint(xi,yi,zi,ptr->col[i]);
                }
            }
        }
}
 

/*==============================*/
/*  Ribbon & Cartoon Functions  */
/*==============================*/

 
static void CalculateVInten( ptr )
    Knot *ptr;
{
    register Real inten;
 
    if( !ptr->vsize )
        ptr->vsize = isqrt( (Long)ptr->vnx*ptr->vnx +
                            (Long)ptr->vny*ptr->vny +
                            (Long)ptr->vnz*ptr->vnz ) + 1;
 
#ifdef ORIGINAL
#ifdef INVERT
    inten = ptr->vnx - ptr->vny + ptr->vnz + ptr->vnz;
#else
    inten = ptr->vnx + ptr->vny + ptr->vnz + ptr->vnz;
#endif
    inten /= ptr->vsize*RootSix;
#else
    inten = (Real)ptr->vnz/ptr->vsize;
#endif
 
    if( ptr->vnz < 0 ) inten = -inten;
 
    if( inten > 0.0 )
    {   ptr->vinten = (char)(ColourMask*inten);
    } else ptr->vinten = 0;
}
 

static void CalculateHInten( ptr )
    Knot *ptr;
{
    register Real inten;
 
    /* The intensity of the sides of a protein cartoon
     * may be calculated using ptr->cx,cy,cz and this
     * should save interpolating ptr->hnx,hny,hnz!
     */
 
    if( !ptr->hsize )
        ptr->hsize = isqrt( (Long)ptr->hnx*ptr->hnx +
                            (Long)ptr->hny*ptr->hny +
                            (Long)ptr->hnz*ptr->hnz ) + 1;
 
#ifdef ORIGINAL
#ifdef INVERT
    inten = ptr->hnx - ptr->hny + ptr->hnz + ptr->hnz;
#else
    inten = ptr->hnx + ptr->hny + ptr->hnz + ptr->hnz;
#endif
    inten /= ptr->hsize*RootSix;
#else
    inten = (Real)ptr->hnz / ptr->hsize;
#endif
    if( ptr->hnz < 0 ) inten = -inten;
 
    if( inten > 0.0 )
    {   ptr->hinten = (char)(ColourMask*inten);
    } else ptr->hinten = 0;
}
 
 
void DisplayRibbon( chain )
    Chain  __far *chain;
{
    register Group __far *group;
    register Atom __far *captr;
    register Atom __far *o1ptr;
    register Atom __far *o2ptr;
    register Atom __far *next;
 
    register int prev,wide;
    register int col1,col2;
    register int bx,by,bz;
    register int dx,dy,dz;
    register int arrow;
    register int size;
 
    static Knot mid1, mid2, mid3;
    static Knot knot1, knot2;
 
    prev = False;
    group = chain->glist;
    if( IsProtein(group->refno) )
    {   captr = FindGroupAtom(group,1);
    } else captr = FindGroupAtom(group,7);
 
    while( group->gnext )
    {   if( IsProtein(group->gnext->refno) )
        {   next = FindGroupAtom(group->gnext,1);
            o1ptr = FindGroupAtom(group,3);
        } else /* Nucleic Acid */
        {   next = FindGroupAtom(group->gnext,7);
            o1ptr = FindGroupAtom(group->gnext,10);
        }
 
        /* When not to have a control point! */
        if( !next || !captr || !o1ptr || (next->flag&BreakFlag) ||
            !((group->flag|group->gnext->flag)&DrawKnotFlag) )
        {   group = group->gnext;
            captr = next;
            prev = False;
            continue;
        }
 
        knot2.tx = next->x - captr->x;
        knot2.ty = next->y - captr->y;
        knot2.tz = next->z - captr->z;
 
        if( IsProtein(group->refno) )
        {   bx = o1ptr->x - captr->x;
            by = o1ptr->y - captr->y;
            bz = o1ptr->z - captr->z;
 
        } else if( !FindGroupAtom(group,17) &&
                   (o2ptr=FindGroupAtom(group,8)) )
        {   /* Deoxyribonucleic Acid */
            o2ptr = FindGroupAtom(group,8);
            bx = (o1ptr->x + o2ptr->x)/2 - captr->x;
            by = (o1ptr->y + o2ptr->y)/2 - captr->y;
            bz = (o1ptr->z + o2ptr->z)/2 - captr->z;
 
        } else /* Ribonucleic Acid */
        {   bx = o1ptr->x - captr->x;
            by = o1ptr->y - captr->y;
            bz = o1ptr->z - captr->z;
        }
 
        knot2.px = (captr->x + next->x)/2;
        knot2.py = (captr->y + next->y)/2;
        knot2.pz = (captr->z + next->z)/2;
 
        /* c := a x b */
        knot2.vnx = knot2.ty*bz - knot2.tz*by;
        knot2.vny = knot2.tz*bx - knot2.tx*bz;
        knot2.vnz = knot2.tx*by - knot2.ty*bx;
 
        if( (group->struc&group->gnext->struc) & HelixFlag )
        {   /* Compensate for narrowing of helices! */
            size = isqrt((Long)knot2.vnx*knot2.vnx +
                         (Long)knot2.vny*knot2.vny +
                         (Long)knot2.vnz*knot2.vnz);
            knot2.vsize = size;
 
            if( size )
            {   /* 1.00 Angstrom Displacement */
                wide = (int)(250*Scale);
#ifdef INVERT
                knot2.px += (int)(((Long)wide*knot2.vnx)/size);
                knot2.py += (int)(((Long)wide*knot2.vny)/size);
                knot2.pz += (int)(((Long)wide*knot2.vnz)/size);
#else
                knot2.px -= (int)(((Long)wide*knot2.vnx)/size);
                knot2.py -= (int)(((Long)wide*knot2.vny)/size);
                knot2.pz -= (int)(((Long)wide*knot2.vnz)/size);
#endif
            }
        } else knot2.vsize = 0;
 
        if( !(group->flag&group->gnext->flag&TraceFlag) )
        {   /* d := c x a */
            dx = (int)(((Long)knot2.vny*knot2.tz -
                        (Long)knot2.vnz*knot2.ty)/96);
            dy = (int)(((Long)knot2.vnz*knot2.tx -
                        (Long)knot2.vnx*knot2.tz)/96);
            dz = (int)(((Long)knot2.vnx*knot2.ty -
                        (Long)knot2.vny*knot2.tx)/96);
 
            knot2.hsize = isqrt((Long)dx*dx + (Long)dy*dy + (Long)dz*dz);
 
            /* Handle Carbonyl Oxygen Flip */
            if( prev && (((Long)knot1.hnx*dx +
                          (Long)knot1.hny*dy +
                          (Long)knot1.hnz*dz)<0) )
            {   knot2.hnx = -dx;   knot2.vnx = -knot2.vnx;
                knot2.hny = -dy;   knot2.vny = -knot2.vny;
                knot2.hnz = -dz;   knot2.vnz = -knot2.vnz;
            } else
            {   knot2.hnx = dx;
                knot2.hny = dy;
                knot2.hnz = dz;
            }
 
            arrow = False;
            if( group->flag&CartoonFlag )
            {   if( DrawBetaArrows && (group->struc&SheetFlag) &&
                    !(group->gnext->struc&SheetFlag) )
                {   wide = (3*group->width)>>1;
                    arrow = True;
                } else wide = group->width;
            } else if( group->flag & WideKnotFlag )
            {   /* Average Ribbon Width */
                if( group->gnext->flag & WideKnotFlag )
                {   wide = (group->width+group->gnext->width)>>1;
                } else if( group->gnext->flag & CartoonFlag )
                {   wide = group->gnext->width;
                } else wide = group->width;
            } else wide = group->gnext->width;
 
            /* Set Ribbon Width */
            wide = (int)(wide*Scale);
 
            if( knot2.hsize && !arrow )
            {   size = knot2.hsize;
                knot2.wx = (int)(((Long)wide*knot2.hnx)/size);
                knot2.wy = (int)(((Long)wide*knot2.hny)/size);
                knot2.wz = (int)(((Long)wide*knot2.hnz)/size);
                knot2.wide = (short)wide;
            } else
            {   knot2.wide = 0;
                knot2.wx = 0;
                knot2.wy = 0;
                knot2.wz = 0;
            }
 
            if( group->flag & CartoonFlag )
                if( prev && (knot1.wide!=wide) && knot1.hsize )
                {   size = knot1.hsize;
                    knot1.wx = (int)(((Long)wide*knot1.hnx)/size);
                    knot1.wy = (int)(((Long)wide*knot1.hny)/size);
                    knot1.wz = (int)(((Long)wide*knot1.hnz)/size);
                }
 
            if( (group->flag|group->gnext->flag)&CartoonFlag )
            {   CalculateVInten( &knot2 );
                CalculateHInten( &knot2 );
 
                size = knot2.vsize;
                wide = (int)(CartoonHeight*Scale);
                knot2.dx = (int)(((Long)wide*knot2.vnx)/size);
                knot2.dy = (int)(((Long)wide*knot2.vny)/size);
                knot2.dz = (int)(((Long)wide*knot2.vnz)/size);
            } else if( (group->flag|group->gnext->flag)&RibbonFlag )
                CalculateVInten( &knot2 );
        }
 
        if( !(col1 = group->col1) )
            col1 = captr->col;
 
        if( prev )
        {   /* Approximate spline segment with plane! */
            /* SolidRibbon( &knot1, &knot2, col1 );   */
 
            /* Calculate Hermite Spline Points */
            mid1.px = (int)(((Long)54*knot1.px + (Long)9*knot1.tx +
                             (Long)10*knot2.px - (Long)3*knot2.tx)/64);
            mid1.py = (int)(((Long)54*knot1.py + (Long)9*knot1.ty +
                             (Long)10*knot2.py - (Long)3*knot2.ty)/64);
            mid1.pz = (int)(((Long)54*knot1.pz + (Long)9*knot1.tz +
                             (Long)10*knot2.pz - (Long)3*knot2.tz)/64);
 
            mid2.px = (int)(((Long)4*knot1.px + knot1.tx +
                             (Long)4*knot2.px - knot2.tx)/8);
            mid2.py = (int)(((Long)4*knot1.py + knot1.ty +
                             (Long)4*knot2.py - knot2.ty)/8);
            mid2.pz = (int)(((Long)4*knot1.pz + knot1.tz +
                             (Long)4*knot2.pz - knot2.tz)/8);
 
            mid3.px = (int)(((Long)10*knot1.px + (Long)3*knot1.tx +
                             (Long)54*knot2.px - (Long)9*knot2.tx)/64);
            mid3.py = (int)(((Long)10*knot1.py + (Long)3*knot1.ty +
                             (Long)54*knot2.py - (Long)9*knot2.ty)/64);
            mid3.pz = (int)(((Long)10*knot1.pz + (Long)3*knot1.tz +
                             (Long)54*knot2.pz - (Long)9*knot2.tz)/64);
 
            if( group->flag & TraceFlag )
            {   wide = (int)(group->width*Scale);
                ClipCylinder( knot1.px, knot1.py, knot1.pz,
                              mid1.px, mid1.py, mid1.pz,
                              col1, col1, wide );
                ClipCylinder( mid1.px, mid1.py, mid1.pz,
                              mid2.px, mid2.py, mid2.pz,
                              col1, col1, wide );
                ClipCylinder( mid2.px, mid2.py, mid2.pz,
                              mid3.px, mid3.py, mid3.pz,
                              col1, col1, wide );
                ClipCylinder( mid3.px, mid3.py, mid3.pz,
                              knot2.px, knot2.py, knot2.pz,
                              col1, col1, wide );
            } else
            {   /* Calculate Hermite Spline Widths */
                mid1.wx = (27*knot1.wx + 5*knot2.wx)/32;
                mid1.wy = (27*knot1.wy + 5*knot2.wy)/32;
                mid1.wz = (27*knot1.wz + 5*knot2.wz)/32;
 
                mid2.wx = (knot1.wx + knot2.wx)/2;
                mid2.wy = (knot1.wy + knot2.wy)/2;
                mid2.wz = (knot1.wz + knot2.wz)/2;
 
                mid3.wx = (5*knot1.wx + 27*knot2.wx)/32;
                mid3.wy = (5*knot1.wy + 27*knot2.wy)/32;
                mid3.wz = (5*knot1.wz + 27*knot2.wz)/32;
 
                /* Draw the Spline Segments */
                if( group->flag & (StrandFlag|DashStrandFlag) )
                {   if( !(col2 = group->col2) )
                        col2 = captr->col;
                    if( group->flag & StrandFlag )
                    {   StrandRibbon( &knot1, &mid1,  col1, col2 );
                        StrandRibbon( &mid1,  &mid2,  col1, col2 );
                        StrandRibbon( &mid2,  &mid3,  col1, col2 );
                        StrandRibbon( &mid3,  &knot2, col1, col2 );
                    } else /* group->flag & DashStrandFlag */
                    {   DashRibbon( &knot1, &mid1,  col1, col2 );
                        DashRibbon( &mid1,  &mid2,  col1, col2 );
                        DashRibbon( &mid2,  &mid3,  col1, col2 );
                        DashRibbon( &mid3,  &knot2, col1, col2 );
                    }
                } else /* Ribbon or Cartoon! */
                {   mid1.vsize = 0;
                    mid1.vnx = (int)(((Long)27*knot1.vnx +
                                      (Long) 5*knot2.vnx)/32);
                    mid1.vny = (int)(((Long)27*knot1.vny +
                                      (Long) 5*knot2.vny)/32);
                    mid1.vnz = (int)(((Long)27*knot1.vnz +
                                      (Long) 5*knot2.vnz)/32);
                    CalculateVInten( &mid1 );
 
                    mid2.vsize = 0;
                    mid2.vnx = (knot1.vnx + knot2.vnx)/2;
                    mid2.vny = (knot1.vny + knot2.vny)/2;
                    mid2.vnz = (knot1.vnz + knot2.vnz)/2;
                    CalculateVInten( &mid2 );
 
                    mid3.vsize = 0;
                    mid3.vnx = (int)(((Long) 5*knot1.vnx +
                                      (Long)27*knot2.vnx)/32);
                    mid3.vny = (int)(((Long) 5*knot1.vny +
                                      (Long)27*knot2.vny)/32);
                    mid3.vnz = (int)(((Long) 5*knot1.vnz +
                                      (Long)27*knot2.vnz)/32);
                    CalculateVInten( &mid3 );
 
                    if( group->flag & RibbonFlag )
                    {   SolidRibbon( &knot1, &mid1,  col1 );
                        SolidRibbon( &mid1,  &mid2,  col1 );
                        SolidRibbon( &mid2,  &mid3,  col1 );
                        SolidRibbon( &mid3,  &knot2, col1 );
                    } else /* Cartoon! */
                    {   /* Calculate Spline Heights */
                        wide = (int)(CartoonHeight*Scale);
 
                        size = mid1.vsize;
                        mid1.dx = (int)(((Long)wide*mid1.vnx)/size);
                        mid1.dy = (int)(((Long)wide*mid1.vny)/size);
                        mid1.dz = (int)(((Long)wide*mid1.vnz)/size);
 
                        size = mid2.vsize;
                        mid2.dx = (int)(((Long)wide*mid2.vnx)/size);
                        mid2.dy = (int)(((Long)wide*mid2.vny)/size);
                        mid2.dz = (int)(((Long)wide*mid2.vnz)/size);
 
                        size = mid3.vsize;
                        mid3.dx = (int)(((Long)wide*mid3.vnx)/size);
                        mid3.dy = (int)(((Long)wide*mid3.vny)/size);
                        mid3.dz = (int)(((Long)wide*mid3.vnz)/size);
 
                        /* Calculate Surface Intensity */
                        mid1.hsize = 0;
                        mid1.hnx = (int)(((Long)27*knot1.hnx +
                                          (Long) 5*knot2.hnx)/32);
                        mid1.hny = (int)(((Long)27*knot1.hny +
                                          (Long) 5*knot2.hny)/32);
                        mid1.hnz = (int)(((Long)27*knot1.hnz +
                                          (Long) 5*knot2.hnz)/32);
                        CalculateHInten( &mid1 );
 
                        mid2.hsize = 0;
                        mid2.hnx = (knot1.hnx + knot2.hnx)/2;
                        mid2.hny = (knot1.hny + knot2.hny)/2;
                        mid2.hnz = (knot1.hnz + knot2.hnz)/2;
                        CalculateHInten( &mid2 );
 
                        mid3.hsize = 0;
                        mid3.hnx = (int)(((Long) 5*knot1.hnx +
                                          (Long)27*knot2.hnx)/32);
                        mid3.hny = (int)(((Long) 5*knot1.hny +
                                          (Long)27*knot2.hny)/32);
                        mid3.hnz = (int)(((Long) 5*knot1.hnz +
                                          (Long)27*knot2.hnz)/32);
                        CalculateHInten( &mid3 );
 
                        RectRibbon( &knot1, &mid1,  col1 );
                        RectRibbon( &mid1,  &mid2,  col1 );
                        RectRibbon( &mid2,  &mid3,  col1 );
                        RectRibbon( &mid3,  &knot2, col1 );
                    }
                }
            }
        } else if( group == chain->glist )
        {   knot1 = knot2;
            knot1.px = captr->x;
            knot1.py = captr->y;
            knot1.pz = captr->z;
 
            if( group->flag & RibbonFlag )
            {   SolidRibbon( &knot1, &knot2, col1 );
            } else if( group->flag & RibbonFlag )
            {   RectRibbon( &knot1, &knot2, col1 );
            } else if( group->flag & StrandFlag )
            {   if( !(col2 = group->col2) )
                    col2 = captr->col;
                StrandRibbon( &knot1,  &knot2, col1, col2 );
            } else if( group->flag & DashStrandFlag )
            {   if( !(col2 = group->col2) )
                    col2 = captr->col;
                DashRibbon( &knot1,  &knot2, col1, col2 );
            } else if( group->flag & TraceFlag )
                ClipCylinder( knot1.px, knot1.py, knot1.pz,
                              knot2.px, knot2.py, knot2.pz,
                              col1, col1, (int)(group->width*Scale) );
            prev = True;
        } else prev = True;
        group = group->gnext;
        captr = next;
 
        knot1 = knot2;
    }
 
    if( prev )
    {   if( !(col1 = group->col1) )
            col1 = captr->col;
 
        if( group->flag & CartoonFlag )
        {   /* Test for arrow head! */
            if( DrawBetaArrows && (group->struc&SheetFlag) )
            {   wide = (3*group->width)>>1;
                knot2.px = captr->x + (knot2.tx/2);
                knot2.py = captr->y + (knot2.ty/2);
                knot2.pz = captr->z + (knot2.tz/2);
 
                arrow = True;
            } else
            {   wide = group->width;
                knot2.px = captr->x;
                knot2.py = captr->y;
                knot2.pz = captr->z;
                arrow = False;
            }
 
            wide = (int)(Scale*wide);
            if( (knot1.wide!=wide) && knot1.hsize )
            {   size = knot1.hsize;
                knot1.wx = (int)(((Long)wide*knot1.hnx)/size);
                knot1.wy = (int)(((Long)wide*knot1.hny)/size);
                knot1.wz = (int)(((Long)wide*knot1.hnz)/size);
 
                if( !arrow )
                {   knot2.wx = knot1.wx;
                    knot2.wy = knot1.wy;
                    knot2.wz = knot1.wz;
                } else
                {   knot2.wx = 0;
                    knot2.wy = 0;
                    knot2.wz = 0;
                }
            } else if( arrow )
            {   knot2.wx = 0;
                knot2.wy = 0;
                knot2.wz = 0;
            }
 
            RectRibbon( &knot1, &knot2, col1 );
        } else /* !Cartoon */
        {   knot2.px = captr->x;
            knot2.py = captr->y;
            knot2.pz = captr->z;
 
            if( group->flag & RibbonFlag )
            {   SolidRibbon( &knot1, &knot2, col1 );
            } else if( group->flag & StrandFlag )
            {   if( !(col2 = group->col2) )
                    col2 = captr->col;
                StrandRibbon( &knot1,  &knot2, col1, col2 );
            } else if( group->flag & DashStrandFlag )
            {   if( !(col2 = group->col2) )
                    col2 = captr->col;
                DashRibbon( &knot1,  &knot2, col1, col2 );
            } else if( group->flag & TraceFlag )
                ClipCylinder( knot1.px, knot1.py, knot1.pz,
                              knot2.px, knot2.py, knot2.pz,
                              col1, col1, (int)(group->width*Scale) );
        }
    }
}
 
 

void ResetRepres()
{
    DeleteSurface();
    DeleteMonitors();
    SolventDots = False;
    ProbeRadius = 0;

    DrawLabels = False;
    ResetLabels();

    DrawMonitDistance = True;
    DrawBetaArrows = True;
    CartoonHeight = 100;
}


void InitialiseRepres()
{
    DotPtr = (DotStruct __far*)0;
    MonitList = (Monitor __far*)0;
    LabelList = (void*)0;

    FreeMonit = (Monitor __far*)0;
    FreeLabel = (void*)0;

    ResetRepres();
}


