/* outfile.c
 * RasMol2 Molecular Graphics
 * Roger Sayle, August 1995
 * Version 2.6
 */
#define OUTFILE
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

#include <stdio.h>
#include <ctype.h>
#include <math.h>

#include "outfile.h"
#include "molecule.h"
#include "command.h"
#include "abstree.h"
#include "transfor.h"
#include "render.h"
#include "repres.h"
#include "graphics.h"
#include "pixutils.h"
#include "script.h"


#ifdef EIGHTBIT
#define RComp(x)   (RLut[LutInv[x]])
#define GComp(x)   (GLut[LutInv[x]])
#define BComp(x)   (BLut[LutInv[x]])
#else
#define RComp(x)   (((x)>>16)&0xff)
#define GComp(x)   (((x)>>8)&0xff)
#define BComp(x)   ((x)&0xff)
#endif


#ifdef INVERT
#define InvertY(y) (y)
#else
#define InvertY(y) (-(y))
#endif


/* Standard A4 size page: 8.267x11.811 inches */
/* U.S. Normal size page: 8.500x11.000 inches */
#define PAGEHIGH  (11.811*72.0)
#define PAGEWIDE  (8.267*72.0)
#define BORDER    0.90


#define Round(x)       ((int)(x))


#define PSBond      0x00
#define PSHBond     0x01
#define PSSSBond    0x02
#define PSAtom      0x03
#define PSRibbon    0x04
#define PSMonit     0x05



typedef void __far* PSItemPtr;








static Byte LutInv[256];
static FILE *OutFile;

static Real LineWidth;
static int VectSolid;
static int VectCol;

/* Macros for commonly used loops */
#define ForEachAtom  for(chain=Database->clist;chain;chain=chain->cnext) \
		     for(group=chain->glist;group;group=group->gnext)    \
		     for(aptr=group->alist;aptr;aptr=aptr->anext)
#define ForEachBond  for(bptr=Database->blist;bptr;bptr=bptr->bnext)
#define ForEachBack  for(chain=Database->clist;chain;chain=chain->cnext) \
		     for(bptr=chain->blist;bptr;bptr=bptr->bnext)



#ifdef APPLEMAC
/* External RasMac Function Declaration! */
void SetFileInfo( char*, OSType, OSType, short );
#endif

static void FatalOutputError( ptr )
    char *ptr;
{
    if( CommandActive ) WriteChar('\n');
    WriteString("Output Error: Unable to create file `");
    WriteString( ptr );  WriteString("'!\n");
    CommandActive = False;
}


static void WriteByte( val )
    int val;
{
    putc( val, OutFile );
}


static void CalcInvColourMap()
{
#ifdef EIGHTBIT
    register int i;

    for( i=0; i<256; i++ )
        if( ULut[i] )
            LutInv[Lut[i]] = i;
#endif
}


#ifdef EIGHTBIT
static int CompactColourMap()
{
    register Pixel __huge *ptr;
    register Long pos, count;
    register int i, cols;

    CalcInvColourMap();
    for( i=0; i<256; i++ )
    {   Buffer[i] = 0;
	Node[i] = 5;
    }

#ifdef IBMPC
    ptr = (Pixel __huge*)GlobalLock(FBufHandle);    
#else
    ptr = FBuffer;
#endif

    cols = 0;
    count = (Long)XRange*YRange;
    for( pos=0; pos<count; pos++ )
    {   i = LutInv[*ptr++];
	if( !Buffer[i] ) 
	{   Node[cols++] = i;
	    Buffer[i] = cols;
	}
    }

    for( i=0; i<256; i++ )
	LutInv[i] = Buffer[LutInv[i]]-1;
#ifdef IBMPC
    GlobalUnlock(FBufHandle);
#endif
    return( cols );
}
#endif




static void WritePPMWord( i )
    int i;
{
    if( i>99 )
    {   putc((i/100)+'0',OutFile); i %= 100;
	putc((i/10) +'0',OutFile); i %= 10;
    } else if( i>9 )
    {   putc((i/10)+'0',OutFile);  i %= 10;
    }
    putc(i+'0',OutFile);
}


int WritePPMFile( name, raw )
    char *name;  int raw;
{
    register Pixel __huge *ptr;
    register int i,col;
    register int x,y;

#if defined(IBMPC) || defined(APPLEMAC)
    OutFile = fopen(name, (raw?"wb":"w") );
#else
    OutFile = fopen(name,"w");
#endif

    if( !OutFile ) 
    {   FatalOutputError(name);
	return( False );
    }

    CalcInvColourMap();
    fprintf(OutFile,"P%c %d %d 255\n",(raw?'6':'3'),XRange,YRange);

#ifdef IBMPC
    FBuffer = (Pixel __huge*)GlobalLock(FBufHandle);
#endif

#ifndef INVERT
    ptr = FBuffer;
#endif

    if( !raw )
    {   col = 0;
	for( y=YRange-1; y>=0; y-- )
	{
#ifdef INVERT
	    ptr = FBuffer + (Long)y*XRange;
#endif
	    for( x=0; x<XRange; x++ )
	    {   i = *ptr++;
		WritePPMWord((int)RComp(i));  WriteByte(' ');
		WritePPMWord((int)GComp(i));  WriteByte(' ');
		WritePPMWord((int)BComp(i));  
		if( ++col == 5 )
		{   WriteByte('\n');
                    col = 0;
		} else WriteByte(' ');
	    }
	}
    } else
	for( y=YRange-1; y>=0; y-- )
	{
#ifdef INVERT
	    ptr = FBuffer + (Long)y*XRange;
#endif
	    for( x=0; x<XRange; x++ )
	    {   i = *ptr++;
		putc((int)RComp(i),OutFile);
		putc((int)GComp(i),OutFile);
		putc((int)BComp(i),OutFile);
	    }
	}

    fclose(OutFile);
#ifdef APPLEMAC
    /* Avoid ANSI trigraph problems! */
    SetFileInfo(name,'\?\?\?\?','\?\?\?\?',134);
#endif
#ifdef IBMPC
    GlobalUnlock(FBufHandle); 
#endif
    return( True );
}


#ifdef FUNCPROTO
static int FindDepth( PSItemPtr, int );
static void DepthSort( PSItemPtr __far*, char __far*, int );
#endif

static int FindDepth( item, type )
     PSItemPtr item;  int type;
{
    register HBond __far *hbond;
    register Atom __far *atom;
    register Bond __far *bond;
    register Monitor *monit;
    register int result;

    switch( type )
    {   case(PSAtom):    atom = (Atom __far*)item;
	                 return( atom->z );

        case(PSBond):    bond = (Bond __far*)item;
	                 result = bond->srcatom->z;
	                 if( result < bond->dstatom->z )
	                     result = bond->dstatom->z;
                         return( result );

        case(PSSSBond):  
        case(PSHBond):   hbond = (HBond __far*)item;
                         if( (type==PSHBond)? HBondMode : SSBondMode )
                         {   result = hbond->srcCA->z;
                             if( result < hbond->dstCA->z )
                                 result = hbond->dstCA->z;
                         } else
                         {   result = hbond->src->z;
                             if( result < hbond->dst->z )
                                 result = hbond->dst->z;
                         }
                         return( result );

        case(PSMonit):   monit = (Monitor*)item;
                         result = monit->src->z;
                         if( result < monit->dst->z )
                             result = monit->dst->z;
                         return( result );
    }
    return( 0 );
}


static void DepthSort( data, type, count )
    PSItemPtr __far *data;
    char __far *type;
    int count;
{
    register char ttmp;
    register void __far *dtmp;
    register int i, j, k;
    register int depth;
    register int temp;

    for( i=1; i<count; i++ )
    {   dtmp = data[i];  
	ttmp = type[i];

	j = i-1;
	depth = FindDepth(dtmp,ttmp);
	temp = FindDepth(data[j],type[j]);
	while( (depth<temp) || ((depth==temp)&&(ttmp<type[j])) )
	    if( j-- ) 
	    {   temp = FindDepth(data[j],type[j]);
	    } else break;
	j++;

	if( j != i )
	{   for( k=i; k>j; k-- )
	    {    data[k] = data[k-1];
		 type[k] = type[k-1];
	    }
	    data[j] = dtmp;
	    type[j] = ttmp;
	}
    }
}

#ifdef FUNCPROTO
static int ClipVectSphere( Atom __far* );
static int ClipVectBond( Atom __far*, Atom __far* );

static void WriteVectSphere( PSItemPtr __far*, char __far*, int );
static void WriteVectStick( Atom __far*, Atom __far*, int, int );
static void WriteVectWire( Atom __far*, Atom __far*, int, int );

static Long CountPSItems();
static void FetchPSItems( PSItemPtr __far*, char __far* );
static void WritePSItems( PSItemPtr __far*, char __far*, int );
#endif


static int ClipVectSphere( ptr )
    Atom __far *ptr;
{
    register int rad;

    rad = ptr->irad;

    if( ptr->x + rad < 0 )  return( True );
    if( ptr->y + rad < 0 )  return( True );
    if( ptr->x - rad >= XRange )  return( True );
    if( ptr->y - rad >= YRange )  return( True );
    return( False );
}


static int ClipVectBond( src, dst )
    Atom __far *src;
    Atom __far *dst;
{
    if( !src || !dst )  return( True );
    if( (src->x<0) && (dst->x<0) )  return( True );
    if( (src->y<0) && (dst->y<0) )  return( True );
    if( (src->x>=XRange) && (dst->x>=XRange) )  return( True );
    if( (src->y>=YRange) && (dst->y>=YRange) )  return( True );
    return( False );
}



static void WriteVectColour( col )
    int col;
{
    if( col != VectCol )
    {   fprintf(OutFile,"%g ",(Real)RLut[col]/255.0);
	fprintf(OutFile,"%g ",(Real)GLut[col]/255.0);
	fprintf(OutFile,"%g ",(Real)BLut[col]/255.0);
	fputs("setrgbcolor\n",OutFile);
	VectCol = col;
    }
}


#define MAXSECT 5
typedef struct {
	/* Ellipse */
	Real ephi,epsi;
	Real etheta;
	Real ex,ey;
	Real erad;

	/* Sphere */
	Real sphi,spsi;
	int sx,sy;
	Real srad;
    } SphereSect;


static int VectClipContain( x, y )
    SphereSect *x; SphereSect *y;
{
    if( x->erad != 0.0 )
    {   if( y->erad != 0.0 )
	    /* Simple segment containment test! */
	    return( ((x->sphi+x->spsi)>=(y->sphi+y->spsi)) &&
		    ((x->sphi-x->spsi)<=(y->sphi-y->spsi)) );
    } else if( y->erad == 0.0 )
	return( x->srad >= y->srad );
    return( False );
}


static void WriteVectSphere( data, type, index )
    PSItemPtr __far*data; 
    char __far *type;
    int index;
{
    register int ecount, count;
    register Atom __far *atm;
    register Atom __far *ptr;
    register Long dist2,dist3;
    register int dx, dy, dz;
    register int i,j,k;

    register Real b,d,f,g,x;
    register Real radf,radb;
    register Real phi1,phi2;
    register Real temp,psi;
    register Real theta;

    register SphereSect *sptr;
    SphereSect sect[MAXSECT];

    ptr = (Atom __far*)data[index];
    radf = ptr->radius*Scale;

    count = 0;
    ecount = 0;
    sptr = sect;
    for( i=index-1; i>=0; i-- )
    {   if( type[i] != PSAtom )
	    continue;

	atm = (Atom __far*)data[i];
	/* Atom can't intersect visibly! */
	if( atm->z + atm->irad < ptr->z )
	    continue;

	dx = atm->x - ptr->x; 
	dy = atm->y - ptr->y; 
	dz = atm->z - ptr->z;

	dist2 = (Long)dx*dx + (Long)dy*dy;
	dist3 = dist2 + dz*dz;

	radb = atm->radius*Scale;  
	temp = radf + radb;

	/* Atoms don't intersect! */
	if( dist3 > temp*temp ) continue;


	d = sqrt( (double)dist3 );
	f = (temp*(radf-radb)+dist3)/(2.0*d);
	theta = -dz/d;

	if( f>0 )
	{   temp = radf*radf;
	    /* Intersection not visible! */
	    if( theta*temp > temp-f*f )
		continue;
	} else if( f < -radf )
	    return;

	x = sqrt( (radf-f)*(radf+f) );

	if( dx || dy )
	{   g = sqrt( (double)dist2 );
	    psi = Rad2Deg*atan2(dy,dx);
	    b = (f*(dz*dz))/(d*g);

	    if( AbsFun(b)>x )
		continue;

	    phi1 = b + (f*g)/d;
	    phi1 = Rad2Deg*acos(phi1/radf);
	    if( phi1!=phi1 ) continue;

	    phi2 = (d*b)/g;
	    if( AbsFun(phi2) < x )
	    {   phi2 = Rad2Deg*acos(phi2/x);
		if( phi2!=phi2 ) continue;
		if( phi2 > 90.0 ) 
		    phi2 = 180.0-phi1;
	    } else phi2 = 90.0;

	    sptr->erad = x;
	    sptr->etheta = -theta;
	    sptr->ephi = psi;
	    sptr->epsi = phi2;

	    temp = f/d;
	    sptr->ex = ptr->x+temp*dx;
	    sptr->ey = ptr->y+temp*dy;

	    sptr->srad = radf;
	    sptr->sphi = psi;
	    sptr->spsi = phi1;
	    sptr->sx = ptr->x;
	    sptr->sy = ptr->y;

	} else
	{   x = sqrt( (radf-g)*(radf+g) );

	    sptr->srad = x;
	    sptr->erad = 0.0;
	    sptr->sx = ptr->x;
	    sptr->sy = ptr->y;
	    sptr->sphi = 180;
	    sptr->spsi = -180;
	}

	/* Optimize Segments */
	j = 0;
	while( j<count )
	    if( VectClipContain(sptr,sect+j) )
	    {   /* Delete Segment sect[j] */
		for( k=j; k<count; k++ )
		    sect[k] = sect[k+1];
		count--;  sptr--;
	    } else if( VectClipContain(sect+j,sptr) )
	    {   break;  /* Exclude Segment */
	    } else j++;
	   

	if( j==count )
	{   count++;  sptr++;
	    if( sptr->erad != 0.0 )
		ecount++;
	    if( count==MAXSECT )
		break;
	}
    }

    if( UseOutLine )
    {   temp = (ptr->z-ZOffset)/ImageSize + 1.0;
	if( temp != LineWidth )
	{   fprintf(OutFile,"%g setlinewidth\n",temp);
	    LineWidth = temp;
	}
    }

    if( !VectSolid )
    {   fputs("[] 0 setdash\n",OutFile);
        VectSolid = True;
    }

    if( count )
    {   fputs("gsave\n",OutFile);
	fprintf(OutFile,"%%%% %d %d\n",count,ecount);

	sptr = sect;
	for( i=0; i<count; i++ )
	{   if( sptr->erad != 0.0 )
	    {   fprintf(OutFile,"%g %g %g %g %g %g ClipEllips\n",
			    sptr->erad,sptr->epsi,sptr->etheta,
			    sptr->ephi,sptr->ex,sptr->ey);
	    }

	    if( (i==count-1) || (sptr->erad==0.0) )
	    {   fprintf(OutFile,"%g %g %g %d %d ClipSphere\n",sptr->srad,
				sptr->sphi+sptr->spsi,sptr->sphi-sptr->spsi,
				sptr->sx, sptr->sy );
	    } else fprintf(OutFile,"%g %g %g %d %d ClipBox\n",
				    sptr->srad+sptr->srad+2,
				    sptr->srad+1, sptr->ephi,
				    sptr->sx, sptr->sy );
	    sptr++;
	}

	i = ptr->col + ColourMask;
	fprintf(OutFile,"%g ",(Real)RLut[i]/255.0);
	fprintf(OutFile,"%g ",(Real)GLut[i]/255.0);
	fprintf(OutFile,"%g ",(Real)BLut[i]/255.0);
	fprintf(OutFile,"%g Shade\n",radf);
	fputs("grestore\n\n",OutFile);
    } else
    {   i = ptr->col + ColourMask;
	fprintf(OutFile,"%g ",(Real)RLut[i]/255.0);
	fprintf(OutFile,"%g ",(Real)GLut[i]/255.0);
	fprintf(OutFile,"%g ",(Real)BLut[i]/255.0);
	fprintf(OutFile,"%g %d %d ",radf,ptr->x,ptr->y);
	fputs("Sphere\n\n",OutFile);
    }
}


static void WriteVectWire( src, dst, col, dash )
    Atom __far *src;
    Atom __far *dst;
    int col, dash;
{
    register Atom __far *tmp;
    register Real radius;
    register Real temp;
    register Real dist;

    register Real midx, midy;
    register Real endx, endy;
    register int col1, col2;
    register int dx, dy, dz;
    register Long dist2;
    register int inten;


    if( src->z > dst->z )
    {   tmp = src;
        src = dst;
	dst = tmp;
    }

    if( !col )
    {   col1 = src->col;
	col2 = dst->col;
    } else col1 = col2 = col;

    if( UseBackFade )
    {   dz = (src->z+dst->z)>>1;
        inten = (ColourDepth*(dz+ImageRadius-ZOffset))/ImageSize;
    } else inten = ColourMask;

    dx = dst->x - src->x;  
    dy = dst->y - src->y;
    dist2 = dx*dx + dy*dy;
    dist = sqrt( (double)dist2 );

    if( dst->flag & SphereFlag )
    {   radius = dst->radius*Scale;
	if( dist <= radius ) return;

	/* Test for second half obscured! */
	if( (col1!=col2) && (0.5*dist < radius) )
	    col2 = col1;
    }

    if( src->flag & SphereFlag )
    {   radius = src->radius*Scale;
	if( dist <= radius ) return;

	/* Test for first half obscured! */
	if( (col1!=col2) && (0.5*dist < radius) )
	    col1 = col2;
    }


    WriteVectColour( col1+inten );

    dz = (src->z+dst->z)>>1;
    temp = (double)(dz-ZOffset)/ImageSize + 1.0;
    if( temp != LineWidth )
    {   fprintf(OutFile,"%g setlinewidth\n",temp);
	LineWidth = temp;
    }

    if( dash )
    {   if( VectSolid )
        {   fputs("[3 3] 0 setdash\n",OutFile);
            VectSolid = False;
        }
    } else
        if( !VectSolid )
        {   fputs("[] 0 setdash\n",OutFile);
            VectSolid = True;
        }


    if( src->flag & SphereFlag )
    {   dz = dst->z - src->z;
	dist = sqrt( (double)(dist2 + dz*dz) );
	endx = src->x + (radius*dx)/dist;
	endy = src->y + (radius*dy)/dist;
	fprintf(OutFile,"%g %g ",endx,endy);
    } else
	fprintf(OutFile,"%d %d ",src->x,src->y);

    if( col1 != col2 )
    {   midx = 0.5*(src->x + dst->x);
	midy = 0.5*(src->y + dst->y);
	fprintf(OutFile,"%g %g Wire\n",midx,midy);

	WriteVectColour( col2+inten );
	fprintf(OutFile,"%g %g ",midx,midy);
    } 
    fprintf(OutFile,"%d %d Wire\n",dst->x,dst->y);
}


static void WriteVectStick( src, dst, col, rad )
    Atom __far *src;  
    Atom __far *dst;
    int col, rad;
{
    register Atom __far *tmp;
    register Real midx, midy;
    register Real relx, rely;
    register Real endx, endy;
    register Real radius, angle;
    register Real dist, dist3;
    register Real temp, ratio;

    register Long dist2;
    register int dx, dy, dz;
    register int col1, col2;
    register int i, inten;

    if( !rad )
    {   WriteVectWire(src,dst,col,False);
	return;
    }

    if( src->z > dst->z )
    {   tmp = src;
        src = dst;
	dst = tmp;
    }

    if( !col )
    {   col1 = src->col;
	col2 = dst->col;
    } else col1 = col2 = col;

    dx = dst->x - src->x;  
    dy = dst->y - src->y;
    dz = dst->z - src->z;
    dist2 = dx*dx + dy*dy;
    dist3 = sqrt( (double)(dist2 + dz*dz) );
    dist = sqrt( (double)dist2 );

    if( dst->flag & SphereFlag )
    {   radius = dst->radius*Scale;
	if( dist <= radius ) return;

	/* Test for nearest half obscured! */
	if( (col1!=col2) && (0.5*dist < radius) )
	    col2 = col1;
    }

    if( src->flag & SphereFlag )
    {   radius = src->radius*Scale;
	if( dist <= radius ) return;

	/* Test for furthest half obscured! */
	if( (col1!=col2) && (0.5*dist < radius) )
	    col1 = col2;
    }

    if( !VectSolid )
    {   fputs("[] 0 setdash\n",OutFile);
        VectSolid = True;
    }

    temp = ((src->z-ZOffset)+(dst->z-ZOffset))/ImageSize + 1.0;
    if( temp != LineWidth )
    {   fprintf(OutFile,"%g setlinewidth\n",temp);
	LineWidth = temp;
    }

    radius = rad*Scale;
    angle = Rad2Deg*atan2((double)dy,(double)dx);
    inten = (int)((dist/dist3)*ColourMask);

    if( col1 != col2 )
    {   midx = 0.5*(src->x + dst->x);
	midy = 0.5*(src->y + dst->y);
	relx = (radius*dx)/dist;
	rely = (radius*dy)/dist;

	fprintf(OutFile,"%g %g moveto\n",midx+rely,midy-relx);
	fprintf(OutFile,"%g %g lineto\n",midx-rely,midy+relx);

	ratio = dz/dist3;

	if( (src->flag&SphereFlag) && (src->radius>rad) )
	{   temp = (Scale*src->radius)/dist3;
	    endx = src->x + temp*dx;
	    endy = src->y + temp*dy;

	    fprintf(OutFile,"%g %g %g ",radius,ratio,angle);
	    fprintf(OutFile,"%g %g StickEnd\n",endx,endy);
	} else
	{   fprintf(OutFile,"%d %d %g ",src->x,src->y,radius);
	    fprintf(OutFile,"%g %g arc\n",angle+90,angle-90);
	}
	fputs("closepath ",OutFile);

	i = col1 + inten;
	fprintf(OutFile,"%g ",(Real)RLut[i]/255.0);
	fprintf(OutFile,"%g ",(Real)GLut[i]/255.0);
	fprintf(OutFile,"%g ",(Real)BLut[i]/255.0);
	fputs("setrgbcolor fill\n",OutFile);

	fprintf(OutFile,"%d %d %g ",dst->x,dst->y,radius);
	fprintf(OutFile,"%g %g arc\n",angle-90,angle+90);
	fprintf(OutFile,"%g %g %g ",radius,ratio,angle);
	fprintf(OutFile,"%g %g StickEnd\n",midx,midy);
	fputs("closepath ",OutFile);

	i = col2 + inten;
	fprintf(OutFile,"%g ",(Real)RLut[i]/255.0);
	fprintf(OutFile,"%g ",(Real)GLut[i]/255.0);
	fprintf(OutFile,"%g ",(Real)BLut[i]/255.0);
	fputs("setrgbcolor fill\n",OutFile);

	if( UseOutLine )
	{   fprintf(OutFile,"%d %d %g ",dst->x,dst->y,radius);
	    fprintf(OutFile,"%g %g arc\n",angle-90,angle+90);
	    if( (src->flag&SphereFlag) && (src->radius>rad) )
	    {   fprintf(OutFile,"%g %g %g ",radius,ratio,angle);
		fprintf(OutFile,"%g %g StickEnd\n",endx,endy);
	    } else
	    {   fprintf(OutFile,"%d %d %g ",src->x,src->y,radius);
		fprintf(OutFile,"%g %g arc\n",angle+90,angle-90);
	    }
	    fputs("closepath 0 setgray stroke\n",OutFile);
	}
    } else /* col1 == col2! */
    {   fprintf(OutFile,"%d %d %g ",dst->x,dst->y,radius);
	fprintf(OutFile,"%g %g arc\n",angle-90,angle+90);

	if( (src->flag&SphereFlag) && (src->radius>rad) )
	{   temp = (Scale*src->radius)/dist3;
	    endx = src->x + temp*dx;
	    endy = src->y + temp*dy;
	    ratio = dz/dist3;

	    fprintf(OutFile,"%g %g %g ",radius,ratio,angle);
	    fprintf(OutFile,"%g %g StickEnd\n",endx,endy);
	} else
	{   fprintf(OutFile,"%d %d %g ",src->x,src->y,radius);
	    fprintf(OutFile,"%g %g arc\n",angle+90,angle-90);
	}

	i = col1 + inten;
	fprintf(OutFile,"%g ",(Real)RLut[i]/255.0);
	fprintf(OutFile,"%g ",(Real)GLut[i]/255.0);
	fprintf(OutFile,"%g ",(Real)BLut[i]/255.0);
	fputs("Stick\n",OutFile);
    }
    VectCol = 0;
}


static void WriteVectDots()
{
    register DotStruct __far *ptr;
    register Real x,y,z;
    register Real xi,yi;
    register int inten;
    register int temp;
    register int zi;
    register int i;


    if( LineWidth != 1.0 )
    {   fputs("1 setlinewidth\n",OutFile);
	LineWidth = 1.0;
    }

    temp = SlabValue - ZOffset;
    for( ptr=DotPtr; ptr; ptr=ptr->next )
	for( i=0; i<ptr->count; i++ )
	{   x = ptr->xpos[i];  
	    y = ptr->ypos[i];  
	    z = ptr->zpos[i];


	    xi = (x*MatX[0]+y*MatX[1]+z*MatX[2]) + XOffset;
	    if( (xi<0.0) || (xi>=XRange) ) continue;
	    yi = (x*MatY[0]+y*MatY[1]+z*MatY[2]) + YOffset;
	    if( (yi<0.0) || (yi>=YRange) ) continue;

	    zi = (int)(x*MatZ[0]+y*MatZ[1]+z*MatZ[2]);
	    if( UseSlabPlane && (zi>=temp) ) continue;

	    inten = (ColourDepth*(zi+ImageRadius))/ImageSize;
	    WriteVectColour( ptr->col[i]+inten );
	    fprintf(OutFile,"%g %g Dot\n",xi,yi);
	}
}


static void WriteVectLabels()
{
    register Chain __far *chain;
    register Group __far *group;
    register Atom __far *aptr;
    register Label *label;
    auto char buffer[80];

    fputs("/Times-Roman",OutFile); /* Courier or Courier-Bold? */
    fprintf(OutFile," findfont %d scalefont setfont\n",FontSize<<1);

    if( UseLabelCol )
    {   if( BackR || BackG || BackB )
	{   fprintf(OutFile,"%g %g %g setrgbcolor\n",
		    LabR/250.0, LabG/250.0, LabB/250.0);
	} else fputs("0 setgray\n",OutFile);
    } else VectCol = 0;

    ForEachAtom
	if( aptr->label )
	{   if( !UseLabelCol && (aptr->col!=VectCol) )
		 WriteVectColour( aptr->col );

	    label = (Label*)aptr->label;
	    FormatLabel(chain,group,aptr,label->label,buffer);
	    fprintf(OutFile,"(%s) %d %d Label\n",buffer,aptr->x,aptr->y);
	}
}


static void WriteVectMonitors()
{
    register Atom __far *s;
    register Atom __far *d;
    register Monitor *ptr;
    register int x,y,col;

    register char *cptr;
    register int dist;
    char buffer[10];
 
    buffer[9] = '\0';
    buffer[6] = '.';

    fputs("/Times-Roman",OutFile); /* Courier or Courier-Bold? */
    fprintf(OutFile," findfont %d scalefont setfont\n",FontSize<<1);

    for( ptr=MonitList; ptr; ptr=ptr->next )
    {   s = ptr->src;
        d = ptr->dst;

        if( ZValid( (s->z+d->z)/2 ) )
        {   x = (s->x+d->x)/2;
            y = (s->y+d->y)/2;
 
            if( !UseLabelCol )
            {   /* Use Source atom colour! */
                if( ptr->col )
                {   col = ptr->col + (ColourMask>>1);
                } else col = s->col + (ColourMask>>1);
            } else col = LabelCol;
            WriteVectColour(col);
 
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
 
	    fprintf(OutFile,"(%s) %d %d Label\n",cptr,x+4,y);
        }
    }
}


static Long CountPSItems()
{
    register Chain __far *chain;
    register Group __far *group;
    register HBond __far *hptr;
    register Bond __far *bptr;
    register Atom __far *aptr;
    register Monitor *mptr;
    register Long result;

    result = 0;
    if( DrawAtoms )
        ForEachAtom 
	    if( aptr->flag&SphereFlag ) 
	        if( !UseClipping || !ClipVectSphere(aptr) )
		    result++;

    if( DrawBonds )
        ForEachBond 
            if( bptr->flag&DrawBondFlag && (!UseClipping ||
	        !ClipVectBond(bptr->srcatom,bptr->dstatom)) )
		    result++;

    ForEachBack 
        if( bptr->flag&DrawBondFlag && (!UseClipping ||
	    !ClipVectBond(bptr->srcatom,bptr->dstatom)) )
		result++;

    for( hptr=Database->hlist; hptr; hptr=hptr->hnext )
        if( hptr->flag&DrawBondFlag )
        {   if( HBondMode )
            {   if( !ClipVectBond(hptr->srcCA,hptr->dstCA) )
                    result++;
            } else if( !ClipVectBond(hptr->src,hptr->dst) )
                result++;
        }

    for( hptr=Database->slist; hptr; hptr=hptr->hnext )
        if( hptr->flag&DrawBondFlag )
        {   if( SSBondMode )
            {   if( !ClipVectBond(hptr->srcCA,hptr->dstCA) )
                    result++;
            } else if( !ClipVectBond(hptr->src,hptr->dst) )
                result++;
        }

    for( mptr=MonitList; mptr; mptr=mptr->next )
        if( !UseClipping || !ClipVectBond(mptr->src,mptr->dst) )
            result++;

    return( result );
}


static void FetchPSItems( data, type )
    PSItemPtr __far *data;
    char __far *type;
{
    register Chain __far *chain;
    register Group __far *group;
    register HBond __far *hptr;
    register Bond __far *bptr;
    register Atom __far *aptr;
    register Monitor *mptr;
    register int i,flag;

    i = 0;
    if( DrawAtoms )
        ForEachAtom
	    if( aptr->flag&SphereFlag )
	        if( !UseClipping || !ClipVectSphere(aptr) )
	        {   type[i] = PSAtom; 
		    data[i++] = aptr;
	        }

    if( DrawBonds )
        ForEachBond
            if( bptr->flag&DrawBondFlag && (!UseClipping ||
	        !ClipVectBond(bptr->srcatom,bptr->dstatom)) )
	    {   type[i] = PSBond;
	        data[i++] = bptr;
	    } 

    ForEachBack
       if( bptr->flag&DrawBondFlag && (!UseClipping ||
           !ClipVectBond(bptr->srcatom,bptr->dstatom)) )
       {   type[i] = PSBond;
           data[i++] = bptr; 
       } 

    for( hptr=Database->hlist; hptr; hptr=hptr->hnext )
        if( hptr->flag&DrawBondFlag )
        {   if( HBondMode )
            {   flag = !ClipVectBond(hptr->srcCA,hptr->dstCA);
            } else flag = !ClipVectBond(hptr->src,hptr->dst);

            if( flag )
            {   type[i] = PSHBond;
                data[i++] = hptr;
            }
        }

    for( hptr=Database->slist; hptr; hptr=hptr->hnext )
        if( hptr->flag&DrawBondFlag )
        {   if( SSBondMode )
            {   flag = !ClipVectBond(hptr->srcCA,hptr->dstCA);
            } else flag = !ClipVectBond(hptr->src,hptr->dst);

            if( flag )
            {   type[i] = PSSSBond;
                data[i++] = hptr;
            }
        }

    for( mptr=MonitList; mptr; mptr=mptr->next )
        if( !UseClipping || !ClipVectBond(mptr->src,mptr->dst) )
        {   type[i] = PSMonit;
            data[i++] = mptr;
        } 
}


static void WritePSItems( data, type, count )
    PSItemPtr __far *data;
    char __far *type;
    int count;
{
    register HBond __far *hbond;
    register Bond __far *bond;
    register Atom __far *src;
    register Atom __far *dst;
    register Monitor *monit;
    register int i;

    for( i=0; i<count; i++ )
        switch( type[i] )
        {   case(PSAtom):   WriteVectSphere(data,type,i);
                            break;

            case(PSBond):   bond = (Bond __far*)data[i];
                            src = bond->srcatom;
                            dst = bond->dstatom;

                            if( bond->flag & WireFlag )
                            {   WriteVectWire(src,dst,bond->col,False);
                            } else if( bond->flag & CylinderFlag )
                            {   WriteVectStick(src,dst,bond->col,bond->radius);
                            } else /* bond->flag & DashFlag */
                                WriteVectWire(src,dst,bond->col,True);
                            break;

            case(PSSSBond): 
            case(PSHBond):  hbond = (HBond __far*)data[i];
                            if( (type[i]==PSHBond)? HBondMode : SSBondMode )
                            {   src = hbond->srcCA;
                                dst = hbond->dstCA;
                            } else
                            {   src = hbond->src;
                                dst = hbond->dst;
                            }

                            if( hbond->flag & WireFlag )
                            {   WriteVectWire(src,dst,hbond->col,True);
                            } else /* bond->flag & CylinderFlag */
                                WriteVectStick(src,dst,hbond->col,
                                                       hbond->radius);
                            break;

            case(PSMonit):  monit = (Monitor*)data[i];
                            WriteVectWire(monit->src,monit->dst,
                                          monit->col,True);
                            break;
        }
}



int WriteVectPSFile( name )
    char *name;
{
    register Real ambi;
    register Real temp, inten;
    register int xsize, ysize;
    register int xpos, ypos;
    register Long count;
    register int i;

    PSItemPtr __far *data;
    char __far *type;

    count = CountPSItems();
    if( !count ) return( True );

#ifdef IBMPC
    if( count > 16383 )
    {   if( CommandActive ) WriteChar('\n');
        WriteString("Output Error: Too many PostScript objects!\n");
        CommandActive = False;
	return( False );
    }
#endif

    /* Allocate arrays for objects! */
    data = (PSItemPtr __far*)_fmalloc((size_t)count*sizeof(PSItemPtr));
    type = (char __far*)_fmalloc((size_t)count*sizeof(char));
    if( !data || !type )
    {   if( CommandActive ) WriteChar('\n');
	WriteString("Output Error: Not enough memory to create PostScript!\n");
	CommandActive = False;

	if( data ) _ffree( data );
	if( type ) _ffree( type );
	return( False );
    }

    OutFile = fopen(name,"w");
    if( !OutFile )
    {   FatalOutputError(name);
	return(False);
    }

    /* Determine the size of the image */
    ysize = (int)(YRange*(BORDER*PAGEWIDE)/XRange);
    if( ysize > (int)(BORDER*PAGEHIGH) )
    {   xsize = (int)(XRange*(BORDER*PAGEHIGH)/YRange);
	ysize = (int)(BORDER*PAGEHIGH);
    } else xsize = (int)(BORDER*PAGEWIDE);

    xpos = (int)(PAGEWIDE-xsize)/2;
    ypos = (int)(PAGEHIGH-ysize)/2;

    fputs("%!PS-Adobe-2.0 EPSF-2.0\n",OutFile);
    fputs("%%Creator: RasMol Version 2.6\n",OutFile);
    fprintf(OutFile,"%%%%Title: %s\n",name);
    fprintf(OutFile,"%%%%BoundingBox: %d %d ",xpos,ypos);
    fprintf(OutFile,"%d %d\n",xpos+xsize,ypos+ysize);

    fputs("%%Pages: 1\n",OutFile);
    fputs("%%EndComments\n",OutFile);
    fputs("%%EndProlog\n",OutFile);
    fputs("%%BeginSetup\n",OutFile);

    fputs("1 setlinecap 1 setlinejoin [] 0 setdash\n",OutFile);
    fputs("1 setlinewidth 0 setgray\n",OutFile);
    fputs("%%EndSetup\n",OutFile);
    fputs("%%Page: 1 1\n",OutFile);

    fputs("gsave\n",OutFile);
    fputs("14 dict begin\n\n",OutFile);
    fputs("/handleerror { showpage } def\n\n",OutFile);
    fputs("/Inten {\n  dup 4 index mul exch\n",OutFile);
    fputs("  dup 4 index mul exch\n",OutFile);
    fputs("  3 index mul setrgbcolor\n} def\n\n",OutFile);

    fputs("/Dot {\n  moveto 0 0 rlineto stroke\n} def\n\n",OutFile);
    fputs("/Wire {\n  moveto lineto stroke\n} def\n\n",OutFile);
#ifdef INVERT
    fputs("/Label {\n  moveto show\n} def\n\n",OutFile);
#else
    fputs("/Label {\n  moveto 1 -1 scale\n",OutFile);
    fputs("  show mtrx setmatrix\n} def\n\n",OutFile);
#endif

    if( UseOutLine )
    {   fputs("/Stick {\n  closepath gsave setrgbcolor fill\n",OutFile);
	fputs("  grestore 0 setgray stroke\n} def\n\n",OutFile);
    } else
	fputs("/Stick {\n  closepath setrgbcolor fill\n} def\n\n",OutFile);

    fputs("/StickEnd {\n  matrix currentmatrix 6 1 roll\n",OutFile);
    fputs("  translate rotate 1 scale\n",OutFile);
    fputs("  0 0 3 2 roll 90 -90 arc\n  setmatrix\n} def\n\n",OutFile);

    if( UseOutLine )
    {   fputs("/Shade {\n  closepath gsave clip\n",OutFile);
    } else fputs("/Shade {\n  closepath clip\n",OutFile);

    if( Ambient < 0.99 )
    {   ambi = 0.5*Ambient;
        fputs("  45 rotate dup -0.81649658092 mul scale\n",OutFile);
        fprintf(OutFile,"  %g Inten fill\n",ambi);
        inten = (1.0-ambi)/31;
        for( i=0; i<31; i++ )
        {   temp = (Real)(i+1)/32;
	    fprintf(OutFile,"  0 %g ",(Real)i/32);
	    fprintf(OutFile,"%g 0 360 arc ",sqrt(1.0-temp*temp));
	    fprintf(OutFile,"%g Inten fill\n",i*inten+ambi);
        }
        if( UseOutLine )
        {   fputs("  grestore 0 setgray stroke",OutFile);
        } else fputc(' ',OutFile);
        fputs(" pop pop pop\n} def\n\n",OutFile);

    } else /* Saturated Colours! */
    {   fputs("  pop setrgbcolor fill\n",OutFile);
        if( UseOutLine )
            fputs("  grestore 0 setgray stroke\n",OutFile);
        fputs("} def\n\n",OutFile);
    }


    fputs("/ClipSphere {\n  translate 0 0 5 2 roll arc\n} def\n\n",OutFile);
    fputs("/ClipBox {\n  translate rotate\n  dup lineto dup neg ",OutFile);
    fputs("dup\n  0 rlineto 0 exch rlineto 0 rlineto closepath\n",OutFile);
    fputs("  clip newpath mtrx setmatrix\n} def\n\n",OutFile);
    fputs("/ClipEllips {\n  translate rotate 1 scale\n",OutFile);
    fputs("  0 0 4 2 roll dup neg arc\n",OutFile);
    fputs("  reversepath mtrx setmatrix\n} def\n\n",OutFile);

    fputs("/Sphere {\n  gsave\n",OutFile);
    fputs("  translate 0 0 2 index 0 360 arc\n",OutFile);
    if( UseOutLine )
    {   fputs("  gsave Shade grestore\n",OutFile);
	fputs("  0 setgray stroke\n",OutFile);
	fputs("  grestore\n} def\n\n",OutFile);
    } else
	fputs("  Shade grestore\n} def\n\n",OutFile);

#ifdef INVERT
    fprintf(OutFile,"%d %d translate\n",xpos,ypos);
    fprintf(OutFile,"%g ",(Real)xsize/XRange);
    fprintf(OutFile,"%g ",(Real)ysize/YRange);
#else
    fprintf(OutFile,"%d %d translate\n",xpos,ypos+ysize);
    fprintf(OutFile,"%g ",(Real)xsize/XRange);
    fprintf(OutFile,"%g ",(Real)-ysize/YRange);
#endif
    fputs("scale\n/mtrx matrix currentmatrix def\n\n",OutFile);

    fputs("newpath 0 0 moveto 0 ",OutFile);
    fprintf(OutFile,"%d rlineto %d 0 rlineto 0 %d",YRange,XRange,-YRange);
    fputs(" rlineto\nclosepath clip ",OutFile);
    if( BackR || BackG || BackB )
    {   fprintf(OutFile,"%g %g %g",BackR/255.0,BackG/255.0,BackB/255.0);
	fputs(" setrgbcolor fill\n\n",OutFile);
    } else fputs("newpath\n\n",OutFile);

    LineWidth = 1.0;
    VectSolid = True;
    VectCol = 0;

    FetchPSItems(data,type);
    if( count>1 )
	DepthSort(data,type,(int)count);

    WritePSItems(data,type,(int)count);
 
    if( !VectSolid )
    {   fputs("[] 0 setdash\n",OutFile);
        VectSolid = True;
    }

    if( DrawDots )
	WriteVectDots();
    if( DrawMonitDistance && MonitList )
        WriteVectMonitors();
    if( DrawLabels )
	WriteVectLabels();

    fputs("newpath 0 0 moveto 0 ",OutFile);
    fprintf(OutFile,"%d rlineto %d 0 rlineto 0 %d",YRange,XRange,-YRange);
    fputs(" rlineto\nclosepath 0 setgray 1 setlinewidth stroke\n",OutFile);
    fputs("end grestore\nshowpage\n",OutFile);
    fputs("%%Trailer\n",OutFile);
    fputs("%%EOF\n",OutFile);

    fclose( OutFile );
#ifdef APPLEMAC
    /* Avoid ANSI trigraph problems! */
    SetFileInfo(name,'vgrd','TEXT',134);
#endif
    _ffree( data );
    _ffree( type );
    return(True);
}



void InitialiseOutFile()
{

    UseTransparent = False;
    KinemageFlag = False;
    UseOutLine = False;
}
