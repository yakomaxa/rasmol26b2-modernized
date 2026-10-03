/* cexio.c
 * RasMol2 Molecular Graphics
 * Roger Sayle, August 1995
 * Version 2.6
 */
#include <stdio.h>
#include "rasmol.h"
#include "molecule.h"

#include "cx.h"
#include "cx_molecule.h"

#define PROP_ATSYM   "atomic symbol"
#define PROP_BONDO   "bond order"

#define PROP_MNAME   "molname"
#define PROP_ALABS   "atomname"
#define PROP_RESID   "resname"
#define PROP_CHAIN   "chain"
#define PROP_RESNO   "resno"
#define PROP_COORD   "coordinates"
#define PROP_BVALU   "bvalue"
#define PROP_CHARG   "charge"

static cx_String prop_alabs;
static cx_String prop_resid;
static cx_String prop_chain;
static cx_String prop_resno;
static cx_String prop_coord;
static cx_String prop_tempr;


static cx_String get_atom_tuple( mol, prop )
    cx_Object mol;  cx_String prop;
{
    cx_Object tupleseq, tuple;
    cx_String result;
    register int i,j;

    if( (tupleseq = cx_prefix2atuples(mol,prop)) )
    {   tuple = cx_next(tupleseq);
    } else tuple = (cx_Object)NULL;
 
    if( !tuple )
    {   /* No such property! */
        cx_destroy( tupleseq );
        return( (cx_String)NULL );
    }
 
    /* Ambiguous atom property prefix! */
    /* cx_next(tupleseq) or !cx_atend(tupleseq) */

    cx_destroy(tupleseq);
    return( cx_atomtuple_name(tuple) );
}
 

static int process_cx_atom( atom, serno )
    cx_Object atom;  int serno;
{
    register Group __far *group;
    register Atom __far *ptr;
    register int refno,resno;
    register int elem,chain;
    cx_String prop;
    float x,y,z,t;

    /* Test for valid 3D co-ordinates */
    if( !(prop=cx_sprop(atom,prop_coord)) )
        return( False );
    if( sscanf(prop,"%g,%g,%g",&x,&y,&z) != 3 )
        return( False );

    /* Determine Chain */
    prop = cx_sprop(atom,prop_chain);
    if( prop && *prop ) 
    {   chain = *prop;
    } else chain = ' ';

    if( CurChain && (CurChain->ident!=chain) )
    {   CurChain = CurMolecule->clist;
        while( CurChain )
            if( CurChain->ident != chain )
            {   CurChain = CurChain->cnext;
            } else break;
    }
    if( !CurChain ) CreateChain( chain );

    /* Determine Residue */
    resno = cx_iprop(atom,prop_resno);
    prop = cx_sprop(atom,prop_resid);
    if( !prop || !*prop ) 
    {   refno = FindResNo("MOL");
    } else refno = FindResNo(prop);

    if( !CurGroup || CurGroup->refno!=refno || CurGroup->serno!=resno )
    {   group = CurChain->glist;
        CurGroup = (Group __far*)0;
        CurAtom = (Atom __far*)0;

        while( group && group->serno <= resno )
        {   CurGroup = group;
            if( (group->serno==resno) && (group->refno==refno) ) break;
            group = group->gnext;
        }

        if( !group || (group->serno>resno) ) 
        {   CreateGroup( 8 );
            CurGroup->serno = resno;
            CurGroup->refno = refno;
            ProcessGroup( False );
        }
    }

    /* Process Atom */
    ptr = CreateAtom();
    ptr->serno = serno;

    if( !(prop=cx_sprop(atom,prop_alabs)) )
    {   prop = cx_sprop(atom,PROP_ATSYM);
        if( !prop ) prop = "Du  ";
        ptr->refno = SimpleAtomType(prop);
    } else ptr->refno = ComplexAtomType(prop);

    t = cx_rprop(atom,prop_tempr);
    ptr->temp = (short)(t*100.0);

    ptr->xorg =  (Long)(x*250.0);
    ptr->yorg =  (Long)(y*250.0);
    ptr->zorg = -(Long)(z*250.0);
    ProcessAtom( ptr );
    return( True );
}


static int process_cx_bond( bond )
    cx_Object bond;
{
    cx_Object atoms,atom;
    register int src,dst;
    register int i,flag;

    i = cx_iprop(bond,PROP_BONDO);
    if( i==2 )
    {   flag = DoubBondFlag;
    } else if( i == 3 )
    {   flag = TripBondFlag;
    } else flag = NormBondFlag;

    atoms = cx_stream(bond,CX_OB_ATOM);
    if( atom = cx_next(atoms) )
        if( src = cx_iprop(atom,"_serno") )
            if( atom = cx_next(atoms) )
                if( dst = cx_iprop(atom,"_serno") )
                    CreateBond(src,dst,flag);
    cx_destroy(atoms);
}


static void load_cx_molecule( mol )
    cx_Object *mol;
{
    cx_String molname;
    cx_Object atoms,atom;
    cx_Object bonds,bond;
    register int serno;
    register int i;

    if( (molname = cx_sprop(mol,PROP_MNAME)) )
    {   for( i=0; molname[i] && (i<62); i++ )
            InfoMoleculeName[i] = molname[i];
        while( i && molname[i-1]==' ' ) i--;
        InfoMoleculeName[i] = '\0';
    }

    /* Molecule must have 3d co-ordinates */
    prop_coord = get_atom_tuple(mol,PROP_COORD);
    if( !prop_coord ) return;

    prop_alabs = get_atom_tuple(mol,PROP_ALABS);
    prop_resid = get_atom_tuple(mol,PROP_RESID);
    if( !(prop_tempr=get_atom_tuple(mol,PROP_BVALU)) )
        prop_tempr = get_atom_tuple(mol,PROP_CHARG);

    prop_chain = get_atom_tuple(mol,PROP_CHAIN);
    prop_resno = get_atom_tuple(mol,PROP_RESNO);

    serno = 1;
    atoms = cx_stream(mol,CX_OB_ATOM);
    while( atom = cx_next(atoms) )
    {   if( process_cx_atom(atom,serno) )
        {   cx_set_iprop(atom,"_serno",serno++);
        } else cx_set_iprop(atom,"_serno",0);
    }
    cx_destroy(atoms);

    bonds = cx_stream(mol,CX_OB_BOND);
    while( bond = cx_next(bonds) )
        process_cx_bond(bond);
    cx_destroy(bonds);
}


int LoadCEXMolecule( fp )
    FILE *fp;
{
    register cx_Object obj;

    cx_molecule_pkg();
    while( (obj = cx_receive(NULL,fp,NULL)) ) 
    {   if( cx_type(obj) == CX_OB_MOLECULE )
        {   load_cx_molecule(obj);
            cx_destroy(obj);
            break;
        } else cx_destroy(obj);
    }
    /* cx_cleanup(); */
    return( True );
}
 
 
int SaveCEXMolecule( filename )
    char *filename;
{
    if( !Database )
        return( False );
    return( True );
}
 

