/* tokens.h
 * RasMol2 Molecular Graphics
 * Roger Sayle, August 1995
 * Version 2.6
 */

/* Lexeme Tokens */
#define IdentTok       256
#define NumberTok      257
#define FloatTok       258
#define StringTok      259

/* Command Tokens */
#define BackboneTok    260
#define CartoonTok     261
#define CentreTok      262
#define ClipboardTok   263
#define ColourTok      264
#define ConnectTok     265
#define DashTok        266
#define DefineTok      267
#define DisplayTok     268
#define EchoTok        269
#define ExitTok        270
#define HelpTok        271
#define LabelTok       272
#define LoadTok        273
#define MonitorTok     274
#define PrintTok       275
#define QuitTok        276
#define RefreshTok     277
#define RenumTok       278
#define ResetTok       279
#define ResizeTok      280
#define RestrictTok    281
#define RotateTok      282
#define SaveTok        283
#define ScriptTok      284
#define SelectTok      285
#define SetTok         286
#define ShowTok        287
#define SlabTok        288
#define SourceTok      289
#define SpacefillTok   290
#define StructureTok   291
#define SymmetryTok    292
#define TraceTok       293
#define TranslateTok   294
#define WaitTok        295
#define WireframeTok   296
#define WriteTok       297
#define ZapTok         298
#define ZoomTok        299

/* Predicate Tokens */
#define IsPredTok(x)   (((x)>=300) && ((x)<=338))
#define PredTokOrd(x)  ((x)-300)
#define PredTokChr(x)  ((x)+300)

#define AlphaTok       300
#define AminoTok       301
#define ATTok          302
#define BondedTok      303
#define CGTok          304
#define CystineTok     305
#define DNATok         306
#define HelixTok       307
#define HeteroTok      308
#define HydrogenTok    309
#define IonTok         310
#define LigandTok      311
#define MainChainTok   312
#define NucleicTok     313
#define ProteinTok     314
#define PurineTok      315
#define PyrimidineTok  316
#define RNATok         317
#define SelectedTok    318
#define SheetTok       319
#define SidechainTok   320
#define SolventTok     321
#define TurnTok        322
#define WaterTok       323

#define AcidicTok      324
#define AcyclicTok     325
#define AliphaticTok   326
#define AromaticTok    327
#define BasicTok       328
#define BuriedTok      329
#define ChargedTok     330
#define CyclicTok      331
#define HydrophobicTok 332
#define LargeTok       333
#define MediumTok      334
#define NeutralTok     335
#define PolarTok       336
#define SmallTok       337
#define SurfaceTok     338


/* Property Tokens */
#define IsPropTok(x)   (((x)>=340) && ((x)<=345))
#define TemperatureTok 340
#define RadiusTok      341
#define AtomNoTok      342
#define ElemNoTok      343
#define ModelTok       344
#define ResNoTok       345

/* File Format Tokens */
/* Warning! Tokens are related to Format values */
#define IsMoleculeFormat(x)  (((x)>=350) && ((x)<=365))

#define PDBTok         350
#define NMRPDBTok      354
#define CIFTok         364

/* Raster Tokens */
#define IsImageFormat(x) (((x)>=366) && ((x)<=378))
#define PPMTok         367
#define VectPSTok      375

/* Feature Tokens */
#define AtomTok        380
#define BondTok        381
#define DotsTok        382
#define HBondTok       383
#define RibbonTok      384
#define SSBondTok      385
#define Ribbon1Tok     386
#define Ribbon2Tok     387

/* Expression Tokens */
#define TrueTok        390
#define FalseTok       391
#define AllTok         392
#define NoneTok        393
#define AndTok         394
#define OrTok          395
#define NotTok         396
#define WithinTok      397
#define XorTok         398

/* Colour Tokens */
#define BlueTok        400
#define BlueTintTok    401
#define BlackTok       402
#define BrownTok       403
#define CyanTok        404
#define GoldTok        405
#define GrayTok        406
#define GreenTok       407
#define GreenblueTok   408
#define GreenTintTok   409
#define HotPinkTok     410
#define MagentaTok     411
#define OrangeTok      412
#define PinkTok        413
#define PinkTintTok    414
#define PurpleTok      415
#define RedTok         416
#define RedorangeTok   417
#define SeaTok         418
#define SkyTok         419
#define VioletTok      420
#define WhiteTok       421
#define YellowTok      422
#define YellowTintTok  423

#define CPKTok         424
#define ShapelyTok     425
#define UserTok        426
#define GroupTok       427
#define ChainTok       428
#define TypeTok        429
#define PotentialTok   430
#define ChargeTok      431

/* Variable Tokens */
#define AmbientTok     440
#define AxesTok        441
#define BackFadeTok    442
#define BackgroundTok  443
#define BondModeTok    444
#define BoundBoxTok    445
#define DepthCueTok    446
#define FontSizeTok    447
#define HourGlassTok   448
#define MenusTok       449
#define MouseTok       450
#define PickingTok     451
#define ShadowTok      452
#define SlabModeTok    453
#define SpecularTok    454
#define SpecPowerTok   455
#define StrandsTok     456
#define TransparentTok 457
#define UnitCellTok    458

/* SlabMode Tokens */
#define RejectTok      460
#define HalfTok        461
#define HollowTok      462
#define SolidTok       463
#define SectionTok     464

/* MouseMode Tokens */
#define RasMolTok      465
#define InsightTok     466
#define QuantaTok      467

/* Information Tokens */
#define InfoTok        470
#define SequenceTok    471
#define VersionTok     472

/* Display Mode Tokens */
#define NormalTok      475
#define StereoTok      476
#define MonoTok        477

/* Axis Tokens */
#define XTok           480
#define YTok           481
#define ZTok           482

/* Picking Tokens */
#define IdentifyTok    485
#define DistanceTok    486
#define AngleTok       487
#define TorsionTok     488

/* Misc Tokens */
#define InLineTok      490


typedef struct {
                char *ident;
                int token;
               } KeywordEntry;


#define MAXKEYLEN 11
static int KeyLen[MAXKEYLEN+1] = {
        0, 3, 7, 21, 48, 82, 128, 160, 186, 201, 204, 210 };

static KeywordEntry Keyword[] = {
            { "X",  XTok },
            { "Y",  YTok },
            { "Z",  ZTok },  /* 3 */

            { "AT", ATTok   },
            { "CG", CGTok   },
            { "ON", TrueTok },
            { "OR", OrTok   },

            { "ALL", AllTok   },
            { "AND", AndTok   },
            { "CIF", CIFTok   },
            { "CPK", CPKTok   },
            { "DNA", DNATok   },
            { "ION", IonTok   },
            { "NOT", NotTok   },
            { "OFF", FalseTok },
            { "PDB", PDBTok   },
            { "PPM", PPMTok   },
            { "RED", RedTok   },
            { "RNA", RNATok   },
            { "SET", SetTok   },
            { "ZAP", ZapTok   }, /* 28 */

            { "ATOM", AtomTok },
            { "AXES", AxesTok },
            { "AXIS", AxesTok },
            { "BLUE", BlueTok },
            { "BOND", BondTok },
            { "CYAN", CyanTok },
            { "DASH", DashTok },
            { "DOTS", DotsTok },
            { "ECHO", EchoTok },
            { "EXIT", ExitTok },
            { "HALF", HalfTok },
            { "HELP", HelpTok },
            { "INFO", InfoTok },
            { "IONS", IonTok  },
            { "LOAD", LoadTok },
            { "MONO", MonoTok },
            { "NONE", NoneTok },
            { "QUIT", QuitTok },
            { "SAVE", SaveTok },
            { "SHOW", ShowTok },
            { "SLAB", SlabTok },
            { "TRUE", TrueTok },
            { "TURN", TurnTok },
            { "TYPE", TypeTok },
            { "USER", UserTok },
            { "WAIT", WaitTok },
            { "ZOOM", ZoomTok }, /* 61 */

            { "ALPHA", AlphaTok    },
            { "AMINO", AminoTok    },
            { "ANGLE", AngleTok    },
            { "ATOMS", AtomTok     },
            { "BASIC", BasicTok    },
            { "BLACK", BlackTok    },
            { "BONDS", BondTok     },
            { "CHAIN", ChainTok    },
            { "COLOR", ColourTok   },
            { "FALSE", FalseTok    },
            { "GREEN", GreenTok    },
            { "GROUP", GroupTok    },
            { "HBOND", HBondTok    },
            { "HELIX", HelixTok    },
            { "IDENT", IdentifyTok },
            { "LABEL", LabelTok    },
            { "LARGE", LargeTok    },
            { "MENUS", MenusTok    },
            { "MODEL", ModelTok    },
            { "MOUSE", MouseTok    },
            { "PAUSE", WaitTok     },
            { "POLAR", PolarTok    },
            { "PRINT", PrintTok    },
            { "RENUM", RenumTok    },
            { "RESET", ResetTok    },
            { "RESNO", ResNoTok    },
            { "SHEET", SheetTok    },
            { "SMALL", SmallTok    },
            { "SOLID", SolidTok    },
            { "TRACE", TraceTok    },
            { "TURNS", TurnTok     },
            { "WATER", WaterTok    },
            { "WHITE", WhiteTok    },
            { "WRITE", WriteTok    },  /* 97 */

            { "ACIDIC", AcidicTok },
            { "ANGLES", AngleTok  },
            { "ATOMNO", AtomNoTok },
            { "BONDED", BondedTok },
            { "BURIED", BuriedTok },
            { "CENTER", CentreTok },
            { "CENTRE", CentreTok },
            { "CHARGE", ChargeTok },
            { "COLORS", ColourTok },
            { "COLOUR", ColourTok },
            { "CYCLIC", CyclicTok },
            { "DASHES", DashTok   },
            { "DEFINE", DefineTok },
            { "ELEMNO", ElemNoTok },
            { "HBONDS", HBondTok  },
            { "HETERO", HeteroTok },
            { "HOLLOW", HollowTok },
            { "INLINE", InLineTok },
            { "LABELS", LabelTok  },
            { "LIGAND", LigandTok },
            { "MEDIUM", MediumTok },
            { "NMRPDB", NMRPDBTok },
            { "NORMAL", NormalTok },
            { "ORANGE", OrangeTok },
            { "PURINE", PurineTok },
            { "PURPLE", PurpleTok },
            { "QUANTA", QuantaTok },
            { "RADIUS", RadiusTok },
            { "RASMOL", RasMolTok },
            { "RASWIN", RasMolTok },
            { "REJECT", RejectTok },
            { "RESIZE", ResizeTok },
            { "RIBBON", RibbonTok },
            { "ROTATE", RotateTok },
            { "SCRIPT", ScriptTok },
            { "SELECT", SelectTok },
            { "SHADOW", ShadowTok },
            { "SHEETS", SheetTok  },
            { "SOURCE", SourceTok },
            { "SSBOND", SSBondTok },
            { "STEREO", StereoTok },
            { "VECTPS", VectPSTok },
            { "VIOLET", VioletTok },
            { "WATERS", WaterTok  },
            { "WITHIN", WithinTok },
            { "YELLOW", YellowTok },  /* 148 */

            { "ACYCLIC", AcyclicTok },
            { "AMBIENT", AmbientTok },
            { "CARTOON", CartoonTok },
            { "CHARGED", ChargedTok },
            { "CHARGES", ChargeTok  },
            { "COLOURS", ColourTok  },
            { "CONNECT", ConnectTok },
            { "CYSTINE", CystineTok },
            { "DISPLAY", DisplayTok },
            { "HELICES", HelixTok   },
            { "INSIGHT", InsightTok },
            { "LIGANDS", LigandTok  },
            { "MAGENTA", MagentaTok },
            { "MONITOR", MonitorTok },
            { "NEUTRAL", NeutralTok },
            { "NUCLEIC", NucleicTok },
            { "PICKING", PickingTok },
            { "PROTEIN", ProteinTok },
            { "PURINES", PurineTok  },
            { "REFRESH", RefreshTok },
            { "RESIDUE", GroupTok   },
            { "RIBBON1", Ribbon1Tok },
            { "RIBBON2", Ribbon2Tok },
            { "RIBBONS", RibbonTok  },
            { "SECTION", SectionTok },
            { "SHADOWS", ShadowTok  },
            { "SHAPELY", ShapelyTok },
            { "SOLVENT", SolventTok },
            { "SSBONDS", SSBondTok  },
            { "STRANDS", StrandsTok },
            { "SURFACE", SurfaceTok },  
            { "TORSION", TorsionTok }, /* 181 */

            { "AROMATIC", AromaticTok },
            { "BACKBONE", BackboneTok },
            { "BACKFADE", BackFadeTok },
            { "BONDMODE", BondModeTok },
            { "BOUNDBOX", BoundBoxTok },
            { "CARTOONS", CartoonTok  },
            { "DEPTHCUE", DepthCueTok },
            { "DISTANCE", DistanceTok },
            { "FONTSIZE", FontSizeTok },
            { "HYDROGEN", HydrogenTok },
            { "IDENTIFY", IdentifyTok },
            { "MONITORS", MonitorTok  },
            { "NEGATIVE", AcidicTok   },
            { "POSITIVE", BasicTok    },
            { "RENUMBER", RenumTok    },
            { "RESTRICT", RestrictTok },
            { "RIBBONS1", Ribbon1Tok  },
            { "RIBBONS2", Ribbon2Tok  },
            { "SELECTED", SelectedTok },
            { "SEQUENCE", SequenceTok },
            { "SLABMODE", SlabModeTok },
            { "SOLVENTS", SolventTok  },
            { "SPECULAR", SpecularTok }, 
            { "SYMMETRY", SymmetryTok },
            { "TORSIONS", TorsionTok  },
            { "UNITCELL", UnitCellTok },  /* 209 */

            { "ALIPHATIC", AliphaticTok },
            { "CLIPBOARD", ClipboardTok },
            { "DISTANCES", DistanceTok  },
            { "GREENBLUE", GreenblueTok },
            { "HOURGLASS", HourGlassTok },
            { "MAINCHAIN", MainChainTok },
            { "MOUSEMODE", MouseTok     },
            { "POTENTIAL", PotentialTok },
            { "REDORANGE", RedorangeTok },
            { "SIDECHAIN", SidechainTok },
            { "SPACEFILL", SpacefillTok },
            { "SPECPOWER", SpecPowerTok },
            { "STRUCTURE", StructureTok },
            { "TRANSLATE", TranslateTok },
            { "WIREFRAME", WireframeTok },  /* 225 */

            { "BACKGROUND", BackgroundTok },
            { "MONOCHROME", MonoTok       },
            { "PYRIMIDINE", PyrimidineTok },  /* 229 */

            { "BOUNDINGBOX", BoundBoxTok    },
            { "HYDROPHOBIC", HydrophobicTok },
            { "INFORMATION", InfoTok        },
            { "PYRIMIDINES", PyrimidineTok, },
            { "TEMPERATURE", TemperatureTok },
            { "TRANSPARENT", TransparentTok }  /* 236 */
                };

