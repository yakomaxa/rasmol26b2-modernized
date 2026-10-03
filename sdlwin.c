/* sdlwin.c
 * RasMol2 Molecular Graphics
 * Minimal SDL2 display backend
 *
 * Replaces x11win.c for modern desktops (Linux and macOS without
 * XQuartz). Implements the same graphics.h entry points that the
 * core renderer calls, but blits the existing software framebuffer
 * straight to an SDL texture instead of drawing a hand-rolled X11
 * menu bar/scrollbar/dials-box GUI. All molecule display, mouse
 * rotate/translate/zoom/slab and atom picking behave the same as
 * before; menus and the hardware dials box are gone in favour of
 * RasMol's command language (typed at the terminal or in scripts),
 * which was always the more complete interface anyway.
 */
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>

#include <SDL.h>

#define GRAPHICS
#include "rasmol.h"
#include "graphics.h"
#include "command.h"

extern int ProcessCommand();

void TransferImage();

static SDL_Window *Window;
static SDL_Renderer *Renderer;
static SDL_Texture *Texture;

static SDL_Cursor *ArrowCursor;
static SDL_Cursor *CrossCursor;
static SDL_Cursor *WaitCursorPtr;

static int InitX, InitY;

/* Determine Mouse Sensitivity! */
#define IsClose(u,v) (((u)>=(v)-1) && ((u)<=(v)+1))


static void FatalGraphicsError( ptr )
    char *ptr;
{
    char buffer[80];

    sprintf(buffer,"Graphics Error: %s!",ptr);
    RasMolFatalExit(buffer);
}


void AllocateColourMap()
{
    /* Lut[] already holds direct 0x00RRGGBB pixel values computed by
     * SetLutEntry() in transfor.c, which exactly matches the layout
     * of SDL_PIXELFORMAT_RGB888. Nothing further to allocate.
     */
}


static void ClampDial( dial, value )
    int dial;  Real value;
{
    register Real temp;

    temp = DialValue[dial] + value;

    if( temp > 1.0 )
    {   DialValue[dial] = 1.0;
    } else if( temp < -1.0 )
    {   DialValue[dial] = -1.0;
    } else DialValue[dial] = temp;
}


static void WrapDial( dial, value )
    int dial;  Real value;
{
    register Real temp;

    temp = DialValue[dial] + value;
    while( temp < -1.0 )  temp += 2.0;
    while( temp > 1.0 )   temp -= 2.0;
    DialValue[dial] = temp;
}


/* Mouse-to-viewpoint mapping, ported from the original X11 backend's
 * "MMRasMol" mode (the default/only mode now that the Insight/Quanta
 * menu toggle is gone):
 *   Left drag                -> rotate about X/Y
 *   Middle/Right drag         -> translate X/Y
 *   Shift + Left drag         -> zoom
 *   Shift + Middle/Right drag -> rotate about Z
 *   Ctrl  + Left drag         -> slab
 */
static void MouseMove( buttons, mods, dx, dy )
    Uint32 buttons; Uint16 mods; int dx, dy;
{
    if( mods & KMOD_SHIFT )
    {   if( buttons & SDL_BUTTON(SDL_BUTTON_LEFT) )
        {   if( dy )
            {   ClampDial( 3, (Real)dy/HRange );
                ReDrawFlag |= RFZoom;
            }
        } else if( buttons & (SDL_BUTTON(SDL_BUTTON_MIDDLE)|SDL_BUTTON(SDL_BUTTON_RIGHT)) )
            if( dx )
            {   WrapDial( 2, (Real)dx/WRange );
                ReDrawFlag |= RFRotateZ;
            }
    } else if( mods & KMOD_CTRL )
    {   if( buttons & SDL_BUTTON(SDL_BUTTON_LEFT) )
            if( dy )
            {   ClampDial( 7, (Real)dy/YRange );
                ReDrawFlag |= RFSlab;
            }
    } else
    {   if( buttons & SDL_BUTTON(SDL_BUTTON_LEFT) )
        {   if( dx )
            {   WrapDial( 1, (Real)dx/WRange );
                ReDrawFlag |= RFRotateY;
            }
            if( dy )
            {   WrapDial( 0, (Real)dy/HRange );
                ReDrawFlag |= RFRotateX;
            }
        } else if( buttons & (SDL_BUTTON(SDL_BUTTON_MIDDLE)|SDL_BUTTON(SDL_BUTTON_RIGHT)) )
        {   if( dx )
            {   ClampDial( 4, (Real)dx/XRange );
                ReDrawFlag |= RFTransX;
            }
            if( dy )
            {   ClampDial( 5, (Real)dy/YRange );
                ReDrawFlag |= RFTransY;
            }
        }
    }
}


static void ResizeWindow( wide, high )
    int wide, high;
{
    if( wide < 40 ) wide = 40;
    if( high < 40 ) high = 40;
    if( (wide==XRange) && (high==YRange) ) return;

    XRange = wide;   WRange = XRange>>1;
    YRange = high;   HRange = YRange>>1;
    Range = MinFun(XRange,YRange);

    ReDrawFlag |= RFReSize;
}


static void ProcessEvent( event )
    SDL_Event *event;
{
    switch( event->type )
    {   case(SDL_QUIT):
            RasMolExit();
            break;

        case(SDL_WINDOWEVENT):
            switch( event->window.event )
            {   case(SDL_WINDOWEVENT_SIZE_CHANGED):
                case(SDL_WINDOWEVENT_RESIZED):
                    ResizeWindow( event->window.data1, event->window.data2 );
                    break;

                case(SDL_WINDOWEVENT_EXPOSED):
                    if( Texture ) TransferImage();
                    break;
            }
            break;

        case(SDL_MOUSEBUTTONDOWN):
            InitX = PointX = event->button.x;
            InitY = PointY = event->button.y;
            break;

        case(SDL_MOUSEMOTION):
            if( event->motion.state )
                if( !IsClose(event->motion.x,InitX) || !IsClose(event->motion.y,InitY) )
                {   MouseMove( event->motion.state, SDL_GetModState(),
                               event->motion.x-PointX, event->motion.y-PointY );
                    PointX = event->motion.x;
                    PointY = event->motion.y;
                }
            break;

        case(SDL_MOUSEBUTTONUP):
            PointX = event->button.x;
            PointY = event->button.y;
            if( IsClose(PointX,InitX) && IsClose(PointY,InitY) )
            {   if( SDL_GetModState() & (KMOD_SHIFT|KMOD_CTRL) )
                {      ReDrawFlag |= RFPoint1;
                } else ReDrawFlag |= RFPoint2;
            }
            break;
    }
}


int FetchEvent( wait )
    int wait;
{
    SDL_Event event;

    if( wait && !ReDrawFlag )
    {   if( SDL_WaitEventTimeout(&event,50) )
            ProcessEvent( &event );
    }

    while( SDL_PollEvent(&event) )
        ProcessEvent( &event );

    /* No menu system, so never a packed menu-selection result */
    return( 0 );
}


int OpenDisplay( x, y )
    int x, y;
{
    register int i;

    MouseMode = MMRasMol;
    UseHourGlass = True;
    DisableMenu = False;
    InitX = InitY = 0;

    for( i=0; i<8; i++ )
        DialValue[i] = 0.0;

    /* Default greyscale ramp, matching SetLutEntry()'s packing */
    Lut[0] = 65793*0;
    Lut[1] = 65793*64;
    Lut[2] = 65793*128;
    Lut[3] = 65793*196;
    Lut[4] = 65793*255;

    XRange = x;  WRange = XRange>>1;
    YRange = y;  HRange = YRange>>1;
    Range = MinFun(XRange,YRange);

    if( !Interactive ) return( False );

    if( SDL_Init(SDL_INIT_VIDEO) )
        return( 0 );

    Window = SDL_CreateWindow( "RasMol Version 2.6",
                               SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED,
                               x, y, SDL_WINDOW_RESIZABLE );
    if( !Window )
    {   SDL_Quit();
        return( 0 );
    }

    Renderer = SDL_CreateRenderer( Window, -1, SDL_RENDERER_ACCELERATED );
    if( !Renderer )
        Renderer = SDL_CreateRenderer( Window, -1, SDL_RENDERER_SOFTWARE );
    if( !Renderer )
    {   SDL_DestroyWindow(Window);
        SDL_Quit();
        return( 0 );
    }

    SDL_SetWindowMinimumSize( Window, 40, 40 );

    ArrowCursor = SDL_CreateSystemCursor( SDL_SYSTEM_CURSOR_ARROW );
    CrossCursor = SDL_CreateSystemCursor( SDL_SYSTEM_CURSOR_CROSSHAIR );
    WaitCursorPtr = SDL_CreateSystemCursor( SDL_SYSTEM_CURSOR_WAIT );
    if( CrossCursor ) SDL_SetCursor( CrossCursor );

    return( True );
}


int CreateImage()
{
    register Long size;

    if( FBuffer ) { free(FBuffer); FBuffer = (Pixel*)NULL; }
    size = (Long)XRange*YRange*sizeof(Pixel);
    FBuffer = (Pixel*)malloc( size+32 );
    if( !FBuffer ) return( False );

    if( !Interactive )
        return( True );

    if( Texture ) SDL_DestroyTexture( Texture );
    Texture = SDL_CreateTexture( Renderer, SDL_PIXELFORMAT_RGB888,
                                 SDL_TEXTUREACCESS_STREAMING, XRange, YRange );
    if( !Texture ) return( False );

    return( True );
}


void TransferImage()
{
    int texW, texH;

    if( !Texture ) return;

    /* A live edge-drag fires EXPOSED (-> here) before CreateImage()
     * resizes Texture/FBuffer to match the already-updated XRange/
     * YRange; blitting at the new pitch into the old-sized buffer
     * overreads FBuffer. Skip stale frames; the pending resize's own
     * RefreshScreen() redraws correctly once buffers are resized.
     */
    if( SDL_QueryTexture(Texture,NULL,NULL,&texW,&texH) ||
        (texW!=XRange) || (texH!=YRange) )
        return;

    SDL_UpdateTexture( Texture, NULL, FBuffer, XRange*sizeof(Pixel) );
    SDL_RenderClear( Renderer );
    SDL_RenderCopy( Renderer, Texture, NULL, NULL );
    SDL_RenderPresent( Renderer );
}


void ClearImage()
{
    if( !Renderer ) return;

    SDL_SetRenderDrawColor( Renderer, 0, 0, 0, 255 );
    SDL_RenderClear( Renderer );
    SDL_RenderPresent( Renderer );
}


int PrintImage()
{
    return( False );
}

int ClipboardImage()
{
    return( False );
}


void SetMouseMode( mode )
    int mode;
{
    MouseMode = mode;
}


void EnableMenus( flag )
    int flag;
{
    DisableMenu = !flag;
}


void UpdateScrollBars()
{
    /* No scrollbars in the minimal backend */
}


/* Small built-in table of common X11/SVG colour names, used as a
 * fallback when a colour name isn't one of RasMol's own predefined
 * colour tokens. Replaces the original XLookupColor() call into the
 * X server's rgb.txt database, which isn't available without X11.
 */
typedef struct {
    char *name;
    Byte r, g, b;
} NamedColour;

static NamedColour ColourTable[] = {
    { "aliceblue",      240,248,255 }, { "antiquewhite",   250,235,215 },
    { "aquamarine",     127,255,212 }, { "azure",          240,255,255 },
    { "beige",          245,245,220 }, { "bisque",         255,228,196 },
    { "blueviolet",     138,43, 226 }, { "brown",          165,42, 42  },
    { "burlywood",      222,184,135 }, { "cadetblue",       95,158,160 },
    { "chartreuse",     127,255,0   }, { "chocolate",      210,105,30  },
    { "coral",          255,127,80  }, { "cornflowerblue", 100,149,237 },
    { "cornsilk",       255,248,220 }, { "crimson",        220,20, 60  },
    { "darkblue",        0,  0,  139}, { "darkcyan",        0,  139,139},
    { "darkgoldenrod",  184,134,11  }, { "darkgray",       169,169,169 },
    { "darkgreen",       0,  100,0  }, { "darkgrey",       169,169,169 },
    { "darkkhaki",      189,183,107 }, { "darkmagenta",    139,0,  139 },
    { "darkolivegreen", 85, 107,47  }, { "darkorange",     255,140,0   },
    { "darkorchid",     153,50, 204 }, { "darkred",        139,0,  0   },
    { "darksalmon",     233,150,122 }, { "darkseagreen",   143,188,143 },
    { "darkslateblue",  72, 61, 139 }, { "darkslategray",  47, 79, 79  },
    { "darkturquoise",   0, 206,209 }, { "darkviolet",     148,0,  211 },
    { "deeppink",       255,20, 147 }, { "deepskyblue",     0, 191,255 },
    { "dimgray",        105,105,105 }, { "dodgerblue",      30,144,255 },
    { "firebrick",      178,34, 34  }, { "forestgreen",     34,139,34  },
    { "gainsboro",      220,220,220 }, { "gold",           255,215,0   },
    { "goldenrod",      218,165,32  }, { "gray",           190,190,190 },
    { "grey",           190,190,190 }, { "greenyellow",    173,255,47  },
    { "honeydew",       240,255,240 }, { "hotpink",        255,105,180 },
    { "indianred",      205,92, 92  }, { "indigo",          75,0,  130 },
    { "ivory",          255,255,240 }, { "khaki",          240,230,140 },
    { "lavender",       230,230,250 }, { "lawngreen",      124,252,0   },
    { "lemonchiffon",   255,250,205 }, { "lightblue",      173,216,230 },
    { "lightcoral",     240,128,128 }, { "lightcyan",      224,255,255 },
    { "lightgray",      211,211,211 }, { "lightgreen",     144,238,144 },
    { "lightgrey",      211,211,211 }, { "lightpink",      255,182,193 },
    { "lightsalmon",    255,160,122 }, { "lightseagreen",   32,178,170 },
    { "lightskyblue",   135,206,250 }, { "lightslategray", 119,136,153 },
    { "lightsteelblue", 176,196,222 }, { "lightyellow",    255,255,224 },
    { "lime",            0, 255,0   }, { "limegreen",       50,205,50  },
    { "linen",          250,240,230 }, { "maroon",         176,48, 96  },
    { "mediumaquamarine",102,205,170}, { "mediumblue",      0,  0,  205},
    { "mediumorchid",   186,85, 211 }, { "mediumpurple",   147,112,219 },
    { "mediumseagreen",  60,179,113 }, { "mediumslateblue",123,104,238 },
    { "mediumspringgreen",0,250,154 }, { "mediumturquoise", 72,209,204 },
    { "mediumvioletred",199,21, 133 }, { "midnightblue",    25,25, 112 },
    { "mintcream",      245,255,250 }, { "mistyrose",      255,228,225 },
    { "moccasin",       255,228,181 }, { "navajowhite",    255,222,173 },
    { "navy",            0,  0,  128}, { "navyblue",        0,  0,  128},
    { "oldlace",        253,245,230 }, { "olive",          128,128,0   },
    { "olivedrab",      107,142,35  }, { "orangered",      255,69, 0   },
    { "orchid",         218,112,214 }, { "palegoldenrod",  238,232,170 },
    { "palegreen",      152,251,152 }, { "paleturquoise",  175,238,238 },
    { "palevioletred",  219,112,147 }, { "papayawhip",     255,239,213 },
    { "peachpuff",      255,218,185 }, { "peru",           205,133,63  },
    { "pink",           255,192,203 }, { "plum",           221,160,221 },
    { "powderblue",     176,224,230 }, { "rosybrown",      188,143,143 },
    { "royalblue",       65,105,225 }, { "saddlebrown",    139,69, 19  },
    { "salmon",         250,128,114 }, { "sandybrown",     244,164,96  },
    { "seagreen",        46,139,87  }, { "seashell",       255,245,238 },
    { "sienna",         160,82, 45  }, { "silver",         192,192,192 },
    { "skyblue",        135,206,235 }, { "slateblue",      106,90, 205 },
    { "slategray",      112,128,144 }, { "snow",           255,250,250 },
    { "springgreen",      0,255,127 }, { "steelblue",       70,130,180 },
    { "tan",            210,180,140 }, { "teal",             0,128,128 },
    { "thistle",        216,191,216 }, { "tomato",         255,99, 71  },
    { "turquoise",       64,224,208 }, { "wheat",          245,222,179 },
    { "whitesmoke",     245,245,245 }, { "yellowgreen",    154,205,50  },
    { NULL, 0,0,0 }
};

int LookUpColour( name, red, grn, blu )
    char *name; int *red, *grn, *blu;
{
    register NamedColour *ptr;
    register char *src, *dst;
    char buffer[64];

    /* Normalise: lower-case, drop spaces (e.g. "Forest Green") */
    dst = buffer;
    for( src=name; *src && (dst-buffer)<63; src++ )
        if( *src != ' ' )
            *dst++ = (char)tolower((unsigned char)*src);
    *dst = '\0';

    for( ptr=ColourTable; ptr->name; ptr++ )
        if( !strcmp(ptr->name,buffer) )
        {   *red = ptr->r;
            *grn = ptr->g;
            *blu = ptr->b;
            return( True );
        }
    return( False );
}


void BeginWait()
{
    if( UseHourGlass && WaitCursorPtr )
        SDL_SetCursor( WaitCursorPtr );
}


void EndWait()
{
    if( UseHourGlass && CrossCursor )
        SDL_SetCursor( CrossCursor );
}


void CloseDisplay()
{
    if( Texture ) { SDL_DestroyTexture(Texture); Texture = NULL; }
    if( Renderer ) { SDL_DestroyRenderer(Renderer); Renderer = NULL; }
    if( Window ) { SDL_DestroyWindow(Window); Window = NULL; }
    if( ArrowCursor ) SDL_FreeCursor(ArrowCursor);
    if( CrossCursor ) SDL_FreeCursor(CrossCursor);
    if( WaitCursorPtr ) SDL_FreeCursor(WaitCursorPtr);
    SDL_Quit();
}
