#include <vector>
#include <functional>
#include <memory>
#include <map>
#include <filesystem>
#include <thread>
#include <mutex>
#include <string>
#include <set>

#include "imgui/examples/libs/gl3w/GL/gl3w.h"
#include "system/system.hpp"
extern bool imGuiMoveWindows;

#if lol == 1
#endif

#include <math.h>
#include <stdarg.h>
#define MIN(__a,__b) std::min(__a,__b) // ((__a)<(__b) ? (__a):(__b))
#define MAX(__a,__b) std::max(__a,__b) // ((__a)>(__b) ? (__a):(__b))

std::string saveDirectory = "C:/!mad/Projekte/Games/newstuff/binaryfeatures/demo/savies/";
std::string trackingFolder = "/Faces";

#define DEBUG(...) printf(__VA_ARGS__);
#define DEBUGTRACE(__name) //printf("%s->%d\n",__name,rand());

//-----------------
// - ESBIFT - Extremely Simple Brightness Invariant Feature Tracker -
//-----------------
#define DEFINED_USE_ESBIFT 1
#if DEFINED_USE_ESBIFT == 1
namespace ESBIFT {

//---------------------------------------
#define ROTATIONALDESCRIPTORS {tunnelPoint, tunnelX, tunnelY, polarCircleQuadratic}
#define NONROTATIONALDESCRIPTORS {polarCircleQuadratic}
//---------------------------------------

//---------------------------------------
#define GUI_DEFAULT_ALGORITHM_CLASS algo_class_resbift
#define GUI_DEFAULT_ALGORITHM algo_distance
//---------------------------------------

//---------------------------------------
#define ADAPTIVE_NONE 0
#define ADAPTIVE_ITERATIVE 1
#define ADAPTIVE_NONITERATIVEFINEADAPTIVE 2
#define ITERATIVECONFIGURATION NONE
#define NONITERATIVECONFIGURATION NONE
#define DEFAULTMIPLEVELSCALING 0.625
#define MULTICHANNEL 0
//---------------------------------------
#define DEFINED_ADAPTIVE_TYPE ADAPTIVE_NONE
//---------------------------------------
#if DEFINED_ADAPTIVE_TYPE == ADAPTIVE_NONE
    #define SAMPLINGTYPE NONITERATIVE
    #define ITERATIVECONFIGURATION FAST
    #define NONITERATIVECONFIGURATION GOOD
#endif

#if DEFINED_ADAPTIVE_TYPE == ADAPTIVE_ITERATIVE
    #define RESAMPLINGTYPE ITERATIVE
    #define SAMPLINGTYPE ITERATIVE
    #define ITERATIVECONFIGURATION FAST
    #define NONITERATIVECONFIGURATION FAST
#endif

#if DEFINED_ADAPTIVE_TYPE == ADAPTIVE_NONITERATIVEFINEADAPTIVE
    #define RESAMPLINGTYPE NONITERATIVE
    #define SAMPLINGTYPE ITERATIVE
    #define ITERATIVECONFIGURATION FAST
    #define NONITERATIVECONFIGURATION FAST
#endif
//---------------------------------------

#define COMPILERBUGGONE 1 // currently the compilers on my computer are modified.. At least that's the best erklärung for a bug...
#define FORCE_SINGLE_THREAD 0

//---------------------------------------
#define DEFINED_FORCESETUP_FOR_ONLYFLAT_TRACKING 1
#define DEFINED_FORCESETUP_FOR_MAXIMUM_PERSPECTIVE 1
//---------------------------------------

//---------------------------------------
#define DEFINED_ADAPTIVEBRIGHTNESSINVARIANCE 1 // some more locally focused brightness invariance (somehow already some sort of weak 180 degree rotations)
//---------------------------------------
#define DEFINED_COLORRANGE 0 // (doesn't work propperly) needs DEFINED_PUREBINARITY == 0, (e.g. for specular) + DEFINED_COLORRANGE_LOW/DEFINED_COLORRANGE_HIGH + needs normalized intensities
#define DEFINED_NORMALIZEDIMAGES 0 //DEFINED_COLORRANGE // feature images are in range 0..1
#define DEFINED_WITHROTATIONINVARIANCE 1 // 360 degree rotationinvariance (not so good looking, have to use sectors or limit to small angles (around 45 degrees))
#define DEFINED_WITHROTATIONINVARIANCE45DEGREES 0 // limit rotation invariance to 45 degrees
#define DEFINED_FULLROTATIONINVARIANCEBYSECTORS 1 // for 360 degrees through sampling of "8 images" else direct 360 degrees
#define DEFINED_ROUNDEDPATCHES 1
#define DEFINED_NORMAL3DPATCHES 1
#define DEFINED_TRIANGULAR_DERIVATIVES 0 // 0 leads to some more precission, and a more "correct" formulation (todo: maybe, add rotating triangular lookup to avoid convergence of features to triangular shape)
#define DEFINED_DESCRIPTORSHAPETOTRIANGLE 0 // somehow multiplies x by y (for more rotation invariance grip)
#define DEFINED_ADVANCEDPATTERNTYPE 0 // creates a more "meaningful" descriptor pattern (actually maybe for help at rotation invariance)
#define DEFINED_DERIVATIVE_ROTATION DEFINED_TRIANGULAR_DERIVATIVES //1 // (may have larger impact on non standard ESBIFT e.g. AIPOP stuff) "rotate" derivatives, is not necessary, TRIANGULAR_DERIVATIVES may benefit from it, maybe quads somehow too (not sure about the maths (didn't check this, yet, please check it if you are unsure about using it.), but it works, somehow)
#define DEFINED_ALSO_RESBIFT 1 // include RESBIFT algorithm classes
//---------------------------------------
#define DEFINED_REMOVEDUST 0 // some denoising if needed
#define DEFINED_COLORRANGE_LOW 0.0
#define DEFINED_COLORRANGE_HIGH (256.0/255.0) // e.g. for specular highlights
#define DEFINED_HASRANDOMMOVEMENT DEFINED_COLORRANGE // descriptors could get get stuck when unusable pixels are encountered, so this is a way to let them do some noise steps in that case
#define DEFINED_CLIPDESCRIPTORLEVELSTOLASTLEVEL 1 // (maybe temporary) disable not needed miplevels
#define DEFINED_CHANNELCOUNT (MULTICHANNEL == 1 ? 4 : 1) // either 4 channels or 1(maybe faster lookups) channel (yet)
#define DEFINED_MAXMULTICORE (FORCE_SINGLE_THREAD ? 1 : 7) // just if you want to limit threads (e.g. for propper youtube playback)
#define DEFINED_DESCRIPTIVITYMAP 1 // if there should be the possibility for building a descriptivity map e.g. for showing em
#define DEFINED_PUREBINARITY 0 // (pure binarity is of course better if descriptors are packed into floats(gpu) or ints or whatever) vundefined is somehow needed here for calculating a propper descriptivity value of the descriptors (selecting best trackable features)
#define DEFINED_BORDERCLAMPING 1 // without it the out of border pixels would be just marked as invalid
#define DEFINED_BORDERCLAMPMIRROR 0 // mirrors the descriptors on the border (DEFINED_BORDERCLAMPING 1 needed)
#define DEFINED_BORDERCLAMPWRAP 1 // mirrors the descriptors on the border (DEFINED_BORDERCLAMPING 1 needed)
#define DEFINED_CORRECTINVALIDFLOATS 1 // !really required! handles denormals, infinites and so on
#define DEFINED_FEATURETRACKING 1
#define DEFINED_VARYINGDESCRIPTORSIZE (DEFINED_ADAPTIVE_TYPE != ADAPTIVE_NONE ? 1 : 0) // e.g. for adaptive modes
#define DEFINED_MAXMIPLEVELLOOPSTEPS 2500 // how many steps to take maximally on one miplevel
#define DEFINED_ONLYNOPROBLEMATIC 1
#define DEFINED_ADVANCED_PERSPECTIVE 1 // tunnel descriptor shape size going wider than 1 (resulting in more on form focussed "tunnel" features)
#define DEFINED_PATCH_RANDOM_DEVIATION 0.5 // 3D patch perspective estimation factor
#define DEFINED_GUI_REMOVE_INVALID_FEATURES_AFTER_ADDINGS 1 // features to be found to far away from point of placing are removed after addition
//---------------------------------------
#define DEFINED_CLUSTER_ATTENUATION_FORMULA ((d*binkDistanceAttenuationFactor)*(d*binkDistanceAttenuationFactor)+0.01) // +0.01 may be useful for division by zero stuff
#define DEFINED_CORNER_BINKS 0 // Bink lookup only uses cornershapes e.g. sectors instead of whole "coordinate cross"
#define DEFINED_WEIGHTED_BINKS 1 // binks use feature importance information to check if they are valuable for binks
#define DEFINED_ITERATIVE_BINKS 0 // binks are recaptured every frame
#define DEFINED_REFINED_BINK_CENTERS 1 // non averaging for centers
//---------------------------------------
#define DEFINED_USE_DIFF37_RESOLVE 1 // if you want to decode the encoded 0x37 files
#define DEFINED_INCLUDE_REDUCED_SHADER_PREPROCESSOR 1 // maybe a preprocessor is needed to produce easily propper shader files later
//---------------------------------------
#define POWER_FEATUREREFINE SIGNED_POWER_FEATUREREFINE // unsigned power features also already perform good
#define FEATURERANDOMSEED 0
#define MAXDESCRIPTORSIZERANDOM (1.0+0.25)
#define MINDESCRIPTORSIZERANDOM (1.0-0.25)
#define DEFINED_POLARCOORDINATES DEFINED_WITHROTATIONINVARIANCE
#define TUBULARPOW 3.0
#define DENOISE_THRESHOLD 0.5
#define DENOISE_KERNEL_RAD 2
#define REMOVE_INVALID_FEATURES_AFTER_ADDINGS_AREA 0.6
#define POW_FOR_BINKS 0
#define UNITTEST 0
#if UNITTEST == 1
#define ActiveConfiguration Unit_Configuration
#define UNITTEST_FUNC unittest();
#else
#define ActiveConfiguration Configuration
#define UNITTEST_FUNC
#endif

#define EXCEPTIONTRY try
#define EXCEPTIONCATCH catch

#if DEFINED_FORCESETUP_FOR_ONLYFLAT_TRACKING == 1
    #define DEFINED_ADAPTIVEBRIGHTNESSINVARIANCE 0
    #define DEFINED_COLORRANGE 0
    #define DEFINED_NORMALIZEDIMAGES 0
#endif

#if DEFINED_FORCESETUP_FOR_MAXIMUM_PERSPECTIVE == 1
    #define DEFINED_PATCH_RANDOM_DEVIATION 1.0
#endif

#define THREADCOUNT std::thread::hardware_concurrency()

#if DEFINED_WITHROTATIONINVARIANCE == 1
    #define DESCRIPTORSHAPES rotationalDescriptors
#else
    #define DESCRIPTORSHAPES nonRotationalDescriptors
#endif

    using DOUBLE = float;
    using ICHANNEL = DOUBLE;
    using FCHANNEL = DOUBLE;
    using SINGLECHANNEL = FCHANNEL;
#define FCHANNELFROM(__v) (__v)

#if DEFINED_PUREBINARITY == 0
    using VECTORBOOL = unsigned char;
    #define vfalse      0
    #define vtrue       1
    #define vundefined  2   // on gpu you can pack the binary descriptor into "floats" (like I have done already for a lost GPU implementation)
    #define vunusable   3
#else
    using VECTORBOOL = bool;
    #define vfalse      false
    #define vtrue       true
    #define vundefined  false // take care here
    #define vunusable   false // take care here
#endif

#define DEFINED_USE_STL_VECTOR 1
#if DEFINED_USE_STL_VECTOR == 1
#define Array(__t) std::vector<__t>
#define ArrayResize(__t,__k) (__t).resize(__k)
#define ArraySize(__t) (__t).size()
#define ArrayAdd(__t,__a) (__t).push_back(__a)
#define ArrayAppend(__t,__a) SMAC(__t,__a,&t,const &a) t.insert(t.end(),a.begin(),a.end()) _SMAC
#define ArrayRemove(__t,__a) (__t).erase((__t).begin() + (__a))
#else
    template <class __T>
    class vectorWrap : public std::vector<__T> {
    public:
        int k = 0;
        vectorWrap() : std::vector<__T>() { ; }
        vectorWrap(const int n) : std::vector<__T>(n) { ; }
        void resize(int n) {
            k = n;
            vector::resize(n);
        }
        const int size() const {
            if (vector::size() != k)
                printf("k:%d l:%d\n", k, vector::size());
            return vector::size();
        }
        void push_back(const __T& c) {
            k++;
            vector::push_back(c);
        }
        void clear() {
            k = 0;
            vector::clear();
        }
    };
#define Array(__t) vectorWrap<__t>
#define ArrayResize(__t,__k) (__t).resize(__k)
#define ArraySize(__t) (__t).size()
#define ArrayAdd(__t,__a) (__t).push_back(__a)
#define ArrayRemove(__t,__a) (__t).erase((__t).begin() + (__a))
#endif

#define SWAP(__t,__a,__b) {__t = __a; __a = __b; __b = __t;}
#define FIXUPDIVBYZERO(__v,__r) ((__v) == 0.0 ? (__r) : (__v))
#define DEGREETORADIANS(__degrees) ((__degrees) * PI2 / 360.0)
#define POSITIONTOVALUE(__x,__y)  ( (__x) / ((__y)+sign(__y)*3.0)*2.0*fabs(__x) - fabs(__y) / fabs((__x)+sign(__x)*1.5)*0.5/fabs(__y) ) // the sign stuff is to prevent division by zero singularity
#define LMAC(__returntype,__p0,__p1,__v0,__v1) [&]() -> __returntype {auto &__v0 = __p0; auto __v1 = __p1; // LAMBDA MACRO
#define LMAC(__returntype) [&]() -> __returntype {
#define LMAC_ [&]() {
#define _LMAC ;}()
//#define SMAC(__p0,__p1,__p2,__v0,__v1,__v2) {auto __v0 = __p0; auto __v1 = __p1; auto __v2 = __p2; // SCOPE MACRO
#define SMAC(__p0,__p1,__v0,__v1) {auto __v0 = __p0; auto __v1 = __p1;
//#define SMAC(__p0,__v0) {auto __v0 = __p0;
#define SMAC_ {
#define _SMAC ;}

#define COMMENTOUT(...)
#define COMMENTIN(...) __VA_ARGS__
#define COMMENT0(...)
#define COMMENT1(...) COMMENTIN(__VA_ARGS__)
#define IFDEFINED(__v, ...) COMMENT(__VA_ARGS__) // seems not possible to implement propper yet? (perhaps an optional "resolve(stuff)" macro name resolving would help)

#define toString(__v) std::to_string(__v)
#define boolToString(__v) std::string(__v ? "yes" : "no")

#if DEFINED_CORRECTINVALIDFLOATS == 1
    #define VALIDFLOAT(__v) (!((std::fpclassify(__v) == FP_NAN) || (std::fpclassify(__v) == FP_INFINITE) || (std::fpclassify(__v) == FP_SUBNORMAL)))
    #define FIXFLOATNUMBER(__v,__r) if (!VALIDFLOAT(__v)) __v = (__r);
    #define FIXFLOATNUMBER2(__v1,__v2,__v3) ((__v1) == (__v2) ? (__v1) : (v__3))
#else
    #define VALIDFLOAT(__v) true
    #define FIXFLOATNUMBER(__v,__r)
    #define FIXFLOATNUMBER2(__v1,__v2,__v3)
#endif

#define SortDescriptorPair descriptorPairNoSort //descriptorPairRightDown

#if DEFINED_CHANNELCOUNT == 1
#define FeatureLayerType 1
#define FeatureLayerDepthFormula(__value,__layer) __value;//(pow(__layer * featureLayerBaseLayerMul + featureLayerBaseBase,__value))
#if FeatureLayerType == 0
#define FeatureLayerFormulaName "Grey"
#define FeatureLayerFormula (red + green + blue)
#endif
#if FeatureLayerType == 1
#define FeatureLayerFormulaName "Luminance"
#define FeatureLayerFormula (luminance(red,green,blue))
#endif
#if FeatureLayerType == 2
#define FeatureLayerFormulaName "Power1"
#define FeatureLayerFormula (pow(pow(luminance(red,green,blue),-red),-green))
#endif
#if FeatureLayerType == 3
#define FeatureLayerFormulaName "Quadratic"
#define FeatureLayerFormula (red*red+green*green+blue*blue)
#endif
#if FeatureLayerType == 4
#define BASE 0.9999999
#define FeatureLayerFormulaName "Expressive"
#define FeatureLayerFormula ((pow(BASE,-red)+pow(BASE,-green)+pow(BASE,-blue)))
#endif
#if FeatureLayerType == 5
#define FeatureLayerFormulaName "Lumacol"
#define FeatureLayerFormula (red + green + blue + red * green * blue)
#endif
#endif

#if DEFINED_CHANNELCOUNT == 4
#define FeatureLayerType 0
#define FeatureLayerDepthFormula(__value,__layer) __value;//(pow(__layer * featureLayerBaseLayerMul + featureLayerBaseBase,__value))
#if FeatureLayerType == 0
#define FeatureLayerFormulaName "PassThru"
#define FeatureLayerFormula channel
#endif
#endif

#define GRADIENTINVERSION + // possible to invert gradient by place + or - here (default is +)
#define DIRECTIONFORMULA bitFormula_forward
#define DIRECTIONFORMULASTRING directionstring_forward
#define DIRECTIONFORMULAOFFSET directionoffset_forward
    const std::string directionstring_forward = "forward";
    const std::string directionstring_invers = "invers";
    const std::string directionstring_alternate = "alternate"; // seems to be a degradation, yet

    //--------------------------------------------------------------------------------------------
    // threading stuff
    //--------------------------------------------------------------------------------------------
#define srandom(__a) srandomWithThreadId(__a)
#define random(__a) randomWithThreadId(__a)
#define randomandback(__a) randomandbackWithThreadId(__a) // debugging
    std::map<std::thread::id, unsigned int> seeds;
    std::map<std::thread::id, std::mutex> randomMutex;
    std::mutex randomMutex_mutex;
    std::mutex loading_mutex;
    std::mutex seeds_mutex;

    inline unsigned char rand8(unsigned char& seedlo, unsigned char& seedhi) {
        const bool k = seedlo & 128;
        seedlo <<= 1;
        if (seedhi & 1)
            seedlo++;
        seedhi >>= 1;
        if (k) {
            seedhi ^= 0xb4;
        }
        return seedlo ^ seedhi;
    }

    //--------------------------------------------------------------------------------------------

    // to be done better (concurrent_map or something)
    inline const double randomWithThreadId(const double randval) {
        const std::thread::id thisThreadId = std::this_thread::get_id();
        randomMutex_mutex.lock(); randomMutex[thisThreadId].lock(); randomMutex_mutex.unlock();
        seeds_mutex.lock(); unsigned int seed = seeds[thisThreadId]; seeds_mutex.unlock();
        int lo = rand8(((unsigned char*)&seed)[0], ((unsigned char*)&seed)[1]);
        int hi = rand8(((unsigned char*)&seed)[2], ((unsigned char*)&seed)[3]);
        seeds_mutex.lock(); seeds[thisThreadId] = seed; seeds_mutex.unlock();
        randomMutex_mutex.lock(); randomMutex[thisThreadId].unlock(); randomMutex_mutex.unlock();
        return (double)(lo + (hi << 8)) / 0x10000 * randval;
    }

    inline const double randomandbackWithThreadId(const double randval) {
        const std::thread::id thisThreadId = std::this_thread::get_id();
        randomMutex_mutex.lock(); randomMutex[thisThreadId].lock(); randomMutex_mutex.unlock();
        seeds_mutex.lock(); unsigned int seed = seeds[thisThreadId]; seeds_mutex.unlock();
        int lo = rand8(((unsigned char*)&seed)[0], ((unsigned char*)&seed)[1]);
        int hi = rand8(((unsigned char*)&seed)[2], ((unsigned char*)&seed)[3]);
        //seeds[thisThreadId] = seed;
        randomMutex_mutex.lock(); randomMutex[thisThreadId].unlock(); randomMutex_mutex.unlock();
        return (double)(lo + (hi << 8)) / 0x10000 * randval;
    }

    void srandomWithThreadId(const int seed) {
        const std::thread::id thisThreadId = std::this_thread::get_id();
        randomMutex_mutex.lock(); randomMutex[thisThreadId].lock(); randomMutex_mutex.unlock();
        unsigned char seedhi1 = 0x6f; unsigned char seedlo1 = 0x3f;
        unsigned char seedhi2 = 0x6f; unsigned char seedlo2 = 0x3f;
        for (int i = 0; i < (seed >> 4) & 15; ++i) rand8(seedlo1, seedhi1);
        for (int i = 0; i < seed & 15; ++i) rand8(seedlo2, seedhi2);
        seeds_mutex.lock(); seeds[thisThreadId] = ((int)seedhi1 << 8) | seedlo1 | ((((int)seedhi2 << 8) | seedlo2) << 16); seeds_mutex.unlock();
        randomMutex_mutex.lock(); randomMutex[thisThreadId].unlock(); randomMutex_mutex.unlock();
    }

    void initThread() {
        srandom(0);
    }

    inline const double randomLike(const int index, const double randval) {
        int b = index ^ (index * 11) ^ (index / 17) ^ (index >> 16) ^ (index * 1877) ^ (index * 8332) ^ (index * 173);
        b = b ^ (b << 8) ^ (b * 23);
        b >>= 3;
        return (double)(b & 0xffff) / 0x10000 * randval;
    }

    // this is the original algorithm
#define ADDITION_FEATUREREFINE(__channel) \
    if ((!bitFound) && (bitRequired)) {\
        movementX -= i0adx - i1adx;\
        movementY -= i0ady - i1ady;\
    }\
    if ((bitFound) && (!bitRequired)) {\
        movementX += i0adx - i1adx;\
        movementY += i0ady - i1ady;\
    }

// this seemed to be more reliable in some contexts however, (3 branches to "down right", 1 to "left up") i have no clue if this is a good idea at all
#define GLIDEADD_FEATUREREFINE(__channel)\
    if ((!bitFound) && bitRequired) {\
        movementX -= i0adx - i1adx;\
        movementY -= i0ady - i1ady;\
    }\
    else {\
        movementX += i0adx - i1adx;\
        movementY += i0ady - i1ady;\
    } 

// this is a more binary one, maybe used for very artificial color spaces (needs more and tinier sub steps)
#define STEP_FEATUREREFINE(__channel)\
    if ((!bitFound) && (bitRequired)) {\
        movementX -= signcmp2(i0adx,i1adx);\
        movementY -= signcmp2(i0ady,i1ady);\
    }\
    if ((bitFound) && (!bitRequired)) {\
        movementX += signcmp2(i0adx,i1adx);\
        movementY += signcmp2(i0ady,i1ady);\
    }

// this seemed to be more reliable in some contexts however, (3 branches to "down right", 1 to "left up") i have no clue if this is a good idea at all
#define GLIDESTEP_FEATUREREFINE(__channel)\
    if ((!bitFound) && bitRequired) {\
        movementX -= signcmp2(i0adx,i1adx);\
        movementY -= signcmp2(i0ady,i1ady);\
    }\
    else {\
        movementX += signcmp2(i0adx,i1adx);\
        movementY += signcmp2(i0ady,i1ady);\
    } 

// this is a version that takes the descriptor shape into account (more coarse features are targetted better by this)
#define DISTANCE_FEATUREREFINE(__channel)\
    DOUBLE s = bitRequired != bitFound ? -1 : 1;\
    if (s < 0) {\
        const DOUBLE dx = pair.x1 - pair.x0;\
        const DOUBLE dy = pair.y1 - pair.y0;\
        const DOUBLE d = sqrt(dx * dx + dy * dy);\
        const DOUBLE gdx = i0adx - i1adx;\
        const DOUBLE gdy = i0ady - i1ady;\
        if ((!bitFound) && (bitRequired)) {\
            movementX -= gdx * d;\
            movementY -= gdy * d;\
        }\
        if ((bitFound) && (!bitRequired)) {\
            movementX += gdx * d;\
            movementY += gdy * d;\
        }\
    }


/// ---------------------------------------------------------------------------------
/// Kinda Sort of Algorithms invented by the Youtube video selection AI or something (I didn't understood these algorithms in the lines below fully) (rotation invariance may be a big problem?)
/// ---------------------------------------------------------------------------------
// maybe I couldn't restore all ( maybe needs some help from AIs again :) )

// step based a version proposed by song lyrics of songs choosen by a year 2023 Artificial Intelligence Algorithm (the AIs or whatever resigned later from these ones a bigger bit, they just wanted to help me on a problem with these. These algorithms worked flawlessly for a long time.)
#define AIPOPADDSTEP_FEATUREREFINE(__channel)\
    if ((!bitFound) && (bitRequired)) {\
        if (((j+kt) & 1)==0) {\
            movementX -= idir1(idiff1(channel(c0adx, n), i0a), idiff1(channel(c1adx, n), i1a));\
            movementY -= idir1(idiff1(channel(c0ady, n), i0a), idiff1(channel(c1ady, n), i1a));\
        } else {\
            movementX -= idir2(idiff1(channel(c0adx, n), i0a), idiff1(channel(c1adx, n), i1a));\
            movementY -= idir2(idiff1(channel(c0ady, n), i0a), idiff1(channel(c1ady, n), i1a));\
        }\
    }\
    if ((bitFound) && (!bitRequired)) {\
        if (((j+kt) & 1)==1) {\
            movementX += idir1(idiff1(channel(c0adx, n), i0a), idiff1(channel(c1adx, n), i1a));\
            movementY += idir1(idiff1(channel(c0ady, n), i0a), idiff1(channel(c1ady, n), i1a));\
        } else {\
            movementX += idir2(idiff1(channel(c0adx, n), i0a), idiff1(channel(c1adx, n), i1a));\
            movementY += idir2(idiff1(channel(c0ady, n), i0a), idiff1(channel(c1ady, n), i1a));\
        }\
    }

// step based a version proposed by song lyrics of songs choosen by a year 2023 Artificial Intelligence Algorithm ( seems to be useful for some more specific structures (lights?))
// seems not working so well like the normal sub version (I suspect no bug, maybe it's just the behaviour, the "sub" versions where always rather "special" (by "exploding" the features instead of searching them)
#define AIPOPSUBSTEP_FEATUREREFINE(__channel) \
    if ((!bitFound) && (bitRequired)) {\
        if (((j+kt) & 1)==0) { \
            movementX += idir1(idiff1(channel(c0adx, n), i0a), idiff1(channel(c1adx, n), i1a));\
            movementY += idir1(idiff1(channel(c0ady, n), i0a), idiff1(channel(c1ady, n), i1a));\
        } else {\
            movementX += idir2(idiff1(channel(c0adx, n), i0a), idiff1(channel(c1adx, n), i1a));\
            movementY += idir2(idiff1(channel(c0ady, n), i0a), idiff1(channel(c1ady, n), i1a));\
        }\
    }\
    if ((bitFound) && (!bitRequired)) {\
        if (((j+kt) & 1)==1) {\
            movementX -= idir1(idiff1(channel(c0adx, n), i0a), idiff1(channel(c1adx, n), i1a));\
            movementY -= idir1(idiff1(channel(c0ady, n), i0a), idiff1(channel(c1ady, n), i1a));\
        } else {\
            movementX -= idir2(idiff1(channel(c0adx, n), i0a), idiff1(channel(c1adx, n), i1a));\
            movementY -= idir2(idiff1(channel(c0ady, n), i0a), idiff1(channel(c1ady, n), i1a));\
        }\
    }

// add version (based on a step based version) proposed by song lyrics of songs choosen by a 2023th Artificial Intelligence Algorithm
#define AIPOPADD_FEATUREREFINE(__channel) \
    if ((!bitFound) && (bitRequired)) {\
        if (((j+kt) & 1)==0) {\
            movementX -= xdir1(i0adx, i1adx);\
            movementY -= xdir1(i0ady, i1ady);\
        } else {\
            movementX -= xdir2(i0adx, i1adx);\
            movementY -= xdir2(i0ady, i1ady);\
        }\
    }\
    if ((bitFound) && (!bitRequired)) {\
        if (((j+kt) & 1)==1) {\
            movementX += xdir1(i0adx, i1adx);\
            movementY += xdir1(i0ady, i1ady);\
        } else {\
            movementX += xdir2(i0adx, i1adx);\
            movementY += xdir2(i0ady, i1ady);\
        }\
    }

// negative add version seems to be useful for some more specific structures (lights?) (based on a step based version) proposed by song lyrics of songs choosen by a 2023th Artificial Intelligence Algorithm
#define AIPOPSUB_FEATUREREFINE(__channel)\
    if ((!bitFound) && (bitRequired)) {\
        if (((j+kt) & 1)==0) {\
            movementX += xdir1(i0adx, i1adx);\
            movementY += xdir1(i0ady, i1ady);\
        } else {\
            movementX += xdir2(i0adx, i1adx);\
            movementY += xdir2(i0ady, i1ady);\
        }\
    }\
    if ((bitFound) && (!bitRequired)) {\
        if (((j+kt) & 1)==1) {\
            movementX -= xdir1(i0adx, i1adx);\
            movementY -= xdir1(i0ady, i1ady);\
        } else {\
            movementX -= xdir2(i0adx, i1adx);\
            movementY -= xdir2(i0ady, i1ady);\
        }\
    }

/// ---------------------------------------------------------------------------------
/// Power Based Functions for local Brightness Invariance
/// ---------------------------------------------------------------------------------
#define SIGNEDPOWER_VAL_2SIDED(__i0ad, __i1ad) ((randomLike(randomLikeBase + j + 123134,1.0) < 0.5) ? sign(__i0ad) * pow(fabs(__i0ad), fabs(__i1ad)) : -sign(__i1ad) * pow(fabs(__i1ad), fabs(__i0ad)))
#define SIGNEDPOWER_VAL(__i0ad, __i1ad) (sign(__i0ad) * pow(fabs(__i0ad), fabs(__i1ad))) // this determines how much influence i0ad should have in a step
#define SIGNEDPOWER SIGNEDPOWER_VAL // single sided version ("developed" by the "youtubeai" (youtube video order) dunno if it makes sense to use the 2sided one

#define SIGNED_POWER_FEATUREREFINE(__channel)\
    if ((!bitFound) && (bitRequired)) {\
        DOUBLE oldMovementX = movementX;\
        DOUBLE oldMovementY = movementY;\
        movementX -= SIGNEDPOWER(i0adx, i1adx);\
        movementY -= SIGNEDPOWER(i0ady, i1ady);\
        FIXFLOATNUMBER(movementX,oldMovementX);\
        FIXFLOATNUMBER(movementY,oldMovementY);\
    }\
    if ((bitFound) && (!bitRequired)) {\
        DOUBLE oldMovementX = movementX;\
        DOUBLE oldMovementY = movementY;\
        movementX += SIGNEDPOWER(i0adx, i1adx);\
        movementY += SIGNEDPOWER(i0ady, i1ady);\
        FIXFLOATNUMBER(movementX,oldMovementX);\
        FIXFLOATNUMBER(movementY,oldMovementY);\
    }

// step is some sort of pow (pow actually is some sort of pow pow pow) :) (here only the biggest value is choosen to be the step (and inverted for i0 to i1))
#define STEPSIGNEDPOWER(__i0ad, __i1ad) (fabs(__i0ad) > fabs(__i1ad) ? sign(__i0ad) : 0) // since it's local brightness invariant the second term is 0, actually -sign(__i1ad) there shows much better results on non brightness invariant features
#define STEP_SIGNED_POWER_FEATUREREFINE(__channel) \
    if ((!bitFound) && (bitRequired)) { \
        movementX -= STEPSIGNEDPOWER(i0adx, i1adx);\
        movementY -= STEPSIGNEDPOWER(i0ady, i1ady);\
    } \
    if ((bitFound) && (!bitRequired)) { \
        movementX += STEPSIGNEDPOWER(i0adx, i1adx);\
        movementY += STEPSIGNEDPOWER(i0ady, i1ady);\
    }

// a combination of local brightness invariance and normal step based one for rather good behaviour in both settings
#define COMBINESTEP_FEATUREREFINE(__channel)\
    if ((!bitFound) && (bitRequired)) {\
        movementX -= signcmp2(i0adx, i1adx) * 0.5;\
        movementY -= signcmp2(i0ady, i1ady) * 0.5;\
        movementX -= STEPSIGNEDPOWER(i0adx, i1adx) * 0.5;\
        movementY -= STEPSIGNEDPOWER(i0ady, i1ady) * 0.5;\
    }\
    if ((bitFound) && (!bitRequired)) {\
        movementX += signcmp2(i0adx,i1adx);\
        movementY += signcmp2(i0ady,i1ady);\
    }

// mix between signed power and aipop algorithm (very strongly by me (but no clue why (using fabs() and so on) anymore, have to check again later)) (it's also local luma dependend there may be even better versions with no 0 -> STEPSIGNEDPOWER(__i0ad, __i1ad))
#define AIPOP_SIGNED_POWER_FEATUREREFINE(__channel)\
    if ((!bitFound) && (bitRequired)) {\
        DOUBLE oldMovementX = movementX;\
        DOUBLE oldMovementY = movementY;\
        movementX -= i1adx < 0 ? fabs(SIGNEDPOWER_VAL(i0adx, i1adx)) : 0;\
        movementY -= i1ady < 0 ? fabs(SIGNEDPOWER_VAL(i0ady, i1ady)) : 0;\
        FIXFLOATNUMBER(movementX,oldMovementX);\
        FIXFLOATNUMBER(movementY,oldMovementY);\
    }\
    if ((bitFound) && (!bitRequired)) {\
        DOUBLE oldMovementX = movementX;\
        DOUBLE oldMovementY = movementY;\
        movementX += i0adx < 0 ? fabs(SIGNEDPOWER_VAL(i1adx, i0adx)) : 0;\
        movementY += i0ady < 0 ? fabs(SIGNEDPOWER_VAL(i1ady, i0ady)) : 0;\
        FIXFLOATNUMBER(movementX,oldMovementX);\
        FIXFLOATNUMBER(movementY,oldMovementY);\
    }

// mix between signed power and !restored! aipop algorithm
#define AIPOP2_SIGNED_POWER_FEATUREREFINE(__channel)\
    if ((!bitFound) && (bitRequired)) {\
        DOUBLE oldMovementX = movementX;\
        DOUBLE oldMovementY = movementY;\
        if (((j+kt) & 1)==0) {\
            movementX -= i1adx < 0 ? fabs(SIGNEDPOWER_VAL(i0adx, i1adx)) : 0;\
            movementY -= i1ady < 0 ? fabs(SIGNEDPOWER_VAL(i0ady, i1ady)) : 0;\
        } else {\
            movementX -= i0adx < 0 ? fabs(SIGNEDPOWER_VAL(i1adx, i0adx)) : 0;\
            movementY -= i0ady < 0 ? fabs(SIGNEDPOWER_VAL(i1ady, i0ady)) : 0;\
        }\
        FIXFLOATNUMBER(movementX,oldMovementX);\
        FIXFLOATNUMBER(movementY,oldMovementY);\
    }\
    if ((bitFound) && (!bitRequired)) {\
        DOUBLE oldMovementX = movementX;\
        DOUBLE oldMovementY = movementY;\
        if (((j+kt) & 1)==1) {\
            movementX += i1adx < 0 ? fabs(SIGNEDPOWER_VAL(i0adx, i1adx)) : 0;\
            movementY += i1ady < 0 ? fabs(SIGNEDPOWER_VAL(i0ady, i1ady)) : 0;\
        } else {\
            movementX += i0adx < 0 ? fabs(SIGNEDPOWER_VAL(i1adx, i0adx)) : 0;\
            movementY += i0ady < 0 ? fabs(SIGNEDPOWER_VAL(i1ady, i0ady)) : 0;\
        }\
        FIXFLOATNUMBER(movementX,oldMovementX);\
        FIXFLOATNUMBER(movementY,oldMovementY);\
    }

#define AIPOP3_SIGNED_POWER_FEATUREREFINE(__channel) \
    if ((!bitFound) && (bitRequired)) {\
        DOUBLE oldMovementX = movementX;\
        DOUBLE oldMovementY = movementY;\
        movementX -= i1adx < 0 ? fabs(SIGNEDPOWER_VAL(i0adx, i1adx)) * 0.5 : 0;\
        movementY -= i1ady < 0 ? fabs(SIGNEDPOWER_VAL(i0ady, i1ady)) * 0.5 : 0;\
        movementX -= i0adx < 0 ? fabs(SIGNEDPOWER_VAL(i1adx, i0adx)) * 0.5 : 0;\
        movementY -= i0ady < 0 ? fabs(SIGNEDPOWER_VAL(i1ady, i0ady)) * 0.5 : 0;\
        FIXFLOATNUMBER(movementX,oldMovementX);\
        FIXFLOATNUMBER(movementY,oldMovementY);\
    }\
    if ((bitFound) && (!bitRequired)) {\
        DOUBLE oldMovementX = movementX;\
        DOUBLE oldMovementY = movementY;\
        movementX += i1adx < 0 ? fabs(SIGNEDPOWER_VAL(i0adx, i1adx)) * 0.5 : 0;\
        movementY += i1ady < 0 ? fabs(SIGNEDPOWER_VAL(i0ady, i1ady)) * 0.5 : 0;\
        movementX += i0adx < 0 ? fabs(SIGNEDPOWER_VAL(i1adx, i0adx)) * 0.5 : 0;\
        movementY += i0ady < 0 ? fabs(SIGNEDPOWER_VAL(i1ady, i0ady)) * 0.5 : 0;\
        FIXFLOATNUMBER(movementX,oldMovementX);\
        FIXFLOATNUMBER(movementY,oldMovementY);\
    }

/// ---------------------------------------------------------------------------------
/// Normal AIPop Power Based Functions
/// ---------------------------------------------------------------------------------
// -fabs is superflous
#define AIPOP4_SIGNED_POWER_FEATUREREFINE(__channel)\
    if ((!bitFound) && (bitRequired)) {\
        DOUBLE oldMovementX = movementX;\
        DOUBLE oldMovementY = movementY;\
        if (((j+kt) & 1)==0) {\
            movementX -= i1adx < 0 ? fabs(SIGNEDPOWER_VAL(i0adx, i1adx)) : -fabs(SIGNEDPOWER_VAL(i1adx, i0adx));\
            movementY -= i1ady < 0 ? fabs(SIGNEDPOWER_VAL(i0ady, i1ady)) : -fabs(SIGNEDPOWER_VAL(i1adx, i0adx));\
        } else {\
            movementX -= i0adx < 0 ? fabs(SIGNEDPOWER_VAL(i1adx, i0adx)) : -fabs(SIGNEDPOWER_VAL(i0adx, i1adx));\
            movementY -= i0ady < 0 ? fabs(SIGNEDPOWER_VAL(i1ady, i0ady)) : -fabs(SIGNEDPOWER_VAL(i0adx, i1adx));\
        }\
        FIXFLOATNUMBER(movementX,oldMovementX);\
        FIXFLOATNUMBER(movementY,oldMovementY);\
    }\
    if ((bitFound) && (!bitRequired)) {\
        DOUBLE oldMovementX = movementX;\
        DOUBLE oldMovementY = movementY;\
        if (((j+kt) & 1)==1) {\
            movementX += i1adx < 0 ? fabs(SIGNEDPOWER_VAL(i0adx, i1adx)) : -fabs(SIGNEDPOWER_VAL(i1adx, i0adx));\
            movementY += i1ady < 0 ? fabs(SIGNEDPOWER_VAL(i0ady, i1ady)) : -fabs(SIGNEDPOWER_VAL(i1adx, i0adx));\
        } else {\
            movementX += i0adx < 0 ? fabs(SIGNEDPOWER_VAL(i1adx, i0adx)) : -fabs(SIGNEDPOWER_VAL(i0adx, i1adx));\
            movementY += i0ady < 0 ? fabs(SIGNEDPOWER_VAL(i1ady, i0ady)) : -fabs(SIGNEDPOWER_VAL(i0adx, i1adx));\
        }\
        FIXFLOATNUMBER(movementX,oldMovementX);\
        FIXFLOATNUMBER(movementY,oldMovementY);\
    }


//#define UNSIGNEDPOWER_VAL(__i0ad, __i1ad) pow(fabs(__i0ad), fabs(__i1ad))
//#define UNSIGNEDPOWER_VAL_EXPSIGN(__k,__i0ad,__i1ad) (((__k & 1)*2-1) * UNSIGNEDPOWER_VAL(__i0ad,__i1ad))
//#define UNSIGNED_POWER_FEATUREREFINE(__channel) \
//    const int r = j;\
//    if ((!bitFound) && (bitRequired)) { \
//        DOUBLE oldMovementX = movementX; \
//        DOUBLE oldMovementY = movementY; \
//        movementX += UNSIGNEDPOWER_VAL_EXPSIGN(r, i0adx, i1adx); \
//        movementY += UNSIGNEDPOWER_VAL_EXPSIGN(r, i0ady, i1ady); \
//        FIXFLOATNUMBER(movementX,oldMovementX);\
//        FIXFLOATNUMBER(movementY,oldMovementY);\
//    } \
//    if ((bitFound) && (!bitRequired)) { \
//        DOUBLE oldMovementX = movementX; \
//        DOUBLE oldMovementY = movementY; \
//        movementX -= UNSIGNEDPOWER_VAL_EXPSIGN(r, i0adx, i1adx); \
//        movementY -= UNSIGNEDPOWER_VAL_EXPSIGN(r, i0ady, i1ady); \
//        FIXFLOATNUMBER(movementX,oldMovementX);\
//        FIXFLOATNUMBER(movementY,oldMovementY);\
//    }

/// ---------------------------------------------------------------------------------
// Some Angels just say that this one will get more important later, however it's rather bad currently
#define MULTIPLICATION_FEATUREREFINE(__channel)\
    if ((!bitFound) && (bitRequired)) {\
        movementX -= i0adx * i1adx;\
        movementY -= i0ady * i1ady;\
    }\
    if ((bitFound) && (!bitRequired)) {\
        movementX += i0adx * i1adx;\
        movementY += i0ady * i1ady;\
    }

#define FILESYSTEM std::filesystem

    namespace Configuration {

        namespace DEFAULTS {
            namespace DEFAULTMIP {
                const double descriptorLevelScale = DEFAULTMIPLEVELSCALING;
                const double stepFactorIncrease = 1.0;
            }
        }

        //-----------------
        namespace ITERATIVE {

            const bool iterative = true;

            //-----------------
            namespace PRECISION {
                const int pairCount = 1024;
                const double refinementSteps = 1024;
                const double mipsToTakeIntoAccountRatio = 0.5;
                using namespace DEFAULTS::DEFAULTMIP;
            }
            //using namespace PRECISION;
            //-----------------

            //-----------------
            namespace GOOD {
                const int pairCount = 100;
                const double refinementSteps = 100;
                const double mipsToTakeIntoAccountRatio = 0.5;
                using namespace DEFAULTS::DEFAULTMIP;
            }
            //using namespace GOOD;
            //-----------------

            //-----------------
            namespace FAST {
                const int pairCount = 32;
                const double refinementSteps = 32;
                const double mipsToTakeIntoAccountRatio = 0.5;
                using namespace DEFAULTS::DEFAULTMIP;
            }
            //using namespace FAST;
            //-----------------

            //-----------------
            namespace COMICTRACE {
                const int pairCount = 200;
                const double refinementSteps = 500;
                const double mipsToTakeIntoAccountRatio = 0.25;
                using namespace DEFAULTS::DEFAULTMIP;
            }
            //using namespace COMICTRACE;
            //-----------------

            //-----------------
            namespace NONE {
                const int pairCount = -1;
                const double refinementSteps = -1;
                const double mipsToTakeIntoAccountRatio = 1.0;
                using namespace DEFAULTS::DEFAULTMIP;
            }
            //using namespace NONE;
            //-----------------

            //-----------------
            using namespace ITERATIVECONFIGURATION;
        }
        //using namespace ITERATIVE;
        //-----------------

        //-----------------
        namespace NONITERATIVE {

            const bool iterative = false;

            //-----------------
            namespace FAST {
                const int pairCount = 32;
                const double refinementSteps = 32;
                const double mipsToTakeIntoAccountRatio = 1.0;
                using namespace DEFAULTS::DEFAULTMIP;
            }
            //using namespace GOOD;
            //-----------------

            //-----------------
            namespace GOOD {
                const int pairCount = 50;
                const double refinementSteps = 50;
                const double mipsToTakeIntoAccountRatio = 1.0;
                using namespace DEFAULTS::DEFAULTMIP;
            }
            //using namespace GOOD;
            //-----------------

            //-----------------
            namespace IMAGETRACKING {
                const int pairCount = 40;
                const double refinementSteps = 400;
                const double mipsToTakeIntoAccountRatio = 1.0;
                using namespace DEFAULTS::DEFAULTMIP;
            }
            //using namespace IMAGETRACKING;
            //-----------------

            //-----------------
            namespace SCALELEVEL { // just for testing
                const int pairCount = 30;
                const double refinementSteps = 400;
                const double mipsToTakeIntoAccountRatio = 1.0;
                using namespace DEFAULTS::DEFAULTMIP;
            }
            //using namespace SCALELEVEL;
            //-----------------

            //-----------------
            namespace VERYGOOD {
                const int pairCount = 100;
                const double refinementSteps = 400;
                const double mipsToTakeIntoAccountRatio = 1.0;
                using namespace DEFAULTS::DEFAULTMIP;
            }
            //using namespace VERYGOOD;
            //-----------------


            //-----------------
            namespace COMICTRACE {
                const int pairCount = 100;
                const double refinementSteps = 500;
                const double mipsToTakeIntoAccountRatio = 1.0;
                using namespace DEFAULTS::DEFAULTMIP;
            }
            //using namespace GOOD;
            //-----------------

            //-----------------
            namespace NONE {
                const int pairCount = -1;
                const double refinementSteps = -1;
                const double mipsToTakeIntoAccountRatio = 1.0;
                using namespace DEFAULTS::DEFAULTMIP;
            }
            //using namespace NONE;
            //-----------------

            //-----------------
            using namespace NONITERATIVECONFIGURATION;
        }
        //using namespace NONITERATIVE;
        //-----------------
        using namespace SAMPLINGTYPE;
        //-----------------

#if DEFINED_ADAPTIVE_TYPE != ADAPTIVE_NONE
        //-----------------
        namespace ADAPTIVERESAMPLING {

            //-----------------
            namespace NONITERATIVE {
                const bool resampleIterative = false;
                const double resampleRefinementSteps = refinementSteps * 4;
                const double resampleStepFactorIncrease = stepFactorIncrease;
                const double resampleMipsToTakeIntoAccountRatio = 1.0;
            }
            //using namespace NONITERATIVE;
            //-----------------

            //-----------------
            namespace ITERATIVE {
                const bool resampleIterative = true;
                const double resampleRefinementSteps = refinementSteps * 4;
                const double resampleStepFactorIncrease = stepFactorIncrease;
                const double resampleMipsToTakeIntoAccountRatio = 0.5;
            }
            //using namespace ITERATIVE;
            //-----------------

            //-----------------
            using namespace RESAMPLINGTYPE;
        }
        using namespace ADAPTIVERESAMPLING;
        //-----------------
#endif

        // algorithm parameters
#define SUBPIXELFACTOR 0.5 // how many pixel one step represents (on each miplevel (without additional scaling)) (steps = config.refinementSteps)
#define STEPBASE (1.0*SUBPIXELFACTOR)
        const DOUBLE stepFactor = STEPBASE * stepFactorIncrease;
#if DEFINED_ADAPTIVE_TYPE != ADAPTIVE_NONE
        const DOUBLE resampleStepFactor = STEPBASE * resampleStepFactorIncrease;
        const DOUBLE resampleRefinementStepPerLevelScaleFactor = 1.0; // multiplication of stepcount (and divisionof stepsize), on every subsequent finer(or signed=coarser) miplevel
#endif
        const DOUBLE refinementStepPerLevelScaleFactor = 1.0; // multiplication of stepcount (and divisionof stepsize), on every subsequent finer(or signed=coarser) miplevel

        // step base (e.g. compare)
        const DOUBLE constantStepRefinementStepsFactor = 0.1; // after which ratio of config.refinementSteps a stepbased(compare) step should be made
        const DOUBLE constantStepMoveFactor = 0.5; // the ratio of a stepbased(compare) step distance(stepFactor) to that of a normal step

        // pruning
        const DOUBLE pruneByZoomLevelRatio = 0.5; // take 0..n..1 most descriptive zoomlevels
        const DOUBLE pruneByAngularRatio = 0.5; // take 0..n..1 most descriptive rotation angles
        const DOUBLE pruneByTranslationalRatio = 0.5; // take 0..n..1 most descriptive distance values

        // debug
        bool animOn = false;
        bool addNewContinuosFeatures = false;

        // gui
        double descriptorThresholdRatioForAdaptiveUpdate = 0.7;
        int currentFrame = 0;
        int successorFrame = 0;
        int trackingFrame = 0;
        int requestedBestFeatureCount = 100;
        int requestedRandomFeatureCount = 100;
        bool useSpice = false;
        DOUBLE spicePercentage = 100.0;
        DOUBLE spiceIntensity = 0.01;
        bool displayFeatureLayer = false;
        int displayFeatureLayerChannel = -1;
        bool useResampling = true;
        int pyramidLayer = 0;
        int descriptivityLayer = 0;
        bool forceLayer = false;
        bool manhattenDistance = true;
        bool allFeatures = false;
        DOUBLE baseZoom = 1.0;
        DOUBLE updatedZoom = 1.0;
        DOUBLE trackingZoom = 1.0;
        double imageFeatureAlpha = 0.6;
        bool viewNormalized = false;
        bool viewAlpha = false;
        bool viewDescriptivenes = true;
        bool viewClassifications = false;
        bool viewFeatures = true;
        bool viewMesh = false;
        int dataSetFirstFrame = 0;
        int dataSetLastFrame = 24;
        int dataSetStepSize = 1;
        bool dataSetReversed = false;
        bool useFirstFrameValue = true;
        bool useLastFrameValue = true;
        bool useOwnWidth = true;
        bool placeMink = false;
        int minkType = 0;
        bool selectMinks = false;
        bool sampleDescriptors = iterative;
        int desiredImageWidth = 640;
        double desiredImageScale = 1.0;
        int requestedCircularFeatureSectors = 16;
        int requestedCircularFeatureRadis = 8;
        bool showTrackingDescriptor = false;
        bool trackingPerFrame = false;
        double pruneThreshold = 0.0;
        int descriptorAlgorithm = 0;
        DOUBLE binkDistanceAttenuationFactor = 0.05;
        DOUBLE binkDistanceAttenuationThreshold = 0.1;
        DOUBLE minkMoveThreshold = 2; // pixels
        DOUBLE binkForceThreshold = 0.05;
        bool circularLanding = false;
        DOUBLE circularLandingFactor = 0.5;
        bool pyramidal = true;
        enum {
            algo_class_esbift = 0,
            algo_class_resbift, // randomized esbift
        };
        enum {
            algo_distance = 0,
            algo_add,
            algo_glideadd,
            algo_step,
            algo_glidestep,
            algo_power,
            algo_powerstep,
            algo_combinestep,
            algo_aipop4poweradd,
            algo_aipopaddstep,
            algo_aipopsubstep,
            algo_multiply,
            algo_aipopadd,
            algo_aipopsub,
            algo_aipoppoweradd,
            algo_aipop2poweradd,
            algo_aipop3poweradd,
        };
        bool RESBIFT = GUI_DEFAULT_ALGORITHM_CLASS == algo_class_resbift;
        DOUBLE resbiftFactor = 2.0;
        int algorithm = GUI_DEFAULT_ALGORITHM;
        const char* algorithmNames[] = {
            "distance",
            "original(add)",
            "glideadd",
            // step based ones require some different performance settings for propper movement
            "step",
            "glidestep",
            // maybe good for local brightness invariance?
            "luma_power",
            "luma_powerstep",
            "luma_combinestep",
            // this one is totally new I didn't check it's properties yet
            "aipop4power",
                // algorithms that only work without rotation invariance or other things
            "aipopaddstep",
            "aipopsubstep",
            // these algorithms got some more problems and where invented mainly for the iterative version (with non existent non iterative mode)
            "multiply",
            "aipopadd",
            "aipopsub",
            // some ai stuff with local luma
            "luma_aipoppower",
            "luma_aipop2power",
            "luma_aipop3power",
        };
        const int algorithmCount = sizeof(algorithmNames) / sizeof(const char*);
        int currentDatasetNr = 0;
        double loopStepsIncreaseMultiplicator = 1.0; // maybe useful some day
        // internal
        std::string dataSetsFolder = "C:/!mad/Daten/Odometry/DataSets/Selection";
        std::string dataSetFileName;
        int descriptorRandomSeed = 0;
        const double quadraticDescriptorGrip = 0.5; // 0.5 seems to be good value
        const double roundedPatchLensFactor = 0.4; // how curvy the patch is
        const double roundedPatchRoundRatio = 1.0; // fade between the linear patch and the curvy patch

#define PATCHES_RANDOM_VALUE(__a) randoma_patches_##__a

#define PATCHVARS \
        DOUBLE PATCHES_RANDOM_VALUE(0) = 0.0;\
        DOUBLE PATCHES_RANDOM_VALUE(1) = 0.0;\
        DOUBLE PATCHES_RANDOM_VALUE(2) = 0.0;\
        DOUBLE PATCHES_RANDOM_VALUE(3) = 0.0;
    }
    //using namespace Configuration;
    //-----------------
#if UNITTEST == 1
#define UNIT_HEADER 1
#include "esbift_unit.hpp"
#endif
    using namespace ActiveConfiguration;

    //-----------------
    namespace Multiprocessor {
        int coreCount = 1;
        Array(std::thread) threadPool;
#define JOIN(...) EXCEPTIONTRY {__VA_ARGS__.join();} EXCEPTIONCATCH(...) {printf("exception1 that should be impossible\n");};
        void threadPoolThread(const std::function<void()>& f) {
#ifdef DEFINED_MAXMULTICORE != 1
            static int k = 0;
            if (threadPool[k % ArraySize(threadPool)].joinable())
                JOIN(threadPool[k % ArraySize(threadPool)]);
            threadPool[k % ArraySize(threadPool)] = std::thread(f);
            k++;
#else
            f();
#endif
        }
        void threadPoolFinish() {
#ifdef DEFINED_MAXMULTICORE != 1
            for (int i = 0; i < ArraySize(threadPool); ++i) {
                if (threadPool[i].joinable()) {
                    EXCEPTIONTRY{ threadPool[i].join(); } EXCEPTIONCATCH(...) {
                        printf("exception2 that should be impossible\n");
                    }
                }
            }
#endif
        }
        void initMultiCore() {
            coreCount = THREADCOUNT;
            coreCount = MIN(coreCount, DEFINED_MAXMULTICORE);
            printf("max used cores:%d\n", coreCount);
            ArrayResize(threadPool, coreCount);
            initThread();
        }

        void rebuildThreadPool() {
            threadPoolFinish();
            coreCount = MIN(coreCount, DEFINED_MAXMULTICORE);
            printf("max used cores:%d\n", coreCount);
            ArrayResize(threadPool, coreCount);
        }
    }
    using namespace Multiprocessor;
    //-----------------

    //-----------------
    namespace Debug {
        int featuresInvolved = 0;
        int sampledDescriptors = 0;
        int adaptiveSampledDescriptors = 0;
        float debugFloat1 = 0;
#define macheweitermitdemloophere continue // continue had some problems in macros and maybe that maybe useful for debugging
#define PROCESSMACHEWEITER 
        //#define macheweitermitdemloophere macheweiter = true;
        //        bool macheweiter = false;
        //        bool lol = false;
        //#define PROCESSMACHEWEITER lol = macheweiter; debugFloat1 += lol ? 1 : 0; macheweiter = false; if (!lol)
    }
    using namespace Debug;
    //-----------------

    // ---------------------------------------------------------------------
#define SINGLECHANNELMINMAX(__l,__mincolor,__maxcolor,__first) \
    {\
        const SINGLECHANNEL &value = __l;\
        if (value < __mincolor || __first) __mincolor = value;\
        if (value > __maxcolor || __first) __maxcolor = value;\
        __first = false;\
    }
#if DEFINED_CHANNELCOUNT == 4
#define CHANNELMINMAX(__l,__mincolor,__maxcolor,__first) \
    {\
        const CHANNELS &value = __l;\
        if (value.r < __mincolor.r || __first) __mincolor.r = value.r;\
        if (value.g < __mincolor.g || __first) __mincolor.g = value.g;\
        if (value.b < __mincolor.b || __first) __mincolor.b = value.b;\
        if (value.r > __maxcolor.r || __first) __maxcolor.r = value.r;\
        if (value.g > __maxcolor.g || __first) __maxcolor.g = value.g;\
        if (value.b > __maxcolor.b || __first) __maxcolor.b = value.b;\
        __first = false;\
    }
#else
#define CHANNELMINMAX(__l,__mincolor,__maxcolor,__first) SINGLECHANNELMINMAX((__l),(__mincolor),(__maxcolor),(__first))
#endif

#define CLAMP(__v, __mn, __mx) ((__v) < (__mn) ? (__mn) : ((__v) > (__mx) ? (__mx) : (__v)))
#define JUSTSETUPIMAGEHULL(__d,__v)\
    {\
        __d.width = (__v).width;\
        __d.height = (__v).height;\
        __d.level = (__v).level;\
        ArrayResize(__d.pixels,ArraySize((__v).pixels));\
    }
    // ---------------------------------------------------------------------

    using UNIQUEID = unsigned int;
    UNIQUEID uniqueId = 0; // some sort of a "uuid"

    //-----------------
    namespace Math {

        inline DOUBLE clamp(const DOUBLE v, const DOUBLE mn, const DOUBLE mx) {
            return CLAMP(v, mn, mx);
        }

        inline DOUBLE sign(DOUBLE f) {
            return f < 0 ? -1 : f > 0 ? 1 : 0;
        }

        inline DOUBLE lerp(const DOUBLE v0, const DOUBLE v1, const DOUBLE f) {
            return (v1 - v0) * f + v0;
        }

        inline DOUBLE pow(const DOUBLE v, const DOUBLE b) {
            return ::pow(v, b);
        }

        const DOUBLE PI = 3.1415927;
        const DOUBLE PI2 = 2.0 * 3.1415927;
    }
    using namespace Math;
    //-----------------

    //-----------------
    namespace IO {

        //-----------------
        namespace BasicImageIO {
#define MAXIMAGEWIDTH 4096
#define MAXIMAGEHEIGHT 4096
            unsigned int pictureS[MAXIMAGEWIDTH * MAXIMAGEHEIGHT];
            unsigned int pictureWriteOut[MAXIMAGEWIDTH * MAXIMAGEHEIGHT];
            int pictureWidth;
            int pictureHeight;
            bool loadJPG(const std::string& jpg);
            bool loadPNG(const char* png);
            void savePNG(const char* file_name);
        }
        using namespace BasicImageIO;
        //-----------------

        //-----------------
        namespace FileFunctions {

            std::string extension(const std::string& fileName) {
                std::string r;
                for (int i = ArraySize(fileName) - 1; i >= 0; --i) {
                    if (fileName[i] == '.') break;
                    r = fileName[i] + r;
                }
                return r;
            }
#define MAXDIRECTORIES 32
#define COLLECTEDDIRECTORIESMAXFILENAMESIZE 8000
            char collectedDirectoriesN[MAXDIRECTORIES][COLLECTEDDIRECTORIESMAXFILENAMESIZE];
            int collectedDirectoriesCount = 0;
            const char* collectedDirectories[MAXDIRECTORIES];

            void collectDirectories(const std::string& folder) {
                collectedDirectoriesCount = 0;
                for (const auto& entry : FILESYSTEM::directory_iterator(folder)) {
                    if (collectedDirectoriesCount >= MAXDIRECTORIES)
                        break;
                    if (entry.is_directory()) {
                        std::string yeah = entry.path().string();
                        memcpy(collectedDirectoriesN[collectedDirectoriesCount], yeah.c_str(), ArraySize(yeah));
                        collectedDirectoriesN[collectedDirectoriesCount][ArraySize(yeah)] = 0;
                        collectedDirectories[collectedDirectoriesCount] = collectedDirectoriesN[collectedDirectoriesCount];
                        collectedDirectoriesCount++;
                    }
                }
            }
        }
        using namespace FileFunctions;
        //-----------------
    }
    using namespace ESBIFT::IO;
    //-----------------

    //-----------------
    namespace Graphical {

        //-----------------
        namespace ColorSpace {

            //-----------------
            namespace RGBA {

                //...
                class PIXEL {
                public:
                    ICHANNEL r = 0, g = 0, b = 0, a = 0;
                };

#define LUMINANCEFACTORR 0.25
#define LUMINANCEFACTORG 0.6
#define LUMINANCEFACTORB 0.15
                inline DOUBLE luminance(DOUBLE r, DOUBLE g, DOUBLE b) {
                    return r * LUMINANCEFACTORR + g * LUMINANCEFACTORG + b * LUMINANCEFACTORB;
                }

                void normalize(PIXEL& v, const PIXEL& mincolor, const PIXEL& maxcolor) {
                    v.r = (v.r - mincolor.r) / FIXUPDIVBYZERO(maxcolor.r - mincolor.r, 1.0);
                    v.g = (v.g - mincolor.g) / FIXUPDIVBYZERO(maxcolor.g - mincolor.g, 1.0);
                    v.b = (v.b - mincolor.b) / FIXUPDIVBYZERO(maxcolor.b - mincolor.b, 1.0);
                    v.a = (v.a - mincolor.a) / FIXUPDIVBYZERO(maxcolor.a - mincolor.a, 1.0);
                }

                void normalizeMax(PIXEL& v, const PIXEL& mincolor, const PIXEL& maxcolor) {
                    const double maxChannel = FIXUPDIVBYZERO(MAX(MAX(MAX(maxcolor.r, maxcolor.g), maxcolor.b), maxcolor.a), 1.0);
                    v.r = v.r / maxChannel;
                    v.g = v.g / maxChannel;
                    v.b = v.b / maxChannel;
                    v.a = v.a / maxChannel;
                }

                void saturate(PIXEL& v) {
                    v.r = Math::clamp(v.r, 0.0, 1.0);
                    v.g = Math::clamp(v.g, 0.0, 1.0);
                    v.b = Math::clamp(v.b, 0.0, 1.0);
                    v.a = Math::clamp(v.a, 0.0, 1.0);
                }

                DOUBLE intensity(const PIXEL& v) {
                    return v.r * 0.3 + v.g * 0.6 + v.b * 0.1;
                }

                void madd(PIXEL& v, const PIXEL& m, const DOUBLE scale) {
                    v.r += m.r * scale;
                    v.g += m.g * scale;
                    v.b += m.b * scale;
                    v.a += m.a * scale;
                }

                void add(PIXEL& v, const PIXEL& m) {
                    v.r += m.r;
                    v.g += m.g;
                    v.b += m.b;
                    v.a += m.a;
                }

                void sub(PIXEL& v, const PIXEL& m) {
                    v.r -= m.r;
                    v.g -= m.g;
                    v.b -= m.b;
                    v.a -= m.a;
                }

                void add(PIXEL& v, const DOUBLE m) {
                    v.r += m;
                    v.g += m;
                    v.b += m;
                    v.a += m;
                }

                void alpha(PIXEL& v, const PIXEL& m) {
                    v.r = (m.r - v.r) * m.a + v.r;
                    v.g = (m.g - v.g) * m.a + v.g;
                    v.b = (m.b - v.b) * m.a + v.b;
                    v.a = (m.a - v.a) * m.a + v.a;
                }

                void div(PIXEL& v, const DOUBLE m) {
                    v.r /= m;
                    v.g /= m;
                    v.b /= m;
                    v.a /= m;
                }

                void div(PIXEL& v, const PIXEL &m) {
                    v.r /= m.r;
                    v.g /= m.g;
                    v.b /= m.b;
                    v.a /= m.a;
                }

                void mul(PIXEL& v, const DOUBLE m) {
                    v.r *= m;
                    v.g *= m;
                    v.b *= m;
                    v.a *= m;
                }

                void pow(PIXEL& v, const DOUBLE m) {
                    v.r = ::pow(v.r, m);
                    v.g = ::pow(v.g, m);
                    v.b = ::pow(v.b, m);
                    v.a = ::pow(v.a, m);
                }

                void base(PIXEL& v, const DOUBLE m) {
                    v.r = ::pow(m, v.r);
                    v.g = ::pow(m, v.g);
                    v.b = ::pow(m, v.b);
                    v.a = ::pow(m, v.a);
                }

                void clamp(PIXEL& v, DOUBLE minval, DOUBLE maxval) {
                    v.r = Math::clamp(v.r, minval, maxval);
                    v.g = Math::clamp(v.g, minval, maxval);
                    v.b = Math::clamp(v.b, minval, maxval);
                    v.a = Math::clamp(v.a, minval, maxval);
                }

                PIXEL lerp(const PIXEL& v0, const PIXEL& v1, const DOUBLE m) {
                    PIXEL v;
                    v.r = Math::lerp(v0.r, v1.r, m);
                    v.g = Math::lerp(v0.g, v1.g, m);
                    v.b = Math::lerp(v0.b, v1.b, m);
                    v.a = Math::lerp(v0.a, v1.a, m);
                    return v;
                }

                unsigned int toRGBA32(const PIXEL& pixel) {
                    unsigned int r;
                    r = ((unsigned int)Math::clamp((pixel.r * 255.0), 0, 255));
                    r |= ((unsigned int)Math::clamp((pixel.g * 255.0), 0, 255)) << 8;
                    r |= ((unsigned int)Math::clamp((pixel.b * 255.0), 0, 255)) << 16;
                    r |= ((unsigned int)Math::clamp((pixel.a * 255.0), 0, 255)) << 24;
                    return r;
                }

                PIXEL fromRGBA32(const unsigned int rgba32) {
                    PIXEL r;
                    r.r = (DOUBLE)(rgba32 & 255) / 255.0;
                    r.g = (DOUBLE)((rgba32 >> 8) & 255) / 255.0;
                    r.b = (DOUBLE)((rgba32 >> 16) & 255) / 255.0;
                    r.a = (DOUBLE)((rgba32 >> 24) & 255) / 255.0;
                    return r;
                }

#define PIXELCOLORMINMAX_RGBA(__l,__mincolor,__maxcolor,__first) {\
                    const PIXEL &value = __l;\
                    if (value.r < __mincolor.r || __first) __mincolor.r = value.r;\
                    if (value.g < __mincolor.g || __first) __mincolor.g = value.g;\
                    if (value.b < __mincolor.b || __first) __mincolor.b = value.b;\
                    if (value.a < __mincolor.a || __first) __mincolor.a = value.a;\
                    if (value.r > __maxcolor.r || __first) __maxcolor.r = value.r;\
                    if (value.g > __maxcolor.g || __first) __maxcolor.g = value.g;\
                    if (value.b > __maxcolor.b || __first) __maxcolor.b = value.b;\
                    if (value.a > __maxcolor.a || __first) __maxcolor.a = value.a;\
                    __first = false;\
                }
            }
            using namespace RGBA;
            //-----------------
#define PIXELCOLORMINMAX PIXELCOLORMINMAX_RGBA

        }
        using namespace ColorSpace;
        //-----------------

        //...
        class PIXELLayer {
        public:
            PIXELLayer() {
                uniqueId++;
            }
            int width, height, level;
            Array(PIXEL) pixels;
            UNIQUEID _id = uniqueId;
        };

        //...
        class ImageConfig {
        public:
            bool normalize = false;
            bool normalizeMax = false;
            bool alphaOnly = false;
            DOUBLE add = 0.0;
            DOUBLE mul = 1.0;
            DOUBLE zoom = 1.0;
            int desiredWidth = -1;
            int borderSize = 0;
            bool forceUpdate = false;
            double imageAlpha = 1.0;
            double imageScale = 1.0;
            double imageRotation = 0.0;
            int onlySingleChannel = -1;
            PIXEL color = fromRGBA32(0xffffffff);
            int frameId = 0;
        };

        //-----------------
        namespace HAL {

            //-----------------
            namespace OpenGL {

#define MAXGLIDCACHESIZE 50
                using HANDLE = unsigned int;

                using HALTEXTURECACHEENTRY = std::pair<UNIQUEID, unsigned int>;
                using HALTEXTURECACHE = std::map<UNIQUEID, HALTEXTURECACHEENTRY>;
                HALTEXTURECACHE halCache;

                int upload(const PIXELLayer& layer, const ImageConfig& config) {
                    PIXEL mincolor, maxcolor; bool first = true;
                    for (int i = 0; i < ArraySize(layer.pixels); ++i) {
                        PIXELCOLORMINMAX(layer.pixels[i], mincolor, maxcolor, first)
                    }
                    unsigned int* rgbaData;
                    rgbaData = new unsigned int[layer.width * layer.height];
                    unsigned int id;
                    glGenTextures(1, &id);
                    glBindTexture(GL_TEXTURE_2D, id);
                    Array(unsigned int) data;
                    for (int y = 0; y < layer.height; ++y) {
                        for (int x = 0; x < layer.width; ++x) {
                            PIXEL a = layer.pixels[x + y * layer.width];
                            if (config.normalize) {
                                if (config.normalizeMax)
                                    normalizeMax(a, mincolor, maxcolor);
                                else
                                    normalize(a, mincolor, maxcolor);
                            }
                            rgbaData[x + y * layer.width] = toRGBA32(a);
                            if (config.alphaOnly) {
                                unsigned int a = (rgbaData[x + y * layer.width] >> 24) & 255;
                                rgbaData[x + y * layer.width] = a | (a << 8) | (a << 16) | 0xff000000;
                            }
                        }
                    }
                    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, layer.width, layer.height, 0, GL_RGBA, GL_UNSIGNED_BYTE, rgbaData);
                    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
                    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
                    glBindTexture(GL_TEXTURE_2D, 0x00);
                    delete[] rgbaData;
                    return id;
                }

                bool isInCache(UNIQUEID _id) {
                    return halCache.find(_id) != halCache.end();
                }

                int cachedUpload(const PIXELLayer& layer, const ImageConfig& config) {
                    HALTEXTURECACHE::iterator f = halCache.find(layer._id);
                    if (f != halCache.end() && (!config.forceUpdate)) {
                        f->second.second++;
                        return f->second.first;
                    }

                    // if cache size exceeded remove one element (with least usage)
                    if (ArraySize(halCache) > MAXGLIDCACHESIZE) {
                        HALTEXTURECACHE::iterator start = halCache.begin();
                        HALTEXTURECACHE::iterator end = halCache.end();
                        HALTEXTURECACHE::iterator found = start;
                        HALTEXTURECACHE::iterator it;
                        int minHits = 0;
                        for (it = start; it != end; ++it) {
                            int hits = it->second.second;
                            if (hits < minHits || it == start) {
                                minHits = hits;
                                found = it;
                            }
                        }
                        unsigned int texId = found->second.first;
                        glDeleteTextures(1, &texId);
                        halCache.erase(found);
                    }

                    unsigned int id = upload(layer, config);
                    static int freshness = 0;
                    halCache[layer._id] = HALTEXTURECACHEENTRY(id, freshness++); // i am no cache specialist
                    return id;
                }

                void clearCache() {
                    HALTEXTURECACHE::iterator start = halCache.begin();
                    HALTEXTURECACHE::iterator end = halCache.end();
                    HALTEXTURECACHE::iterator it;
                    for (it = start; it != end; ++it) {
                        unsigned int texId = it->second.first;
                        glDeleteTextures(1, &texId);
                    }
                    halCache.clear();
                }

                //...
                class DISPLAYLayer {
                public:
                    int width, height, level;
                    HANDLE handle;
                    int uniqueIdOfBaseObject;
                };
            }
            using namespace OpenGL;
            //-----------------

        }
        using namespace HAL;
        //-----------------

        //-----------------
        namespace Processing {
            PIXELLayer add(const PIXELLayer& v, const PIXELLayer& m, DOUBLE scale) {
                PIXELLayer r = v;
                for (int i = 0; i < v.width * v.height; ++i) {
                    ColorSpace::madd(r.pixels[i], m.pixels[i], scale);
                }
                return r;
            }

            PIXELLayer mul(const PIXELLayer& v, DOUBLE scale) {
                PIXELLayer r = v;
                for (int i = 0; i < v.width * v.height; ++i) {
                    mul(r.pixels[i], scale);
                }
                return r;
            }

            PIXELLayer pow(const PIXELLayer& v, DOUBLE power) {
                PIXELLayer r = v;
                for (int i = 0; i < v.width * v.height; ++i) {
                    pow(r.pixels[i], power);
                }
                return r;
            }

            PIXELLayer base(const PIXELLayer& v, DOUBLE basep) {
                PIXELLayer r = v;
                for (int i = 0; i < v.width * v.height; ++i) {
                    base(r.pixels[i], basep);
                }
                return r;
            }

            PIXELLayer addNoise(const PIXELLayer& v, DOUBLE noiseIntensity) {
                PIXELLayer r = v;
                for (int i = 0; i < v.width * v.height; ++i) {
                    r.pixels[i].r += random(noiseIntensity);
                    r.pixels[i].g += random(noiseIntensity);
                    r.pixels[i].b += random(noiseIntensity);
                    r.pixels[i].a += random(noiseIntensity);
                }
                return r;
            }

            PIXELLayer modification(const PIXELLayer& v, const std::function<PIXEL(const PIXEL&)> &t) {
                PIXELLayer r = v;
                for (int i = 0; i < v.width * v.height; ++i) {
                    r.pixels[i] = t(v.pixels[i]);
                }
                return r;
            }

            PIXELLayer normalize(const PIXELLayer& v) {
                PIXELLayer r = v;
                PIXEL mincolor, maxcolor;
                bool first = true;
                for (int i = 0; i < ArraySize(v.pixels); ++i) {
                    PIXELCOLORMINMAX(v.pixels[i], mincolor, maxcolor, first)
                }
                for (int i = 0; i < v.width * v.height; ++i) {
                    ColorSpace::normalize(r.pixels[i], mincolor, maxcolor);
                }
                return r;
            }

            PIXELLayer boxFilter(const PIXELLayer& v, int boxWidth, int boxHeight) {
                PIXELLayer r = v;
                for (int y = 0; y < v.height; ++y) {
                    for (int x = 0; x < v.width; ++x) {
                        PIXEL c;
                        DOUBLE w = 0.0;
                        for (int y2 = -boxHeight; y2 <= boxHeight; ++y2) {
                            for (int x2 = -boxWidth; x2 <= boxWidth; ++x2) {
                                const int x3 = x2 + x;
                                const int y3 = y2 + y;
                                if ((unsigned int)x3 < v.width && (unsigned int)y3 < v.height) {
                                    PIXEL d = v.pixels[x3 + y3 * v.width];
                                    ColorSpace::add(c, d);
                                    w += 1.0;
                                }
                            }
                        }
                        if (w != 0.0)
                            div(c, w);
                        r.pixels[x + y * r.width] = c;
                    }
                }
                return r;
            }

            void doAveragingFilter(const PIXELLayer& source, PIXELLayer& dest) {
                if (dest.height <= 0 || dest.width <= 0)
                    return;
                if (dest.width > source.width || dest.height > source.height) {
                    dest = source;
                    return;
                }
                std::fill(dest.pixels.begin(), dest.pixels.end(), PIXEL());
                for (int y = 0; y < source.height; ++y) {
                    const int ny = y * dest.height / source.height;
                    for (int x = 0; x < source.width; ++x) {
                        const int nx = x * dest.width / source.width;
                        PIXEL b = source.pixels[x + y * source.width];
                        b.a = 1.0;
                        ColorSpace::add(dest.pixels[nx + ny * dest.width], b);
                    }
                }

                for (int i = 0; i < dest.width * dest.height; ++i) {
                    div(dest.pixels[i], dest.pixels[i].a);
                }
            }
        }
        using namespace Processing;
        //-----------------

        //-----------------
        namespace IO {
            ImageConfig loadConfig(int borderSize = 0, int imageWidth = -1, DOUBLE brightness = 1.0, DOUBLE contrast = 1.0, DOUBLE scale = 1.0) {
                ImageConfig r;
                r.mul = contrast;
                r.add = brightness;
                r.desiredWidth = imageWidth;
                r.borderSize = borderSize;
                r.imageScale = scale;
                return r;
            }

            PIXELLayer pixelLayerFromData(const ImageConfig& config) {
                const int borderSize = config.borderSize;
                const DOUBLE brightness = config.add;
                const DOUBLE contrast = config.mul;
                int desiredWidth = config.desiredWidth;
                PIXELLayer l;
                if (desiredWidth > 1 || config.imageScale != 1.0) {
                    if (desiredWidth < 1) // just image scale
                        desiredWidth = pictureWidth;
                    PIXELLayer l1, l2;
                    l1.width = pictureWidth;
                    l1.height = pictureHeight;
                    l2.width = desiredWidth * config.imageScale;
                    l2.height = desiredWidth * pictureHeight / pictureWidth * config.imageScale;
                    ArrayResize(l1.pixels, l1.width * l1.height);
                    ArrayResize(l2.pixels, l2.width * l2.height);
                    for (int y = 0; y < pictureHeight; ++y) {
                        for (int x = 0; x < pictureWidth; ++x) {
                            l1.pixels[x + y * l1.width] = fromRGBA32(pictureS[x + y * pictureWidth]);
                        }
                    }
                    doAveragingFilter(l1, l2);
                    pictureWidth = l2.width;
                    pictureHeight = l2.height;
                    for (int y = 0; y < l2.height; ++y) {
                        for (int x = 0; x < l2.width; ++x) {
                            pictureS[x + y * pictureWidth] = toRGBA32(l2.pixels[x + y * l2.width]);
                        }
                    }
                }
                l.width = pictureWidth + config.borderSize * 2;
                l.height = pictureHeight + config.borderSize * 2;
                ArrayResize(l.pixels,l.width * l.height);
                for (int y = 0; y < pictureHeight; ++y) {
                    for (int x = 0; x < pictureWidth; ++x) {
                        l.pixels[x + y * l.width + borderSize + borderSize * l.width] = fromRGBA32(pictureS[x + y * pictureWidth]);
                        mul(l.pixels[x + y * l.width + borderSize + borderSize * l.width], contrast);
                        add(l.pixels[x + y * l.width + borderSize + borderSize * l.width], brightness);
                        clamp(l.pixels[x + y * l.width + borderSize + borderSize * l.width], 0.0, 1.0);
                    }
                }
                return l;
            }

            PIXELLayer pixelLayerFromJPG(const std::string& jpgFileName, const ImageConfig &config = ImageConfig()) {
                loadJPG(jpgFileName);
                return pixelLayerFromData(config);
            }

            PIXELLayer pixelLayerFromPNG(const std::string& pngFileName, const ImageConfig& config = ImageConfig()) {
                loadPNG(pngFileName.c_str());
                return pixelLayerFromData(config);
            }

            void savePixelLayer(const PIXELLayer& layer, const std::string& pngFileName) {
                PIXELLayer b = layer;
                pictureWidth = b.width;
                pictureHeight = b.height;
                if (pictureWidth > 0 && pictureHeight > 0) {
                    for (int y = 0; y < pictureHeight; ++y) {
                        for (int x = 0; x < pictureWidth; ++x) {
                            pictureWriteOut[x + y * pictureWidth] = toRGBA32(b.pixels[x + y * b.width]);
                        }
                    }
                    savePNG(pngFileName.c_str());
                }
            }
        }
        using namespace IO;
        //-----------------
    }
    using namespace Graphical;
    //-----------------
    
    //-----------------
    namespace Perceptional {

#if DEFINED_CHANNELCOUNT != 1
        using CHANNELS = struct { SINGLECHANNEL channel[DEFINED_CHANNELCOUNT]; };
#else
        using CHANNELS = SINGLECHANNEL;
#endif

#if DEFINED_CHANNELCOUNT == 4
        inline int channelElementCount() {return 4;}
        inline DOUBLE channel(const CHANNELS& channels, const int c) {
                return channels.channel[c];
        }

        DOUBLE intensity(const CHANNELS& channels) {
            return luminance(channels.channel[0], channels.channel[1], channels.channel[2]);
        }

        CHANNELS lerp(const CHANNELS& v0, const CHANNELS& v1, const DOUBLE m) {
            CHANNELS v;
            v.channel[0] = Math::lerp(v0.channel[0], v1.channel[0], m);
            v.channel[1] = Math::lerp(v0.channel[1], v1.channel[1], m);
            v.channel[2] = Math::lerp(v0.channel[2], v1.channel[2], m);
            v.channel[3] = Math::lerp(v0.channel[3], v1.channel[3], m);
            return v;
        }

        void standardSrcToDst(PIXEL& dst, const CHANNELS& src) {
            dst.r = src.channel[0];
            dst.g = src.channel[1];
            dst.b = src.channel[2];
            dst.a = src.channel[3];
        }
#else
        inline int channelElementCount() {
            return 1;
        }

        inline DOUBLE channel(const CHANNELS& channels, const int c) {
            return channels;
        }

        DOUBLE intensity(const CHANNELS& channels) {
            return channels;
        }

        void standardSrcToDst(PIXEL& dst, const CHANNELS& src) {
            dst.r = src;
            dst.g = src;
            dst.b = src;
            dst.a = src;
        }
#endif

        //-----------------
        namespace Descriptors {

            //...
            class DescriptorPair {
            public:
                DOUBLE x0, y0;
                DOUBLE x1, y1;
            };

            //...
            class DescriptorDescLevel {
            public:
                Array(DescriptorPair) pairs; // vector because of possibly more pairs on important points
            };

            //...
            class DescriptorDesc {
            public:
                Array(DescriptorDescLevel) levels;
            };

            //...
            class Descriptor {
            public:
                Array(Array(VECTORBOOL)) bits;
            };

            //-----------------
            namespace Shapes {
#define SHAPE(__x0,__y0,__x1,__y1) {pair.x0 = (__x0); pair.y0 = (__y0); pair.x1 = (__x1); pair.y1 = (__y1);}
#define SHAPEFUNCTION(__v) void __v(DescriptorPair& pair)
#define POLARANGLEDISTANCE1DISTANCE2 \
                pair.x0 = sin(angle) * distance1; \
                pair.y0 = cos(angle) * distance1; \
                pair.x1 = sin(angle) * distance2; \
                pair.y1 = cos(angle) * distance2;
#define POLARANGLEANDANGLEADDDISTANCE \
                pair.x0 = sin(angle) * distance; \
                pair.y0 = cos(angle) * distance; \
                pair.x1 = sin(angle + angleAdd) * distance; \
                pair.y1 = cos(angle + angleAdd) * distance;
#define TUNNELNEARPLANE (DEFINED_ADVANCED_PERSPECTIVE ? 0.01 : 1.0)
#define TUNNELFARPLANE 16.0
#define nDist (random(2.0)-1.0)
#define rSign (random(1.0)<0.5 ? -1 : 1)
#define randOne (random(1.0))
#define randPI2 (random(PI2))
#define defaultPersp (1.0 / (random(TUNNELFARPLANE-TUNNELNEARPLANE)+TUNNELNEARPLANE))
                SHAPEFUNCTION(tunnelPoint) {const double distance1 = defaultPersp; const double distance2 = random(distance1); const double angle = randPI2; const double twirl = random(0.1); const double angle1 = angle + twirl * distance1; const double angle2 = angle + twirl * distance2; pair.x0 = sin(angle1) * distance1; pair.y0 = cos(angle1) * distance1; pair.x1 = sin(angle2) * distance2; pair.y1 = cos(angle2) * distance2;}
                SHAPEFUNCTION(tunnelX) {const double x1 = defaultPersp * rSign; const double x2 = random(fabs(x1))*sign(x1); SHAPE(x1,nDist,x2, nDist)}
                SHAPEFUNCTION(tunnelY) {const double y1 = defaultPersp * rSign; const double y2 = random(fabs(y1))*sign(y1); SHAPE(nDist,y1, nDist,y2)}
                SHAPEFUNCTION(normQuad) SHAPE(nDist, nDist, nDist, nDist)
                SHAPEFUNCTION(normCircle) {while (true) {pair.x0 = nDist; pair.y0 = nDist; const double d0 = sqrt(pair.x0 * pair.x0 + pair.y0 * pair.y0); if (d0 <= 1.0) break;} while (true) {pair.x1 = nDist; pair.y1 = nDist; const double d1 = sqrt(pair.x1 * pair.x1 + pair.y1 * pair.y1); if (d1 <= 1.0) break;}}
                SHAPEFUNCTION(polarCircle) { double d0 = randOne; double d1 = randOne; double p0 = randPI2; double p1 = randPI2; SHAPE(sin(p0) * d0, cos(p0) * d0, sin(p1) * d1, cos(p1) * d1)}
                SHAPEFUNCTION(tubularX) SHAPE(Math::pow(randOne, TUBULARPOW) * rSign, Math::pow(randOne, TUBULARPOW) * rSign, nDist, nDist)
                SHAPEFUNCTION(tubularY) SHAPE(nDist, nDist, Math::pow(random(1.0), TUBULARPOW) * rSign, Math::pow(random(1.0), TUBULARPOW) * rSign)
                SHAPEFUNCTION(polarCircleQuadratic) {double d0 = randOne; double d1 = randOne; d0 *= d0; d1 *= d1; double p0 = randPI2; double p1 = randPI2; SHAPE(sin(p0) * d0, cos(p0) * d0, sin(p1) * d1, cos(p1) * d1)}
                SHAPEFUNCTION(manhattenSimple) { {const double d0 = nDist; const double d1 = randOne; pair.x0 = d0 * (1.0-d1); pair.y0 = d1;} {const double d0 = nDist; const double d1 = randOne; pair.x1 = d0 * (1.0-d1); pair.y1 = -d1;}}
                SHAPEFUNCTION(manhattenQuadratic) { {double d0 = nDist; double d1 = randOne; d0 *= fabs(d0); d1 *= d1; pair.x0 = d0 * (1.0-d1); pair.y0 = d1;} {double d0 = nDist; double d1 = randOne; d0 *= fabs(d0); d1 *= d1; pair.x1 = d0 * (1.0-d1); pair.y1 = -d1;}}
                SHAPEFUNCTION(manhattenSquareRoot) { {double d0 = nDist; double d1 = randOne; d0 = sqrt(fabs(d0)) * sign(d0); d1 = sqrt(d1); pair.x0 = d0 * (1.0-d1); pair.y0 = d1;} {double d0 = nDist; double d1 = randOne; d0 = sqrt(fabs(d0)) * sign(d0); d1 = sqrt(d1); pair.x1 = d0 * (1.0-d1); pair.y1 = -d1;}}
                SHAPEFUNCTION(distanceStar) {const double angle = randPI2; const double distance1 = randOne; const double distance2 = random(distance1); POLARANGLEDISTANCE1DISTANCE2}
                SHAPEFUNCTION(distanceStarLens) {const double angle = randPI2; const double distance1 = sqrt(randOne); const double distance2 = random(distance1); POLARANGLEDISTANCE1DISTANCE2}
                SHAPEFUNCTION(distanceStarQuadratic) {const double angle = randPI2; const double distance1 = pow(randOne,2.0); const double distance2 = random(distance1); POLARANGLEDISTANCE1DISTANCE2}
                SHAPEFUNCTION(distanceStarTunnel) {const double angle = randPI2; const double distance1 = 1.0 / random(TUNNELFARPLANE); const double distance2 = random(distance1); POLARANGLEDISTANCE1DISTANCE2}
                SHAPEFUNCTION(angleCircle) {const double angle = randPI2; const double distance = randOne; const double angleAdd = random(PI*0.5); POLARANGLEANDANGLEADDDISTANCE}
                SHAPEFUNCTION(angleLens) {const double angle = randPI2; const double distance = sqrt(randOne); const double angleAdd = random(PI*0.5); POLARANGLEANDANGLEADDDISTANCE}
                SHAPEFUNCTION(angleTunnel) {const double angle = randPI2; const double distance = defaultPersp; const double angleAdd = random(PI*0.5); POLARANGLEANDANGLEADDDISTANCE}
                SHAPEFUNCTION(descriptorPairRightDown) { if (pair.x0 > pair.x1) { const double tx = pair.x0; pair.x0 = pair.x1; pair.x1 = tx;} if (pair.y0 > pair.y1) { const double ty = pair.y0; pair.y0 = pair.y1; pair.y1 = ty;}}
                SHAPEFUNCTION(descriptorPairNoSort) {};
                SHAPEFUNCTION(toPolar) {double d0 = sqrt(pair.x0 * pair.x0 + pair.y0 * pair.y0); double d1 = sqrt(pair.x1 * pair.x1 + pair.y1 * pair.y1); double a0 = atan2(pair.x0, pair.y0); double a1 = atan2(pair.x1, pair.y1); SHAPE(a0,d0,a1,d1)}
                SHAPEFUNCTION(fromPolar) { DescriptorPair p = pair; SHAPE(sin(p.x0) * p.y0, cos(p.x0) * p.y0, sin(p.x1) * p.y1, cos(p.x1) * p.y1);}
                SHAPEFUNCTION(polarCircleSquareRoot) { double d0 = randOne; double d1 = randOne; d0 = sqrt(d0); d1 = sqrt(d1); const double p0 = randPI2; const double p1 = randPI2; SHAPE(sin(p0) * d0, cos(p0) * d0, sin(p1) * d1, cos(p1) * d1)}
                std::vector<std::function<SHAPEFUNCTION()>> nonRotationalDescriptors = NONROTATIONALDESCRIPTORS;
                std::vector<std::function<SHAPEFUNCTION()>> rotationalDescriptors = ROTATIONALDESCRIPTORS;
            }
            using namespace Shapes;
            //-----------------
        }
        using namespace Descriptors;
        //-----------------

        //-----------------
        namespace Features {

            //...
            class Feature {
            public:
                Feature() {;}
                Feature(const DOUBLE x, const DOUBLE y, const Feature &b) {
                    this->x = x;
                    this->y = y;
                    // :(          V
                    this->strength = b.strength;
                    this->scale = b.scale;
                    this->rotation = b.rotation;
                    this->force = b.force;
                    this->frame = b.frame;
                    this->id = b.id;
                    // :(          A
                }
                Feature(const DOUBLE x, const DOUBLE y)  {
                    this->x = x;
                    this->y = y;
                    this->force = 1; // notifying a working feature (hopefully)
                    this->id = uniqueId;
                    uniqueId++; // todo: thread safe
                }
                DOUBLE x = 0, y = 0;
                DOUBLE strength = 0; // used for removing
                DOUBLE scale = 0; // not used yet
                DOUBLE rotation = 0; // not used yet
                DOUBLE force = 0; // strength of gradients on that point
                int frame = 0; // not used yet
                int id = 0;
            };

            using FeatureMap = std::map<int, Feature>;

            //...
            class FEATURELayer {
            public:
                FEATURELayer() {
                    uniqueId++;
                }
                int width, height, level;
                double fwidth, fheight;
                Array(CHANNELS) pixels;
                int _id = uniqueId;
            };
        }
        using namespace Features;
        //-----------------

        //...
        class FeatureConfig {
        public:
            Array(DescriptorDesc) descriptorPatterns;
            int deepestLevel = 0;
            int lastLevel = 0;
            int sampleLevel = 0;
            bool scaleInvariance = true;
            int coordReferenceLevel = 0;
            double refinementSteps = Configuration::refinementSteps;
            bool iterative = Configuration::iterative;
            double stepFactor = 1.0;
            DOUBLE refinementStepPerLevelScaleFactor = Configuration::refinementStepPerLevelScaleFactor;
            DOUBLE descriptorSizeCoarse = 0.5;
            DOUBLE descriptorSizeFine = 8.0 / 640.0 / 2.0; // n pixels
            double descriptorLevelScale = Configuration::descriptorLevelScale;
            double sourceToDestScaleRatio = 1.0;
            int scaleLevelChecksCount = 50;
            double descriptorShapeSize = 1.0;
        };

        //-----------------
        namespace Implementation {
            FeatureConfig baseConfig;
            FeatureConfig resampleConfig;
        }
        using namespace Implementation;
        //-----------------

        //-----------------
        namespace Algorithms {

            const double distance(const Feature &a, const Feature &b) {
                const double dx = b.x - a.x;
                const double dy = b.y - a.y;
                return sqrt(dx*dx + dy*dy);
            }

            double descriptorRatio(const Descriptor& v1, const Descriptor& v2) {
                if (ArraySize(v1.bits) != ArraySize(v2.bits)) return -1;
                int k = 0;
                int v = 0;
                for (int i = 0; i < ArraySize(v1.bits); ++i) {
                    if (ArraySize(v1.bits[i]) != ArraySize(v2.bits[i])) 
                        return -1;
                    for (int j = 0; j < ArraySize(v1.bits[i]); ++j) {
                        if (v1.bits[i][j] == v2.bits[i][j]) v++;
                        k++;
                    }
                }
                return (double)v / FIXUPDIVBYZERO(k,1.0);
            }

            inline const DOUBLE descriptorSizeForLevel(double level, const FeatureConfig& config) {
                const double coarseLevel = config.deepestLevel;
                const double base = 1.0 / config.descriptorLevelScale;
                const double a1 = log(config.descriptorSizeFine) / log(base);
                const double a2 = log(config.descriptorSizeCoarse) / log(base);
                const double easedValue = pow(level / (coarseLevel - 1), quadraticDescriptorGrip);
                const double toMipLevel = pow(base, Math::lerp(a1, a2, easedValue));
                return toMipLevel;
            }

#define DESCRIPTORSIZE(__imageLevel) (descriptorSizeForLevel(__imageLevel, config) * config.descriptorShapeSize)
#define DESCRIPTORLEVELCOUNT (config.deepestLevel + 1)
#define RANDOMSCALEINVARIANCE(__baseSize) ((randomLike(randomLikeBase+3333,maxDescriptorSizeRandom - minDescriptorSizeRandom) + minDescriptorSizeRandom) * (__baseSize))

            PIXELLayer toPixelLayer(const FEATURELayer& v, const ImageConfig &config) {
                PIXELLayer r;
                JUSTSETUPIMAGEHULL(r, v);
                for (int i = 0; i < ArraySize(v.pixels); ++i) {
#if DEFINED_CHANNELCOUNT != 1
                    if (config.onlySingleChannel != -1) {
                        CHANNELS b;
                        b.channel[0] = v.pixels[i].channel[config.onlySingleChannel];
                        b.channel[1] = v.pixels[i].channel[config.onlySingleChannel];
                        b.channel[2] = v.pixels[i].channel[config.onlySingleChannel];
                        b.channel[3] = 1.0;
                        standardSrcToDst(r.pixels[i], b);
                    } else
                        standardSrcToDst(r.pixels[i], v.pixels[i]);
#else
                    standardSrcToDst(r.pixels[i], v.pixels[i]);
#endif
                }
                return r;
            }

            //-----------------
            namespace Access {

#define MARKREMOVE(__f) __f.strength = -1;
#define SHOULDREMOVE(__f) (__f.strength == -1)
                void pruneByGradientForce(Array(Feature) &features, double threshold) {
                    for (int i = features.size() - 1; i >= 0; --i) {
                        if (features[i].force < threshold) {
                            MARKREMOVE(features[i]);
                        }
                    }
                }

                FeatureMap featureMap(const std::vector<Feature> &features) {
                    FeatureMap r;
                    for (int i = 0; i < ArraySize(features); ++i) {
                        r.insert(std::make_pair(features[i].id, features[i]));
                    }
                    return r;
                }

                const Feature get(const FeatureMap &m, const int id, const Feature &defaultV) {
                    FeatureMap::const_iterator t = m.find(id);
                    if (t != m.end()) return t->second;
                    return defaultV;
                }

                //-----------------
                namespace FeatureNeighborhood {
                    Array(int) sorted[2];
                    Array(int) featureIdToSorted[2];
                    void createNeighborHood(const Array(Feature) &features) {
                        for (int k = 0; k < 2; ++k) {
                            ArrayResize(sorted[k], ArraySize(features));
                            ArrayResize(featureIdToSorted[k], ArraySize(features));
                            for (int i = 0; i < ArraySize(features); ++i) {
                                sorted[k][i] = i;
                            }
                            std::function<bool(const int,const int)> sort;

                            switch(k) {
                                case 0: {
                                    sort = [=](const int a, const int b) -> bool {
                                        return POSITIONTOVALUE(features[a].x, features[a].y) >
                                               POSITIONTOVALUE(features[b].x, features[b].y);
                                    };
                                } break;
                                case 1: {
                                    sort = [=](const int a, const int b) -> bool {
                                        return POSITIONTOVALUE(features[a].y, features[a].x) >
                                               POSITIONTOVALUE(features[b].y, features[b].x);
                                    };
                                } break;
                                case 2: { // not used, yet (may lead to problems in Eineindeutigkeit
                                    sort = [=](const int a, const int b) -> bool {
                                        return POSITIONTOVALUE(features[a].y, features[a].x) >
                                               POSITIONTOVALUE(features[b].x, features[b].y);
                                    };
                                } break;
                            }
                            std::sort(sorted[k].begin(), sorted[k].end(), sort);
                            for (int i = 0; i < ArraySize(features); ++i) {
                                featureIdToSorted[k][sorted[k][i]] = i;
                            }
                        }
                    }

                    int nearestNeighbor(const Array(Feature) &features, int index, const std::set<int> &excludeList, const int maxDelta = 8) {
#define SORTI(__id) if (excludeList.find(id) == excludeList.end()) {\
                            if (distance(i0, features[__id]) < lastD || lastD < 0) {\
                                lastD = distance(i0, features[__id]);\
                                r = id;\
                            }\
                        }
                        const Feature i0 = features[index];
                        int r = 0; double lastD = -1;
                        for (int i = 1; i < maxDelta; ++i) {
                            for (int j = 0; j < 2; ++j) {
                                const int base = featureIdToSorted[j][index];
                                {
                                    const int id = sorted[j][CLAMP(base + i, 0, ArraySize(features) - 1)];
                                    SORTI(id)
                                }
                                {
                                    const int id = sorted[j][CLAMP(base - i, 0, ArraySize(features) - 1)];
                                    SORTI(id)
                                }
                            }
                        }
                        return r;
                    }
                }
                using namespace FeatureNeighborhood;
                //-----------------

                inline CHANNELS sample_featurepixel_direct(const FEATURELayer& layer, int x, int y) {
                    return layer.pixels[x + y * layer.width];
                }

                inline CHANNELS sample_featurepixel_direct_clamped(const FEATURELayer& layer, DOUBLE x, DOUBLE y) {
                    const int ix = (int)clamp(trunc(x), 0, layer.width - 1);
                    const int iy = (int)clamp(trunc(y), 0, layer.height - 1);
                    return layer.pixels[ix + iy * layer.width];
                }

                inline CHANNELS sample_featurepixel_linear_2d(const FEATURELayer& layer, DOUBLE x, DOUBLE y) {
                    const int ix0 = (int)trunc(x);
                    const int iy0 = (int)trunc(y);
                    const int ix1 = ix0 + 1;
                    const int iy1 = iy0 + 1;
                    const DOUBLE nx = x - ix0;
                    const DOUBLE ny = y - iy0;
                    const CHANNELS p00 = sample_featurepixel_direct(layer, ix0, iy0);
                    const CHANNELS p10 = sample_featurepixel_direct(layer, ix1, iy0);
                    const CHANNELS p01 = sample_featurepixel_direct(layer, ix0, iy1);
                    const CHANNELS p11 = sample_featurepixel_direct(layer, ix1, iy1);
                    const CHANNELS mx0 = lerp(p00, p10, nx);
                    const CHANNELS mx1 = lerp(p01, p11, nx);
                    return lerp(mx0, mx1, ny);
                }

                inline void lookupPairFiltered(DOUBLE x0, DOUBLE y0, DOUBLE x1, DOUBLE y1, const FEATURELayer& layer, CHANNELS& c1, CHANNELS& c2) {
                    c1 = sample_featurepixel_linear_2d(layer, x0, y0);
                    c2 = sample_featurepixel_linear_2d(layer, x1, y1);
                }

                inline void lookupPairDirect(DOUBLE x0, DOUBLE y0, DOUBLE x1, DOUBLE y1, const FEATURELayer& layer, CHANNELS& c1, CHANNELS& c2) {
                    c1 = sample_featurepixel_direct(layer, x0, y0);
                    c2 = sample_featurepixel_direct(layer, x1, y1);
                }

                inline void lookupPairDirectClamped(DOUBLE x0, DOUBLE y0, DOUBLE x1, DOUBLE y1, const FEATURELayer& layer, CHANNELS& c1, CHANNELS& c2) {
                    c1 = sample_featurepixel_direct_clamped(layer, x0, y0);
                    c2 = sample_featurepixel_direct_clamped(layer, x1, y1);
                }

#if DEFINED_BORDERCLAMPING == 1
    #if DEFINED_BORDERCLAMPMIRROR == 1
            #define BORDERCLAMP(__d, __v, __l, __h) \
            DOUBLE __d = (__v);\
            FIXFLOATNUMBER(__d, 0.0)\
            if ((__d) < (__l)) __d = (__l) + ((__l) - (__d));\
            if ((__d) > (__h)) __d = (__h) - ((__d) - (__h));\
            if ((__d) < (__l)) __d = (__l); \
            if ((__d) > (__h)) __d = (__h);
    #else
            #if DEFINED_BORDERCLAMPWRAP == 1
                #define BORDERCLAMP(__d, __v, __l, __h) \
                DOUBLE __d = (__v);\
                FIXFLOATNUMBER(__d, 0.0)\
                if ((__d) < (__l)) __d = (__h) - ((__d) - (__l));\
                if ((__d) > (__h)) __d = (__l) + ((__d) - (__h));\
                if ((__d) < (__l)) __d = (__l); \
                if ((__d) > (__h)) __d = (__h);
            #else
                #define BORDERCLAMP(__d, __v, __l, __h) \
                DOUBLE __d = (__v);\
                FIXFLOATNUMBER(__d, 0.0)\
                if ((__d) < (__l)) __d = (__l);\
                if ((__d) > (__h)) __d = (__h);
            #endif
    #endif
#else
    #define BORDERCLAMP(__d, __v, __l, __h) 
            DOUBLE __d = (__v);\
            FIXFLOATNUMBER(__d, 0.0)\
            if ((__d < (__l)) || (__d > (__h))) return false; // not used anymore (not the best idea, yet) would return false in lookupDescriptor
#endif

#define ANGLEMOD(__angle) (__angle)
#if DEFINED_WITHROTATIONINVARIANCE == 1
    #if DEFINED_FULLROTATIONINVARIANCEBYSECTORS == 1
        #define ANGLEMOD(__angle) fmod(__angle,45.0)
    #endif
#endif

#define ROTATE(__x, __y, __radians) (DOUBLE)(cos(__radians) * (__x) + sin(__radians) * (__y)), (DOUBLE)(-sin(__radians) * (__x) + cos(__radians) * (__y))

#if DEFINED_NORMAL3DPATCHES == 1
#define PREPAREPATCHPOLY\
    PATCHES_RANDOM_VALUE(0) = randomLike(randomLikeBase + 1234, DEFINED_PATCH_RANDOM_DEVIATION * 2.0) - DEFINED_PATCH_RANDOM_DEVIATION;\
    PATCHES_RANDOM_VALUE(1) = randomLike(randomLikeBase + 2345, DEFINED_PATCH_RANDOM_DEVIATION * 2.0) - DEFINED_PATCH_RANDOM_DEVIATION;\
    PATCHES_RANDOM_VALUE(2) = randomLike(randomLikeBase + 3456, DEFINED_PATCH_RANDOM_DEVIATION * 2.0) - DEFINED_PATCH_RANDOM_DEVIATION;\
    PATCHES_RANDOM_VALUE(3) = randomLike(randomLikeBase + 4567, DEFINED_PATCH_RANDOM_DEVIATION * 2.0) - DEFINED_PATCH_RANDOM_DEVIATION;
#else
#define PREPAREPATCHPOLY
#endif

#if COMPILERBUGGONE == 1
#define LOOKUPDESCRIPTORPART_NORMAL3DPATCHES\
                {const DescriptorPair inputPair = __pair__;\
                __pair__.x0 *= 1.0 + PATCHES_RANDOM_VALUE(0) * inputPair.y0;\
                __pair__.y0 *= 1.0 + PATCHES_RANDOM_VALUE(1) * inputPair.x0;\
                __pair__.x1 *= 1.0 + PATCHES_RANDOM_VALUE(2) * inputPair.y1;\
                __pair__.y1 *= 1.0 + PATCHES_RANDOM_VALUE(3) * inputPair.x1;}
#else
#define LOOKUPDESCRIPTORPART_NORMAL3DPATCHES\
                {const DescriptorPair inputPair = __pair__;\
                __pair__.x0 *= 1.0 + (randomLike(randomLikeBase + j + 12234423,DEFINED_PATCH_RANDOM_DEVIATION * 2.0) - DEFINED_PATCH_RANDOM_DEVIATION) * inputPair.y0;\
                __pair__.y0 *= 1.0 + (randomLike(randomLikeBase + j + 122375,DEFINED_PATCH_RANDOM_DEVIATION * 2.0) - DEFINED_PATCH_RANDOM_DEVIATION) * inputPair.x0;\
                __pair__.x1 *= 1.0 + (randomLike(randomLikeBase + j + 12277,DEFINED_PATCH_RANDOM_DEVIATION * 2.0) - DEFINED_PATCH_RANDOM_DEVIATION) * inputPair.y1;\
                __pair__.y1 *= 1.0 + (randomLike(randomLikeBase + j + 1223759,DEFINED_PATCH_RANDOM_DEVIATION * 2.0) - DEFINED_PATCH_RANDOM_DEVIATION) * inputPair.x1;}
#endif

#define LOOKUPDESCRIPTORPART_ROUNDEDPATCHES\
                {const DOUBLE maxhypotenuse = sqrt(2.0);\
                const DOUBLE d0 = 1.0/(pow(normdistance0/maxhypotenuse,roundedPatchLensFactor));\
                const DOUBLE d1 = 1.0/(pow(normdistance1/maxhypotenuse,roundedPatchLensFactor));\
                __pair__.x0 = Math::lerp(__pair__.x0, __pair__.x0 * d0, ActiveConfiguration::roundedPatchRoundRatio);\
                __pair__.y0 = Math::lerp(__pair__.y0, __pair__.y0 * d0, ActiveConfiguration::roundedPatchRoundRatio);\
                __pair__.x1 = Math::lerp(__pair__.x1, __pair__.x1 * d1, ActiveConfiguration::roundedPatchRoundRatio);\
                __pair__.y1 = Math::lerp(__pair__.y1, __pair__.y1 * d1, ActiveConfiguration::roundedPatchRoundRatio);}
#define LOOKUPDESCRIPTORPART_FROMPOLAR(__angleRadians)\
                {const DescriptorPair inputPair = __pair__;\
                __pair__.x0 = sin(inputPair.x0 + (__angleRadians)) * inputPair.y0;\
                __pair__.y0 = cos(inputPair.x0 + (__angleRadians)) * inputPair.y0;\
                __pair__.x1 = sin(inputPair.x1 + (__angleRadians)) * inputPair.y1;\
                __pair__.y1 = cos(inputPair.x1 + (__angleRadians)) * inputPair.y1;}
#define LOOKUPDESCRIPTORPART_CLAMP(__x, __y, __layer, __descriptorSize, __ddx, __ddy, __angleRadians)\
                const DOUBLE gradients[2] = {ROTATE((__ddx),(__ddy),(__angleRadians))};\
                BORDERCLAMP(x0, ((__pair__.x0 * (__descriptorSize) + (__x)) * (__layer).fwidth + gradients[0]), 1, (__layer).fwidth - 2);\
                BORDERCLAMP(y0, ((__pair__.y0 * (__descriptorSize) + (__y)) * (__layer).fheight + gradients[1]), 1, (__layer).fheight - 2);\
                BORDERCLAMP(x1, ((__pair__.x1 * (__descriptorSize) + (__x)) * (__layer).fwidth + gradients[0]), 1, (__layer).fwidth - 2);\
                BORDERCLAMP(y1, ((__pair__.y1 * (__descriptorSize) + (__y)) * (__layer).fheight + gradients[1]), 1, (__layer).fheight - 2);
#define LOOKUPDESCRIPTORPART_LOOKUP(__layer, __c1, __c2)\
                lookupPairFiltered(x0, y0, x1, y1, (__layer), (__c1), (__c2));

#if DEFINED_POLARCOORDINATES == 0
#define LOOKUPDESCRIPTORPART_FROMPOLAR(__angleRadians)
#endif
#if DEFINED_NORMAL3DPATCHES == 0
#define LOOKUPDESCRIPTORPART_NORMAL3DPATCHES
#endif
#if DEFINED_ROUNDEDPATCHES == 0
#define LOOKUPDESCRIPTORPART_ROUNDEDPATCHES
#endif

#define LOOKUPDESCRIPTORPART_PAIRPOLARDISTANCE const double normdistance0 = sqrt(__pair__.x0*__pair__.x0 + __pair__.y0*__pair__.y0);\
                            const double normdistance1 = sqrt(__pair__.x1*__pair__.x1 + __pair__.y1*__pair__.y1);


// as macro
#define lookupDescriptorMac(__returnValue, __x, __y, __layer, __inputPair, __descriptorSize, __c1, __c2, __ddx, __ddy, __angleRadians) \
    {DescriptorPair __pair__ = __inputPair;\
        LOOKUPDESCRIPTORPART_FROMPOLAR(__angleRadians)\
        LOOKUPDESCRIPTORPART_PAIRPOLARDISTANCE\
        LOOKUPDESCRIPTORPART_NORMAL3DPATCHES\
        LOOKUPDESCRIPTORPART_ROUNDEDPATCHES\
        LOOKUPDESCRIPTORPART_CLAMP(__x, __y, __layer, __descriptorSize, __ddx, __ddy, __angleRadians)\
        LOOKUPDESCRIPTORPART_LOOKUP(__layer, __c1, __c2);\
        __returnValue = true;\
    }
#define lookupRawDescriptorMac(__returnValue, __x, __y, __layer, __inputPair, __descriptorSize, __c1, __c2, __ddx, __ddy, __angleRadians) \
    {DescriptorPair __pair__ = __inputPair;\
        LOOKUPDESCRIPTORPART_FROMPOLAR(__angleRadians)\
        LOOKUPDESCRIPTORPART_CLAMP(__x, __y, __layer, __descriptorSize, __ddx, __ddy, __angleRadians)\
        LOOKUPDESCRIPTORPART_LOOKUP(__layer, __c1, __c2);\
        __returnValue = true;\
    }
// as function
    //inline void lookupDescriptorInline(bool &returnValue, DOUBLE x, DOUBLE y, const FEATURELayer& layer, const DescriptorPair& inputPair, DOUBLE descriptorSize, CHANNELS& c1, CHANNELS& c2, const DOUBLE ddx = 0, const DOUBLE ddy = 0, const DOUBLE angleRadians = 0) {
    //    DescriptorPair __pair__ = inputPair;
    //    LOOKUPDESCRIPTORPART_FROMPOLAR(angleRadians)
    //    LOOKUPDESCRIPTORPART_PAIRPOLARDISTANCE
    //    LOOKUPDESCRIPTORPART_NORMAL3DPATCHES
    //    LOOKUPDESCRIPTORPART_ROUNDEDPATCHES
    //    LOOKUPDESCRIPTORPART_CLAMP(x, y, layer, descriptorSize, ddx, ddy, angleRadians)
    //    LOOKUPDESCRIPTORPART_LOOKUP(layer, c1, c2)
    //    returnValue = true;
    //}
    //inline void lookupRawDescriptorInline(bool &returnValue, DOUBLE x, DOUBLE y, const FEATURELayer& layer, const DescriptorPair& inputPair, DOUBLE descriptorSize, CHANNELS& c1, CHANNELS& c2, const DOUBLE ddx = 0, const DOUBLE ddy = 0, const DOUBLE angleRadians = 0) {
    //    DescriptorPair __pair__ = inputPair;
    //    LOOKUPDESCRIPTORPART_FROMPOLAR(angleRadians)
    //    LOOKUPDESCRIPTORPART_CLAMP(x, y, layer, descriptorSize, ddx, ddy, angleRadians)
    //    LOOKUPDESCRIPTORPART_LOOKUP(layer, c1, c2)
    //    returnValue = true;
    //}

#define lookupDescriptorV lookupDescriptorMac
#define lookupRawDescriptorV lookupRawDescriptorMac
            }
            using namespace Access;
            //-----------------
        }
        using namespace Algorithms;
        //-----------------

    }
    using namespace Perceptional;
    //-----------------

    //-----------------
    namespace Descriptivity {

        //...
        class DESCRIPTIVITYLayer {
        public:
            DESCRIPTIVITYLayer() {
                uniqueId++;
            }
            int width, height, level;
            Array(DOUBLE) pixels;
            int _id = uniqueId;
        };

        PIXEL descriptivityColor(const DOUBLE d) {
            PIXEL r;
            r.r = sin(CLAMP(d * 0.25, 0, 1) * PI2);
            r.g = sin(CLAMP(d * 0.5, 0, 1) * PI2);
            r.b = sin(CLAMP(d * 1.0, 0, 1) * PI2);
            r.a = CLAMP(d*2.0,0,1);
            return r;
        }

        PIXELLayer toPixelLayer(const DESCRIPTIVITYLayer& v, const ImageConfig &config) {
            PIXELLayer r;
            JUSTSETUPIMAGEHULL(r, v);
            DOUBLE mincolor, maxcolor; bool first = true;
            for (int i = 0; i < ArraySize(v.pixels); ++i) {
                SINGLECHANNELMINMAX(v.pixels[i], mincolor, maxcolor, first)
            }
            for (int i = 0; i < ArraySize(v.pixels); ++i) {
                DOUBLE k = v.pixels[i];
                k = (k - mincolor) / FIXUPDIVBYZERO((maxcolor - mincolor),1.0);
                r.pixels[i] = descriptivityColor(k);
            }
            return r;
        }
    }
    using namespace Descriptivity;
    //-----------------

    //...
    class Image {
    public:
        Array(PIXELLayer) pixelLayers;
        Array(FEATURELayer) featureLayers;
#if DEFINED_DESCRIPTIVITYMAP == 1
        DESCRIPTIVITYLayer descriptivityLayer;
#endif
        std::string name;
    };

    DISPLAYLayer toDisplay(const PIXELLayer& v, const ImageConfig &config) {
        DISPLAYLayer r;
        r.width = v.width;
        r.height = v.height;
        r.uniqueIdOfBaseObject = v._id;
        r.handle = cachedUpload(v, config);
        return r;
    }

    DISPLAYLayer toDisplay(const FEATURELayer& v, const ImageConfig& config) {
        PIXELLayer b = toPixelLayer(v, config);
        b._id = v._id;
        return toDisplay(b, config);
    }

    DISPLAYLayer toDisplay(const DESCRIPTIVITYLayer& v, const ImageConfig& config) {
        PIXELLayer b = toPixelLayer(v, config);
        b._id = v._id;
        return toDisplay(b, config);
    }

    //-----------------
    namespace Generation {
#define RETURNLOADEDIMAGE(__t, __n, ...) \
                            {\
                                ImageConfig config = loadConfig(__VA_ARGS__);\
                                PIXELLayer l = pixelLayerFrom##__t(__n,config);\
                                Image i;\
                                ArrayAdd(i.pixelLayers, l);\
                                i.name = __n;\
                                return i;\
                            }

        Image imageFromFileNameExtension(const std::string& imageFileName, int borderSize = 0, int imageWidth = -1, DOUBLE brightness = 1.0, DOUBLE contrast = 1.0, DOUBLE scale = 1.0) {
            std::string ext = extension(imageFileName);
            if (ext == "jpg") RETURNLOADEDIMAGE(JPG, imageFileName, borderSize, imageWidth, brightness, contrast, scale);
            if (ext == "png") RETURNLOADEDIMAGE(PNG, imageFileName, borderSize, imageWidth, brightness, contrast, scale);
            return Image();
        }

        FEATURELayer toFeatureLayer(const PIXELLayer& v, const std::function<CHANNELS(const PIXEL&)> &t, double fWidth = -1, double fHeight = -1) {
            FEATURELayer r;
            JUSTSETUPIMAGEHULL(r, v);
            r.fheight = fHeight < 0 ? r.height : fHeight;
            r.fwidth = fWidth < 0 ? r.width : fWidth;
            for (int i = 0; i < v.width * v.height; ++i) {
                r.pixels[i] = t(v.pixels[i]);
            }
            return r;
        }
    
#if DEFINED_CHANNELCOUNT == 4
#define FEATURELAYERFORMULA(__c,__f) \
                        {const ICHANNEL channel = red; __c.channel[0] = FCHANNELFROM(__f);} \
                        {const ICHANNEL channel = green; __c.channel[1] = FCHANNELFROM(__f);} \
                        {const ICHANNEL channel = blue; __c.channel[2] = FCHANNELFROM(__f);} \
                        {const ICHANNEL channel = alpha; __c.channel[3] = FCHANNELFROM(__f);}
#else
#define FEATURELAYERFORMULA(__c,__f) \
                        __c = FCHANNELFROM(__f);
#endif

        FEATURELayer buildFeatureLayer(const PIXELLayer& v, int layer, double fWidth = -1, double fHeight = -1) {
#if DEFINED_CHANNELCOUNT == 1
            return toFeatureLayer(v, [=](const PIXEL& v)->CHANNELS {
                const DOUBLE red = FeatureLayerDepthFormula(v.r, layer);
                const DOUBLE green = FeatureLayerDepthFormula(v.g, layer);
                const DOUBLE blue = FeatureLayerDepthFormula(v.b, layer);
                const DOUBLE alpha = FeatureLayerDepthFormula(v.a, layer);
                CHANNELS r;
                FEATURELAYERFORMULA(r, FeatureLayerFormula)
                return r;
                }, fWidth, fHeight);
#endif
#if DEFINED_CHANNELCOUNT == 4
            return toFeatureLayer(v, [=](const PIXEL& v)->CHANNELS {
                const DOUBLE red = FeatureLayerDepthFormula(v.r, layer);
                const DOUBLE green = FeatureLayerDepthFormula(v.g, layer);
                const DOUBLE blue = FeatureLayerDepthFormula(v.b, layer);
                const DOUBLE alpha = FeatureLayerDepthFormula(v.a, layer);
                CHANNELS r;
                FEATURELAYERFORMULA(r, FeatureLayerFormula);
                return r;
            }, fWidth, fHeight);
#endif
        }

        void denoiseRGB(PIXELLayer &v, int kx, int ky, double threshhold) {
            std::map<unsigned int, int> most;
            const PIXELLayer v2 = v; // copy
            for (int y = 0; y < v.height; ++y) {
                for (int x = 0; x < v.width; ++x) {
                    PIXEL maxP = fromRGBA32(0x00000000);
                    PIXEL minP = fromRGBA32(0xffffffff);
                    most.clear();
                    for (int y2 = -ky; y2 <= ky; ++y2) {
                        for (int x2 = -kx; x2 <= kx; ++x2) {
                            const int x3 = x2 + x;
                            const int y3 = y2 + y;
                            if ((unsigned int)x3 < v.width && (unsigned int)y3 < v.height) {
                                const PIXEL &c = v2.pixels[x3 + y3 * v.width];
                                const unsigned int hash = toRGBA32(c);
                                most[hash]++;
                                if (c.r > maxP.r) maxP.r = c.r; if (c.g > maxP.g) maxP.g = c.g; if (c.b > maxP.b) maxP.b = c.b;
                                if (c.r < minP.r) minP.r = c.r; if (c.g < minP.g) minP.g = c.g; if (c.b < minP.b) minP.b = c.b;
                            }
                        }
                    }
                    const double thresh = ((maxP.r - minP.r) + (maxP.g - minP.g) + (maxP.b - minP.b)) / 3.0 * threshhold;
                    if (((maxP.r - minP.r) < thresh) && ((maxP.g - minP.g) < thresh) && ((maxP.b - minP.b) < thresh)) {
                        int mostFound = 0;
                        PIXEL mostC = v2.pixels[x + y * v.width];
                        for (const auto &p : most) {
                            if (p.second > mostFound) {
                                mostC = fromRGBA32(p.first);
                                mostFound = p.second;
                            }
                        }
                        v.pixels[x + y * v.width] = mostC;
                    }
                }
            }
        }

        void preprocessImage(Image &image) {
            if (DEFINED_NORMALIZEDIMAGES == 1) {
                PIXELLayer &v = image.pixelLayers[0];
                PIXEL minColor, maxColor;
                bool first = true;
                for (int i = 0; i < ArraySize(v.pixels); ++i) {
                    PIXELCOLORMINMAX(v.pixels[i], minColor, maxColor, first)
                }
                image.pixelLayers[0] = modification(v, [&](const PIXEL &p) -> PIXEL {
                   PIXEL r = p;
                   normalizeMax(r, minColor, maxColor);
                   return r;
                });
            }
            if (DEFINED_REMOVEDUST == 1) {
                denoiseRGB(image.pixelLayers[0],DENOISE_KERNEL_RAD,DENOISE_KERNEL_RAD, DENOISE_THRESHOLD);
            }
        }

        int setUpImagePyramid(Image& image) {
            ArrayResize(image.featureLayers, 1);
            ArrayResize(image.pixelLayers, 1);
            image.featureLayers[0] = buildFeatureLayer(image.pixelLayers[0], 0);
            double fImageWidth = image.pixelLayers[0].width;
            double fImageHeight = image.pixelLayers[0].height;
            int i = 1;
            while (true) {
                fImageWidth *= descriptorLevelScale;
                fImageHeight *= descriptorLevelScale;
                if (fImageWidth < 2 || fImageHeight < 2)
                    return i;
                PIXELLayer layer;
                layer.width = (int)fImageWidth;
                layer.height = (int)fImageHeight;
                layer.level = ArraySize(image.pixelLayers);
                ArrayResize(layer.pixels, layer.width * layer.height);
                doAveragingFilter(image.pixelLayers[i - 1], layer);
                ArrayAdd(image.pixelLayers, layer);
                ArrayAdd(image.featureLayers, buildFeatureLayer(layer, i, fImageWidth, fImageHeight));
                i++;
            }
        }

        void prepareDescriptorPattern(int levels, int pairs, DescriptorDesc &pattern, const std::function<SHAPEFUNCTION()> &shapeFunction) {
            double scaling = 1.0;
#if DEFINED_ADVANCEDPATTERNTYPE == 1
            scaling = 1.0 / (sqrt(2.0)*sqrt(2.0))*1.75; // aspect ratio yet missing here
            pairs /= 2;
#endif
            ArrayResize(pattern.levels, levels+1);
            srandom(descriptorRandomSeed);
            for (int i = 0; i < levels; ++i) {
                DescriptorDescLevel level;
                for (int j = 0; j < pairs; ++j) {
                    DescriptorPair pair;
                    shapeFunction(pair); pair.x0 *= scaling; pair.y0 *= scaling; pair.x1 *= scaling; pair.y1 *= scaling;
                    SortDescriptorPair(pair);
#if DEFINED_POLARCOORDINATES == 1
                    toPolar(pair);
#endif
                    ArrayAdd(level.pairs, pair);
#if DEFINED_ADVANCEDPATTERNTYPE == 1
                    {
#if DEFINED_POLARCOORDINATES == 1
                        fromPolar(pair);
#endif
                        const double nx = (pair.y1 - pair.y0)*0.5;
                        const double ny = -(pair.x1 - pair.x0)*0.5;
                        pair.x0 = pair.x0 + nx;
                        pair.y0 = pair.y0 + ny;
#if DEFINED_POLARCOORDINATES == 1
                        toPolar(pair);
#endif
                        ArrayAdd(level.pairs, pair);
                    }
#endif
                }
#if DEFINED_DESCRIPTORSHAPETOTRIANGLE == 1
                for (int j = 0; j < ArraySize(level.pairs); ++j) {
                    DescriptorPair& pair = level.pairs[j];
                    pair.x0 *= CLAMP((pair.y0 + 1.0) * 0.5 + 0.5, 0.0, 1.0);
                    pair.x1 *= CLAMP((pair.y1 + 1.0) * 0.5 + 0.5, 0.0, 1.0);
                }
#endif
                pattern.levels[i] = level;
            }
        }

        FeatureConfig createFeatureConfig(int levels, int pairs, bool iterative, DOUBLE stepFactor, int steps, double lastLevelRatio = 1.0, bool resampling = false) {
            FeatureConfig base;
            base.iterative = iterative;
            base.stepFactor = stepFactor;
            base.refinementSteps = steps;
            base.refinementStepPerLevelScaleFactor = refinementStepPerLevelScaleFactor;
            base.deepestLevel = std::max(levels-1,0);
            base.lastLevel = std::max((int)(levels * lastLevelRatio-1),0);
#if DEFINED_ADAPTIVE_TYPE != ADAPTIVE_NONE
            if (resampling) base.refinementStepPerLevelScaleFactor = resampleRefinementStepPerLevelScaleFactor;
#endif
            base.descriptorPatterns.resize(ArraySize(DESCRIPTORSHAPES));
            for (int i = 0; i < base.descriptorPatterns.size(); ++i) {
                prepareDescriptorPattern(
                        (DEFINED_CLIPDESCRIPTORLEVELSTOLASTLEVEL == 1 && DEFINED_VARYINGDESCRIPTORSIZE == 0) ? (
                                base.lastLevel + 1) : (base.deepestLevel + 1), pairs,
                        base.descriptorPatterns[i], DESCRIPTORSHAPES[i]);
            }
            return base;
        }

        FeatureConfig createDefaultFeatureConfig(int levels) {
            FeatureConfig config;
            config = createFeatureConfig(levels, ActiveConfiguration::pairCount, ActiveConfiguration::iterative, ActiveConfiguration::stepFactor, ActiveConfiguration::refinementSteps, ActiveConfiguration::mipsToTakeIntoAccountRatio);
            return config;
        }

#if DEFINED_ADAPTIVE_TYPE != ADAPTIVE_NONE
        FeatureConfig createResampleFeatureConfig(int levels) {
            FeatureConfig config;
            config = createFeatureConfig(levels, pairCount, resampleIterative, resampleStepFactor, resampleRefinementSteps, resampleMipsToTakeIntoAccountRatio, true);
            return config;
        }
#endif

    }
    using namespace Generation;
    //-----------------

    //-----------------
    namespace Hash {
        int hash(const DescriptorDescLevel& level) {
            int k = 0;
            for (int i = 0; i < ArraySize(level.pairs); ++i) {
                k ^= *((int*)(&level.pairs[i].x0));
                k ^= *((int*)(&level.pairs[i].y0));
                k ^= *((int*)(&level.pairs[i].x1));
                k ^= *((int*)(&level.pairs[i].y1));
            }
            return k;
        }
        int hash(const Descriptor& desc) {
            int k = 0;
            for (int i = 0; i < ArraySize(desc.bits); ++i) {
                for (int j = 0; j < ArraySize(desc.bits[i]); ++j) {
                    k ^= desc.bits[i][j];
                    k *= 2;
                    if (k > (1 << 15))
                        k += 1;
                }
            }
            return k;
        }
    }
    using namespace Hash;
    //-----------------

    namespace Functions {

        void globalState() {
            srandom(0);
        }

        void paintFeature(PIXELLayer& layer, const Feature& f, DOUBLE featureOuterRadius = 4.f, DOUBLE featureInnerRadius = 3.f, const PIXEL& featureColor = fromRGBA32(0xffffffff)) {
            PIXEL c = featureColor;
            for (DOUBLE y = f.y - featureOuterRadius; y <= f.y + featureOuterRadius; y += 1.0) {
                for (DOUBLE x = f.x - featureOuterRadius; x <= f.x + featureOuterRadius; x += 1.0) {
                    DOUBLE cx = clamp(x, 0, layer.width - 1);
                    DOUBLE cy = clamp(y, 0, layer.height - 1);
                    DOUBLE dx = (cx - f.x);
                    DOUBLE dy = (cy - f.y);
                    DOUBLE d = sqrt(dx * dx + dy * dy);
                    if (d <= featureOuterRadius && d > featureInnerRadius) {
                        d -= featureInnerRadius;
                        d /= (featureOuterRadius - featureInnerRadius);
                        d = sin(d * 3.1415927);
                        int ix = (int)cx;
                        int iy = (int)cy;
                        c.a = CLAMP(d * 1.5, 0.0, 1.0);
                        alpha(layer.pixels[ix + iy * layer.width], c);
                    }
                }
            }
        }
#define COLORRANGECHECK(__c0, __c1, __channel)
#define BITUNUSABLECHECK
#if DEFINED_COLORRANGE == 1 && DEFINED_PUREBINARITY == 0
    #define COLORRANGECHECK(__c0, __c1, __channel) if ( \
        (channel((__c0), (__channel)) >= DEFINED_COLORRANGE_LOW) && \
        (channel((__c1), (__channel)) >= DEFINED_COLORRANGE_LOW) && \
        (channel((__c0), (__channel)) < DEFINED_COLORRANGE_HIGH) && \
        (channel((__c1), (__channel)) < DEFINED_COLORRANGE_HIGH))
    #define BITUNUSABLECHECK if (bitRequired == vtrue || bitRequired == vfalse) // should check for real binary ones here
#endif

        Descriptor sampleDescriptorSimple(const Image& image, const Feature& feature, const FeatureConfig& config) {
            CHANNELS c0a, c1a;
            Descriptor descriptor;
            globalState();

            initThread(); // for random seed reset

            const double x = feature.x / image.featureLayers[config.coordReferenceLevel].width;
            const double y = feature.y / image.featureLayers[config.coordReferenceLevel].height;
            const double descriptorLevelAdd = log(config.sourceToDestScaleRatio) / log(1.0/config.descriptorLevelScale);

            const int angleCount = ((DEFINED_WITHROTATIONINVARIANCE == 1) && (DEFINED_FULLROTATIONINVARIANCEBYSECTORS == 1)) ? 8 : 1;
            const int descriptorLevelCount = DESCRIPTORLEVELCOUNT;
            const int patternCount = config.descriptorPatterns.size();
            ArrayResize(descriptor.bits, descriptorLevelCount * angleCount * patternCount);

            for (int shapeId = 0; shapeId < patternCount; ++shapeId) {
                const DescriptorDesc &pattern = config.descriptorPatterns[shapeId];
                for (int angle45Id = 0; angle45Id < angleCount; ++angle45Id)
                {
                    for (int i = 0; i < descriptorLevelCount; ++i) {
                        const int imageLevel = i;
                        const int descriptorLevel = i;
                        if (descriptorLevel > config.lastLevel)
                            continue;
                        const DescriptorDescLevel &level = pattern.levels[descriptorLevel];
                        const DOUBLE descriptorSize = DESCRIPTORSIZE(imageLevel - descriptorLevelAdd) /
                                                      config.sourceToDestScaleRatio; // maybe better due to non linearity (and clipping (to be fixed))
                        int o = 0;
                        int nonValids = 0;
                        if (imageLevel < ArraySize(image.featureLayers)) {
                            const int descriptorChannels = channelElementCount();
                            Array(VECTORBOOL) &desc = descriptor.bits[descriptorLevel + (angle45Id * descriptorLevelCount) + (shapeId * descriptorLevelCount * angleCount)];
                            ArrayResize(desc, ArraySize(level.pairs) * descriptorChannels);
                            for (int j = 0; j < ArraySize(level.pairs); ++j) {
                                const DescriptorPair &pair = level.pairs[j];
                                bool isValid = true;
                                lookupRawDescriptorV(isValid, x, y, image.featureLayers[imageLevel], pair,
                                                     descriptorSize, c0a, c1a, 0, 0,
                                                     DEGREETORADIANS(-angle45Id * 45.0));
                                for (int n = 0; n < descriptorChannels; ++n) {
                                    VECTORBOOL b = vunusable;
                                    // already raw bits non inverted, just works forward this way, but this may even be a feature later to do it exactly this way (resbift would be a problem if done differently here)
                                    COLORRANGECHECK(c0a, c1a, n) b = (channel(c0a, n) < channel(c1a, n)) ? vtrue
                                                                                                         : vfalse;
                                    desc[o++] = isValid ? b : vundefined;
                                }
                            }
                        }
                    }
                }
                sampledDescriptors++;
            }
            return descriptor;
        }

        Descriptor sampleDescriptor(const Image& image, const Feature& feature, const FeatureConfig& config) {
            return sampleDescriptorSimple(image, feature, config);
        }

#define SINGLEFEATURESTART\
            const int imageLevelBegin = config.sampleLevel;\
            const int imageLevelEnd = config.sampleLevel;\
            DOUBLE minDescriptorSizeRandom = MINDESCRIPTORSIZERANDOM;\
            DOUBLE maxDescriptorSizeRandom = MAXDESCRIPTORSIZERANDOM;
#define PYRAMIDALFEATURESTART\
            const int imageLevelBegin = 0;\
            const int imageLevelEnd = config.lastLevel;\
            DOUBLE minDescriptorSizeRandom = MINDESCRIPTORSIZERANDOM;\
            DOUBLE maxDescriptorSizeRandom = MAXDESCRIPTORSIZERANDOM;

        // only for ascending/descending position x0->x1 y0->y1 sorted descriptors
        // somehow some sort of "linear descriptivity" (but seems ok).
        // quadratic would make no difference, yet.
        DOUBLE descriptivityOfDescriptor(const Descriptor& descriptor, const FeatureConfig& config) {
            if (descriptor.bits.size() <= config.lastLevel)
                return 0.0;
            PYRAMIDALFEATURESTART
            double onBits = 0;
            double allD = 0;
            const DescriptorDesc &pattern = config.descriptorPatterns[0]; // todo: loop here later or something
            for (int imageLevel = imageLevelEnd; imageLevel >= imageLevelBegin; --imageLevel) {
                const DOUBLE refinementScaleFactor = config.refinementStepPerLevelScaleFactor;
                const DOUBLE perLevelScale = pow(fabs(refinementScaleFactor), refinementScaleFactor < 0 ? (imageLevel - imageLevelBegin) : (imageLevelEnd - imageLevel));
                const int featureRefinementSteps = FIXUPDIVBYZERO((int)(config.refinementSteps * perLevelScale), 1.0);
                const int descriptorLevel = imageLevel;
                if (descriptorLevel > config.lastLevel)
                    continue;
                for (int j = 0; j < ArraySize(descriptor.bits[descriptorLevel]); ++j) {
                    DescriptorPair pair = pattern.levels[descriptorLevel].pairs[j / channelElementCount()];
#if DEFINED_POLARCOORDINATES == 1
                    fromPolar(pair);
#endif
                    const VECTORBOOL b = descriptor.bits[descriptorLevel][j];
                    if (b != vundefined)
                    {
                        const double dx = (pair.x1 - pair.x0);
                        const double dy = (pair.y1 - pair.y0);
                        double d = sqrt(dx * dx + dy * dy);
                        d *= 1.0 / featureRefinementSteps; // bigger steps suggest more descriptivity by that (featureRefinementSteps is a stepcount)
                        if (dx < 0) d *= -1;
                        if (dy < 0) d *= -1;
                        onBits += b ? d : -d;
                        allD += fabs(d);
                    }
                }
            }
            onBits /= FIXUPDIVBYZERO(allD, 1.0);
            onBits = fabs(onBits);
            double scale = 1.0;
            //const double clipTooPerfectOnes = 0.005;
            //const double clipTooBadOnes = 0.3;
            //const double clipBoundaryPow = 10.0;
            //if (onBits < clipTooPerfectOnes) {
            //    scale *= pow(onBits / clipTooPerfectOnes, 1.0 / clipBoundaryPow);
            //}
            //if (onBits > (1.0 - clipTooBadOnes)) {
            //    scale *= pow((1.0 - (onBits - (1.0 - clipTooBadOnes)) / clipTooBadOnes), clipBoundaryPow);
            //}
            const double expressiveness = 40.0;
            const double attenuation = 1.0;
            const double humanReadable = 100.0;
            return pow(scale * (1.0 / FIXUPDIVBYZERO((double)(onBits * attenuation + 1.0), 1.0)), expressiveness) * humanReadable;
        }
 
        void savePixelLayerWithFeatures(const PIXELLayer& layer, const Array(Feature)& features, const std::string& pngFileName, int borderSize = 0) {
            PIXELLayer b = layer;
            pictureWidth = b.width - borderSize * 2;
            pictureHeight = b.height - borderSize * 2;
            if (pictureWidth > 0 && pictureHeight > 0) {
                for (int i = 0; i < ArraySize(features); ++i) {
                    if (features[i].force >= Configuration::pruneThreshold) {
                        paintFeature(b, features[i], 8.f, 6.f, fromRGBA32(0x00ffffff));
                    }
                    paintFeature(b, features[i], 5.f, 3.f, fromRGBA32(0x000000ff));
                }
                for (int y = 0; y < pictureHeight; ++y) {
                    for (int x = 0; x < pictureWidth; ++x) {
                        pictureWriteOut[x + y * pictureWidth] = toRGBA32(b.pixels[x + y * b.width + borderSize + borderSize * b.width]);
                    }
                }
                savePNG(pngFileName.c_str());
            }
        }
    }
    using namespace Functions;
    //-----------------

    //-----------------
    namespace Implementation {
            Array(Image) frames;
            Array(Image) trackingFrames;
            Array(Feature) baseFeatures;
            Array(Feature) updatedFeatures;
            Array(Descriptor) originalDescriptors;
            Array(Descriptor) baseDescriptors;
            Array(Descriptor) updatedDescriptors;
    }
    using namespace ESBIFT::Implementation;
    //-----------------


//--------------------------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------------------------

//-------------------------------------------------------------
// BINARY KLT MACROS
//-------------------------------------------------------------
#define PREPAREROTATIONINVARIANCE
    #define PREPAREROTATIONINVARIANCE const DOUBLE patchAngle = 0; const DOUBLE sectorAngle = 0;
#if DEFINED_WITHROTATIONINVARIANCE == 1
        #define PREPAREROTATIONINVARIANCE const DOUBLE patchAngle = randomLike(randomLikeBase-7777,DEFINED_WITHROTATIONINVARIANCE45DEGREES ? 45.0 : 360.0); const DOUBLE sectorAngle = 0;
        #if DEFINED_FULLROTATIONINVARIANCEBYSECTORS == 1
            #define PREPAREROTATIONINVARIANCE const DOUBLE patchAngle = randomLike(randomLikeBase-7777,360.0); const DOUBLE sectorAngle = (((int)(patchAngle/45.0)) & 7) * 45.0; angleDescriptorLevelAdd = DESCRIPTORLEVELCOUNT * (int)(sectorAngle/45.0);
        #endif
#endif

#define STEPMODULOCHECK  kt++; if ((kt % stepModulo) == 0)

#define COLORVARS \
        CHANNELS c0a, c1a;\
        CHANNELS c0adx, c1adx;\
        CHANNELS c0ady, c1ady;\
        CHANNELS d0adx, d1adx;\
        CHANNELS d0ady, d1ady;

#define ADDSPICE(__x, __y, __l, __iteration, __iterationCount, __spicePercentage, __spiceIntensity) \
        if ((__iteration) < (__iterationCount) * (__spicePercentage) / 100.0) {\
            __x += (DOUBLE)(randomLike(randomLikeBase + 5671, 2.0) - 1.0) * (__spiceIntensity) / (__l).width;\
            __y += (DOUBLE)(randomLike(randomLikeBase + 3781, 2.0) - 1.0) * (__spiceIntensity) / (__l).height;\
        }

#define SELECTPATTERN const int currentPattern = (int)floor(randomLike(randomLikeBase+3378,config.descriptorPatterns.size()));\
                      const DescriptorDesc &pattern = config.descriptorPatterns[currentPattern];

#define CONVERGEFEATURESBEGIN \
        globalState();\
        srandom(FEATURERANDOMSEED);\
        int debugInt = 0;\
        int kt = 0;\
        int angleDescriptorLevelAdd = 0;\
        COLORVARS\
        PATCHVARS\
        DOUBLE baseX = lastPosition.x / image.featureLayers[0].width;\
        DOUBLE baseY = lastPosition.y / image.featureLayers[0].height;\
        int stayingLevel = 0;\
        DOUBLE featureStrength = 0.0;\
        DOUBLE featureStrengthHelper = 0.0;\
        for (int imageLevel = imageLevelEnd; imageLevel >= imageLevelBegin; --imageLevel) {\
                featureStrength *= 2.0; featureStrengthHelper *= 2.0;\
                const DOUBLE descriptorLevel = imageLevel;\
                if (descriptorLevel > config.lastLevel)\
                    continue; \
                const DOUBLE refinementScaleFactor = config.refinementStepPerLevelScaleFactor;\
                const DOUBLE perLevelScale = pow(fabs(refinementScaleFactor) , refinementScaleFactor < 0 ? (imageLevel - imageLevelBegin) : (imageLevelEnd - imageLevel));\
                const int featureRefinementSteps = (int)CLAMP(config.refinementSteps * perLevelScale * loopStepsIncreaseMultiplicator,0,DEFINED_MAXMIPLEVELLOOPSTEPS) / config.sourceToDestScaleRatio;\
                const int stepModulo = std::max((int)(featureRefinementSteps * constantStepRefinementStepsFactor),1);\
                const DOUBLE widthHere = image.featureLayers[imageLevel].width;\
                const DOUBLE heightHere = image.featureLayers[imageLevel].height;\
                if (image.featureLayers[imageLevel].width > 1 && image.featureLayers[imageLevel].height > 1) {\
                    const DOUBLE descriptorSizeInLevel = DESCRIPTORSIZE(imageLevel) * config.sourceToDestScaleRatio;\
                    for (int k = 0; k < featureRefinementSteps; ++k) { \
                        const DOUBLE stepFactorScale = circularLanding ? ((DOUBLE)(featureRefinementSteps-k*circularLandingFactor) / featureRefinementSteps) : 1.0;\
                        const DOUBLE featureStepFactor = (DOUBLE)config.stepFactor / ((DOUBLE)featureRefinementSteps / (config.refinementSteps * perLevelScale)) * stepFactorScale;\
                        const int randomLikeBase = k ^ (imageLevel*55);\
                        SELECTPATTERN\
                        const DescriptorDescLevel& level = pattern.levels[(int)descriptorLevel];\
                        PREPAREROTATIONINVARIANCE\
                        PREPAREPATCHPOLY\
                        const DOUBLE descriptorSize = (config.scaleInvariance ? RANDOMSCALEINVARIANCE(descriptorSizeInLevel) : descriptorSizeInLevel);\
                        DOUBLE movementX = 0.0;\
                        DOUBLE movementY = 0.0;\
                        const int pairCount = ArraySize(level.pairs);\
                        if (useSpice) ADDSPICE(movementX, movementY, image.featureLayers[imageLevel], k, config.refinementSteps, spicePercentage, spiceIntensity);

#if DEFINED_HASRANDOMMOVEMENT == 0
    #define bitFormula_forward(__i0,__i1) ((__i0) < (__i1))
    #define bitFormula_invers(__i0,__i1) ((__i0) >= (__i1))
    #define bitFormula_alternate(__i0,__i1) ((j & 1) ? ((__i0) < (__i1) : ((__i0) >= (__i1)))

    bool resbift_forward(const int randomBase, const DOUBLE i0, const DOUBLE i1) {
        const DOUBLE maxRand = (i1 - i0) * resbiftFactor;
        const DOUBLE newI0 = RESBIFT ? (i0 + randomLike(randomBase,maxRand)) : i0;
        return bitFormula_forward(newI0,i1);
    }

#if DEFINED_ALSO_RESBIFT == 1
    #define bitFormula_forward(__i0,__i1) resbift_forward(j,__i0,__i1)
    #define bitFormula_invers(__i0,__i1) (!resbift_forward(j,__i0,__i1))
    #define bitFormula_alternate(__i0,__i1) (resbift_forward(j,__i0,__i1) ^ ((j & 1) ? true : false))
#endif

#else
    resbift missing here
    #define bitFormula_forward(__i0,__i1) (isUsable ? ((__i0) < (__i1)) : bitRequired ^ (randomLike(j+randomLikeBase,1.0) < 0.5))
    #define bitFormula_invers(__i0,__i1) (isUsable ? ((__i0) >= (__i1)) : bitRequired ^ (randomLike(j+randomLikeBase,1.0) < 0.5))
    #define bitFormula_alternate(__i0,__i1) (isUsable ? ((j & 1) ? ((__i0) < (__i1) : ((__i0) >= (__i1)))) : bitRequired ^ (randomLike(j+randomLikeBase,1.0) < 0.5))
#endif

#if DEFINED_CHANNELCOUNT != 1
// gradient bit required by descriptor
// direct gradient bit
#define DESCRIPTORCHANNELSLOOPHEAD const int descriptorChannels = channelElementCount();\
            const int n = (int)randomLike(randomLikeBase + j,descriptorChannels); { bool isUsable = false; COLORRANGECHECK(c0a, c1a, n) isUsable = true;
#define BITFOUND_BITREQUIRED \
            const bool bitRequired = descriptor.bits[descriptorLevel + angleDescriptorLevelAdd][j * descriptorChannels + n];\
            const bool bitFound = DIRECTIONFORMULA(i0a,i1a);\
            BITUNUSABLECHECK {
#else
#define DESCRIPTORCHANNELSLOOPHEAD const int n = 0; bool isUsable = false; COLORRANGECHECK(c0a, c1a, n) isUsable = true; {
#define BITFOUND_BITREQUIRED \
            const VECTORBOOL bitRequired = descriptor.bits[descriptorLevel + angleDescriptorLevelAdd][j];\
            const bool bitFound = DIRECTIONFORMULA(i0a,i1a);\
            BITUNUSABLECHECK {
#endif

#define DETERMINEFEATURESCALE(__f) __f.scale = getFeatureScale(f.x, f.y, descriptor, image, config);
#define CHECKFORFEATUREIMPORTANCE(__f) __f.force = getFeatureForce(f.x, f.y, descriptor, image, config);

#define CONVERGEFEATURESEND \
                }\
            }\
        }\
        Feature f(baseX * image.featureLayers[0].width, baseY * image.featureLayers[0].height, lastPosition);\
        f.strength = (featureStrengthHelper != 0.0) ? (featureStrength / featureStrengthHelper) : 0.0;\
        CHECKFORFEATUREIMPORTANCE(f)\
        DETERMINEFEATURESCALE(f)\
        return f;

#define DISTANCECHECK(__varname) \
        const DOUBLE __varname = ActiveConfiguration::manhattenDistance ? (fabs(movementX) + fabs(movementY)) : sqrt(movementX * movementX + movementY * movementY);\
        if (__varname != 0.0 && VALIDFLOAT(__varname))

#if DEFINED_ADAPTIVEBRIGHTNESSINVARIANCE == 1
#define PAIRLOOKUP(__channel) \
        const bool iforward = randomLike(randomLikeBase + 92712 + j,1.0) < 0.5;\
        const DOUBLE i0a = iforward ? channel(c0a, __channel) : channel(c1a, __channel);\
        const DOUBLE i1a = iforward ? channel(c1a, __channel) : channel(c0a, __channel);
#else
#define PAIRLOOKUP(__channel) \
        const DOUBLE i0a = channel(c0a, __channel);\
        const DOUBLE i1a = channel(c1a, __channel);
#endif

#if DEFINED_TRIANGULAR_DERIVATIVES == 0
#define GRADIENTLOOKUP(__channel) \
    const DOUBLE i0adx = GRADIENTINVERSION(channel(d0adx, __channel) - channel(c0adx, __channel));\
    const DOUBLE i1adx = GRADIENTINVERSION(channel(d1adx, __channel) - channel(c1adx, __channel));\
    const DOUBLE i0ady = GRADIENTINVERSION(channel(d0ady, __channel) - channel(c0ady, __channel));\
    const DOUBLE i1ady = GRADIENTINVERSION(channel(d1ady, __channel) - channel(c1ady, __channel));
#else
#define GRADIENTLOOKUP(__channel) \
    const DOUBLE i0adx = GRADIENTINVERSION(d0adx - channel(c0adx, __channel));\
    const DOUBLE i1adx = GRADIENTINVERSION(d0ady - channel(c1adx, __channel));\
    const DOUBLE i0ady = GRADIENTINVERSION(d0adx - channel(c0ady, __channel));\
    const DOUBLE i1ady = GRADIENTINVERSION(d0ady - channel(c1ady, __channel));
#endif

    // rather suboptimal debugging functionality but rather convenient..
#ifdef DEBUG
#define DLOOP(__m) {debugValues.clear(); debugMarkers[__COUNTER__+1].first = __m; debugMarkers[__COUNTER__+1].second++; }
#define DFOR(...) for(__VA_ARGS__,debugLoop(__COUNTER__,"_forAdvance",1))
#define LD(__v) ,(debugValues[###__v].first++, debugValues[###__v].second = toString(__v))
#define LN(__v) ,(debugValues[###__v].first++, debugValues[###__v].second = "")
#define D(__v) (debugValues[###__v].first++, debugValues[###__v].second = toString(__v))
#define N(__v) (debugValues[###__v].first++, debugValues[###__v].second = "")
#else
#define DLOOP(__m)
#define DFOR(...) for(__VA_ARGS__)
#define LD(__v) 
#define LN(__v) 
#define D(__v) 
#define N(__v)
#endif
    std::map<unsigned int, std::pair<std::string,int>> debugMarkers;
    std::map<std::string, std::pair<int,std::string>> debugValues;
    void log(const std::pair<std::string, std::pair<int, std::string>> &v) {
        if (v.second.second.empty()) 
            printf(":%s<%d>", v.first.c_str(), v.second.first);
        else
            printf(":%s=%s<%d>", v.first.c_str(), v.second.second.c_str(), v.second.first);
    }
    void logDebug() {
        for (auto& v : debugValues) {
            log(v);
        }
        printf("\n");
    }
    void debugLoop(int debugId, const std::string &suffix = "", int add = 0) {
        if (debugMarkers.find(debugId) != debugMarkers.end()) {
            debugMarkers[debugId].second++;
            std::string debugMarker = debugMarkers[debugId].first + "<" + toString(debugMarkers[debugId].second) + ">" + suffix;
            printf("loop:%s->", debugMarker.c_str());
            logDebug();
        }
    }

    inline int signcmp2(DOUBLE f1, DOUBLE f2) {
        return f1 > f2 ? 1 : -1;
    }

    inline int signcmp(DOUBLE f1, DOUBLE f2) {
        return f1 < f2 ? 1 : -1;
    }

    inline int dircmp(int f1, int f2) {
        return f1 != f2 ? -1 : 1;
    }

    inline int idir1(int a1, int a2) {
        return (a1 == -1) ? -a2 : a2;
    }
    
    inline int idir2(int a1, int a2) {
        return (a2 == 1) ? a1 : -a1;
    }

    inline int idiff1(DOUBLE f1, DOUBLE f2) {
        return f1 < f2 ? 1 : -1;
    }

    inline double xdir1(DOUBLE a1, DOUBLE a2) {
        return (a1 < 0) ? a2 : 0;
    }

    inline double xdir2(DOUBLE a1, DOUBLE a2) {
        return (a2 < 0) ? a1 : 0;
    }

    std::string currentAlgorithmShortDescription = "";
#define PREPAREDESCRIPTION(__a) currentAlgorithmShortDescription = (__a);

#define DESCRIPTORLOOP(__DESCRIPTORFUNCTION__)\
            DESCRIPTORCHANNELSLOOPHEAD\
                descriptorCountAll++;\
                PAIRLOOKUP(n)\
                GRADIENTLOOKUP(n)\
                BITFOUND_BITREQUIRED\
                        if (bitRequired == bitFound) descriptorCountMatching++;\
                        if ((!bitFound)) descriptorCountFlat++;\
                        __DESCRIPTORFUNCTION__(n)\
                }\
            }

#define STILLMOVEMENT\
            if (movementX == 0.0 && movementY == 0.0) { \
                movementX = randomLike(randomLikeBase,3.0)-1.5;\
                movementY = randomLike(randomLikeBase+2222,3.0)-1.5;\
            }

#define FLOATINGUPDATE\
            STILLMOVEMENT DISTANCECHECK(distance) {\
                movementX /= distance;\
                movementY /= distance;\
                baseX += movementX * featureStepFactor / widthHere;\
                baseY += movementY * featureStepFactor / heightHere;\
                baseX = CLAMP(baseX, 0.0, 1.0);\
                baseY = CLAMP(baseY, 0.0, 1.0);\
                movementX = movementY = 0;\
            }

#define STEPUPDATE \
            STILLMOVEMENT DISTANCECHECK(distance) {\
                movementX /= distance;\
                movementY /= distance;\
                baseX += constantStepMoveFactor * movementX * featureStepFactor / widthHere;\
                baseY += constantStepMoveFactor * movementY * featureStepFactor / heightHere;\
                baseX = CLAMP(baseX, 0.0, 1.0);\
                baseY = CLAMP(baseY, 0.0, 1.0);\
                movementX = movementY = 0;\
            }

#define LOOKUPDESCRIPTOR(__destc1, __destc2, __gdx, __gdy, __radians) lookupDescriptorV(isValid, baseX, baseY, image.featureLayers[imageLevel], pair, descriptorSize, (__destc1), (__destc2), (__gdx), (__gdy), (__radians)); if (!isValid) continue;

#if DEFINED_TRIANGULAR_DERIVATIVES == 0
#define FEATURELOOKUP_BASE LOOKUPDESCRIPTOR(c0a, c1a, 0, 0, angleModRadians)
    #if DEFINED_DERIVATIVE_ROTATION == 0
        #define FEATURELOOKUP_DERIVATIVES \
            LOOKUPDESCRIPTOR(c0adx, c1adx, 1, 0, angleModRadians)\
            LOOKUPDESCRIPTOR(c0ady, c1ady, 0, 1, angleModRadians)\
            LOOKUPDESCRIPTOR(d0adx, d1adx, -1, 0, angleModRadians)\
            LOOKUPDESCRIPTOR(d0ady, d1ady, 0, -1, angleModRadians)
        #else
        #define FEATURELOOKUP_DERIVATIVES\
            if (j & 1) {\
                LOOKUPDESCRIPTOR(c0adx, c1adx, 1, 0, angleModRadians)\
                LOOKUPDESCRIPTOR(c0ady, c1ady, 0, 1, angleModRadians)\
            }\
            else {\
                LOOKUPDESCRIPTOR(c0adx, c1ady, 0, 1, angleModRadians)\
                LOOKUPDESCRIPTOR(c0ady, c1adx, 1, 0, angleModRadians)\
            }\
            if (j & 2) {\
                LOOKUPDESCRIPTOR(d0adx, d1adx, -1, 0, angleModRadians)\
                LOOKUPDESCRIPTOR(d0ady, d1ady, 0, -1, angleModRadians)\
            }\
            else {\
                LOOKUPDESCRIPTOR(d0adx, d1ady, 0, -1, angleModRadians)\
                LOOKUPDESCRIPTOR(d0ady, d1adx, -1, 0, angleModRadians)\
            }
    #endif
#else
    #define FEATURELOOKUP_BASE LOOKUPDESCRIPTOR(c0a, c1a, 0, 0, angleModRadians)
    #if DEFINED_DERIVATIVE_ROTATION == 0
        #define FEATURELOOKUP_DERIVATIVES\
            LOOKUPDESCRIPTOR(d0adx, d0ady, -1 + 0.25, -1 + 0.25, angleModRadians)\
            LOOKUPDESCRIPTOR(c0adx, c1adx, 1 + 0.25, 0 + 0.25, angleModRadians)\
            LOOKUPDESCRIPTOR(c0ady, c1ady, 0 + 0.25, 1 + 0.25, angleModRadians)
    #else
        #define FEATURELOOKUP_DERIVATIVES\
            if (j & 1) {\
                LOOKUPDESCRIPTOR(d0adx, d0ady, -1 + 0.25, -1 + 0.25, angleModRadians)\
                LOOKUPDESCRIPTOR(c0adx, c1adx, 1 + 0.25, 0 + 0.25, angleModRadians)\
                LOOKUPDESCRIPTOR(c0ady, c1ady, 0 + 0.25, 1 + 0.25, angleModRadians)\
            }\
            else {\
                LOOKUPDESCRIPTOR(d0ady, d0adx, 1 - 0.25, 1 - 0.25, angleModRadians)\
                LOOKUPDESCRIPTOR(c1adx, c0adx, -1 - 0.25, 0 - 0.25, angleModRadians)\
                LOOKUPDESCRIPTOR(c1ady, c0ady, 0 - 0.25, -1 - 0.25, angleModRadians)\
            }
    #endif
#endif

#define FEATUREINTENSITYLOOKUP const std::vector<VECTORBOOL> &bits = descriptor.bits[descriptorLevel + angleDescriptorLevelAdd];\
    if (j >= bits.size() || j >= level.pairs.size())\
        continue;\
    const DescriptorPair& pair = level.pairs[j];\
    const double angleModRadians = DEGREETORADIANS(ANGLEMOD(patchAngle));\
    const double sectorAngleRadians = DEGREETORADIANS(sectorAngle);\
    bool isValid = true;\
    FEATURELOOKUP_BASE\
    FEATURELOOKUP_DERIVATIVES

//--------------------------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------------------------
#define FEATURESTRENGTHSTUFF featureStrength += (1.0 - descriptorFlatRatio); featureStrengthHelper += 1.0;\
    if (descriptorFlatRatio > 0.75) stayingLevel = imageLevel; // currently not useful in any way

#define PAIRLOOPHEADER int descriptorCountMatching = 0; int descriptorCountAll = 0; int descriptorCountFlat = 0; for(int j = 0; j < pairCount; ++j) {
#define PAIRLOOPDONE } \
    const double descriptorMatchingRatio = descriptorCountAll != 0 ? (double)descriptorCountMatching / descriptorCountAll : 0;\
    const double descriptorFlatRatio = descriptorCountAll != 0 ? (double)descriptorCountFlat / descriptorCountAll : 0;\
    FEATURESTRENGTHSTUFF

#define __NOFUNCTION(...)
#define __NOSTEP
#define __STEP STEPMODULOCHECK {STEPUPDATE};
#define __FLOATINGSTEP STEPMODULOCHECK {FLOATINGUPDATE};

#define PAIRLOOPFORFUNCTION(__ESBIFTFUNCTION, __STEPFUNCTION) PAIRLOOPHEADER\
                            FEATUREINTENSITYLOOKUP\
                            DESCRIPTORLOOP(__ESBIFTFUNCTION)\
                            __STEPFUNCTION\
                            PAIRLOOPDONE

DOUBLE getFeatureScale(const double x, const double y, const Descriptor &descriptor, const Image &image, const FeatureConfig &config) {
    //const DescriptorDesc &pattern = config.descriptorPatterns[0]; // todo: loop or something here
    //const DOUBLE baseX = x; const DOUBLE baseY = y;
    //COLORVARS;
    //PATCHVARS;
    //int descriptorLevel = 0;
    //int descriptorLevelAdd = 0;
    //int angleDescriptorLevelAdd = 0;
    //double imageLevel = 0;
    //DOUBLE patchAngle = 0;
    //DOUBLE sectorAngle = 0;
    //double scale = 0.0;
    //double weight = 0.0;
    //PYRAMIDALFEATURESTART
    //const int randomSamples = config.scaleLevelChecksCount;
    //DOUBLE featureStrength = 0.0;
    //DOUBLE featureStrengthHelper = 0.0;
    //int stayingLevel = 0;
    //for (int i = 0; i < randomSamples; ++i) { imageLevel = random(imageLevelEnd - imageLevelBegin) + imageLevelBegin;
    //    const int randomLikeBase = imageLevel * 55;
    //    const DOUBLE descriptorLevel = imageLevel;
    //    const DOUBLE descriptorSize = DESCRIPTORSIZE(imageLevel);
    //    const DescriptorDescLevel& level = pattern.levels[descriptorLevel];
    //    PAIRLOOPFORFUNCTION(__NOFUNCTION, __NOSTEP)
    //    double scalev = pow(1.0/DEFAULTMIPLEVELSCALING, imageLevel);
    //    scale += scalev * descriptorMatchingRatio;
    //    weight += scalev;
    //}
    //scale /= weight;
    //scale *= 4; // no clue about this factor, yet (i thinks somehow the technique in the loop above, may work or will work if improved some tiny bit)
    //return scale;
    return 1;
}

DOUBLE getFeatureForce(const double x, const double y, const Descriptor &descriptor, const Image &image, const FeatureConfig &config) { // calculating gradient strengths on these feature points
    PYRAMIDALFEATURESTART
    double l1 = 0; double l2 = 0;
    for (int i = 0; i < config.descriptorPatterns.size(); ++i) {
        const DescriptorDesc &pattern = config.descriptorPatterns[i];
        for (int j = imageLevelBegin; j <= imageLevelEnd; ++j) {
            const FEATURELayer &layer = image.featureLayers[j];
            DOUBLE dSize = DESCRIPTORSIZE(j) * image.featureLayers[j].fwidth; // it's wrong but looks better (j)
            const DescriptorDescLevel &level = pattern.levels[j];
            const double sx = layer.fwidth / image.featureLayers[0].fwidth;
            const double sy = layer.fheight / image.featureLayers[0].fheight;
            for (int k = 0; k < level.pairs.size(); ++k) {
                const DescriptorPair &pair = level.pairs[k];
                CHANNELS c1,c2; lookupPairDirectClamped((pair.x0 * dSize + x) * sx, (pair.y0 * dSize + y) * sy, (pair.x1 * dSize + x) * sx, (pair.y1 * dSize + y) * sy, layer, c1, c2);
                const double i1 = fabs(intensity(c1));
                const double i2 = fabs(intensity(c2));
                l1 += fabs(i1 - i2);
                l2 += fabs(i1 + i2) * 0.5;
            }
        }
    }
    const double flatRatio = (l2 == 0.0) ? 0 : (double)l1 / l2;
    return flatRatio;
}

//--------------------------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------------
// all ESBIFT functions
//----------------------------------------------------------------------
#define ESBIFTS\
    ESBIFT(algo_add, Addition, "Addition", PAIRLOOPFORFUNCTION(ADDITION_FEATUREREFINE, __NOSTEP) FLOATINGUPDATE)\
    ESBIFT(algo_step, Step, "Step", PAIRLOOPFORFUNCTION(STEP_FEATUREREFINE, __STEP))\
    ESBIFT(algo_glideadd, GlideAdd, "GlideAdd", PAIRLOOPFORFUNCTION(GLIDEADD_FEATUREREFINE, __NOSTEP) FLOATINGUPDATE)\
    ESBIFT(algo_glidestep, GlideStep, "GlideStep", PAIRLOOPFORFUNCTION(GLIDESTEP_FEATUREREFINE, __STEP))\
    ESBIFT(algo_distance, Distance, "Distance", PAIRLOOPFORFUNCTION(DISTANCE_FEATUREREFINE, __NOSTEP) FLOATINGUPDATE)\
    ESBIFT(algo_multiply, Multiply, "Multiply", PAIRLOOPFORFUNCTION(MULTIPLICATION_FEATUREREFINE, __NOSTEP) FLOATINGUPDATE)\
    ESBIFT(algo_power, Power, "Power", PAIRLOOPFORFUNCTION(POWER_FEATUREREFINE, __NOSTEP) FLOATINGUPDATE)\
    ESBIFT(algo_powerstep, PowerStep, "PowerStep", PAIRLOOPFORFUNCTION(STEP_SIGNED_POWER_FEATUREREFINE, __STEP))\
    ESBIFT(algo_combinestep, CombineStep, "CombineStep", PAIRLOOPFORFUNCTION(COMBINESTEP_FEATUREREFINE, __STEP))\
    ESBIFT(algo_aipopaddstep, AIPopAddStep, "AIPopAddStep", PAIRLOOPFORFUNCTION(AIPOPADDSTEP_FEATUREREFINE, __STEP))\
    ESBIFT(algo_aipopsubstep, AIPopSubStep, "AIPopSubStep", PAIRLOOPFORFUNCTION(AIPOPSUBSTEP_FEATUREREFINE, __STEP))\
    ESBIFT(algo_aipopadd, AIPopAdd, "AIPopAdd", PAIRLOOPFORFUNCTION(AIPOPADD_FEATUREREFINE, __FLOATINGSTEP))\
    ESBIFT(algo_aipopsub, AIPopSub, "AIPopSub", PAIRLOOPFORFUNCTION(AIPOPSUB_FEATUREREFINE, __FLOATINGSTEP))\
    ESBIFT(algo_aipoppoweradd, AIPopPowerAdd, "AIPopPowerAdd", PAIRLOOPFORFUNCTION(AIPOP_SIGNED_POWER_FEATUREREFINE, __FLOATINGSTEP))\
    ESBIFT(algo_aipop2poweradd, AIPop2PowerAdd,"AIPop2PowerAdd", PAIRLOOPFORFUNCTION(AIPOP2_SIGNED_POWER_FEATUREREFINE, __FLOATINGSTEP))\
    ESBIFT(algo_aipop3poweradd, AIPop3PowerAdd,"AIPop3PowerAdd", PAIRLOOPFORFUNCTION(AIPOP3_SIGNED_POWER_FEATUREREFINE, __FLOATINGSTEP))\
    ESBIFT(algo_aipop4poweradd, AIPop4PowerAdd,"AIPop4PowerAdd", PAIRLOOPFORFUNCTION(AIPOP4_SIGNED_POWER_FEATUREREFINE, __FLOATINGSTEP))\


//----------------------------------------------------------------------
// CPU ESBIFTS
//----------------------------------------------------------------------
#define ESBIFT(__case, __algorithm, __string, ...)\
Feature ESBIFT_FeatureLookup_SingleLayer_##__algorithm(const Image& image, const Descriptor& descriptor, Feature lastPosition, const FeatureConfig& config) {\
    PREPAREDESCRIPTION("SingleLayer_" + std::string(__string))\
    SINGLEFEATURESTART\
    CONVERGEFEATURESBEGIN\
        __VA_ARGS__\
    CONVERGEFEATURESEND\
}\
Feature ESBIFT_FeatureLookup_Pyramidal_##__algorithm(const Image& image, const Descriptor& descriptor, Feature lastPosition, const FeatureConfig& config) {\
    PREPAREDESCRIPTION("Pyramidal_" + std::string(__string))\
    PYRAMIDALFEATURESTART\
    CONVERGEFEATURESBEGIN\
        __VA_ARGS__\
    CONVERGEFEATURESEND\
}
ESBIFTS // ------- actual ESBIFTS
//----------------------------------------------------------------------

//#define SHADERBEGIN(...) "void main() {"
//#define SHADEREND(...) "}"
//
//void shader() {
//    SHADERBEGIN()
//    SHADEREND()
//}


//--------------------------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------------------------

#define SELECTBINARYKLTFEATURELOOKUPVERSION(...) \
        if (pyramidal) \
            return ESBIFT_FeatureLookup_Pyramidal_##__VA_ARGS__(image, descriptor, lastPosition, config); \
        else \
            return ESBIFT_FeatureLookup_SingleLayer_##__VA_ARGS__(image, descriptor, lastPosition, config);
#define ALGOCASE(__case,__formula) case __case: {\
    SELECTBINARYKLTFEATURELOOKUPVERSION(__formula)\
}

Feature lookupFeatureType(const Image& image, const Descriptor& descriptor, Feature lastPosition, const FeatureConfig& config) {
        featuresInvolved++;
#define ESBIFT(__case, __algorithm, __string, ...) ALGOCASE(__case, __algorithm)
       switch (algorithm) {
            ESBIFTS
        };
        return Feature();
}

Feature lookupFeature(const Image& image, const Descriptor& descriptor, Feature lastPosition, const FeatureConfig& config) {
        return lookupFeatureType(image, descriptor, lastPosition, config);
}

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------


//-----------------
namespace Minks {

    //...
    class BinkEdge {
    public:
        int fId0 = -1, fId1 = -1;
        double wn = 0,wt = 0; // normal weight and tangent weight (for choosing of quadrant and influence)
        double dn = 0,dt = 0; // normal distance and tangent distance for positioning
        int randomBase = - 1;
    };
    using BinkEdgeMap = std::map<std::pair<int,int>,BinkEdge>;

    BinkEdgeMap edgeMap(const std::vector<BinkEdge> &e) {
        BinkEdgeMap r;
        for (int i = 0; i < e.size(); ++i) {
            r.insert(std::make_pair(std::make_pair(e[i].fId0,e[i].fId1),e[i]));
        }
        return r;
    }

    BinkEdge getEdge(const BinkEdgeMap &map, const int fId0, const int fId1) {
        BinkEdgeMap::const_iterator t = map.find(std::make_pair(fId0,fId1));
        if (t != map.end())
            return t->second;
        return BinkEdge();
    }

    //...
    class Bink {
    public:
        Bink() {;}
        Bink(DOUBLE x, DOUBLE y) {
            this->x = x;
            this->y = y;
            this->px = x;
            this->py = y;
            this->featureId = uniqueId;
            uniqueId++; // todo: thread safe
        }
        std::vector<Feature> minkies;
        std::vector<BinkEdge> edgs;
        DOUBLE x = -1, y = -1;
        DOUBLE px = 0, py = 0;
        int id = -1;
        int frameId = -1;
        bool selected = false;
        bool isFloating = false;
        bool isMetaBink = false;
        int featureId;
    };

#define INVALIDBINK(__v) ((__v).px == 0.0 && (__v).py == 0.0)

    const double binkForce(const Bink &b) {
        double f = 0;
        for (int i = 0; i < ArraySize(b.minkies); ++i) f += b.minkies[i].force;
        return f;
    }

    Feature feature(const Bink &b) {
        Feature r(b.x, b.y);
        r.id = b.featureId;
        r.frame = b.frameId;
        r.force = binkForce(b);
        return r;
    }

    Feature calculatedFeature(const Bink &b) {
        Feature r(b.px, b.py);
        r.id = b.featureId;
        r.frame = b.frameId;
        r.force = binkForce(b);
        return r;
    }

    double distance(const double x, const double &y, const Bink &b) { // compiler put here it seems
        const double dx = b.px - x;
        const double dy = b.py - y;
        return sqrt(dx * dx + dy * dy);
    }

    Feature getMinkCenter(const Bink &b) {
        double px = 0; double py = 0;
        for (int i = 0; i < ArraySize(b.minkies); ++i) {
            py += b.minkies[i].y;
            px += b.minkies[i].x;
        }
        if (ArraySize(b.minkies) != 0) {
            px /= ArraySize(b.minkies);
            py /= ArraySize(b.minkies);
        }
        Feature f(px, py);
        return f;
    }

    Bink getBink(const std::vector<Feature> &features, DOUBLE x, DOUBLE y, const int id, int frameId) {
        Bink r(x,y);
        r.id = id;
        r.frameId = frameId;
        for (int i = 0; i < features.size(); ++i) {
            r.minkies.push_back(features[i]);
        }
        return r;
    }

    Bink getMetaBink(const std::vector<Bink> &binks, DOUBLE x, DOUBLE y, const int id, int frameId) {
        Bink r(x,y);
        r.id = id;
        r.frameId = frameId;
        r.isMetaBink = true;
        for (int i = 0; i < binks.size(); ++i) {
            r.minkies.push_back(feature(binks[i]));
        }
        return r;
    }

    void updateCenterBink(const Image &image, const Array(Feature) &features, Bink &b, const ImageConfig &config, bool isMeta = false) {
        const FeatureMap fMap = featureMap(features);
        double cx = 0;
        double cy = 0;
        // center
        for (int i = 0; i < ArraySize(b.minkies); ++i) {
            cx += get(fMap, b.minkies[i].id, b.minkies[i]).x;
            cy += get(fMap, b.minkies[i].id, b.minkies[i]).y;
        }
        cx /= ArraySize(b.minkies) == 0 ? 1 : ArraySize(b.minkies);
        cy /= ArraySize(b.minkies) == 0 ? 1 : ArraySize(b.minkies);

        // updated center
#if DEFINED_REFINED_BINK_CENTERS == 1
        double cx2 = 0;
        double cy2 = 0;
        double w = 0;
        for (int i = 0; i < ArraySize(b.minkies); ++i) {
            const double x = get(fMap, b.minkies[i].id, b.minkies[i]).x;
            const double y = get(fMap, b.minkies[i].id, b.minkies[i]).y;
            const double dx = cx - x;
            const double dy = cy - y;
            const double d = 1.0 / ((dx * dx + dy * dy) + 0.01); // 0.01 to avoid division by zero
            cx2 += x * d;
            cy2 += y * d;
            w += d;
        }
        cx2 /= w == 0.0 ? 1 : w;
        cy2 /= w == 0.0 ? 1 : w;
        cx = cx2;
        cy = cy2;
#endif
        b.px = cx;
        b.py = cy;
    }

    void updateBink(const Image &image, const Array(Feature) &features, Bink &b, const ImageConfig &config, bool isMeta = false) {
        if (b.minkies.size() < 2)
            return;

        const FeatureMap fMap = featureMap(features);

        // shuffle minkies
        std::vector<int> strip;
        strip.resize(b.minkies.size());
        for (int i = 0; i < b.minkies.size(); ++i) {
            strip[i] = i;
        }
        // sorting source features after there mostly angula positions (means the snapshot taken at addBlinks, can be done at addBlinks already)
#define SORTBYANGLE(__x, __y)  ( (__x) / ((__y)+sign(__y)*3.0)*2.0*fabs(__x) - fabs(__y) / fabs((__x)+sign(__x)*1.5)*0.5/fabs(__y) ) // position to value formula
        //return SORTBYANGLE(f1.x, f1.y) - SORTBYANGLE(f2.x, f2.y) < 0;
        std::sort(strip.begin(), strip.end(), [&](const int v1, const int v2) -> bool {
            const Feature &f1 = b.minkies[v1];
            const Feature &f2 = b.minkies[v2];
            return f1.id < f2.id;
        });

        double xp1 = 0;
        double yp1 = 0;
        double xp2 = 0;
        double yp2 = 0;
        double weight1 = 0;
        double weight2 = 0;

        const double minLineWidth = image.featureLayers[0].width * (double) 10.0 / 640.0;

        ArrayResize(b.edgs,ArraySize(b.minkies));

        for (int i = 0; i < ArraySize(b.minkies); ++i) {

            const Feature &firstMinkie = b.minkies[strip[i]];
            const Feature &secondMinkie = b.minkies[strip[(i + 1) % strip.size()]]; // ensures same pairs are more common through before sorting this array

            BinkEdge &edg = b.edgs[i];
            edg.fId0 = firstMinkie.id;
            edg.fId1 = secondMinkie.id;
            edg.wt = 0; edg.wn = 0;
            edg.dn = 0; edg.wn = 0;
            edg.randomBase = i + 234234;

            const Feature firstMinkieNow = get(fMap, edg.fId0, firstMinkie);
            const Feature secondMinkieNow = get(fMap, edg.fId1, secondMinkie);

            const int quadrant = (int) floor(randomLike(edg.randomBase, 4));

            const double px0_a = firstMinkie.x;
            const double py0_a = firstMinkie.y;
            const double px1_a = secondMinkie.x;
            const double py1_a = secondMinkie.y;
            const double px0_b = firstMinkieNow.x;
            const double py0_b = firstMinkieNow.y;
            const double px1_b = secondMinkieNow.x;
            const double py1_b = secondMinkieNow.y;

            const double xmove1 = px0_a - px0_b;
            const double xmove2 = px1_a - px1_b;
            const double ymove1 = py0_a - py0_b;
            const double ymove2 = py1_a - py1_b;
            const double m1 = sqrt( xmove1 * xmove1 + ymove1 * ymove1 );
            const double m2 = sqrt( xmove2 * xmove2 + ymove2 * ymove2 );

            if ((!isMeta) && (m1 < minkMoveThreshold || m2 < minkMoveThreshold))
                continue;

            const double line1dx = px1_a - px0_a;
            const double line1dy = py1_a - py0_a;
            const double line2dx = px1_b - px0_b;
            const double line2dy = py1_b - py0_b;
            const double line1d = sqrt(line1dx * line1dx + line1dy * line1dy);
            const double line2d = sqrt(line2dx * line2dx + line2dy * line2dy);
            if (line1d < minLineWidth || line2d < minLineWidth)
                continue;

            // barycentric lookup with coord cross and quadrants
            const double tx_a = px1_a - px0_a;
            const double ty_a = py1_a - py0_a;
            const double tx_b = px1_b - px0_b;
            const double ty_b = py1_b - py0_b;
            const double nx_a = ty_a;
            const double ny_a = -tx_a;
            const double nx_b = ty_b;
            const double ny_b = -tx_b;

            // the weighted distance in reference frame of the original picture / minkie
            const double d_n = (b.x - px0_a) * nx_a + (b.y - py0_a) * ny_a;
            const double d_t = (b.x - px0_a) * tx_a + (b.y - py0_a) * ty_a;
            const double ln = sqrt(nx_a * nx_a + ny_a * ny_a) * sqrt(nx_b * nx_b + ny_b * ny_b); // normalizing weights
            const double lt = sqrt(tx_a * tx_a + ty_a * ty_a) * sqrt(tx_b * tx_b + ty_b * ty_b); // normalizing weights

#if DEFINED_CORNER_BINKS == 1
            const bool useNQuadrant = (fabs(sign(d_n) - ((quadrant & 1) * 2 - 1)) < 0.001);
            const bool useTQuadrant = (fabs(sign(d_t) - ((quadrant & 2) * 1 - 1)) < 0.001);
#else
            const bool useNQuadrant = true;
            const bool useTQuadrant = true;
#endif

            // weighted add by strength of found features (->hopefully)
#if DEFINED_WEIGHTED_BINKS == 1
            double w = firstMinkieNow.force * secondMinkieNow.force * firstMinkie.force * secondMinkie.force;
#else
            double w = 1.0;
#endif
            w *= pow( (b.x - px0_a) * (b.x - px0_a) + (b.y - py0_a) * (b.y - py0_a), POW_FOR_BINKS );
            if (useNQuadrant) {
                edg.dn = d_n * w;
                edg.wn = w * ln * 0.5;
                weight1 += edg.wn;
                xp1 += edg.dn * nx_b;
                yp1 += edg.dn * ny_b;
                xp1 += edg.wn * px0_b;
                yp1 += edg.wn * py0_b;
            }
            if (useTQuadrant) {
                edg.dt = d_t * w;
                edg.wt = w * lt * 0.5;
                weight2 += edg.wt;
                xp2 += edg.dt * tx_b;
                yp2 += edg.dt * ty_b;
                xp2 += edg.wt * px0_b;
                yp2 += edg.wt * py0_b;
            }
        }
        if (weight1 == 0.0 && weight2 == 0.0)
            return;
        double x = (weight1 != 0) ? xp1 / weight1 * 0.5 : 0;
        double y = (weight1 != 0) ? yp1 / weight1 * 0.5 : 0;
        x += (weight2 != 0) ? xp2 / weight2 * 0.5 : 0;
        y += (weight2 != 0) ? yp2 / weight2 * 0.5 : 0;
        b.px = x;
        b.py = y;
    }

    void updateMetaBink(const Image &image, const Array(Bink) &binks, Bink &b, const ImageConfig &config) {
        Array(Feature) features; ArrayResize(features,binks.size());
        for (int i = 0; i < binks.size(); ++i) features[i] = calculatedFeature(binks[i]);
        updateBink(image, features, b, config, true);
    }

    Bink getFloatingBink(const Image &image, const Array(Bink) &binks, const std::vector<Feature> &features, int binkId, const ImageConfig &config) {

        Array(BinkEdgeMap) edgeMaps;

        for (int i = 0; i < binks.size(); ++i) {
            const Bink &b = binks[i];
            if (b.id != binkId)
                continue;
            edgeMaps.push_back(edgeMap(b.edgs));
        }

        if (edgeMaps.empty())
            return Bink();

        // .....

        const FeatureMap fMap = featureMap(features);
        const BinkEdgeMap &baseEdgeMap = edgeMaps[0]; // it's just needed that there is the same pairs in all edgemaps, so it's mutual exclusive

        double xp1 = 0;
        double yp1 = 0;
        double xp2 = 0;
        double yp2 = 0;
        double weight1 = 0;
        double weight2 = 0;

        for (const auto &e : baseEdgeMap) {

            const BinkEdge &edg = e.second;
            const Feature firstMinkieNow = get(fMap, edg.fId0, Feature());
            const Feature secondMinkieNow = get(fMap, edg.fId1, Feature());
            const double px0_b = firstMinkieNow.x;
            const double py0_b = firstMinkieNow.y;
            const double px1_b = secondMinkieNow.x;
            const double py1_b = secondMinkieNow.y;
            const double tx_b = px1_b - px0_b;
            const double ty_b = py1_b - py0_b;
            const double nx_b = ty_b;
            const double ny_b = -tx_b;

            double wn = 0;
            double wt = 0;
            double dn = 0;
            double dt = 0;
            for (int i = 0; i < ArraySize(edgeMaps); ++i) {
                const BinkEdgeMap &m = edgeMaps[i];
                const BinkEdge &v = getEdge(m, edg.fId0, edg.fId1);
#define INVALIDBINKEDGE(__v) ((__v).wn == 0.0 && (__v).wt == 0.0)
                if (INVALIDBINKEDGE(v)) {
                    wn = 0.0; wt = 0.0; dn = 0.0; wn = 0.0;
                    break;
                }
                wn += v.wn;
                wt += v.wt;
                dn += v.dn;
                dt += v.dt;
            }
            wn /= ArraySize(edgeMaps);
            wt /= ArraySize(edgeMaps);
            dn /= ArraySize(edgeMaps);
            dt /= ArraySize(edgeMaps);
            if (wn == 0.0 && wt == 0.0 && dn == 0.0 && dt == 0.0)
                continue;

            weight1 += wn;
            xp1 += dn * nx_b;
            yp1 += dn * ny_b;
            xp1 += wn * px0_b;
            yp1 += wn * py0_b;

            weight2 += wt;
            xp2 += dt * tx_b;
            yp2 += dt * ty_b;
            xp2 += wt * px0_b;
            yp2 += wt * py0_b;
        }
        Bink r((weight1 != 0) ? xp1 / weight1 * 0.5 : 0, (weight1 != 0) ? yp1 / weight1 * 0.5 : 0);
        r.isFloating = true;
        r.frameId = config.frameId;
        r.x += (weight2 != 0) ? xp2 / weight2 * 0.5 : 0;
        r.y += (weight2 != 0) ? yp2 / weight2 * 0.5 : 0;
        r.px += (weight2 != 0) ? xp2 / weight2 * 0.5 : 0;
        r.py += (weight2 != 0) ? yp2 / weight2 * 0.5 : 0;
        return r;
    }

    Bink &getNearestBink(Array(Bink) &binks, const double x, const double y,
                         const ImageConfig &config) { // todo better (no references please)
        double currentNearest = 20;
        static Bink tempBink;
        Bink *nearestBink = &tempBink;
        for (int i = 0; i < binks.size(); ++i) {
            if (distance(x, y, binks[i]) < currentNearest) {
                currentNearest = distance(x, y, binks[i]);
                nearestBink = &binks[i];
            }
        }
        return *nearestBink;
    }

    void clearFrameBinks(Array(Bink) &binks, const int baseFrameId) {
        for (int i = binks.size() - 1; i >= 0; --i) {
            if (binks[i].frameId == baseFrameId)
                ArrayRemove(binks, i);
        }
    }

    const int getLastBinkIdInFrame(const Array(Bink) &binks, const int baseFrameId) {
        int l = 0;
        for (int i = binks.size() - 1; i >= 0; --i) {
            if (binks[i].frameId == baseFrameId)
                if (binks[i].id > l)
                    l = binks[i].id;
        }
        return l;
    }

    const int getLastBinkId(const Array(Bink) &binks) {
        int l = 0;
        for (int i = binks.size() - 1; i >= 0; --i) {
                if (binks[i].id > l)
                    l = binks[i].id;
        }
        return l;
    }

    Array(Bink) getRecapturedBinksOfFrame(const Image &image, Array(Bink) &binks, const std::vector<Feature> &features, const ImageConfig &config) {
        Array(Bink) r;
        for (int i = 0; i < binks.size(); ++i) {
            const Bink &b = binks[i];
            if (b.frameId == config.frameId) {
                Bink n = b;
                // maybe look for id?
                n.x = n.px;
                n.y = n.py;
                updateBink(image,features,n,config);
                ArrayAdd(r, n);
            }
        }
        return r;
    }

    void replaceFrameId(Array(Bink) &binks, const int sourceFrameId, const int destFrameId) {
        for (int i = binks.size() - 1; i >= 0; --i) {
            if (binks[i].frameId == sourceFrameId) {
                binks[i].frameId = destFrameId;
            }
        }
    }

#define IDENTICALBINK(__b1,__b2) ( ((__b1).frameId == (__b2).frameId) && ((__b1).id == (__b2).id) )
    void replaceAddBinks(Array(Bink) &binks, const Array(Bink) &recapturedBinks) {
        for (int j = recapturedBinks.size() - 1; j >= 0; --j) {
            for (int i = binks.size() - 1; i >= 0; --i) {
                if (IDENTICALBINK(binks[i],recapturedBinks[j])) {
                    ArrayRemove(binks,i);
                }
            }
        }
        ArrayAppend(binks,recapturedBinks);
    }

    Array(Bink) getClusterBinks(const Array(Feature) &features, double forceThreshold, const ImageConfig &config) {
        Array(Bink) r;
        Array(int) featureList; ArrayResize(featureList, ArraySize(features)); int i = 0; for (int j = 0; j < features.size(); ++j) {
            if (features[j].force >= forceThreshold) {
                featureList[i] = j;
                i++;
            }
        }
        ArrayResize(featureList, i);

            // accelrator structures not uesed, yet (too few features)
            //Array(int) xOrdered;
            //Array(int) yOrdered;
            //Array(int) toXOrdered;
            //Array(int) toYOrdered;
            //ArrayResize(xOrdered, ArraySize(features));
            //ArrayResize(yOrdered, ArraySize(features));
            //ArrayResize(toXOrdered, ArraySize(features));
            //ArrayResize(toYOrdered, ArraySize(features));
            //for (int j = 0; j < features.size(); ++j) {
            //    xOrdered[j] = j;
            //    yOrdered[j] = j;
            //    toXOrdered[j] = j;
            //    toYOrdered[j] = j;
            //}
            //std::sort(xOrdered.begin(), xOrdered.end(), [&](const int a, const int b) -> bool {
            //    return features[a].x < features[b].x;
            //});
            //std::sort(yOrdered.begin(), yOrdered.end(), [&](const int a, const int b) -> bool {
            //    return features[a].y < features[b].y;
            //});
            //for (int i = 0; i < ArraySize(features); ++i) {
            //    toXOrdered[xOrdered[i]] = i;
            //    toYOrdered[yOrdered[i]] = i;
            //}

        int k = 0;
        while (!featureList.empty()) {
            const int randomFeatureListIndex = (int)floor(randomLike(124234+k*19, featureList.size())); k++;
            const Feature &f1 = features[featureList[randomFeatureListIndex]];
            Bink b(f1.x, f1.y);
            double weight = 0;
            for (int i = 0; i < ArraySize(featureList); ++i) {
                const Feature &f2 = features[featureList[i]];
                const double d = distance(f2.x, f2.y, b);
                weight += DEFINED_CLUSTER_ATTENUATION_FORMULA; // am anfang riesenblob because very large number here
            }
            weight /= featureList.size(); // add all cones together to get sort of "bigger" (cone) blob shape
            weight = weight == 0.0 ? 0.0 : 1.0 / weight; // riesenblob is very small number here then
            bool someFound = false;
            for (int i = 0; i < ArraySize(featureList); ++i) {
                const Feature &f2 = features[featureList[i]];
                double d =  distance(f2.x, f2.y, b);
                double a = 1.0 / DEFINED_CLUSTER_ATTENUATION_FORMULA; // small number means wider away
                if (a - weight > binkDistanceAttenuationThreshold) { // wider away should be bigger than riesenblob to suffice (weight is to floor height and a is the feature height)
                    someFound = true;
                    ArrayAdd(b.minkies, f2);
                    ArrayRemove(featureList,i);
                    --i;
                }
                k++;
            }
            if (!someFound)
                break;
            Feature minkCenter = getMinkCenter(b);
            b.x = minkCenter.x; b.px = b.x;
            b.y = minkCenter.y; b.py = b.y;
            b.frameId = config.frameId;
            b.id = getLastBinkIdInFrame(r,config.frameId) + 1;
            r.push_back(b);
        }
        printf("%d, %d\n", r.size(), ArraySize(featureList));
        return r;
    }

    //-----------------
    namespace Implementation {
        Array(Bink) binks;
        Array(Bink) metaBinks;
    }
    using namespace Minks::Implementation;
    //-----------------
}
using namespace Minks;
//-----------------


//-----------------
namespace ImageLookup {

    //...
    class ImageFeature {
    public:
        Feature position;
    };

    ImageFeature tracked;

    namespace IMAGELOOKUP {
        const int pairCount = 200;
        const double stepFactor = 0.5;
        const double refinementSteps = 1000;
        const double mipsToTakeIntoAccountRatio = 1.0;
    };

    FeatureConfig createImageLookupFeatureConfig(int levels) {
        FeatureConfig config = createFeatureConfig(levels, IMAGELOOKUP::pairCount, false, IMAGELOOKUP::stepFactor, IMAGELOOKUP::refinementSteps, IMAGELOOKUP::mipsToTakeIntoAccountRatio);
        return config;
    }

    void scaleDescriptorPattern(DescriptorDesc &pattern, const double scale = 1.0) {
        for (int i = 0; i < pattern.levels.size(); ++i) {
            for (int j = 0; j < pattern.levels[i].pairs.size(); ++j) {
                DescriptorPair &pair = pattern.levels[i].pairs[j];
                pair.x0 *= scale;
                pair.y0 *= scale;
                pair.x1 *= scale;
                pair.y1 *= scale;
            }
        }
    }

    ImageFeature lookupImagePatch(const Image &patchImage, const Image &baseImage, const FeatureConfig &config) {
        FeatureConfig lookupConfig = config;
#if false //this stuff should be removed actually its not needed since we do have scale invariance
        // can be used as sorta hint for scale invariance but not needed
        lookupConfig.sourceToDestScaleRatio = (double)patchImage.featureLayers[0].fwidth / baseImage.featureLayers[0].fwidth;
#endif
        Feature f1, f2;
        ImageFeature r;
        Descriptor desc;
        f1.x = patchImage.featureLayers[0].width * 0.5;
        f1.y = patchImage.featureLayers[0].height * 0.5;
        desc = sampleDescriptor(patchImage, f1, lookupConfig);
        f1.x = baseImage.featureLayers[0].width * 0.5;
        f1.y = baseImage.featureLayers[0].height * 0.5;
        f2 = lookupFeature(baseImage, desc, f1, lookupConfig);
        r.position = f2;
        return r;
    }

}
using namespace ImageLookup;
//-----------------

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------

//-----------------
namespace Tracking {

    //-----------------
    namespace TRACKING {
        const int pairCount = pairCount;
        const double refinementSteps = refinementSteps;
        const double stepFactorIncrease = stepFactorIncrease;
        const int definedPyramidCoarseLevel = definedPyramidCoarseLevel;
        const bool iterative = iterative;
        const double refinementStepPerLevelScaleFactor = refinementStepPerLevelScaleFactor;
        const double descriptorLevelScale = descriptorLevelScale;
        const double mipsToTakeIntoAccountRatio = mipsToTakeIntoAccountRatio;
    }
    //using namespace TRACKING;
    //-----------------

    Array(Feature) trackFeatures(const Array(Feature) &features, const Image &frame1, const Image &frame2, const FeatureConfig &config) {
        Array(Feature) r;
        for (int i = 0; i < ArraySize(features); ++i) {
            const Feature& f = features[i];
            Descriptor descriptor = sampleDescriptor(frame1, f, config);
            Feature feature = lookupFeature(frame2, descriptor, f, config);
            ArrayAdd(r, feature);
        }
        return r;
    }
}
using namespace Tracking;
//-----------------

//-----------------
namespace Classification {

    //-----------------
    namespace Correspondence {

        //...
        class FeatureCorrespondence {
        public:
            int id[2];
        };

        //-----------------
        namespace Iterative {
            Array(FeatureCorrespondence) buildCorrespondence(const Array(Feature)& frame1, const Array(Feature)& frame2) {
                Array(FeatureCorrespondence) r;
                for (int i = 0; i < ArraySize(frame1); ++i) {
                    FeatureCorrespondence v;
                    v.id[0] = i;
                    v.id[1] = i;
                    ArrayAdd(r, v);
                }
                return r;
            }
        }
        using namespace Iterative;
        //-----------------
    }
    using namespace Correspondence;
    //-----------------

    //-----------------
    namespace Maths {

        //...
        class Position2D {
        public:
            double x,y;
        };

        //...
        class Matrix2D {
        public:
            double rotation[4];
            Position2D translation;
        };

        //...
        class Rectangle {
        public:
            double x0 = 0, y0 = 0, x1 = 0, y1 = 0;
        };

        bool hasArea(const Rectangle &r) {
            return (r.x0 < r.x1) && (r.y0 < r.y1);
        }
    }
    using namespace Maths;
    //-----------------

    //-----------------
    namespace Simple {

        //...
        class MovementClass {
        public:
            double angularMovement = 0;
            double deltaX = 0, deltaY = 0;
            double speed = 0;
            double angularClassification = 0;
            double translationalClassification = 0;
            double zoomClassification = 0;
            double classificationValue = 0;
            double classificationStrength = 0;
            double mat2d[4];
        };

        MovementClass classifyMovement(const Feature& f1, const Feature& f2) {
            static int k = 0; k++;
            MovementClass r;
            r.deltaX = f2.x - f1.x;
            r.deltaY = f2.y - f1.y;
            r.angularMovement = atan2(r.deltaX, r.deltaY);
            r.speed = sqrt(r.deltaX * r.deltaX + r.deltaY * r.deltaY);
            return r;
        }

        //...
        class FeatureClass {
        public:
            int byFeature;
            int byFrame;
            bool byOrder;
            MovementClass byMovement;
            Feature sourceFeature;
            Feature destFeature;
            Matrix2D transformation; // put later somewhere else
        };

        FeatureClass classify(const Feature& f1, const Feature& f2, bool forward = true) {
            FeatureClass r;
            r.sourceFeature = f1;
            r.destFeature = f2;
            r.byFrame = f1.frame;
            r.byFeature = f1.id;
            r.byOrder = forward;
            r.byMovement = classifyMovement(f1, f2);
            return r;
        }

        Array(FeatureClass) cleanFeatures(const Array(FeatureClass)& v) {
            Array(FeatureClass) r;
            for (int i = 0; i < ArraySize(v); ++i) {
                const MovementClass& m = v[i].byMovement;
                if (VALIDFLOAT(m.deltaX) && VALIDFLOAT(m.deltaY))
                    ArrayAdd(r, v[i]);
            }
            return r;
        }

        Array(FeatureClass) angularClassification(const Array(FeatureClass)& v) {
            Array(FeatureClass) r = v;
            for (int i = 0; i < ArraySize(r); ++i) {
                FeatureClass& c = r[i];
                const double kx = c.byMovement.deltaX / c.byMovement.speed;
                const double ky = c.byMovement.deltaY / c.byMovement.speed;
// somehow provides a almost linear factor between two directly succeeding angles on a circle (found just by trying and thinking whilst trying)
// even working on different zoom scales, so seems to be completely ok with all i tested
//#define POSITIONTOVALUE(__x,__y)  ( (__x) / ((__y)+sign(__y)*3.0)*2.0*fabs(__x) - fabs(__y) / fabs((__x)+sign(__x)*1.5)*0.5/fabs(__y) ) // the sign stuff is to prevent division by zero singularity
// huiuiui!!!!!!! (at least it's enough for the use case in this ESBIFT classifier here)
//                for (int i = 0; i < 100; ++i) {
//                    const double ri = i / 100.0;
//                    const double kx2 = (random(2.0) - 1.0) * ri;
//                    const double ky2 = (random(2.0) - 1.0) * ri;
//                    double r1, r2;
//                    {
//                        const double kx = kx2;
//                        const double ky = ky2;
//                        r1 = LINEARANGLE;
//                    }
//                    {
//                        const double kx = kx2 * 0.9;
//                        const double ky = ky2 * 0.9;
//                        r2 = LINEARANGLE;
//                    }
//                    printf("%f %f %f\n", r1 / r2, r1, r2);
//                }

                c.byMovement.classificationValue = POSITIONTOVALUE(kx,ky);
            }
            // sort by ascending "angle like" values
            std::sort(r.begin(), r.end(), [&](const FeatureClass& a, const FeatureClass& b) -> bool {
                return fabs(b.byMovement.classificationValue - a.byMovement.classificationValue);});
            for (int i = 1; i < ArraySize(r); ++i) {
                r[i-1].byMovement.classificationStrength = fabs(r[i].byMovement.classificationValue - r[i-1].byMovement.classificationValue);
            }
            if (ArraySize(r) > 2) r[ArraySize(r) - 1] = r[ArraySize(r) - 2];
            // sort by nearest pairs (of the ascending) "angle like" values
            std::sort(r.begin(), r.end(), [&](const FeatureClass& a, const FeatureClass& b) -> bool {
                return a.byMovement.classificationStrength < b.byMovement.classificationStrength;
                });
            double d = 0.0;
            for (int i = 0; i < ArraySize(r); ++i) {
                r[i].byMovement.angularClassification = d;
                d += r[i].byMovement.classificationStrength;
            }
            return r;
        }

        Array(FeatureClass) translationalClassification(const Array(FeatureClass)& v) {
            Array(FeatureClass) r = v;
            for (int i = 0; i < ArraySize(r); ++i) {
                FeatureClass& c = r[i];
                c.byMovement.classificationValue = sqrt(c.byMovement.deltaX * c.byMovement.deltaX + c.byMovement.deltaY * c.byMovement.deltaY);
            }
            // sort by ascending distance values
            std::sort(r.begin(), r.end(), [&](const FeatureClass& a, const FeatureClass& b) -> bool {
                const double d = fabs(b.byMovement.classificationValue - a.byMovement.classificationValue);
                return d;
                });
            for (int i = 1; i < ArraySize(r); ++i) {
                FeatureClass& c0 = r[i - 1];
                const FeatureClass& c1 = r[i];
                double delta = c1.byMovement.classificationValue - c0.byMovement.classificationValue;
                c0.byMovement.classificationStrength = fabs(delta);
            }
            if (ArraySize(r) > 2) r[ArraySize(r) - 1] = r[ArraySize(r) - 2];
            // sort by nearest pairs (of the ascending) distance values
            std::sort(r.begin(), r.end(), [&](const FeatureClass& a, const FeatureClass& b) -> bool {
                return a.byMovement.classificationStrength < b.byMovement.classificationStrength;
                });
            double d = 0.0;
            for (int i = 0; i < ArraySize(r); ++i) {
                r[i].byMovement.translationalClassification = d;
                d += r[i].byMovement.classificationValue;
            }
            return r;
        }

        Array(FeatureClass) zoomClassification(const Array(FeatureClass)& v) {
            Array(FeatureClass) r = v;
            for (int i = 0; i < ArraySize(r); ++i) {
                FeatureClass& c = r[i];
                const double kx = c.byMovement.deltaX;
                const double ky = c.byMovement.deltaY;
                c.byMovement.classificationValue = POSITIONTOVALUE(kx,ky);
            }
            // sort by ascending distance values
            std::sort(r.begin(), r.end(), [&](const FeatureClass& a, const FeatureClass& b) -> bool {
                // the more negative the value the more zoom
                double d = b.byMovement.classificationValue * a.byMovement.classificationValue;
                return d;
                });
            for (int i = 1; i < ArraySize(r); ++i) {
                FeatureClass& c0 = r[i-1];
                const FeatureClass& c1 = r[i];
                double delta = c1.byMovement.classificationValue - c0.byMovement.classificationValue;
                c0.byMovement.classificationStrength = fabs(delta);
            }
            if (ArraySize(r) > 2) r[ArraySize(r) - 1] = r[ArraySize(r) - 2];
            // sort by nearest pairs (of the ascending) distance values
            std::sort(r.begin(), r.end(), [&](const FeatureClass& a, const FeatureClass& b) -> bool {
                return a.byMovement.classificationStrength < b.byMovement.classificationStrength;
                });
            double d = 0.0;
            for (int i = 0; i < ArraySize(r); ++i) {
                r[i].byMovement.zoomClassification = d;
                d += r[i].byMovement.classificationValue;
            }
            return r;
        }

        Position2D normalizedDeltas(const Array(FeatureClass)& r) {
            double dx = 0, dy = 0;
            for (int i = 0; i < ArraySize(r); ++i) {
                dx += r[i].byMovement.deltaX;
                dy += r[i].byMovement.deltaY;
            }
            const double l = sqrt(dx * dx + dy * dy);
            if (l > 0.0) {
                dx /= l; dy /= l;
            }
            return { dx, dy };
        }

        Array(FeatureClass) objectClassification(const Array(FeatureClass)& r) {
            Array(FeatureClass) v = r;

            // sort by movement classification (more separate feature movements in the back of the array and more connected ones in front)
            std::sort(v.begin(), v.end(), [&](const FeatureClass& a, const FeatureClass& b) -> bool {
                return a.byMovement.zoomClassification < b.byMovement.zoomClassification;
                });
            ArrayResize(v, ArraySize(v) * 0.5);

            // sort by angular classification
            std::sort(v.begin(), v.end(), [&](const FeatureClass& a, const FeatureClass& b) -> bool {
                return a.byMovement.angularClassification < b.byMovement.angularClassification;
                });
            // just take good angular values into account
            ArrayResize(v, ArraySize(v) * 0.375);

            // sort by movement classification (more separate feature movements in the back of the array and more connected ones in front)
            std::sort(v.begin(), v.end(), [&](const FeatureClass& a, const FeatureClass& b) -> bool {
                return a.byMovement.translationalClassification < b.byMovement.translationalClassification;
                });

            Position2D rot = normalizedDeltas(v);
            for (int i = 0; i < ArraySize(v); ++i) {
                Matrix2D mat = { rot.x,rot.y,rot.y,-rot.x, {100,100} };
                v[i].transformation = mat;
            }
            return v;
        }

        Array(Feature) prune(const Array(FeatureClass)& r) {
            Array(FeatureClass) v = r;

            // sort by movement classification (more separate feature movements in the back of the array and more connected ones in front)
            std::sort(v.begin(), v.end(), [&](const FeatureClass& a, const FeatureClass& b) -> bool {
                return a.byMovement.zoomClassification < b.byMovement.zoomClassification;
                });
            ArrayResize(v, ArraySize(v) * pruneByZoomLevelRatio);
            
            // sort by angular classification
            std::sort(v.begin(), v.end(), [&](const FeatureClass& a, const FeatureClass& b) -> bool {
                return a.byMovement.angularClassification < b.byMovement.angularClassification;
                });
            // just take good angular values into account
            ArrayResize(v, ArraySize(v) * pruneByAngularRatio);

            // sort by movement classification
            std::sort(v.begin(), v.end(), [&](const FeatureClass& a, const FeatureClass& b) -> bool {
                return a.byMovement.translationalClassification < b.byMovement.translationalClassification;
                });
            // just take good angular values into account
            ArrayResize(v, ArraySize(v) * pruneByTranslationalRatio);

            Array(Feature) f;
            for (int i = 0; i < ArraySize(v); ++i) {
                if (v[i].byOrder) {
                    ArrayAdd(f, v[i].sourceFeature);
                }
            }
            return f;
        }

        Array(FeatureClass) classify(const Array(Feature)& features1, const Array(Feature)& features2, const Array(FeatureCorrespondence)& correspondence) {
            Array(FeatureClass) r;
            for (int i = 0; i < ArraySize(correspondence); ++i) {
                const FeatureCorrespondence& c = correspondence[i];
                const Feature& f1 = features1[c.id[0]];
                const Feature& f2 = features2[c.id[1]];
                FeatureClass cl;
                cl = classify(f1, f2, true); ArrayAdd(r, cl);
                cl = classify(f2, f1, false); ArrayAdd(r, cl);
            }
            r = cleanFeatures(r);
            r = angularClassification(r); // (cleanest describable rotations in front) features with angles which are morly available in the scene are in front
            r = translationalClassification(r); // (cleanest describable movements in front) features with tranlsation which are more available in the scene in front
            r = zoomClassification(r); // (cleanest describable zooms in front) (seems like objects may be classified through their zoom values)
            r = objectClassification(r);
            return r;
        }
    }
    using namespace Simple;
    //-----------------

    Array(FeatureClass) classify(const Image& frame1, const Image& frame2, const Array(Feature) &features1, const FeatureConfig& config) {
        Array(Feature) features2 = trackFeatures(features1, frame1, frame2, config);
        Array(FeatureCorrespondence) correspondence = buildCorrespondence(features1, features2);
        return classify(features1, features2, correspondence);
    }


    Array(Feature) prune(const Array(Feature)& features1, const Array(Feature)& features2, const Array(FeatureCorrespondence)& correspondence) {
        Array(FeatureClass) r;
        for (int i = 0; i < ArraySize(correspondence); ++i) {
            const FeatureCorrespondence& c = correspondence[i];
            const Feature& f1 = features1[c.id[0]];
            const Feature& f2 = features2[c.id[1]];
            FeatureClass cl = classify(f1, f2); 
            ArrayAdd(r, cl);
        }
        r = cleanFeatures(r);
        r = angularClassification(r); // (cleanest describable rotations in front) features with angles which are morly available in the scene are in front
        r = translationalClassification(r); // (cleanest describable movements in front) features with tranlsation which are more available in the scene in front
        r = zoomClassification(r); // (cleanest describable zooms in front) (seems like objects may be classified through their zoom values)
        return prune(r);
    }

    //-----------------
    namespace Implementation {
        FeatureConfig trackConfig;
        Array(FeatureClass) updatedFeatureClasses;

        FeatureConfig createTrackFeatureConfig(const int levels) {
            FeatureConfig config;
            config = createFeatureConfig(levels, TRACKING::pairCount, TRACKING::iterative, stepFactor, TRACKING::refinementSteps, TRACKING::mipsToTakeIntoAccountRatio);
            return config;
        }
    }
    using namespace Classification::Implementation;
    //-----------------
}
using namespace Classification;
//-----------------

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------

int baseImageCountForSaving = 0;
DOUBLE windowMouseStrength = 0.0;
int imagePyramidSize = -1;

DOUBLE dataSetContrastRandom = 0.0;
int datasetLoadingRandomSeed = 0;

void loadTrackingSet(const std::string& folder) {
    loading_mutex.lock();
    DEBUG("Load TrackingSet start\n")
    trackingFrames.clear();
    EXCEPTIONTRY {
        Array(std::string) fileNames;
        for (const auto &entry: FILESYSTEM::directory_iterator(folder)) {
            if (entry.is_regular_file()) {
                const std::string fileName = entry.path().string();
                ArrayAdd(fileNames, fileName);
            }
        }
        std::sort(fileNames.begin(), fileNames.end(), [&](std::string &a, std::string &b) -> bool {
            return a < b;
        });
        for (int i = 0; i < fileNames.size(); ++i) {
            ArrayAdd(trackingFrames, imageFromFileNameExtension(fileNames[i], 0, -1, 0.0, 1.0, desiredImageScale));
            DEBUG("loaded frame:%s\n", fileNames[i].c_str());
        }
        for (int i = 0; i < ArraySize(trackingFrames); ++i) {
            preprocessImage(trackingFrames[i]);
            setUpImagePyramid(trackingFrames[i]);
        }
    } EXCEPTIONCATCH(...) {
    }
    DEBUG("Load TrackingSet end\n")
    loading_mutex.unlock();
}

void loadDataSet(const std::string& folder, int startFrameNr = -1, int endFrameNr = -1, int stepSize = 1, bool timeReversed = false) {
    loading_mutex.lock();
    DEBUG("Load DataSet start\n")
    // mutex missing (of course) :)
    frames.clear();
    Array(std::string) fileNames;
    for (const auto& entry : FILESYSTEM::directory_iterator(folder)) {
        if (entry.is_regular_file()) {
            std::string fileName = entry.path().string();
            ArrayAdd(fileNames, fileName);
        }
    }
    std::sort(fileNames.begin(), fileNames.end(), [&](std::string& a, std::string& b)->bool {
        if (timeReversed)
            return b < a;
        return a < b;
        });
    startFrameNr = CLAMP(startFrameNr, 0, ArraySize(fileNames));
    endFrameNr = CLAMP(endFrameNr, 0, ArraySize(fileNames));
    if (startFrameNr > endFrameNr) startFrameNr = endFrameNr;
    srandom(datasetLoadingRandomSeed);
    for (int i = startFrameNr; i < endFrameNr; i += stepSize) {
        DOUBLE brightness = random(1.0) * 0.0;
        DOUBLE contrast = 1.0 - random(dataSetContrastRandom);
        if (i == startFrameNr) {contrast = 1; brightness = 0;}
        ArrayAdd(frames, imageFromFileNameExtension(fileNames[i], 0, useOwnWidth ? desiredImageWidth : -1, brightness, contrast, desiredImageScale));
        DEBUG("loaded frame:%s\n", fileNames[i].c_str());
    }
    for (int i = 0; i < ArraySize(frames); ++i) {
        preprocessImage(frames[i]);
        imagePyramidSize = setUpImagePyramid(frames[i]);
    }
    descriptivityLayer = 2;
    DEBUG("Load DataSet end\n")
    loading_mutex.unlock();
}

#if DEFINED_DESCRIPTIVITYMAP == 1
void buildDescriptivityMap(Image& image, int level = 0) {
    //DESCRIPTIVITYLayer layer = &image.descriptivityLayer.pixels;
    DESCRIPTIVITYLayer layer;
    ArrayResize(layer.pixels,0);
    JUSTSETUPIMAGEHULL(layer, image.featureLayers[level]);
    for (int y = 0; y < layer.height; ++y) {
        for (int x = 0; x < layer.width; ++x) {
            Feature f;
            f.x = (double)x * image.featureLayers[0].width / layer.width;
            f.y = (double)y * image.featureLayers[0].height / layer.height;
            threadPoolThread([&, f, x, y]() {
                initThread();
                Descriptor descriptor = sampleDescriptor(image, f, baseConfig);
                const DOUBLE strength = descriptivityOfDescriptor(descriptor, baseConfig);
                layer.pixels[x + y * layer.width] = strength;
                });
        }
        printf(":%02d", y * 100 / layer.height);
    }
    JUSTSETUPIMAGEHULL(image.descriptivityLayer, layer);
    for (int y = 1; y < layer.height - 1; ++y) {
        for (int x = 1; x < layer.width - 1; ++x) {
            const double kx0 = layer.pixels[(x - 1) + y * layer.width];
            const double kx1 = layer.pixels[(x + 1) + y * layer.width];
            const double ky0 = layer.pixels[x + (y - 1) * layer.width];
            const double ky1 = layer.pixels[x + (y + 1) * layer.width];
            const double dx = kx1 - kx0;
            const double dy = ky1 - ky0;
            const double d = sqrt(dx * dx + dy * dy);
            const double gradient = d;
            const double angle = atan2(dx, dy);
            const double color = fabs(dx);
            image.descriptivityLayer.pixels[x + (y + 1) * image.descriptivityLayer.width] = color;
        }
    }
    threadPoolFinish();
}
#endif

Array(Feature) placeCircularFeatures(const Image& image, const int radiusCount, const int sectorCount, const DOUBLE innerRadius = 0.15, const DOUBLE outerRadius = 0.95) {
    Array(Feature) allFeatures;
    for (int y = 0; y < radiusCount; ++y) {
        DOUBLE py = (DOUBLE)y / (radiusCount-1);
        DOUBLE r = pow(Math::lerp(innerRadius, outerRadius, py),0.75) * 0.5;
        for (int x = 0; x < sectorCount; ++x) {
            DOUBLE px = (DOUBLE)x / sectorCount + py;
            DOUBLE cx = (sin(px * 2.f * 3.1415927f) * r + 0.5) * image.featureLayers[0].width;
            DOUBLE cy = (cos(px * 2.f * 3.1415927f) * r + 0.5) * image.featureLayers[0].height;
            Feature k(cx, cy);
            ArrayAdd(allFeatures, k);
        }
    }
    return allFeatures;
}

Array(Feature) placeRandomFeatures(const Image& image, const int featureCount, const Rectangle &rectangle = Rectangle()) {
    Array(Feature) allFeatures;
    for (int i = 0; i < featureCount; ++i) {
        Feature k(random(image.featureLayers[0].width), random(image.featureLayers[0].height));
        if (rectangle.x0 < rectangle.x1 && rectangle.y0 < rectangle.y1) {
            k.x = random(rectangle.x1 - rectangle.x0) + rectangle.x0;
            k.y = random(rectangle.y1 - rectangle.y0) + rectangle.y0;
        }
        ArrayAdd(allFeatures, k);
    }
    return allFeatures;
}

}
using namespace ESBIFT;
//-----------------

using UPDATEREQUEST = unsigned int;

const UPDATEREQUEST TEXTURESUPDATE = 1;
const UPDATEREQUEST FEATURERECALC = 2;
const UPDATEREQUEST LASTFRAMEUPDATE = 4;
const UPDATEREQUEST SAVECURRENTFRAME = 8;
const UPDATEREQUEST RELOAD = 16;
const UPDATEREQUEST WINDOWFREEZE = 32;
const UPDATEREQUEST READDFINISHEDFEATURES = 64;
const UPDATEREQUEST FULLFEATURERECALC = 128;
const UPDATEREQUEST REDOCLASSIFICATION = 256;
const UPDATEREQUEST DOPRUNING = 512;
const UPDATEREQUEST PLACEFEATURES = 1024;
const UPDATEREQUEST DESCRIBEFEATURES = 2048;
const UPDATEREQUEST PLACEFEATURESRESAMPLE = 4096;
const UPDATEREQUEST RECTANGULARSELECTION = 8192;
const UPDATEREQUEST LOCATEPATCH = 16384;
const UPDATEREQUEST JUSTLASTFEATURE = 32768;
const UPDATEREQUEST EXPLICITUPDATEORIGINAL = 0x10000;
const UPDATEREQUEST REMOVEUNSUFFICIENTFEATURES = 0x20000;
const UPDATEREQUEST BINKSNEEDMOD = 0x40000;

UPDATEREQUEST requestedUpdates = 0;

void requestUpdate(const UPDATEREQUEST updateType) {
    requestedUpdates |= updateType;
}

void clearUpdate(const UPDATEREQUEST updateType) {
    requestedUpdates &= ~(updateType);
}

bool shouldUpdate(const UPDATEREQUEST updateType) {
    return (requestedUpdates & updateType) != 0;
}

namespace g {
 
    bool checkbox(UPDATEREQUEST uprq, const char* label, bool* v) {
        if (Checkbox(label, v)) {
            requestUpdate(uprq);
            return true;
        }
        return false;
    }

    bool inputInt(UPDATEREQUEST uprq, const char* label, int* v) {
        if (InputInt(label, v)) {
            requestUpdate(uprq);
            return true;
        }
        return false;
    }

    bool inputDOUBLE(UPDATEREQUEST uprq, const char* label, DOUBLE* v, DOUBLE step = 0.0) {
        if (InputFloat(label, v, step)) {
            requestUpdate(uprq);
            return true;
        }
        return false;
    }

    bool button(UPDATEREQUEST uprq, const char* label) {
        if (Button(label)) {
            requestUpdate(uprq);
            return true;
        }
        return false;
    }

    void description(const DISPLAYLayer& layer) {
        g::Text("width:%d height:%d", layer.width, layer.height);
    }

    void imageLayer(const DISPLAYLayer& layer, const ImageConfig &config) {
        int id = layer.handle;
        g::Image((void*)id, ImVec2(layer.width * config.zoom, layer.height * config.zoom));
    }

    void features(const Array(Feature)& features, const ImageConfig& config) {
        ImVec2 w = g::GetWindowPos();
        w.x -= g::GetScrollX();
        w.y -= g::GetScrollY();
        DOUBLE r = 10;
        PIXEL c0 = fromRGBA32(0xff000000);
        PIXEL c1 = fromRGBA32(0xffffffff);
        PIXEL c2 = fromRGBA32(0xff0000ff);
        const double z = config.zoom;
        for (const auto& f : features) {
            g::GetWindowDrawList()->AddCircle(ImVec2((f.x * z + w.x), (f.y * z + w.y + 2)), r * z, toRGBA32(c0));
            if (f.force >= Configuration::pruneThreshold) {
                g::GetWindowDrawList()->AddCircle(ImVec2((f.x * z + w.x), (f.y * z + w.y)), r * z, toRGBA32(c1));
                g::GetWindowDrawList()->AddCircle(ImVec2((f.x * z + w.x), (f.y * z + w.y)), r * f.force * z, toRGBA32(c2));
            } else {
                double rz = r * z * 0.5;
                g::GetWindowDrawList()->AddLine(ImVec2((f.x * z + w.x + rz), (f.y * z + w.y - rz)), ImVec2((f.x * z + w.x - rz), (f.y * z + w.y + rz)), toRGBA32(c2));
                g::GetWindowDrawList()->AddLine(ImVec2((f.x * z + w.x - rz), (f.y * z + w.y - rz)), ImVec2((f.x * z + w.x + rz), (f.y * z + w.y + rz)), toRGBA32(c2));
            }
        }
    }

    void mink(const Feature &f, int minkId, const std::string &minkType, double radius, const ImageConfig& config) {
        ImVec2 w = g::GetWindowPos();
        w.x -= g::GetScrollX();
        w.y -= g::GetScrollY();
        unsigned int color = toRGBA32(config.color);
        PIXEL c1 = fromRGBA32(color & 0xff000000);
        const double z = config.zoom;
        for (int i = 0; i < 2; ++i) {
            for (double r = 1; r < radius; r *= 1.3) {
                g::GetWindowDrawList()->AddCircle(ImVec2((f.x * z + w.x), (f.y * z + w.y)), (r+i) * z, toRGBA32(c1));
            }
            c1 = config.color;
        }
        g::SetCursorPos(ImVec2((f.x * z + z * radius * 0.5),(f.y * z + z * radius * 0.5 + 1)));
        g::PushStyleColor(ImGuiCol_Text,(color & 0xff000000) | 0x60000000);
        g::Text("%04d", minkId);
        g::PopStyleColor();
        g::SetCursorPos(ImVec2((f.x * z + z * radius * 0.5),(f.y * z + z * radius * 0.5)));
        g::PushStyleColor(ImGuiCol_Text,(color & 0xff000000) | 0x608080ff);
        g::Text("%04d", minkId);
        g::PopStyleColor();
        g::SetCursorPos(ImVec2(0,0));
    }

    void bink(const Bink &b, const ImageConfig& config) {
        if (INVALIDBINK(b))
            return;
        Feature f(b.px + g::GetCursorPosX(), b.py + g::GetCursorPosY());
        ImageConfig bc = config;
        bc.color = b.selected ? fromRGBA32(0xff00000ff) : (b.frameId == config.frameId ? fromRGBA32(0xffffffff) : fromRGBA32(0x30600000));
        double radius = 16;
        if (b.isFloating) {
            bc.color = fromRGBA32(0xff80ff20);
            radius = 32;
        }
        if (b.isMetaBink) {
            bc.color.r = 0.5;
            bc.color.g = 1.0;
            bc.color.b = 2.0;
            bc.color.a = 1.0;
        }
        mink(f, b.id, "", radius,bc);
    }

    void feature(const Feature& feature, const double r, const ImageConfig& config) {
        ImVec2 w = g::GetWindowPos();
        w.x -= g::GetScrollX();
        w.y -= g::GetScrollY();
        PIXEL c1 = fromRGBA32(0xffffffff);
        PIXEL c2 = fromRGBA32(0xff0000ff);
        const double z = config.zoom;
        g::GetWindowDrawList()->AddCircle(ImVec2((feature.x * z + w.x), (feature.y * z + w.y)), r * z, toRGBA32(c1));
        g::GetWindowDrawList()->AddCircle(ImVec2((feature.x * z + w.x), (feature.y * z + w.y)), r * 0.9 * z, toRGBA32(c2));
    }

    void withRotationMatrix(const Position2D& pos, const Matrix2D& matrix, const double scale, const ImageConfig& config) {
        const double z = config.zoom;
        ImVec2 w = g::GetWindowPos();
        w.x -= g::GetScrollX();
        w.y -= g::GetScrollY();
        const double m00 = matrix.rotation[0] * scale * z * config.imageScale;
        const double m01 = matrix.rotation[1] * scale * z * config.imageScale;
        const double m10 = matrix.rotation[2] * scale * z * config.imageScale;
        const double m11 = matrix.rotation[3] * scale * z * config.imageScale;
        PIXEL c = fromRGBA32(0xffffffff);
        const double mx = pos.x * z + w.x;
        const double my = pos.y * z + w.y;
        const double lx = mx - m00 * 0.5 - m10 * 0.5;
        const double ly = my - m01 * 0.5 - m11 * 0.5;
        g::GetWindowDrawList()->AddLine(ImVec2(mx + m00 * 0.5, my + m01 * 0.5), ImVec2(mx, my), toRGBA32(c), 3.f);
        g::GetWindowDrawList()->AddLine(ImVec2(lx, ly), ImVec2(lx + m00, ly + m01), toRGBA32(c), 3.f);
        g::GetWindowDrawList()->AddLine(ImVec2(lx, ly), ImVec2(lx + m10, ly + m11), toRGBA32(c), 3.f);
        g::GetWindowDrawList()->AddLine(ImVec2(lx + m00, ly + m01), ImVec2(lx + m00 + m10, ly + m11 + m01), toRGBA32(c), 3.f);
        g::GetWindowDrawList()->AddLine(ImVec2(lx + m10, ly + m11), ImVec2(lx + m00 + m10, ly + m11 + m01), toRGBA32(c), 3.f);
    }

    void imageFeature(const Feature& feature, const double scale, const ESBIFT::Image &image, const ImageConfig& config) {
        DISPLAYLayer layer = toDisplay(image.pixelLayers[0], config);
        const int id = layer.handle;
        const double z = config.zoom;
        const double px = (feature.x - layer.width * 0.5 * config.imageScale) * z;
        const double py = (feature.y - layer.height * 0.5 * config.imageScale) * z;
        g::SetCursorPos(ImVec2(px,py));
        g::PushStyleVar(ImGuiStyleVar_Alpha,config.imageAlpha);
        g::Image((void*)id, ImVec2(layer.width * z * config.imageScale, layer.height * z * config.imageScale));
        g::SetCursorPos(ImVec2(px,py));
        g::Text("%02.02f", feature.scale);
        g::SetCursorPos(ImVec2(0,0));
        g::PopStyleVar();
    }


    void withRotationMatrix(const Feature& feature, const Matrix2D& matrix, const double scale, const ImageConfig& config) {
        Position2D p;
        p.x = feature.x;
        p.y = feature.y;
        withRotationMatrix(p, matrix, scale, config);
    }

    void imageFeature(const Feature& feature, const ESBIFT::Image &image, const ImageConfig& config) {
        Matrix2D matrix;
        matrix.rotation[0] = image.featureLayers[0].fwidth;
        matrix.rotation[1] = 0;
        matrix.rotation[2] = 0;
        matrix.rotation[3] = image.featureLayers[0].fheight;
        ImageConfig config2 = config;
        config2.imageAlpha = Configuration::imageFeatureAlpha;
        // config2.forceUpdate = true; //not needed
        config2.imageScale = feature.scale;
        withRotationMatrix(feature, matrix, 1.0, config2);
        imageFeature(feature, 1.0, image, config2);
    }

    void classifications(const Array(FeatureClass)& classes, const ImageConfig& config) {
        ImVec2 w = g::GetWindowPos();
        w.x -= g::GetScrollX();
        w.y -= g::GetScrollY();
        PIXEL c1 = fromRGBA32(0xffffffff);
        PIXEL c1b = fromRGBA32(0x80ffffff);
        PIXEL c2 = fromRGBA32(0xff0000ff);
        PIXEL cb = fromRGBA32(0xff000000);
        const double z = config.zoom;
        for (const auto& f : classes) {
            DOUBLE r = 7;
            DOUBLE ra = r * 2.0;
            PIXEL c1 = fromRGBA32(0xffffffff);
            if (!f.byOrder) { // currently its from dest frame to source frame
                g::GetWindowDrawList()->AddCircle(ImVec2((f.sourceFeature.x * z + w.x), (f.sourceFeature.y * z + w.y)), r * z, toRGBA32(c1b));
                g::GetWindowDrawList()->AddCircle(ImVec2((f.destFeature.x * z + w.x), (f.destFeature.y * z + w.y)), r * z, toRGBA32(c1));
                double angularX = f.byMovement.deltaX;
                double angularY = f.byMovement.deltaY;
                double l = sqrt(angularX * angularX + angularY * angularY);
                angularX /= l;
                angularY /= l;
                angularX *= 24;
                angularY *= 24;
                g::GetWindowDrawList()->AddLine(
                    ImVec2((f.sourceFeature.x * z + w.x), (f.sourceFeature.y * z + w.y - 1)),
                    ImVec2((f.sourceFeature.x * z + w.x + angularX), (f.sourceFeature.y * z + w.y - 1) + angularY),
                    toRGBA32(cb), 3.f);
                g::GetWindowDrawList()->AddLine(
                    ImVec2((f.sourceFeature.x * z + w.x), (f.sourceFeature.y * z + w.y)),
                    ImVec2((f.sourceFeature.x * z + w.x + angularX), (f.sourceFeature.y * z + w.y) + angularY),
                    toRGBA32(c2), 3.f);
            }
        }
        if (!classes.empty()) {
            Feature f = classes[0].sourceFeature;
            Matrix2D m = classes[0].transformation;
            Position2D p = m.translation;
            withRotationMatrix(p, m, 20.0, config);
        }
    }

    bool listBox(const char* label, int* v, const Array(std::string)& elements) {
        Array(const char*) nativeElements; // not sure if this is really a "good" idea (it's at home and no risk for anyone! :))
        // would have to write a propper wrapper for that later...
        for (int i = 0; i < ArraySize(elements); ++i) {
            ArrayAdd(nativeElements,(&(elements[i][0])));
        }
        *v = CLAMP(*v, 0, (ArraySize(nativeElements) - 1));
        return g::ListBox(label, v, &nativeElements[0], ArraySize(nativeElements));
    }

    bool multiListBox(int mid, const char *label, int *v, const Array(Array(std::string)) elements) {
        return g::listBox(label, v, elements[CLAMP(mid, 0, ArraySize(elements) - 1)]);
    }

    void featureDescs(const Array(Feature)& features) {
        bool f2 = true;
        for (const auto& f : features) {
            if (!f2) g::SameLine();
            f2 = false;
            g::Text("[%.03f:%.03f:%.03f]", f.x, f.y, f.strength);

        }
    }

    void pairs(const DISPLAYLayer &layer, const int shapeId, const FeatureConfig &featureConfig, const ImageConfig& config) {
        ImVec2 w = g::GetWindowPos();
        w.x -= g::GetScrollX();
        w.y -= g::GetScrollY();
        w.x += g::GetCursorPosX();
        w.y += g::GetCursorPosY();
        const PIXEL c = fromRGBA32(0xff0000ff);
        const double z = config.zoom;
        const double wx = layer.width * 0.5;
        const double wy = layer.height * 0.5;
        w.x += wx * z;
        w.y += wy * z;
        const DescriptorDesc &pattern = featureConfig.descriptorPatterns[shapeId];
        for (int j = 0; j < pattern.levels.size(); ++j) {
            for (int i = 0; i < pattern.levels[j].pairs.size(); ++i) {
                DescriptorPair p = pattern.levels[j].pairs[i];
#if DEFINED_POLARCOORDINATES == 1
                fromPolar(p);
#endif
                g::GetWindowDrawList()->AddLine(ImVec2((p.x0 * wx * z + w.x), (p.y0 * wy * z + w.y)), ImVec2((p.x1 * wx * z + w.x), (p.y1 * wy * z + w.y)), toRGBA32(c));
            }
        }
    }

    void overlayPoint(DOUBLE mx, DOUBLE my) {
        DOUBLE r = 10;
        PIXEL c = fromRGBA32(0xffffffff);
        g::GetForegroundDrawList()->AddCircle(ImVec2(mx, my), r, toRGBA32(c));
    }

    //...
    class Poly {
    public: int id[3];
    };

    void featureMesh(const Array(Feature) &features, const ImageConfig& config) {
        Array(Poly) polies;
        ArrayResize(polies, ArraySize(features));
        for (int i = 0; i < ArraySize(features); ++i) {

        }
        ImVec2 w = g::GetWindowPos();
        w.x -= g::GetScrollX();
        w.y -= g::GetScrollY();
        const double z = config.zoom;
        createNeighborHood(features);
        PIXEL c = fromRGBA32(0xc0ffffff);
        std::set<int> alreadyFoundOnes;
        for (int i = 0; i < ArraySize(features); ++i) {
            const Feature here = features[i];
            const int n = nearestNeighbor(features,i,alreadyFoundOnes);
            const Feature next = features[n];
            alreadyFoundOnes.insert(i);
            g::GetWindowDrawList()->AddLine(ImVec2((here.x * z + w.x), (here.y * z + w.y)), ImVec2((next.x * z + w.x), (next.y * z + w.y)), toRGBA32(c), 3);
        }
    }
    
#define SELECTION_AWAIT 1
#define SELECTION_SELECTING 2
#define SELECTION_FINISHED 3
#define SELECTION_DONE 4
    int doRectangularSelection(UPDATEREQUEST uprq, Maths::Rectangle& rect, int &state, const PIXEL &color = fromRGBA32(0xf0ff00ff), int button = 0) {
        // just for fun..
        ImVec2 wp = g::GetWindowPos();
        wp.x -= g::GetScrollX();
        wp.y -= g::GetScrollY();
        if (g::IsWindowHovered()) {
            switch (state) {
            case SELECTION_AWAIT: {
                if (g::IsMouseDown(button)) {
                    requestUpdate(WINDOWFREEZE);
                    rect.x0 = g::GetMousePos().x - wp.x;
                    rect.y0 = g::GetMousePos().y - wp.y;
                    rect.x1 = rect.x0;
                    rect.y1 = rect.y0;
                    state = SELECTION_SELECTING;
                }
            } break;
            case SELECTION_SELECTING: {
                if (g::IsMouseDown(button)) {
                    requestUpdate(WINDOWFREEZE);
                    rect.x1 = g::GetMousePos().x - wp.x;
                    rect.y1 = g::GetMousePos().y - wp.y;
                } else {
                    clearUpdate(WINDOWFREEZE);
                    state = SELECTION_FINISHED;
                }
            } break;
            case SELECTION_FINISHED: {
                if (g::IsMouseDown(button)) {
                    state = SELECTION_DONE;
                }
            } break;
            }
        }

        switch (state)
        {
        case SELECTION_FINISHED:
        case SELECTION_SELECTING: {
            const unsigned int col = toRGBA32(color);
            g::GetWindowDrawList()->AddLine(ImVec2(rect.x0 + wp.x, rect.y0 + wp.y), ImVec2(rect.x1 + wp.x, rect.y0 + wp.y), col, 2.f);
            g::GetWindowDrawList()->AddLine(ImVec2(rect.x0 + wp.x, rect.y0 + wp.y), ImVec2(rect.x0 + wp.x, rect.y1 + wp.y), col, 2.f);
            g::GetWindowDrawList()->AddLine(ImVec2(rect.x1 + wp.x, rect.y1 + wp.y), ImVec2(rect.x1 + wp.x, rect.y0 + wp.y), col, 2.f);
            g::GetWindowDrawList()->AddLine(ImVec2(rect.x1 + wp.x, rect.y1 + wp.y), ImVec2(rect.x0 + wp.x, rect.y1 + wp.y), col, 2.f);
        } break;
        }
        if (state == SELECTION_FINISHED) {
            requestUpdate(uprq);
            return true;
        }
        return false;
    }

    void samplingPoints(int count,std::function<void(double &x, double &y, double &r, PIXEL &c)> samplingFunction) {
        for (int i = 0; i < count; ++i) {
            double x,y,r;
            PIXEL c;
            samplingFunction(x,y,r,c);
            ImVec2 w = g::GetWindowPos();
            w.x -= g::GetScrollX();
            w.y -= g::GetScrollY();
            w.x += g::GetCursorPosX();
            w.y += g::GetCursorPosY();
            g::GetForegroundDrawList()->AddCircle(ImVec2(x+w.x, y+w.y), r, toRGBA32(c));
        }
    }
}

ImageConfig viewConfig(DOUBLE zoom, bool normalize, bool alphaOnly, bool forceUpdate) {
    ImageConfig r;
    r.zoom = zoom;
    r.normalize = normalize;
    r.alphaOnly = alphaOnly;
    r.forceUpdate = forceUpdate;
    r.onlySingleChannel = Configuration::displayFeatureLayerChannel;
    return r;
}

std::string statusWindowTitleString() {
    std::string k = "Channels["+toString(DEFINED_CHANNELCOUNT)+"] " + std::string(FeatureLayerFormulaName) + " " + currentAlgorithmShortDescription + " Threads:" + toString(ArraySize(threadPool)) + " " + DIRECTIONFORMULASTRING;
    if (DEFINED_TRIANGULAR_DERIVATIVES == 1) {
        k += " tri";
    }
    else {
        k += " quad";
    }
    if (DEFINED_DERIVATIVE_ROTATION == 1) {
        k += ":rot";
    }
    if (baseConfig.iterative) {
        k += ":iter";
    }
    else {
        k += ":noiter";
    }
    if (DEFINED_ADAPTIVE_TYPE != ADAPTIVE_NONE) {
        k += " adaptive";
    }
    return k;
}

std::string substr(const std::string& str, int start = 0, int end = -1) {
    if (start < 0) start += str.size();
    if (start >= str.size())
        return "";
    if (end < 0) end = str.size();
    if (end > str.size()) end = str.size();
    if (start >= end)
        return "";
    return str.substr(start, end);
}

bool shouldReAdd(const FEATURELayer &layer, const Feature &f) {
    DOUBLE borderRatio = 0.975;
    if (f.x > layer.width * borderRatio)
        return true;
    if (f.y > layer.height * borderRatio)
        return true;
    if (f.x < layer.width * (1.0-borderRatio))
        return true;
    if (f.y < layer.height * (1.0-borderRatio))
        return true;
    return false;
}

bool isValid(const Descriptor& v) { // not needed later
    for (int i = 0; i < ArraySize(v.bits); ++i)
        if (ArraySize(v.bits[i]) != ActiveConfiguration::pairCount)
            return false;
    return true;
}

void createFeatureConfigs(const int layers) {
    baseConfig = createDefaultFeatureConfig(layers);
#if DEFINED_ADAPTIVE_TYPE != ADAPTIVE_NONE
    resampleConfig = createResampleFeatureConfig(layers);
#endif
#if DEFINED_FEATURETRACKING == 1
    trackConfig = createTrackFeatureConfig(layers);
#endif
}

void displayStuff() {
    DEBUGTRACE("displayStuffStart");
    static Maths::Rectangle baseFrameRect;
    static int baseFrameSelectionState = SELECTION_DONE;

    bool clearFeatures = false;
    if (shouldUpdate(RELOAD)) {
        DEBUG("Dataset Reload begin\n")
        clearUpdate(RELOAD);
        loadDataSet(collectedDirectories[Configuration::currentDatasetNr], Configuration::useFirstFrameValue ? Configuration::dataSetFirstFrame : -1,
                    Configuration::useLastFrameValue ? Configuration::dataSetLastFrame : -1, Configuration::dataSetStepSize, Configuration::dataSetReversed);
        loadTrackingSet(collectedDirectories[Configuration::currentDatasetNr] + trackingFolder);
        createFeatureConfigs(imagePyramidSize);
        clearFeatures = true;
        requestUpdate(FULLFEATURERECALC);
        DEBUG("Dataset Reload finished\n")
    }

    loading_mutex.lock();

    g::Begin("DataSet");
    {
        g::Separator();
        g::Text("%s", Configuration::dataSetFileName.c_str());
        g::Checkbox("reversed dataset", &Configuration::dataSetReversed);
        g::Checkbox("use####aaa1", &Configuration::useFirstFrameValue); if (Configuration::useFirstFrameValue) { g::SameLine(); g::SetNextItemWidth(100); g::InputInt("firstFrame", &Configuration::dataSetFirstFrame); }
        g::Checkbox("use####aaa2", &Configuration::useLastFrameValue); if (Configuration::useLastFrameValue) { g::SameLine(); g::SetNextItemWidth(100); g::InputInt("lastFrame", &Configuration::dataSetLastFrame); }
        g::SetNextItemWidth(100); g::InputInt("stepSize", &Configuration::dataSetStepSize);
        g::Separator();
        g::Checkbox("use####aaa3", &Configuration::useOwnWidth); if (Configuration::useOwnWidth) {
            g::SameLine(); g::SetNextItemWidth(100); g::InputInt("desiredWidth", &Configuration::desiredImageWidth);
        } else {
            g::SameLine(); g::SetNextItemWidth(100); g::InputDouble("desiredImageScale", &Configuration::desiredImageScale);
        }
        g::SetNextItemWidth(100); g::inputDOUBLE(0, "contrastRandom", &dataSetContrastRandom);
        g::Separator();
        g::button(RELOAD | TEXTURESUPDATE | FEATURERECALC, "load");
        g::Separator();
        g::SetNextItemWidth(1000000); g::ListBox("Datasets", &Configuration::currentDatasetNr, collectedDirectories, collectedDirectoriesCount);
        g::Separator();
        g::Text("frame count:%d", ArraySize(frames));
        for (int i = 0; i < frames.size(); ++i) {
            g::Text("%08d: %s", i, frames[i].name.c_str());
        }
    }
    g::End();

    DEBUGTRACE("datasetEnd");

    if (!frames.empty()) {

        const DOUBLE mx = g::GetMousePos().x;
        const DOUBLE my = g::GetMousePos().y;
        g::overlayPoint(mx, my);
        DOUBLE windowMouseX = mx;
        DOUBLE windowMouseY = my;

        const int baseFrameId = Configuration::currentFrame % ArraySize(frames);
        const int updatedFrameId = Configuration::successorFrame % ArraySize(frames);
        Image& baseFrame = frames[baseFrameId];
        Image& updatedFrame = frames[updatedFrameId];
        const std::string baseFrameName = substr(baseFrame.name, -30);
        const std::string updatedFrameName = substr(updatedFrame.name, -30);
        const int pyramidLayerId = Configuration::pyramidLayer % ArraySize(baseFrame.featureLayers);
        const DOUBLE zoomBaseImage = Configuration::baseZoom * baseFrame.featureLayers[0].width / baseFrame.featureLayers[pyramidLayerId].width;
        const DOUBLE zoomUpdatedImage = Configuration::updatedZoom * updatedFrame.featureLayers[0].width / updatedFrame.featureLayers[pyramidLayerId].width;
        ImageConfig configBaseImage = viewConfig(zoomBaseImage, Configuration::viewNormalized, Configuration::viewAlpha, shouldUpdate(TEXTURESUPDATE));
        ImageConfig configUpdatedImage = viewConfig(zoomUpdatedImage, Configuration::viewNormalized, Configuration::viewAlpha, shouldUpdate(TEXTURESUPDATE));
        configBaseImage.frameId = baseFrameId; configUpdatedImage.frameId = updatedFrameId;

        CHANNELS cBase, cUpdated;

        const int windowFlags = 0;

        if (!trackingFrames.empty()) {
            Configuration::trackingFrame %= trackingFrames.size();
            if (shouldUpdate(LOCATEPATCH)) {
                tracked = lookupImagePatch(trackingFrames[Configuration::trackingFrame], baseFrame, baseConfig);
                clearUpdate(LOCATEPATCH);
            }

            g::Begin("Tracking", nullptr, windowFlags);
            {
                ImageConfig configTrackingImage = viewConfig(Configuration::trackingZoom, false, false, shouldUpdate(TEXTURESUPDATE));
                DISPLAYLayer trackingLayer = toDisplay(trackingFrames[Configuration::trackingFrame].pixelLayers[0], configTrackingImage);
                if (Configuration::displayFeatureLayer)
                    trackingLayer = toDisplay(trackingFrames[Configuration::trackingFrame].featureLayers[0], configTrackingImage);
                else
                    trackingLayer = toDisplay(trackingFrames[Configuration::trackingFrame].pixelLayers[0], configTrackingImage);
                g::description(trackingLayer);
                g::SetNextItemWidth(100); g::inputDOUBLE(0, "zoom", &Configuration::trackingZoom, 0.25);
                g::SetNextItemWidth(100); g::inputInt(TEXTURESUPDATE, "frame", &Configuration::trackingFrame); g::SameLine(); g::Text("frameCount:%d", trackingFrames.size());
                g::Checkbox("showDescriptor", &Configuration::showTrackingDescriptor);
                Configuration::trackingFrame = CLAMP(Configuration::trackingFrame,0,trackingFrames.size()-1);

                ImVec2 cursor = g::GetCursorScreenPos();
                g::imageLayer(trackingLayer, configTrackingImage);
                ImVec2 cursor2 = g::GetCursorScreenPos();
                if (Configuration::showTrackingDescriptor) {
                    g::SetCursorScreenPos(cursor);
                    g::pairs( trackingLayer, 0, baseConfig, configTrackingImage); // todo: show propper descriptor shape here (not 0)
                    g::SetCursorScreenPos(cursor2);
                }
                g::button(LOCATEPATCH, "Locate");
                g::checkbox(LOCATEPATCH,"per frame locate", &Configuration::trackingPerFrame);
            }
            g::End();
        }

        DEBUGTRACE("trackingEnd");


        g::Begin((std::string("UPDATED:" + toString(updatedFrameId) + " " + updatedFrameName + "####ImageWindow1").c_str()), nullptr, windowFlags);
        {
            DISPLAYLayer updatedLayer;
            if (Configuration::displayFeatureLayer)
                updatedLayer = toDisplay(updatedFrame.featureLayers[pyramidLayerId], configUpdatedImage);
            else
                updatedLayer = toDisplay(updatedFrame.pixelLayers[pyramidLayerId], configUpdatedImage);
            g::description(updatedLayer);
            g::featureDescs(updatedFeatures);
            g::SetNextItemWidth(100); g::inputDOUBLE(0, "zoom", &Configuration::updatedZoom, 0.25);
#if DEFINED_DESCRIPTIVITYMAP == 1
            g::SameLine();
            g::checkbox(TEXTURESUPDATE, "descriptiveMaps##dsc2", &Configuration::viewDescriptivenes);
#endif
            g::BeginChild((ImGuiID)1000, ImVec2(0, 0), false, ImGuiWindowFlags_AlwaysVerticalScrollbar | ImGuiWindowFlags_AlwaysHorizontalScrollbar);
            ImVec2 cursor = g::GetCursorScreenPos();
            cUpdated = sample_featurepixel_direct_clamped(updatedFrame.featureLayers[pyramidLayerId], (mx - cursor.x) / zoomUpdatedImage, (my - cursor.y) / zoomUpdatedImage);
            g::SetCursorScreenPos(cursor);
            if (g::IsWindowHovered()) {
                windowMouseX = mx - cursor.x;
                windowMouseY = my - cursor.y;
            }
            g::imageLayer(updatedLayer, configUpdatedImage);
#if DEFINED_DESCRIPTIVITYMAP == 1
            const DESCRIPTIVITYLayer& descriptivityLayer = updatedFrame.descriptivityLayer;
            if (!descriptivityLayer.pixels.empty() && Configuration::viewDescriptivenes) {
                const DOUBLE zoomDescriptivity = Configuration::updatedZoom * updatedFrame.featureLayers[0].width / updatedFrame.descriptivityLayer.width;
                ImageConfig configDescriptivityLayer = viewConfig(zoomDescriptivity, 1.0, false, shouldUpdate(TEXTURESUPDATE));
                cursor.x -= zoomUpdatedImage * baseFrame.featureLayers[0].width / baseFrame.descriptivityLayer.width * 0.5;
                cursor.y -= zoomUpdatedImage * baseFrame.featureLayers[0].width / baseFrame.descriptivityLayer.width * 0.5;
                g::SetCursorScreenPos(cursor);
                g::imageLayer(toDisplay(descriptivityLayer, configDescriptivityLayer), configDescriptivityLayer);
            }
#endif
            configUpdatedImage.zoom = Configuration::updatedZoom; // sorry..
            if (Configuration::viewFeatures)
                g::features(updatedFeatures, configUpdatedImage);
            //g::descriptorLevel(0,sourceImage, sourceFeatures, sourceDescriptors, baseDesc);
            g::EndChild();
        }
        g::End();

        DEBUGTRACE("UpdateEnd");

        g::Begin((std::string("BASE:" + toString(baseFrameId) + " " + baseFrameName + "####ImageWindow2").c_str()), nullptr, windowFlags);
        {
            DISPLAYLayer baseLayer;
            if (Configuration::displayFeatureLayer)
                baseLayer = toDisplay(baseFrame.featureLayers[pyramidLayerId], configBaseImage);
            else
                baseLayer = toDisplay(baseFrame.pixelLayers[pyramidLayerId], configBaseImage);
            g::description(baseLayer);
            g::Text("Tracked[x%f,y%f,s%f,r%f]\n", tracked.position.x, tracked.position.y, tracked.position.scale, tracked.position.rotation);
            g::featureDescs(baseFeatures);
            g::SetNextItemWidth(100); g::inputDOUBLE(0, "zoom", &Configuration::baseZoom, 0.25);
#if DEFINED_DESCRIPTIVITYMAP == 1
            g::SameLine();
            g::checkbox(TEXTURESUPDATE, "descriptiveMaps##dsc3", &Configuration::viewDescriptivenes);
#endif
            g::SameLine();
            g::SetNextItemWidth(100); g::InputDouble("featureAlpha", &Configuration::imageFeatureAlpha, 0.1);
            if (g::Button("New Binks")) {
                binks.clear();
                metaBinks.clear();
                ArrayAppend(binks, getClusterBinks(baseFeatures, Configuration::binkForceThreshold, configBaseImage));
            }
            g::SameLine();g::SetNextItemWidth(100.0);g::inputDOUBLE(0, "attenuation", &Configuration::binkDistanceAttenuationThreshold);
            g::SameLine();g::SetNextItemWidth(100.0);g::inputDOUBLE(0, "threshold", &Configuration::binkDistanceAttenuationFactor);
            g::SameLine();g::SetNextItemWidth(100.0);g::inputDOUBLE(0, "forceThreshold", &Configuration::binkForceThreshold);

            if (!binks.empty()) {
                g::button(BINKSNEEDMOD, "Clear Binks");
                g::Checkbox("Place MetaBink", &Configuration::placeMink); g::SameLine(); //g::SetNextItemWidth(100); g::InputInt("Mink Type", &minkType);
                //g::Checkbox("Select Minks", &selectMinks);
                g::SameLine(); if (g::Button("Clear Current Frame MetaBinks")) {
                    clearFrameBinks(metaBinks, baseFrameId);
                }
                for (int i = 0; i < metaBinks.size(); ++i) {
                    g::SameLine(); g::Text("[%f,%f]", metaBinks[i].px, metaBinks[i].py);
                }
            }
            g::BeginChild((ImGuiID)1010, ImVec2(0, 0), false, ImGuiWindowFlags_AlwaysVerticalScrollbar | ImGuiWindowFlags_AlwaysHorizontalScrollbar);
            ImVec2 cursor = g::GetCursorScreenPos();
            cBase = sample_featurepixel_direct_clamped(baseFrame.featureLayers[pyramidLayerId], (mx - cursor.x) / zoomBaseImage, (my - cursor.y) / zoomBaseImage);
            g::SetCursorScreenPos(cursor);
            if (g::IsWindowHovered()) {
                windowMouseX = mx - cursor.x;
                windowMouseY = my - cursor.y;
                Descriptor descriptor = sampleDescriptor(baseFrame, { windowMouseX / zoomBaseImage, windowMouseY / zoomBaseImage}, baseConfig);
                windowMouseStrength = descriptivityOfDescriptor(descriptor, baseConfig);
                if (g::IsMouseClicked(0) && (!Configuration::placeMink)) {
#define isInselection(__a) ((__a) != SELECTION_DONE)
                    if (!isInselection(baseFrameSelectionState)) {
                        Feature f(windowMouseX, windowMouseY);
                        ArrayAdd(baseFeatures, f);
                        requestUpdate(JUSTLASTFEATURE | DESCRIBEFEATURES | PLACEFEATURES | EXPLICITUPDATEORIGINAL | REMOVEUNSUFFICIENTFEATURES);
                    }
                }
            }
            g::imageLayer(baseLayer, configBaseImage);
#if DEFINED_DESCRIPTIVITYMAP == 1
            const DESCRIPTIVITYLayer& descriptivityLayer = baseFrame.descriptivityLayer;
            if (!descriptivityLayer.pixels.empty() && Configuration::viewDescriptivenes) {
                const DOUBLE zoomDescriptivity = Configuration::baseZoom * baseFrame.featureLayers[0].width / baseFrame.descriptivityLayer.width;
                ImageConfig configDescriptivityLayer = viewConfig(zoomDescriptivity, 1.0, false, shouldUpdate(TEXTURESUPDATE));
                cursor.x -= zoomBaseImage * baseFrame.featureLayers[0].width / baseFrame.descriptivityLayer.width * 0.5;
                cursor.y -= zoomBaseImage * baseFrame.featureLayers[0].width / baseFrame.descriptivityLayer.width * 0.5;
                g::SetCursorScreenPos(cursor);
                g::imageLayer(toDisplay(descriptivityLayer, configDescriptivityLayer), configDescriptivityLayer);
            }
#endif
            configBaseImage.zoom = Configuration::baseZoom; // sorry..
            if (Configuration::viewFeatures)
                g::features(baseFeatures, configBaseImage);
            if (!trackingFrames.empty())
                g::imageFeature(tracked.position, trackingFrames[Configuration::trackingFrame], configBaseImage );
            if (Configuration::viewMesh)
                g::featureMesh(baseFeatures, configBaseImage);
            if (Configuration::viewClassifications)
                g::classifications(updatedFeatureClasses, configBaseImage);
            g::SetCursorScreenPos(cursor);
            if (g::IsWindowHovered()) {
                ImVec2 cursor = g::GetCursorScreenPos();
                const double fx = (mx - cursor.x) / zoomBaseImage;
                const double fy = (my - cursor.y) / zoomBaseImage;
                Feature f(fx, fy);
                if (Configuration::selectMinks) {
                    if (g::IsMouseClicked(1)) {
                        getNearestBink(metaBinks, f.x,f.y, configBaseImage).selected ^= true;
                    }
                }
                if (Configuration::placeMink) {
                    g::mink(f, getLastBinkIdInFrame(metaBinks, baseFrameId) + 1, "", 24,configBaseImage);
                    if (g::IsMouseClicked(0)) {
                        metaBinks.push_back(getMetaBink(binks, f.x, f.y, getLastBinkIdInFrame(metaBinks, baseFrameId) + 1, baseFrameId));
                        Configuration::placeMink = false;
                    }
                }
            }

            for (int i = 0; i < binks.size(); ++i) {
                updateCenterBink(baseFrame, baseFeatures, binks[i], configBaseImage);
            }

            for (int i = 0; i < metaBinks.size(); ++i) {
                updateMetaBink(baseFrame, binks, metaBinks[i], configBaseImage);
            }

            //for (int i = 0; i <= getLastBinkId(binks); ++i) {
            //    Bink floating = getFloatingBink(baseFrame, binks, baseFeatures, i, configBaseImage);
            //    g::bink(floating, configBaseImage);
            //}

            for (int i = 0; i < binks.size(); ++i) {
                g::bink(binks[i], configBaseImage);
            }

            for (int i = 0; i < metaBinks.size(); ++i) {
                g::bink(metaBinks[i], configBaseImage);
            }

            g::doRectangularSelection(0, baseFrameRect, baseFrameSelectionState);
            g::EndChild();
        }
        g::End();

        DEBUGTRACE("BaseEnd");

        clearUpdate(TEXTURESUPDATE);

        g::Begin((statusWindowTitleString() + "####StatusWindow").c_str());
        g::Separator();
        if (g::inputInt(REDOCLASSIFICATION | TEXTURESUPDATE | FEATURERECALC | LASTFRAMEUPDATE | SAVECURRENTFRAME, "Next Frame", &Configuration::successorFrame)) {
            if (Configuration::trackingPerFrame)
                requestUpdate(LOCATEPATCH);
        }
#if DEFINED_ADAPTIVE_TYPE != ADAPTIVE_NONE
        g::SetNextItemWidth(100); g::InputDouble("adaptiveThreshold", &descriptorThresholdRatioForAdaptiveUpdate, 0.1);
        g::Text("adaptiveSampledDescriptors: %d->%s of %d", adaptiveSampledDescriptors, resampleIterative ? "iter" : "noiter", sampledDescriptors);
        g::SameLine();
        g::Checkbox("resampling", &useResampling);
#endif
        g::Separator();
        if (g::button(TEXTURESUPDATE | FULLFEATURERECALC | BINKSNEEDMOD, "Reset")) {
            Configuration::currentFrame = 0;
        }
        g::SameLine();
        if (g::button(TEXTURESUPDATE | FULLFEATURERECALC | BINKSNEEDMOD, "Clear Features") || clearFeatures) {
            updatedFeatureClasses.clear();
            baseDescriptors.clear();
            originalDescriptors.clear();
            baseFeatures.clear();
            updatedDescriptors.clear();
            updatedFeatures.clear();
        }
        g::SameLine();
        g::button(DOPRUNING, "prune"); g::SameLine(); g::SetNextItemWidth(100);g::InputDouble("pruneThreshold", &Configuration::pruneThreshold);


        if (g::button(TEXTURESUPDATE | FULLFEATURERECALC | BINKSNEEDMOD, "Pop Feature")) {
            baseDescriptors.pop_back();
            originalDescriptors.clear();
            baseFeatures.pop_back();
            updatedDescriptors.pop_back();
            updatedFeatures.pop_back();
        }
        g::Separator();
#if DEFINED_DESCRIPTIVITYMAP == 1
        if (g::button(TEXTURESUPDATE, "Show Descriptivity Map")) {
            buildDescriptivityMap(baseFrame, Configuration::descriptivityLayer);
        }
        g::SameLine(); g::SetNextItemWidth(100); g::InputInt("Layer###DLayer", &Configuration::descriptivityLayer);
#endif
        g::SameLine(); g::SetNextItemWidth(100); g::InputInt("Count", &Configuration::requestedBestFeatureCount);
        g::SameLine(); g::SetNextItemWidth(100); g::InputInt("Layer", &Configuration::descriptivityLayer);
        if (g::button(TEXTURESUPDATE | DESCRIBEFEATURES | BINKSNEEDMOD, "CircleFeatures")) {
            baseDescriptors.clear();
            originalDescriptors.clear();
            baseFeatures.clear();
            updatedFeatures.clear();
            updatedDescriptors.clear();
            baseFeatures = placeCircularFeatures(baseFrame, Configuration::requestedCircularFeatureRadis , Configuration::requestedCircularFeatureSectors);
        }
        g::SameLine(); g::SetNextItemWidth(100); g::InputInt("R", &Configuration::requestedCircularFeatureSectors); g::SameLine(); g::SetNextItemWidth(100); g::InputInt("D", &Configuration::requestedCircularFeatureRadis);
        g::PushStyleVar(ImGuiStyleVar_Alpha, isInselection(baseFrameSelectionState) ? (baseFrameSelectionState == SELECTION_FINISHED ? 0.7 : 0.5) : 1.0);
        g::PushStyleColor(ImGuiCol_Button, (baseFrameSelectionState == SELECTION_FINISHED) ? 0xffff00ff : g::GetColorU32(ImGuiCol_Button));
        if (g::Button((baseFrameSelectionState == SELECTION_FINISHED) ? "Place##RandomFeaturesInRectangle" : "RandomFeaturesInRectangle##RandomFeaturesInRectangle")) {
            switch (baseFrameSelectionState) {
            case SELECTION_FINISHED: {
                requestUpdate(RECTANGULARSELECTION);
            } break;
            case SELECTION_AWAIT: {
                baseFrameSelectionState = SELECTION_DONE;
            } break;
            default: {
                baseFrameSelectionState = SELECTION_AWAIT;
            } break;
            }
        }
        g::PopStyleColor();
        g::PopStyleVar();
        //g::SameLine();
        //if (g::button(TEXTURESUPDATE | DESCRIBEFEATURES, "RandomFeatures")) {
        //    baseDescriptors.clear();
        //    originalDescriptors.clear();
        //    baseFeatures.clear();
        //    updatedFeatures.clear();
        //    updatedDescriptors.clear();
        //    baseFeatures = placeRandomFeatures(baseFrame, requestedRandomFeatureCount);
        //}
        g::SameLine(); g::SetNextItemWidth(100); g::InputInt("Count##RCount", &Configuration::requestedRandomFeatureCount);
        if (shouldUpdate(RECTANGULARSELECTION)) {
            clearUpdate(RECTANGULARSELECTION);
            baseFrameSelectionState = SELECTION_DONE;
            Array(Feature) features = placeRandomFeatures(baseFrame, Configuration::requestedRandomFeatureCount, baseFrameRect);
            ArrayAppend(baseFeatures, features);
            requestUpdate(TEXTURESUPDATE | DESCRIBEFEATURES | PLACEFEATURES | REMOVEUNSUFFICIENTFEATURES | BINKSNEEDMOD);
        }
        g::button(SAVECURRENTFRAME, "SaveImage");
        std::string labelContinuosRender = (Configuration::animOn  ? "ContinuosRender:" + toString(baseImageCountForSaving) : "ContinuosRender") + "##ContinuosRender";
        if (g::Checkbox(labelContinuosRender.c_str(), &Configuration::animOn)) {
            baseImageCountForSaving = 0;
        }
        if (Configuration::animOn) {
            g::Checkbox("AddNewFeatures", &Configuration::addNewContinuosFeatures);
        }

        g::button(TEXTURESUPDATE | DESCRIBEFEATURES, "describeFeatures");
        g::SameLine();
        if (g::button(TEXTURESUPDATE | PLACEFEATURES, "place"));
#if DEFINED_ADAPTIVE_TYPE == 1
        if (useResampling) {
            g::SameLine();
            if (g::button(TEXTURESUPDATE | PLACEFEATURESRESAMPLE, "placeResample"));
        }
#endif
        g::Separator();
        
        if (Configuration::iterative)
            Configuration::sampleDescriptors = true;
        else
            if (g::Checkbox("SampleDescriptors", &Configuration::sampleDescriptors));
        g::Checkbox("Pyramidal", &Configuration::pyramidal);
        if (Configuration::pyramidal)
        {
            g::SetNextItemWidth(100);
            g::InputInt("lastLevel", &baseConfig.lastLevel); g::SameLine();
            g::Text("deepest:%d array:%d", baseConfig.deepestLevel, baseConfig.descriptorPatterns[0].levels.size());
        } else {
            g::SetNextItemWidth(100);
            g::InputInt("level", &baseConfig.sampleLevel); g::SameLine();
            g::Text("deepest:%d array:%d", baseConfig.deepestLevel, baseConfig.descriptorPatterns[0].levels.size());
        }
#if DEFINED_ALSO_RESBIFT == 1
        g::checkbox(PLACEFEATURES, "RESBIFT Algorithms", &Configuration::RESBIFT);
        if (Configuration::RESBIFT) {
            g::SameLine();
            g::SetNextItemWidth(100.0);
            g::inputDOUBLE(PLACEFEATURES, "resbift intensity", &Configuration::resbiftFactor);
        }
#endif
        g::checkbox(PLACEFEATURES, "Circular Landing", &Configuration::circularLanding);
        if (Configuration::circularLanding) {
            g::SameLine();
            g::SetNextItemWidth(100.0);
            g::inputDOUBLE(PLACEFEATURES, "landing factor", &Configuration::circularLandingFactor);
        }

        if (g::ListBox((std::string("Algorithms:") + toString(Configuration::algorithmCount) + "##Algorithms").c_str() , &Configuration::algorithm, Configuration::algorithmNames, Configuration::algorithmCount)) {
            requestUpdate(PLACEFEATURES);
        }
        if (g::Button("Rebuild Descriptor Shape")) {
            createFeatureConfigs(imagePyramidSize);
        }
        g::Separator();
#if DEFINED_DESCRIPTIVITYMAP == 1
        g::checkbox(TEXTURESUPDATE, "viewDescriptiveMaps", &Configuration::viewDescriptivenes);
#endif
#if DEFINED_CHANNELCOUNT == 1
        g::checkbox(TEXTURESUPDATE, "viewFeatureIntensities", &Configuration::displayFeatureLayer);
#else
        g::checkbox(TEXTURESUPDATE, "viewFeatureChannel", &displayFeatureLayer);
        g::inputInt(TEXTURESUPDATE, "featureChannelToView", &displayFeatureLayerChannel);
        displayFeatureLayerChannel = CLAMP(displayFeatureLayerChannel,-1,channelElementCount()-1);
#endif
        g::checkbox(TEXTURESUPDATE, "viewNormalized", &Configuration::viewNormalized);
#define SELECT(__b, __v1, __v2) (__b ? __v1 : __v2)
        //if (ArraySize(SELECT(displayFeatureLayer, channelNamesFeatureSpace, channelNamesPixelSpace)) > 1)
        //    g::checkbox(TEXTURESUPDATE, "viewAlpha", &viewAlpha);
        int chMode = Configuration::displayFeatureLayer ? 1 : 0;
        //g::multiListBox(displayFeatureLayer ? 0 : 1, "Channel", &channelFormat, { channelNamesFeatureSpace, channelNamesPixelSpace });
        g::SetNextItemWidth(100); if (g::inputInt(TEXTURESUPDATE, "viewPyramidLayer", &Configuration::pyramidLayer))
            Configuration::forceLayer = true;
        if (!Configuration::pyramidal)
            g::Checkbox("useLayer", &Configuration::forceLayer);
        g::Separator();
        g::Checkbox("manhattenDistance", &Configuration::manhattenDistance);
        g::Checkbox("viewMesh", &Configuration::viewMesh);
        g::Checkbox("useSpice", &Configuration::useSpice);
        if (Configuration::useSpice) {
            g::SetNextItemWidth(100);
            g::inputDOUBLE(0, "spiceIterPercent", &Configuration::spicePercentage, 1);
            g::SetNextItemWidth(100);
            g::inputDOUBLE(0, "spiceIntensity", &Configuration::spiceIntensity, 0.1);
            g::SetNextItemWidth(100);
            g::inputDOUBLE(0, "perLevelFeatureReduceFactor", &baseConfig.refinementStepPerLevelScaleFactor, 0.1);
        }
        g::Checkbox("viewFeatures", &Configuration::viewFeatures);
        g::Checkbox("viewClassification", &Configuration::viewClassifications);
        if (g::InputInt("coreCount", &coreCount)) {
            rebuildThreadPool();
        }
        g::Separator();
        g::Text("MousePos %.03f:%.03f", mx, my);
        g::Text("WindowMousePos %.03f:%.03f", windowMouseX, windowMouseY);
        g::Text("DescriptivityBase: %.03f", windowMouseStrength);
        g::Text("MouseIntensityBase: %.03f", intensity(cBase));
        g::Text("MouseIntensityUpdated: %.03f", intensity(cUpdated));
        g::Text("DescriptorLevels: %d", baseConfig.deepestLevel);
        g::Text("DescriptorPairs: %d", Configuration::pairCount);
        g::Text("RefinementSteps: %d", Configuration::refinementSteps);
        g::Text("SampledFeatures: %d", featuresInvolved);
        g::Text("sampledDescriptors: %d->%s", sampledDescriptors, Configuration::iterative ? "iter" : "noiter");
        if (g::Button("ClearTextureCache")) {
            Graphical::clearCache();
        }
        g::Separator();
        g::Text("-- Configuration");
        g::Separator();
        g::Text("Feature Image Channel Count:%d", DEFINED_CHANNELCOUNT);
        g::Text("Iterative (baseConfig):%s", boolToString(baseConfig.iterative).c_str());
        g::Text("Adaptive Mode:%d", DEFINED_ADAPTIVE_TYPE);
        g::Text("Brightness Invariance:%s", boolToString(DEFINED_ADAPTIVEBRIGHTNESSINVARIANCE).c_str());
        g::Text("3D Patches:%s", boolToString(DEFINED_NORMAL3DPATCHES).c_str());
        g::Text("Lens Patches:%s", boolToString(DEFINED_ROUNDEDPATCHES).c_str());
        g::Text("Normalized Images:%s", boolToString(DEFINED_NORMALIZEDIMAGES).c_str());
        g::Text("Denoised Images:%s", boolToString(DEFINED_REMOVEDUST).c_str());
        g::Text("Rotation Invariance:%s", boolToString(DEFINED_WITHROTATIONINVARIANCE).c_str());
        g::Text("Only 45 Degree Invariance:%s", boolToString(DEFINED_WITHROTATIONINVARIANCE45DEGREES).c_str());
        g::Text("Sector Based 360 Degree Invariance:%s", boolToString(DEFINED_FULLROTATIONINVARIANCEBYSECTORS).c_str());
        g::Text("Just Pure Binary Descriptors:%s", boolToString(DEFINED_PUREBINARITY).c_str());
        g::Text("Triangular Derivatives:%s", boolToString(DEFINED_TRIANGULAR_DERIVATIVES).c_str());
        g::Text("Descriptor Shape with Rotation Grip:%s", boolToString(DEFINED_DESCRIPTORSHAPETOTRIANGLE).c_str());
        g::Text("More Rotational Descriptor Shape:%s", boolToString(DEFINED_ADVANCEDPATTERNTYPE).c_str());
        g::Text("Randomly Rotated Derivative Lookups:%s", boolToString(DEFINED_DERIVATIVE_ROTATION).c_str());
        g::Text("Random Movement for Unusable Pixels:%s", boolToString(DEFINED_HASRANDOMMOVEMENT).c_str());
        g::Text("Thread Pool Size:%d", threadPool.size());
        g::Text("Border Clamping:%s", boolToString(DEFINED_BORDERCLAMPING).c_str());
        g::Text("Border Clamp Mirror:%s", boolToString(DEFINED_BORDERCLAMPMIRROR).c_str());
        g::Text("Border Clamp Wrap:%s", boolToString(DEFINED_BORDERCLAMPWRAP).c_str());
        g::Text("Max Refinement Steps per Level:%d", DEFINED_MAXMIPLEVELLOOPSTEPS);
        g::Text("Bit Formula:%s", DIRECTIONFORMULASTRING.c_str());
        g::Separator();
        g::Text("-- baseConfig:");
        g::Text("Iterative:%s", boolToString(baseConfig.iterative).c_str());
        g::Text("Deepest Mip Map:%d", baseConfig.deepestLevel);
        g::Text("Last Looked Up Mip Map:%d", baseConfig.lastLevel);
        g::Text("Scale Invariance:%s", boolToString(baseConfig.scaleInvariance).c_str());
        g::Text("Refinement Steps per Level:%f", baseConfig.refinementSteps);
        g::Text("Step Factor:%f", baseConfig.stepFactor);
        g::Text("Refinement Step Count Mul per Level:%f", baseConfig.refinementStepPerLevelScaleFactor);
        g::Text("Descriptor Size Coarse:%f", baseConfig.descriptorSizeCoarse);
        g::Text("Descriptor Size Fine:%f", baseConfig.descriptorSizeFine);
        g::Text("Source Image To Dest Image Scale Ratio Hint:%f", baseConfig.sourceToDestScaleRatio);
        g::Text("'Featur.scale' Calculation Steps:%d", baseConfig.scaleLevelChecksCount);
        g::Text("Number of different Descriptor Shapes:%d", baseConfig.descriptorPatterns.size());
        g::Text("Descriptor Pattern Count:%d", baseConfig.descriptorPatterns.size());
        g::Text("Descriptor Pattern0 Depth:%d", baseConfig.descriptorPatterns[0].levels.size());
        g::Text("Descriptor Pattern0 Level0 Pairs:%d", baseConfig.descriptorPatterns[0].levels[0].pairs.size());
        g::Separator();
        g::End();

        if (Configuration::animOn) {
            static int ticki = 0;
            ticki++;
            if (!(ticki & 15)) {
                requestUpdate(TEXTURESUPDATE | FEATURERECALC | LASTFRAMEUPDATE | SAVECURRENTFRAME | READDFINISHEDFEATURES);
                if (Configuration::trackingPerFrame)
                    requestUpdate(LOCATEPATCH);
                Configuration::successorFrame++;
            }
        }

        if (shouldUpdate(DOPRUNING)) {
            pruneByGradientForce(baseFeatures, Configuration::pruneThreshold);
            clearUpdate(DOPRUNING);
        }

        if (shouldUpdate(SAVECURRENTFRAME)) {
            clearUpdate(SAVECURRENTFRAME);
            char buffer[1000];
            sprintf(buffer, (saveDirectory + "yo%04d.png").c_str(), baseImageCountForSaving);
            savePixelLayerWithFeatures(frames[baseFrameId].pixelLayers[0], baseFeatures, buffer, 0);
            baseImageCountForSaving++;
        }

        if (shouldUpdate(DESCRIBEFEATURES)) {
            DEBUG("Describe Features begin\n");
            ArrayResize(baseDescriptors, ArraySize(baseFeatures));
            ArrayResize(updatedFeatures, ArraySize(baseFeatures));
            ArrayResize(updatedDescriptors, ArraySize(baseDescriptors));
            ArrayResize(originalDescriptors,  ArraySize(baseDescriptors));
            for (int i = 0; i < ArraySize(baseDescriptors); ++i) threadPoolThread([&, i]() {
                initThread();
                baseDescriptors[i] = sampleDescriptor(baseFrame, baseFeatures[i], baseConfig);
                updatedDescriptors[i] = baseDescriptors[i];
#if DEFINED_ADAPTIVE_TYPE != ADAPTIVE_NONE
                originalDescriptors[i] = sampleDescriptor(baseFrame, baseFeatures[i], resampleConfig);
#endif
                });
            threadPoolFinish();
            clearUpdate(DESCRIBEFEATURES);
            DEBUG("Describe Features end\n");
        }

        if (shouldUpdate(PLACEFEATURES)) {
            DEBUG("Place Features begin\n");
            const int baseFrameId = Configuration::currentFrame % ArraySize(frames);
            const Image& baseFrame = frames[baseFrameId];
            ArrayResize(baseDescriptors, ArraySize(baseFeatures));
            ArrayResize(updatedFeatures, ArraySize(baseFeatures));
            ArrayResize(updatedDescriptors, ArraySize(baseDescriptors));
            ArrayResize(originalDescriptors,  ArraySize(baseDescriptors));
            int firstFeatureToUpdate = 0;
            int lastFeatureToUpdate = ArraySize(baseDescriptors) - 1;
            if (shouldUpdate(JUSTLASTFEATURE)) {
                firstFeatureToUpdate = ArraySize(baseDescriptors) - 1;
                lastFeatureToUpdate = firstFeatureToUpdate;
            }
            for (int i = firstFeatureToUpdate; i <= lastFeatureToUpdate; ++i) threadPoolThread([&, i]() {
                    initThread();
                    const Feature f1 = Feature(baseFrame.featureLayers[0].width * 0.5,baseFrame.featureLayers[0].height * 0.5, baseFeatures[i]);
                    const Feature f2 = baseFeatures[i];
                    Feature feature = lookupFeature(baseFrame, baseDescriptors[i], baseConfig.iterative ? f2 : f1, baseConfig);
                    baseFeatures[i] = feature;
                    updatedFeatures[i] = feature;
                    if (shouldUpdate(REMOVEUNSUFFICIENTFEATURES)) {
                        const double dx = feature.x - f2.x;
                        const double dy = feature.y - f2.y;
                        const double d = sqrt(dx * dx + dy * dy);
#if DEFINED_GUI_REMOVE_INVALID_FEATURES_AFTER_ADDINGS == 1
                        const double tooFar = baseFrame.featureLayers[0].width * REMOVE_INVALID_FEATURES_AFTER_ADDINGS_AREA;
                        if (d > tooFar) {
                            MARKREMOVE(baseFeatures[i]);
                        }
#endif
                    }
                });
            threadPoolFinish();
            clearUpdate(PLACEFEATURES | JUSTLASTFEATURE | REMOVEUNSUFFICIENTFEATURES);
            DEBUG("Place Features end\n");
        }

        if (shouldUpdate(FULLFEATURERECALC | FEATURERECALC)) {
#define ADAPTIVE_RESAMPLEORIGINALS (DEFINED_ADAPTIVE_TYPE == ADAPTIVE_ITERATIVE)

            DEBUG("FULLFEATURERECALC | FEATURERECALC begin\n")
            const int baseFrameId = Configuration::currentFrame % ArraySize(frames);
            const int updatedFrameId = Configuration::successorFrame % ArraySize(frames);
            Image& baseFrame = frames[baseFrameId];
            Image& updatedFrame = frames[updatedFrameId];
            ArrayResize(baseDescriptors, ArraySize(baseFeatures));
            ArrayResize(originalDescriptors, ArraySize(baseFeatures));
            ArrayResize(updatedFeatures, ArraySize(baseFeatures));
            ArrayResize(updatedDescriptors, ArraySize(baseDescriptors));
            int firstFeatureToUpdate = 0;
            int lastFeatureToUpdate = ArraySize(baseDescriptors) - 1;
            if (shouldUpdate(JUSTLASTFEATURE)) {
                firstFeatureToUpdate = ArraySize(baseDescriptors) - 1;
                lastFeatureToUpdate = firstFeatureToUpdate;
                requestUpdate(PLACEFEATURES); // just last feature stays
            }
            for (int i = firstFeatureToUpdate; i <= lastFeatureToUpdate; ++i) threadPoolThread([&, i]() {
                initThread();
                bool updateOriginal = false;
#if DEFINED_ADAPTIVE_TYPE != ADAPTIVE_NONE
                updateOriginal = descriptorRatio(baseDescriptors[i], originalDescriptors[i]) < descriptorThresholdRatioForAdaptiveUpdate ? true : false;
                updateOriginal &= useResampling;
                updateOriginal |= shouldUpdate(EXPLICITUPDATEORIGINAL);
#endif
                const Feature f1 = Feature(baseFrame.featureLayers[0].width * 0.5,baseFrame.featureLayers[0].height * 0.5, baseFeatures[i]);
                const Feature f2 = baseFeatures[i];
                if ((!shouldUpdate(FULLFEATURERECALC)) && updateOriginal) {
                    const Feature f = lookupFeature(baseFrame, originalDescriptors[i], resampleConfig.iterative ? f2 : f1, resampleConfig);
                    if (ADAPTIVE_RESAMPLEORIGINALS || shouldUpdate(EXPLICITUPDATEORIGINAL)) {
                        originalDescriptors[i] = sampleDescriptor(baseFrame, f, resampleConfig);
                    }
                    baseDescriptors[i] = sampleDescriptor(baseFrame, f, baseConfig);
                    updatedDescriptors[i] = baseDescriptors[i]; // currently not used further than here
                    adaptiveSampledDescriptors++;

                    Feature feature = lookupFeature(baseFrame, originalDescriptors[i], resampleConfig.iterative ? f : f1, resampleConfig);
                    updatedFeatures[i] = feature;
                }
                else
                {
                    if (Configuration::sampleDescriptors || shouldUpdate(FULLFEATURERECALC)) {
                        baseDescriptors[i] = sampleDescriptor(baseFrame, baseFeatures[i], baseConfig);
#if ADAPTIVE_RESAMPLEORIGINALS
                            if (shouldUpdate(FULLFEATURERECALC)) originalDescriptors[i] = sampleDescriptor(baseFrame,
                                                                                                           {baseFeatures[i].x,
                                                                                                            baseFeatures[i].y},
                                                                                                           resampleConfig);
#endif
                    }
                    Feature feature = lookupFeature(updatedFrame, baseDescriptors[i], baseConfig.iterative ? f2 : f1, baseConfig);
                    updatedDescriptors[i] = baseDescriptors[i]; // currently not used further than here
                    updatedFeatures[i] = feature;
                }
            });
            threadPoolFinish();
            clearUpdate(FULLFEATURERECALC | FEATURERECALC | EXPLICITUPDATEORIGINAL);
            DEBUG("FULLFEATURERECALC | FEATURERECALC end\n")
        }

#if DEFINED_ONLYNOPROBLEMATIC == 0
        if (shouldUpdate(REDOCLASSIFICATION)) {
            updatedFeatureClasses = classify(updatedFeatures, baseFeatures, buildCorrespondence(updatedFeatures, baseFeatures));
            clearUpdate(REDOCLASSIFICATION);
        }
#endif

        if (shouldUpdate(LASTFRAMEUPDATE)) {
            DEBUG("LASTFRAMEUPDATE begin\n")
            const int baseFrameId = Configuration::currentFrame % ArraySize(frames);
            const int updatedFrameId = Configuration::successorFrame % ArraySize(frames);
            baseFeatures = updatedFeatures;
            baseDescriptors = updatedDescriptors;
            updatedDescriptors.clear();
            updatedFeatures.clear();
            if (shouldUpdate(READDFINISHEDFEATURES)) {
                clearUpdate(READDFINISHEDFEATURES);
                for (int i = 0; i < ArraySize(baseFeatures); ++i) {
                    if (Configuration::addNewContinuosFeatures) {
                        if (shouldReAdd(baseFrame.featureLayers[0], baseFeatures[i])) {
                            Feature k;
                            DOUBLE lastStrength = 0;
                            // check for possible good strength of new added features
                            for (int i = 0; i < 5; ++i) {
                                Feature b(((random(1.0) - 0.5) * 0.2 + 0.5) * baseFrame.featureLayers[0].width, ((random(1.0) - 0.5) * 0.2 + 0.5) * baseFrame.featureLayers[0].height);
                                Descriptor descriptor = sampleDescriptor(baseFrame, b, baseConfig);
                                const DOUBLE strength = descriptivityOfDescriptor(descriptor, baseConfig);
                                if (i == 0 || strength > lastStrength) {
                                    lastStrength = strength;
                                    k = b;
                                }
                            }
                            baseFeatures[i] = k;
                        }
                    }
                }
            }

#if DEFINED_ITERATIVE_BINKS == 1
            ImageConfig config1 = configBaseImage;
            ImageConfig config2 = configUpdatedImage;
            config1.frameId = baseFrameId;
            config2.frameId = successorFrame;

            Array(Bink) recatpuredBinks = getRecapturedBinksOfFrame(baseFrame, binks, baseFeatures, config1);
            replaceFrameId(recatpuredBinks, config1.frameId, config2.frameId);
            replaceAddBinks(binks, recatpuredBinks);
#endif
            clearUpdate(LASTFRAMEUPDATE);
            Configuration::currentFrame = Configuration::successorFrame;
            DEBUG("LASTFRAMEUPDATE end\n")
        }

        // remove features
        for (int i = baseFeatures.size() - 1; i >= 0; --i) {
            if (SHOULDREMOVE(baseFeatures[i])) {
                ArrayRemove(baseFeatures, i);
                ArrayRemove(updatedFeatures, i);
                ArrayRemove(baseDescriptors, i);
                ArrayRemove(updatedDescriptors, i);
                ArrayRemove(originalDescriptors, i);
            }
        }

        if (shouldUpdate(BINKSNEEDMOD)) {
            binks.clear();
            metaBinks.clear();
            clearUpdate(BINKSNEEDMOD);
        }

        // normalizing feature strengths (dunno if that's needed?)
        {
            double minStrength = 1.0;
            double maxStrength = 0.0;
            for (int i = baseFeatures.size() - 1; i >= 0; --i) {
                if (baseFeatures[i].strength < minStrength) minStrength = baseFeatures[i].strength;
                if (baseFeatures[i].strength > maxStrength) maxStrength = baseFeatures[i].strength;
            }
            if (minStrength != maxStrength)
                for (int i = baseFeatures.size() - 1; i >= 0; --i) {
                    baseFeatures[i].strength = (baseFeatures[i].strength - minStrength) / (maxStrength - minStrength);
                }
        }

        imGuiMoveWindows = !shouldUpdate(WINDOWFREEZE);
    }

    loading_mutex.unlock();
    DEBUGTRACE("displayStuffEnd");
}

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------

//-----------------
namespace ESBIFT {
    //-----------------
    namespace Evaluation {
        //-----------------
        namespace Disparity { // not fully implemented, yet

            FeatureConfig disparityConfig;

            FeatureConfig createDisparityFeatureConfig(int imageWidth, double descriptorWidthPixels = 1.0) {
                FeatureConfig config;
                const int levels = 1;
                const int pairs = 100;
                const int steps = 100;
                config.refinementSteps = steps;
                config.descriptorPatterns.resize(1);
                DescriptorDesc &pattern = config.descriptorPatterns[0];
                ArrayResize(pattern.levels,levels+1);
                srandom(descriptorRandomSeed);
                for (int i = 0; i < levels; ++i) {
                    DescriptorDescLevel level;
                    for (int j = 0; j < pairs; ++j) {
                        DescriptorPair pair;
                        polarCircleQuadratic(pair);
                        SortDescriptorPair(pair);
#if DEFINED_POLARCOORDINATES == 1
                        toPolar(pair);
#endif
                        level.pairs[j] = pair;
                    }
                    pattern.levels[i] = level;
                }
                return config;
            }

            double disparityFromColor(const PIXEL& pixel) {
                return pixel.r;
            }

#define DONTDOTHAT -10000.0
            double lookupDisparity(DOUBLE x, DOUBLE y, const Image& left, const Image& right, const FeatureConfig& config) {
                // somehow a bit strange way to do that, but whatever (just to somehow compare feature tracker performances)
                // I don't want to implement klt or other feature trackers (here) / dont't want to use opencv or whatever, yet..
                Descriptor descriptor = sampleDescriptor(left, { x, y }, config);
                Feature r1 = lookupFeature(right, descriptor, { x, y }, config);
                r1.y = y;
                Feature r2 = lookupFeature(left, descriptor, r1, config);
                if (fabs(r2.x - x) > 0.5)
                    return DONTDOTHAT;
                return r1.x - r2.x;
            }

            const PIXEL disparityToColor(double disparity) {
                PIXEL r;
                r.a = disparity == DONTDOTHAT ? 0.0 : 1.0;
                if (disparity < 0)
                    r.a *= -1;
                r.r = disparity * r.a;
                r.g = disparity * r.a;
                r.b = disparity * r.a;
                return r;
            }

            void createArtificialDisparityImage(Image& result, const Image& left, const Image& right, const FeatureConfig &config) {
                result = left;
                const int width = result.pixelLayers[0].width;
                const int height = result.pixelLayers[0].height;
                int i = 0;
                for (int y = 0; y < height; ++y) {
                    for (int x = 0; x < width; ++x) {
                        threadPoolThread([&,x,y](){
                            initThread();
                            const double disparity = lookupDisparity(x, y, left, right, config); // positive means right is more right than left
                            result.pixelLayers[0].pixels[x + y * width] = disparityToColor(disparity);
                            });
                    }
                    printf(":%02d", y * 100 / height);
                }
                threadPoolFinish();
            }

            Image disparityImage(const std::string &leftFileName, const std::string& rightFileName) {
                Image artificialDisparity;
                Image left, right;
                ArrayAdd(left.pixelLayers, pixelLayerFromPNG(leftFileName));
                ArrayAdd(right.pixelLayers, pixelLayerFromPNG(rightFileName));
                ArrayAdd(left.featureLayers, buildFeatureLayer(left.pixelLayers[0], 0));
                ArrayAdd(right.featureLayers, buildFeatureLayer(right.pixelLayers[0], 0));
                disparityConfig = createDisparityFeatureConfig(left.pixelLayers[0].width,32.0);
                createArtificialDisparityImage(artificialDisparity, left, right, disparityConfig);
                return artificialDisparity;
            }

            void alphaMul(PIXELLayer& v) {
                for (int i = 0; i < v.width * v.height; ++i) {
                    PIXEL& p = v.pixels[i];
                    double k = p.a;
                    p.r *= k;
                    p.g *= k;
                    p.b *= k;
                    p.a *= 1.0;
                }
            }
        }
        using namespace Disparity;
        //-----------------
    }
    using namespace Evaluation;
    //-----------------
}
using namespace ESBIFT;
//-----------------

extern std::string imGuiWindowName;

void encode(const std::string &fileName) {
    FILE *in = fopen(fileName.c_str(),"rb");
    std::vector<unsigned char> data;
    fseek(in,0,SEEK_END);
    int fileLength = ftell(in);
    data.resize(fileLength);
    fseek(in,0,SEEK_SET);
    fread(&(data[0]),1,fileLength,in);
    fclose(in);
    std::string output = std::string("c:/save/esbiftxor0x37_") + toString(rand()) + std::string(".txt");
    for (int i = 0; i < fileLength; ++i) data[i] ^= 0x37;
    FILE *out = fopen(output.c_str(),"wb");
    fwrite(&(data[0]),1,fileLength,out);
    fclose(out);
}

#define LFCHAR 0x02 // attention here

// c kram, no wstring etc..
std::vector<unsigned char> withoutExoticStuff(const std::vector<unsigned char> &c) {
    std::vector<unsigned char> r;
    r.resize(c.size()); int i = 0;
    for (int k = 0; k < c.size(); ++k) {
        if (c[k] != 0x0d) {
            r[i++] = (c[k] == 0x0a) ? LFCHAR : c[k];
        }
    }
    return r;
}

std::vector<unsigned char> exoticStuffAgain(const std::vector<unsigned char> &c) {
    std::vector<unsigned char> r;
    r.resize(c.size());
    for (int k = 0; k < c.size(); ++k) {
        unsigned char l = c[k];
        if (l == LFCHAR) l = '\n';
        r[k] = l;
    }
    return r;
}

std::string getNoLineFeeds(const std::string &c) {
    std::string r;
    r.resize(c.size()); int i = 0;
    for (int k = 0; k < c.size(); ++k) {
        if (c[k] != 0x0d && c[k] != 0x0a) {
            r[i++] = c[k];
        }
    }
    return r;
}

#if DEFINED_USE_DIFF37_RESOLVE == 1
#include <time.h>
#include <sys/stat.h>
#define CREATIONTIME st_ctime // just for fun and preliminary
#define MODIFICATIONTIME st_mtime // just for fun and preliminary
//#define APPENDINGTIME st_atime // just for fun and preliminary
#define FILETIMESTRING(__file, __timeType) LMAC(std::string) struct _stat t; _fstat(fileno(__file),&t); time_t c = t.__timeType; return getNoLineFeeds(asctime(gmtime(&c))) _LMAC

const int diff37(const std::string &inputFolder, const std::string &outputFileName) {
    printf("diff37\n");
    std::vector<std::pair<FILESYSTEM::file_time_type, std::string>> data;
    for (const auto &entry: FILESYSTEM::directory_iterator(inputFolder)) {
        if (entry.is_regular_file()) {
            data.push_back(std::make_pair(entry.last_write_time(), entry.path().string()));
        }
    }
    printf("files:%d\n", data.size());
    std::sort(data.begin(), data.end(), [&](const auto &a, const auto &b)->bool {
        return a.first < b.first;
    });

    FILE *out = fopen(outputFileName.c_str(), "w");
    if (out == NULL) {
        printf("error: [%s] not creatable.\n", outputFileName.c_str());
        return -1;
    }

#define isLineFeed(__c) ((__c)=='\n')
    for (int i = 1; i < data.size(); ++i) {
        std::vector<unsigned char> data1;
        std::vector<unsigned char> data2;
        const std::string &fileName1 = data[i-1].second;
        const std::string &fileName2 = data[i].second;
        std::string timeString1;
        std::string timeString2;
        {

            const std::string &fileName = fileName1;
            FILE *in = fopen(fileName.c_str(),"rb");
            timeString1 = FILETIMESTRING(in,CREATIONTIME);
            std::vector<unsigned char> data;
            fseek(in,0,SEEK_END);
            int fileLength = ftell(in);
            data.resize(fileLength);
            fseek(in,0,SEEK_SET);
            fread(&(data[0]),1,fileLength,in);
            fclose(in);
            for (int i = 0; i < data.size(); ++i) data[i] ^= 0x37;
            data1 = withoutExoticStuff(data);
        }
        {
            const std::string &fileName = fileName2;
            FILE *in = fopen(fileName.c_str(),"rb");
            timeString2 = FILETIMESTRING(in,CREATIONTIME);
            std::vector<unsigned char> data;
            fseek(in,0,SEEK_END);
            int fileLength = ftell(in);
            data.resize(fileLength);
            fseek(in,0,SEEK_SET);
            fread(&(data[0]),1,fileLength,in);
            fclose(in);
            for (int i = 0; i < data.size(); ++i) data[i] ^= 0x37;
            data2 = withoutExoticStuff(data);
        }
        fprintf(out, ";//##-----------------------------------\n");
        fprintf(out, ";//##--- file1:[%s][%s] -- file2:[%s][%s] ---\n", fileName1.c_str(), timeString1.c_str(), fileName2.c_str(), timeString2.c_str());
        fprintf(out, ";//##---------------------------\n");
        int i1 = 0; int i2 = 0;
        int lineId1 = 0; int lineId2 = 0;
        while (i1 < data1.size() && i2 < data2.size()) {
            const char c1 = data2[i1]; const char c2 = data1[i2];
            if (c1 == c2) {
                i1++;
                i2++;
            } else {
                const int lineId1_b = lineId1;
                const int lineId2_b = lineId2;
                int lineId1_c = lineId1;
                int lineId2_c = lineId2;
                std::string s1, s2;
                int k1 = i1, k2 = i2;
                for (; k1 < data1.size(); ++k1) {
                    const char c1 = data2[k1];
                    const char c2 = data1[i2];
                    if (c1 == c2)
                        break;
                    s1 += c1;
                    if (isLineFeed(c1)) lineId1_c++;
                }
                for (; k2 < data2.size(); ++k2) {
                    const char c1 = data2[i1];
                    const char c2 = data1[k2];
                    if (c1 == c2)
                        break;
                    s2 += c2;
                    if (isLineFeed(c2)) lineId2_c++;
                }
#define LINEID(type) ("#=" + toString(lineId1_b) + "=:" + std::string(type) + ":=" + toString(lineId2_b) + "=>#").c_str()
                if (k2 < k1) {
                    i2 = k2; fprintf(out,"%s>%s", LINEID("FILE2"), s2.c_str()); lineId2 = lineId2_c;
                } else {
                    i1 = k1; fprintf(out,"%s>%s", LINEID("FILE1"), s1.c_str()); lineId1 = lineId1_c;
                }
            }
            if (isLineFeed(c1)) lineId1++;
            if (isLineFeed(c2)) lineId2++;
        }
        fprintf(out, "\n");
    }
    fclose(out);

    return 1;
}
#endif

#endif // DEFINED_USE_ESBIFT

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------

#if DEFINED_INCLUDE_REDUCED_SHADER_PREPROCESSOR == 1
enum {
    TYPETOLOOKFOR_NONE = 0,
    TYPETOLOOKFOR_DEFINEHEADER_START,
    TYPETOLOOKFOR_DEFINEHEADER_BASE,
    TYPETOLOOKFOR_DEFINEHEADER_PARAMETERS,
    TYPETOLOOKFOR_DEFINE_CONTENT_START,
    TYPETOLOOKFOR_DEFINE_CONTENT,
};

using DEFINE = std::pair<std::pair<std::string, std::string>, std::string>;

std::string trim(const std::string &s) { int i = 0; for (i = 0; i < s.length(); ++i) if (s[i] != ' ' && s[i] != '\t' && s[i] != LFCHAR) break; int start = i; for (i = s.length()-1; i >= start ; --i) if (s[i] != ' ' && s[i] != '\t' && s[i] != LFCHAR) break; return s.substr(start,i-start+1); }

const std::string getDefineString(const DEFINE &define) { const bool h = !define.first.second.empty(); return define.first.first + define.first.second + " " + define.second; }

void outputDefines(const std::string &fileName, const std::vector<DEFINE> &defines) {
    FILE *out = fopen(fileName.c_str(), "w");
    for (int i = 0; i < defines.size(); ++i) {
        fprintf(out, "//;------------------==>>>;\n");
        fprintf(out, "#define %s\n", getDefineString(defines[i]).c_str());
    }
    fclose(out);
}

std::vector<std::string> getParameters(const std::string &parameters) {
    std::vector<std::string> r;
    int pos = 0;
    int closure = 0;
    std::string word;
    pos++; // "(" at front and ")" at end
    while(pos < parameters.size()-1) {
        char c = parameters[pos];
        if (c == '(') closure++;  if (c == ')') closure--;
        if (c == '{') closure++;  if (c == '}') closure--;
        if (c == '[') closure++;  if (c == ']') closure--;
        if (c == ',' && closure == 0) {
            r.push_back(trim(word));
            word = "";
        } else
            word += c;
        pos++;
    }
    r.push_back(trim(word));
    return r;
}

std::string replaceMacroParameters(const std::string &macroContent, const std::vector<std::string> &sourceParams, const std::vector<std::string> &destParams) {
    if (sourceParams.size() != destParams.size())
        return macroContent;
    std::string currentReplaceState = macroContent;

    std::vector<int> parameterOccurences;
    parameterOccurences.resize(macroContent.size()); // needs 0 init
    for (int i = 0; i < sourceParams.size(); ++i) {
        const std::string &name = sourceParams[i];
        if (macroContent.length() >= name.length()) {
            for (int k = 0; k < (int) macroContent.length() - name.length(); ++k) {
                int found = k;
                for (int j = 0; j < name.length(); ++j) {
                    if (macroContent[k + j] != name[j]) {
                        found = -1;
                        break;
                    }
                }
                if (found != -1) {
                    parameterOccurences[k] = i + 1;
                }
            }
        }
    }

    int suffixPos = 0;
    for (int i = 0; i < macroContent.length(); ++i) {
        int here = parameterOccurences[i];
        if (here == 0)
            continue;
        here--;
        const std::string &name = sourceParams[here];
        if (currentReplaceState.length() >= name.length()) {
            for (int k = suffixPos; k <= (int) currentReplaceState.length() - name.length(); ++k) {
                int found = k;
                for (int j = 0; j < name.length(); ++j) {
                    if (currentReplaceState[k + j] != name[j]) {
                        found = -1;
                        break;
                    }
                }
                if (found != -1) {
                    const std::string prefix = currentReplaceState.substr(0, k);
                    const int endPos = k + name.length();
                    const std::string suffix = currentReplaceState.substr(endPos, currentReplaceState.length() - endPos);
                    currentReplaceState = prefix + destParams[here] + suffix;
                    suffixPos = currentReplaceState.length() - suffix.length();
                    break;
                }
            }
        }
    }
    return currentReplaceState;
}

void combineVAARGS(std::vector<std::string> &macroParameters, std::vector<std::string> &closureParameters) {
    if (macroParameters[macroParameters.size()-1] == "...") {
        macroParameters[macroParameters.size()-1] = "__VA_ARGS__";
        for (int i = macroParameters.size(); i < closureParameters.size(); ++i) {
            closureParameters[macroParameters.size()-1] += " , " + closureParameters[i];
        }
        closureParameters.resize(macroParameters.size());
    }
}

const std::string resolveBasicMacro(const DEFINE &macro) {
    return macro.second;
}

const std::string resolveClosureMacro(const DEFINE &macro, const std::string &closure) {
    std::vector<std::string> macroParameters = getParameters("(" + macro.first.second);
    std::vector<std::string> closureParameters = getParameters("(" + closure);
    combineVAARGS(macroParameters, closureParameters);
    return replaceMacroParameters(macro.second, macroParameters, closureParameters);
}

const std::string resolved(const DEFINE &macro, const std::string &closure, const std::vector<DEFINE> &defines, int depth = 0) {
    const bool isBasicMacro = macro.first.second.empty();
    if (!isBasicMacro && closure.empty()) printf("macro with parameters but no parameters provided!");
    std::string r = isBasicMacro ? resolveBasicMacro(macro) : resolveClosureMacro(macro, closure);
    return r;
}

const std::string resolve(const std::string &word, const std::vector<DEFINE> &defines) {
    if (word.empty())
        return "";
    std::string currentReplaceState = word;
    for (int k = 0; k < currentReplaceState.length(); ++k) {
        for (int here = 0; here < defines.size(); ++here) {
            const std::string &name = defines[here].first.first;
            if (currentReplaceState.length() >= name.length() && !(name.empty())) {
                int found = k;
                for (int j = 0; j < name.length(); ++j) {
                    if ((k + j) >= currentReplaceState.length() || currentReplaceState[k + j] != name[j]) {
                        found = -1;
                        break;
                    }
                }
                if (found != -1) {
                    const std::string prefix = currentReplaceState.substr(0, k);
                    std::string closure = "";
                    int endPos = k + name.length();
                    if (!defines[here].first.second.empty()) {
                        int closureCount = 1;
                        int t = endPos;
                        while (closureCount > 0 && endPos < currentReplaceState.length()) {
                            if (currentReplaceState[endPos] == '(') closureCount++; if (currentReplaceState[endPos] == ')') closureCount--;
                            if (currentReplaceState[endPos] == '[') closureCount++; if (currentReplaceState[endPos] == ']') closureCount--;
                            if (currentReplaceState[endPos] == '{') closureCount++; if (currentReplaceState[endPos] == '}') closureCount--;
                            // comments may? be some tiny problem here, yet..
                            endPos++;
                        }
                        closure = currentReplaceState.substr(t, endPos - t);
                    }
                    const std::string suffix = currentReplaceState.substr(endPos,
                                                                          currentReplaceState.length() -
                                                                          endPos);
                    const std::string resl = resolved(defines[here], closure, defines);
                    const std::string newReplaceState = prefix + resl;
                    const bool same = newReplaceState == currentReplaceState;
                    currentReplaceState = newReplaceState;
                    k = currentReplaceState.length();
                    currentReplaceState = currentReplaceState + suffix;
                    break;
                }
            }
        }
    }
    return currentReplaceState;
}

std::string stringFromArray(const std::vector<unsigned char> &data) { std::string r; for (int i = 0; i < data.size(); ++i) { r += data[i];} return r;}
std::vector<unsigned char> arrayFromString(const std::string &data) { std::vector<unsigned char> r; for (int i = 0; i < data.size(); ++i) { r.push_back(data[i]);} return r;}

const bool evaluateExpression(const std::string &s, const std::vector<DEFINE> &defines) {
    return true; // to be implemented
}

void preprocessFile(const std::string &inputFileName, const std::string &outputFileName) {
    std::vector<unsigned char> data;
    {
        FILE *in = fopen(inputFileName.c_str(), "rb");
        fseek(in, 0, SEEK_END);
        int fileLength = ftell(in);
        data.resize(fileLength);
        fseek(in, 0, SEEK_SET);
        fread(&(data[0]), 1, fileLength, in);
        data.push_back(LFCHAR);
        fclose(in);
        data = withoutExoticStuff(data);
    }
    std::string nextPass;
    std::vector<DEFINE> defines;
    const int passes = 1;
    for (int pass = 0; pass < passes; ++pass) {
        defines.clear();
        printf("pass:%d [lastPassSize:%d]\n", pass, nextPass.length());
        std::string lastPass = nextPass;
        nextPass.clear();
        int readPos = 0;
        const int fileEnd = data.size();
        bool inPreprocessor = false;
        bool inDefine = false;
        bool lookForDefineHeader = false;
        int nextTypeToLookFor = TYPETOLOOKFOR_NONE;
        int inIf = 0;
        int inElse = -1;
        std::string word;
        std::string currentDefineString;
        DEFINE currentDefine;
        char lc = 0;
        int mathClosureLevel = 0;
        bool additionalIfComments = false;
        std::string additionalComment;

        int activationIndentLevel = 0;
        bool ifDefineSection = false;
        std::set<int> commentedOutIndents;
        std::map<int, std::string> expressionsForIndent;
        std::string expression;

#define READPOSUPDATE lc = c; readPos++;
        int lineNr = 0;
        while (readPos < fileEnd) {

            char c = data[readPos];
            const char nc = ((readPos + 1) < fileEnd) ? data[readPos + 1] : 0;
            const bool macroSharp = (c == '#' && nc == '#');
            if (macroSharp) { READPOSUPDATE; READPOSUPDATE; continue; }
            if (c == '/' && nc == '/') { for (; readPos < fileEnd; ++readPos) {if (data[readPos] == LFCHAR) {break;} additionalComment += data[readPos]; } readPos--; c = lc; READPOSUPDATE; continue; }
            if (c == '/' && nc == '*') { for (; readPos < fileEnd-1; ++readPos) {if (data[readPos] == '*' && data[readPos+1] == '/') {break;}; additionalComment += data[readPos]; }  readPos += 2; readPos--; c = lc; READPOSUPDATE; continue; }
            if (c == '\n') printf("normal linefeed encoutnered error\n");

            if (c == '#') { inPreprocessor = true;}
            if (c == LFCHAR && lc != '\\') { inPreprocessor = false; inDefine = false; lookForDefineHeader = false; }
            if (c == LFCHAR) lineNr++;

#define ADDADDITIONALCOMMENT { if (!additionalComment.empty()) { nextPass += additionalComment; additionalComment = ""; } }

#define COMMENTEDOUTSECTION ( commentedOutIndents.find(activationIndentLevel) != commentedOutIndents.end() )
#define COMMENTINSECTION commentedOutIndents.erase(activationIndentLevel);
#define COMMENTOUTSECTION commentedOutIndents.insert(activationIndentLevel);
#define INVERTSECTION if (COMMENTEDOUTSECTION) COMMENTINSECTION else COMMENTOUTSECTION
#define EVALUATEEXPRESSION(__e) evaluateExpression(__e,defines)
#define MAYBEEVALUATEIFSECTION {\
if (ifDefineSection) {\
if (c == LFCHAR && lc != '\\') { ifDefineSection = false;\
    if (pass == 0) {\
        expressionsForIndent[activationIndentLevel] = expression; additionalComment = "// " + expression;\
    }\
    if(!EVALUATEEXPRESSION(expression)) {COMMENTOUTSECTION; nextPass += "/*";}                           \
    expression = "";            \
    }       \
    else { expression += c;}\
}}
            MAYBEEVALUATEIFSECTION

            if (nextTypeToLookFor != TYPETOLOOKFOR_NONE) {
                ADDADDITIONALCOMMENT
                nextPass += c;
                while (1) {
                    bool redo = true;
                    const int typeHere = nextTypeToLookFor;
                    switch (typeHere) {
                        case TYPETOLOOKFOR_DEFINEHEADER_START: {
                            if (c != ' ' && c != '\t') nextTypeToLookFor = TYPETOLOOKFOR_DEFINEHEADER_BASE;
                            if (c == LFCHAR && lc != '\\') nextTypeToLookFor = TYPETOLOOKFOR_NONE;
                            break;
                        }
                        case TYPETOLOOKFOR_DEFINEHEADER_BASE: {
                            word += c;
                            if (c == '(') {
                                currentDefine.first.first = trim(word);
                                nextTypeToLookFor = TYPETOLOOKFOR_DEFINEHEADER_PARAMETERS;
                                redo = false;
                                word = "";
                                break;
                            }
                            if (c == LFCHAR && lc == '\\') {
                                currentDefine.first.first = trim(word.substr(0, word.length() - 2));
                                word = "";
                                nextTypeToLookFor = TYPETOLOOKFOR_DEFINE_CONTENT_START;
                                redo = false;
                                break;
                            }
                            if (c == ' ' || c == '\t' || c == LFCHAR) {
                                currentDefine.first.first = trim(word);
                                word = "";
                                nextTypeToLookFor = TYPETOLOOKFOR_DEFINE_CONTENT_START;
                            }
                            break;
                        }
                        case TYPETOLOOKFOR_DEFINEHEADER_PARAMETERS: {
                            word += c;
                            if (c == ')') {
                                currentDefine.first.second = word;
                                word = "";
                                nextTypeToLookFor = TYPETOLOOKFOR_DEFINE_CONTENT_START;
                                redo = false;
                                break;
                            }
                            break;
                        }
                        case TYPETOLOOKFOR_DEFINE_CONTENT_START: {
                            if (c != ' ' && c != '\t') nextTypeToLookFor = TYPETOLOOKFOR_DEFINE_CONTENT;
                            if (c == LFCHAR) nextTypeToLookFor = TYPETOLOOKFOR_DEFINE_CONTENT; // not perfect but whatever
                            break;
                        }
                        case TYPETOLOOKFOR_DEFINE_CONTENT: {
                            if (c == LFCHAR && lc == '\\') {
                                word = word.substr(0, word.length() - 1) + c;
                                redo = false;
                                break;
                            }
                            if (c == LFCHAR && lc != '\\') {
                                currentDefine.second = trim(word);
                                word = "";
                                bool found = false;
                                for (int i = 0; i < defines.size(); ++i) {
                                    if (defines[i].first == currentDefine.first) {
                                        defines[i] = currentDefine;
                                        found = true;
                                    }
                                }
                                if (!found)
                                    defines.push_back(currentDefine);
                                std::string defineString = getDefineString(currentDefine);
                                //fprintf(out, "/* DEFINE %s */\n", defineString.c_str());
                                nextTypeToLookFor = TYPETOLOOKFOR_NONE;
                                break;
                            }
                            word += c;
                            break;
                        }
                    }
                    if ((!redo) || (redo && (typeHere == nextTypeToLookFor)))
                        break;
                }
            } else {
                if (c == '(') mathClosureLevel++;
                if (c == ')') mathClosureLevel--;
                if (macroSharp || ((c == '\t' || c == ' ') && mathClosureLevel == 0) || c == LFCHAR) {
                        word = resolve(word, defines);
                        nextPass += word;
                        ADDADDITIONALCOMMENT
                        nextPass += c;
                        word = "";
                } else {
                    word += c;
                }
            }

//----------------------------
#define DEFINED_MAC_PREPROCESSOR_DEFINES 1
#define DEFINED_MAC_PREPROCESSOR_IFS 1
//----------------------------
#if DEFINED_MAC_PREPROCESSOR_DEFINES == 1
            if (inPreprocessor && word == "#define") {
                inDefine = true;
                nextPass += word;
                word = "";
                currentDefine = DEFINE();
                nextTypeToLookFor = TYPETOLOOKFOR_DEFINEHEADER_START;
            }
#endif
#if DEFINED_MAC_PREPROCESSOR_IFS == 1
            if (inPreprocessor && word == "#if" && nc != 'd') { ifDefineSection = true; activationIndentLevel++; inIf++; nextPass += word; word = ""; expression = ""; additionalIfComments = true; }
            if (inPreprocessor && word == "#if" && nc == 'd') { ifDefineSection = true; activationIndentLevel++; inIf++; nextPass += word; word = ""; expression = ""; additionalIfComments = true; }
            if (inPreprocessor && word == "#else") { nextPass += (COMMENTEDOUTSECTION ? "*/" : ""); INVERTSECTION; nextPass += word; nextPass += (COMMENTEDOUTSECTION ? "/*" : ""); word = ""; additionalComment = "// " + expressionsForIndent[activationIndentLevel]; }
            if (inPreprocessor && word == "#endif") { inIf--; if (inIf < 0) printf("macro Error too much #endifs\n"); nextPass += (COMMENTEDOUTSECTION ? "*/" : "") + word; word = ""; COMMENTINSECTION; ; additionalComment = "// " + expressionsForIndent[activationIndentLevel]; activationIndentLevel--; }
#endif
//----------------------------
            READPOSUPDATE
        }

        data = arrayFromString(nextPass);
        FILE *out = fopen(outputFileName.c_str(), "w");
        std::vector<unsigned char> data2 = exoticStuffAgain(data);
        fwrite(&(data2[0]),1,data2.size(),out);
        fclose(out);
        if (nextPass.size() == lastPass.size()) {
            bool found = true;
            for (int i = 0; i < lastPass.size(); ++i) {
                if (lastPass[i] != nextPass[i]) {
                    found = false;
                    break;
                }
            }
            if (found)
                break;
        }
    }
    outputDefines(outputFileName + "_defines" + toString(0) + ".txt", defines);
}
#endif // DEFINED_INCLUDE_REDUCED_SHADER_PREPROCESSOR

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------

#if UNITTEST == 1
#define UNIT_HEADER 0
#include "esbift_unit.hpp"
#endif


int main(int argc, char** argv)
{
    time_t t; time(&t); srand(t);

    encode("c:/!mad/Projekte/Games/newstuff/binaryfeatures/source/esbift.cpp");
    COMMENTOUT(diff37("c:/save", "c:/save_all/esbifti" + std::to_string(rand()) + "_all.txt"));
    preprocessFile("c:/!mad/Projekte/Games/newstuff/binaryfeatures/source/esbift.cpp", "c:/!mad/Projekte/Games/newstuff/binaryfeatures/source/somemore_esbift.cpp");

#if DEFINED_USE_ESBIFT == 1
    initMultiCore();
    imGuiWindowName = "ESBIFT - Extremely Simple Brightness Invariant Feature Tracking";  initImGui();

    UNITTEST_FUNC

    requestUpdate(RELOAD);
    collectDirectories(Configuration::dataSetsFolder);

    if (false) {
        std::string disparityFolder = "data/disparity/1/";
        Image artificialDisparity = disparityImage(disparityFolder + "leftright0001.png", disparityFolder + "leftright0002.png");
        artificialDisparity.pixelLayers[0] = normalize(artificialDisparity.pixelLayers[0]);
        alphaMul(artificialDisparity.pixelLayers[0]);
        savePixelLayer(artificialDisparity.pixelLayers[0], disparityFolder + "disparityImage.png");
    }

    while (mainLoopImGui());
    closeImGui();
#endif

    return 0;
}

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------


//-----------------
namespace ESBIFT {
    //-----------------
    namespace IO {
        //-----------------
        namespace BasicImageIO {

#define JPG_INTERFACE 1
#include "jpeg/jpeg3.cpp"

#include "png/png.h"

            bool loadJPG(const std::string& jpg) {
                FILE* f = fopen(jpg.c_str(), "rb");
                if (!f) {
                    return false;
                }
                fseek(f, 0, SEEK_END);
                int size = (int)ftell(f);
                char* buf = (char*)malloc(size);
                fseek(f, 0, SEEK_SET);
                size = (int)fread(buf, 1, size, f);
                fclose(f);

                njInit();
                if (njDecode(buf, size)) {
                    free((void*)buf);
                    return false;
                }
                free((void*)buf);

                pictureWidth = njGetWidth();
                pictureHeight = njGetHeight();
                unsigned char* p = njGetImage();
                for (int y = 0; y < pictureHeight; ++y) {
                    for (int x = 0; x < pictureWidth; ++x) {
                        pictureS[x + y * pictureWidth] = (((int)p[0]) << 0) | (((int)p[1]) << 8) | (((int)p[2]) << 16) | 0xff000000;
                        p += 3;
                    }
                }
                njDone();
                return true;
            }
            bool loadPNG(const char* png) // I - File to read
            {
                int		i;			// Looping var
                FILE* fp;			// File pointer
                int		channels;		// Number of color channels
                png_structp	pp;			// PNG read pointer
                png_infop	info;			// PNG info pointers
                png_bytep* rows;			// PNG row pointers

                memset(pictureS, 0, sizeof(pictureS));

                // Open the PNG file...
                if ((fp = fopen(png, "rb")) == NULL) return false;

                // Setup the PNG data structures...
                pp = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
                info = png_create_info_struct(pp);

                if (setjmp(pp->jmpbuf))
                {
                    int debug = 1;
                    return false;
                }

                // Initialize the PNG read "engine"...
                png_init_io(pp, fp);

                // Get the image dimensions and convert to grayscale or RGB...
                png_read_info(pp, info);

                if (info->color_type == PNG_COLOR_TYPE_PALETTE)
                    png_set_expand(pp);

                if (info->color_type & PNG_COLOR_MASK_COLOR)
                    channels = 3;
                else
                    channels = 1;

                if ((info->color_type & PNG_COLOR_MASK_ALPHA) || info->num_trans)
                    channels++;

                int w = (int)(info->width);
                int h = (int)(info->height);
                int d = channels;
                pictureWidth = w;
                pictureHeight = h;

                if (info->bit_depth < 8)
                {
                    png_set_packing(pp);
                    png_set_expand(pp);
                }
                else if (info->bit_depth == 16)
                    png_set_strip_16(pp);

#  if defined(HAVE_PNG_GET_VALID) && defined(HAVE_PNG_SET_TRNS_TO_ALPHA)
                // Handle transparency...
                if (png_get_valid(pp, info, PNG_INFO_tRNS))
                    png_set_tRNS_to_alpha(pp);
#  endif // HAVE_PNG_GET_VALID && HAVE_PNG_SET_TRNS_TO_ALPHA

                unsigned char* array = (unsigned char*)pictureS;

                // Allocate pointers...
                rows = new png_bytep[h];

                for (i = 0; i < h; i++)
                    rows[i] = (png_bytep)(array + i * w * d); // we flip it

                    // Read the image, handling interlacing as needed...
                for (i = png_set_interlace_handling(pp); i > 0; i--)
                    png_read_rows(pp, rows, NULL, h);

#ifdef WIN32
                // Some Windows graphics drivers don't honor transparency when RGB == white
                if (channels == 4)
                {
                    // Convert RGB to 0 when alpha == 0...
                    unsigned char* ptr = (unsigned char*)array;
                    for (i = w * h; i > 0; i--, ptr += 4)
                        if (!ptr[3]) ptr[0] = ptr[1] = ptr[2] = 0;
                }
#endif // WIN32

                if (channels == 3)
                {
                    unsigned char* array2 = new unsigned char[pictureWidth * pictureHeight * 4];
                    for (int i = w * h - 1; i >= 0; --i)
                    {
                        array2[i * 4 + 0] = array[i * 3 + 0];
                        array2[i * 4 + 1] = array[i * 3 + 1];
                        array2[i * 4 + 2] = array[i * 3 + 2];
                        array2[i * 4 + 3] = 255;
                    }
                    memcpy(array, array2, w * h * 4);
                    delete[] array2;
                }

                // Free memory and return...
                delete[] rows;

                png_read_end(pp, info);
                png_destroy_read_struct(&pp, &info, NULL);

                fclose(fp);
                return true;
            }

            void abort_(const char* s, ...)
            {
                va_list args;
                va_start(args, s);
                vfprintf(stderr, s, args);
                fprintf(stderr, "\n");
                va_end(args);
                abort();
            }

            void savePNG(const char* file_name)
            {
                png_structp png_ptr;
                png_infop info_ptr;

                // create file
                FILE* fp = fopen(file_name, "wb");
                if (!fp)
                    abort_("[write_png_file] File %s could not be opened for writing", file_name);


                // initialize stuff
                png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);

                if (!png_ptr) {
                    fclose(fp);
                    abort_("[write_png_file] png_create_write_struct failed");
                }

                info_ptr = png_create_info_struct(png_ptr);
                if (!info_ptr) {
                    fclose(fp);
                    abort_("[write_png_file] png_create_info_struct failed");
                }

                if (setjmp(png_jmpbuf(png_ptr))) {
                    fclose(fp);
                    abort_("[write_png_file] Error during init_io");
                }

                png_init_io(png_ptr, fp);

                // write header
                if (setjmp(png_jmpbuf(png_ptr))) {
                    fclose(fp);
                    abort_("[write_png_file] Error during writing header");
                }

                png_set_IHDR(png_ptr, info_ptr, pictureWidth, pictureHeight,
                    8, PNG_COLOR_TYPE_RGBA, PNG_INTERLACE_NONE,
                    PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);

                png_write_info(png_ptr, info_ptr);


                // write bytes
                if (setjmp(png_jmpbuf(png_ptr))) {
                    fclose(fp);
                    abort_("[write_png_file] Error during writing bytes");
                }

                png_bytep* row_pointers;

                unsigned char* array = new unsigned char[pictureWidth * pictureHeight * 4 + 4];
                for (int i = pictureWidth * pictureHeight; i >= 0; --i)
                {
                    array[i * 4 + 0] = ((unsigned char*)pictureWriteOut)[i * 4 + 0];
                    array[i * 4 + 1] = ((unsigned char*)pictureWriteOut)[i * 4 + 1];
                    array[i * 4 + 2] = ((unsigned char*)pictureWriteOut)[i * 4 + 2];
                    array[i * 4 + 3] = ((unsigned char*)pictureWriteOut)[i * 4 + 3]*0+255;
                }


                // Allocate pointers...
                row_pointers = new png_bytep[pictureHeight];

                for (int i = 0; i < pictureHeight; i++)
                    row_pointers[i] = (png_bytep)(array + (i)*pictureWidth * 4); // we flip it

                png_write_image(png_ptr, row_pointers);

                // end write
                if (setjmp(png_jmpbuf(png_ptr)))
                    abort_("[write_png_file] Error during end of write");

                png_write_end(png_ptr, NULL);

                // cleanup heap allocation
                delete[] row_pointers;

                fclose(fp);
            }
        }
    }
}