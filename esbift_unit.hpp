#if UNIT_HEADER == 1
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
#define ArrayResize(__t, __k) (__t).resize(__k)
#define ArraySize(__t) (__t).size()
#define ArrayAdd(__t, __a) (__t).push_back(__a)
#define ArrayAppend(__t, __a) SMAC(__t,__a,&t,const &a) t.insert(t.end(),a.begin(),a.end()) _SMAC
#define ArrayRemove(__t, __a) (__t).erase((__t).begin() + (__a))
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

#define SWAP(__t, __a, __b) {__t = __a; __a = __b; __b = __t;}
#define FIXUPDIVBYZERO(__v, __r) ((__v) == 0.0 ? (__r) : (__v))
#define DEGREETORADIANS(__degrees) ((__degrees) * PI2 / 360.0)
#define POSITIONTOVALUE(__x, __y)  ( (__x) / ((__y)+sign(__y)*3.0)*2.0*fabs(__x) - fabs(__y) / fabs((__x)+sign(__x)*1.5)*0.5/fabs(__y) ) // the sign stuff is to prevent division by zero singularity
#define LMAC(__returntype, __p0, __p1, __v0, __v1) [&]() -> __returntype {auto &__v0 = __p0; auto __v1 = __p1; // LAMBDA MACRO
#define LMAC(__returntype) [&]() -> __returntype {
#define LMAC_ [&]() {
#define _LMAC ;}()
//#define SMAC(__p0,__p1,__p2,__v0,__v1,__v2) {auto __v0 = __p0; auto __v1 = __p1; auto __v2 = __p2; // SCOPE MACRO
#define SMAC(__p0, __p1, __v0, __v1) {auto __v0 = __p0; auto __v1 = __p1;
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
#define FIXFLOATNUMBER(__v, __r) if (!VALIDFLOAT(__v)) __v = (__r);
#define FIXFLOATNUMBER2(__v1, __v2, __v3) ((__v1) == (__v2) ? (__v1) : (v__3))
#else
#define VALIDFLOAT(__v) true
#define FIXFLOATNUMBER(__v,__r)
#define FIXFLOATNUMBER2(__v1,__v2,__v3)
#endif

#define SortDescriptorPair descriptorPairNoSort //descriptorPairRightDown

#if DEFINED_CHANNELCOUNT == 1
#define FeatureLayerType 1
#define FeatureLayerDepthFormula(__value, __layer) __value;//(pow(__layer * featureLayerBaseLayerMul + featureLayerBaseBase,__value))
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

//--------------------------------------------------------------------------------------------
// threading stuff
//--------------------------------------------------------------------------------------------
#define srandom(__a) srandomWithThreadId(__a)
#define random(__a) randomWithThreadId(__a)

    namespace Unit_Configuration {

        namespace DEFAULTS {
            namespace DEFAULTMIP {
                const double descriptorLevelScale = DEFAULTMIPLEVELSCALING;
                const double stepFactorIncrease = 1.0;
            }
        }

        //-----------------
        namespace ITERATIVE {

            bool iterative = true;

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

            bool iterative = false;

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
    //using namespace Unit_Configuration;
    //-----------------
#endif


#if UNIT_HEADER == 0
    //-----------------
    namespace Unit {

        using namespace Unit_Configuration;

        // no multi threading, yet
        int unitTestStartLine = 0;
        int unitTestLineIdAdder = 0;
        Array(int) unitTestStartLines;
        int unitTick = 0;
#define uprintf(...) {unitTick++; if ((unitTick % 10) == 0) printf("."); if ((unitTick % 100) == 0) printf("%d", unitTick); FILE *out = fopen((configName + ".txt").c_str(),"a+"); fprintf(out, "[0x%x]", (__LINE__ - unitTestStartLine) + unitTestLineIdAdder); fprintf(out, __VA_ARGS__); fclose(out); }
#define vprintf(...) {const int __rb = unitTestStartLine; unitTestStartLine = __LINE__; vcprintf(configName, __VA_ARGS__); unitTestStartLine = __rb;}
#define upush { unitTestStartLine = __LINE__; unitTestStartLines.push_back(unitTestStartLine); unitTestLineIdAdder = 0;} // todo: implement later for (int i = 0; i < ArraySize(unitTestStartLines); ++i) unitTestLineIdAdder += unitTestStartLines[i];}
#define upop {unitTestStartLine = unitTestStartLines.back(); unitTestStartLines.pop_back();}
#define uCompareWith(__v) vcCompareWith(configName,__v)

        Array(unsigned char) removeLineIds(const Array(unsigned char) &data) {
            Array(unsigned char) r;
            bool scanHere = true;
            for (int i = 0 ; i < ArraySize(data); ++i) {
                unsigned char c = data[i];
                if (c == '\n') {
                    scanHere = true;
                }
                if (scanHere) {
                    if (c == ']') scanHere = false;
                }
                else
                    ArrayAdd(r, c);
            }
            return r;
        }

        std::string pathComponent(const std::string &fileName) {
            for (int i = fileName.size()-1; i >= 0; --i)
                if (fileName[i] == '/' || fileName[i] == '\\')
                    return fileName.substr(0, i);
            return fileName;
        }

        std::string fileExtension(const std::string &fileName) {
            for (int i = fileName.size()-1; i >= 0; --i)
                if (fileName[i] == '.')
                    return fileName.substr(i + 1, fileName.size() - (i + 1));
            return fileName;
        }

        std::string nameComponent(const std::string &fileName) {
            for (int i = fileName.size()-1; i >= 0; --i)
                if (fileName[i] == '.') {
                    const int ei = i;
                    for (; i >= 0; --i) {
                        if (fileName[i] == '/' || fileName[i] == '\\')
                            return fileName.substr(i + 1, ei - (i + 1));
                    }
                }
            return fileName;
        }

        const bool rgbCompare(unsigned int a, unsigned int b) {
            return a == b;
        }

        const bool comparePNGS(const std::string &configName, const std::string &secondConfigName) {

            const std::string folder1 = pathComponent(configName);
            const std::string folder2 = pathComponent(secondConfigName);

            Array(std::string) fileNames1;
            for (const auto &entry: FILESYSTEM::directory_iterator(folder1)) {
                if (entry.is_regular_file()) {
                    const std::string fileName = entry.path().string();
                    if (fileExtension(fileName) == "png") {
                        printf("unittest png file:[%s]\n", fileName.c_str());
                        ArrayAdd(fileNames1, fileName);
                    }
                }
            }

            bool allgood = true;
            for (int i = 0; i < ArraySize(fileNames1); ++i) {
                const std::string &fileName1 = fileNames1[i];
                const std::string fileName2 = pathComponent(secondConfigName) + "/" + nameComponent(fileName1) + "." + fileExtension(fileName1);
                printf("comparing pngs: [%s] [%s]", fileName1.c_str(), fileName2.c_str());
                if (!loadPNG(fileName1.c_str())) {
                    printf(" file not found:[%s]", fileName1.c_str());
                }
                Array(unsigned int) r; ArrayResize(r, pictureWidth * pictureHeight);
                for (int i = 0; i < pictureWidth * pictureHeight; ++i) {
                    r[i] = pictureS[i];
                }
                if (!loadPNG(fileName2.c_str())) {
                    printf(" file not found:[%s]", fileName2.c_str());
                }
                bool alright = true;
                for (int i = 0; i < pictureWidth * pictureHeight; ++i) {
                    if (!rgbCompare(r[i],pictureS[i])) {
                        alright = false;
                        break;
                    }
                }
                printf(alright ? " succeeded\n" : " failed\n" );
                allgood &= alright;
            }

            printf("!!! png unittest %s !!!\n", allgood ? "succeeded" : "failed");
            return allgood;
        }

        bool vcCompareWith(const std::string &configName, const std::string &secondConfigName, const bool withLineIds = false) {
            printf("\n");
            std::vector<unsigned char> data1;
            std::vector<unsigned char> data2;
            {
                FILE *in = fopen((configName + ".txt").c_str(),"rb");
                fseek(in,0,SEEK_END);
                int fileSize = ftell(in);
                fseek(in,0,SEEK_SET);
                ArrayResize(data1,fileSize);
                fread(&(data1[0]),1,fileSize,in);
                fclose(in);
            }
            {
                FILE *in = fopen((secondConfigName + ".txt").c_str(),"rb");
                fseek(in,0,SEEK_END);
                int fileSize = ftell(in);
                fseek(in,0,SEEK_SET);
                ArrayResize(data2,fileSize);
                fread(&(data2[0]),1,fileSize,in);
                fclose(in);
            }
            if (!withLineIds) {data1 = removeLineIds(data1); data2 = removeLineIds(data2);}
            if (ArraySize(data1) != ArraySize(data2)) {
                printf("!!! unittest failed !!!\n");
                printf("error: files differ in size %s:[%s[%d]] [%s[%d]]\n", withLineIds ? "" : "after removing line ids", (configName + ".txt").c_str(), data1.size(), (secondConfigName + ".txt").c_str(), data2.size() );
                return false;
            }
            std::string currentLine;
            int lnr = 0;
            for (int i = 0; i < ArraySize(data1); ++i) {
                const unsigned char c = data1[i];
                currentLine += c;
                if (data1[i] != data2[i]) {
                    for (i ; i < ArraySize(data1); ++i) {
                        if (data1[i] == '\n')
                            break;
                        currentLine += data1[i];
                    }
                    printf("!!! unittest failed !!!\n");
                    printf("error in file:%s line:%d content:%s\n", (configName + ".txt").c_str(), lnr, currentLine.c_str());
                    return false;
                }
                if (c == '\n') {
                    currentLine = "";
                    lnr++;
                }
            }
            printf("!!! data unittest succeeded !!!\n");
            bool additionalTests = comparePNGS(configName, secondConfigName);
            printf("---------------------------------------------\n");
            printf("---------------------------------------------\n");
            printf("---- !!!   UNITTEST !!%s!!   !!! ----\n", additionalTests ? "SUCEEDED" : "FAILED");
            printf("---------------------------------------------\n");
            printf("---------------------------------------------\n");
        }

        void vcprintf(const std::string &configName, const Image &image) {
            upush;
            uprintf("image[\"%s\"] featureLayers:%d pixelLayers:%d\n", image.name.c_str(), ArraySize(image.featureLayers), ArraySize(image.pixelLayers) );
            for (int i = 0; i < ArraySize(image.featureLayers); ++i) uprintf("->featureLayer: width: %d height: %d - fwidth: %f fheight: %f\n",
                                                                     image.featureLayers[i].width, image.featureLayers[i].height,
                                                                     image.featureLayers[i].fwidth, image.featureLayers[i].fheight);
            for (int i = 0; i < ArraySize(image.pixelLayers); ++i) uprintf("->pixelLayer: width: %d height: %d\n",
                                                                     image.featureLayers[i].width, image.featureLayers[i].height);
            upop;
        }

        void vcprintf(const std::string &configName, const char *s, double v) {
            upush;
            uprintf("%s = %f\n", s, v);
            upop;
        }

        void vcprintf(const std::string &configName, const char *s, const DescriptorPair &pair) {
            upush;
            uprintf("%s->%f,%f,%f,%f\n", s, pair.x0, pair.y0, pair.x1, pair.y1);
            upop;
        }

        void vcprintf(const std::string &configName, const char *s, const DescriptorDescLevel &level) {
            upush;
            for (int i = 0; i < ArraySize(level.pairs); ++i) {
                vprintf((std::string(s) + "->pair[" + toString(i) + "]").c_str(), level.pairs[i]);
            }
            upop;
        }

        void vcprintf(const std::string &configName, const char *s, const DescriptorDesc &desc) {
            upush;
            for (int i = 0; i < ArraySize(desc.levels); ++i) {
                vprintf((std::string(s) + "->level[" + toString(i) + "]").c_str(), desc.levels[i]);
            }
            upop;
        }

        void vcprintf(const std::string &configName, const char *s, const Array(DescriptorDesc) &descriptorPatterns) {
            upush;
            for (int i = 0; i < ArraySize(descriptorPatterns); ++i) {
                vprintf((std::string("->pattern[") + toString(i)+"]").c_str(),descriptorPatterns[i]);
            }
            upop;
        }

        void vcprintf(const std::string &configName, const char *s, const Descriptor &descriptor) {
            upush;
            uprintf("%s\n", s);
            uprintf("levels:%d\n", ArraySize(descriptor.bits));
            for (int r = 0; r < ArraySize(descriptor.bits); ++r) {
                std::string v;
                for (int j = 0; j < ArraySize(descriptor.bits[r]); ++j) {
                    v += toString(descriptor.bits[r][j]);
                }
                const int angleCount = ((DEFINED_WITHROTATIONINVARIANCE == 1) && (DEFINED_FULLROTATIONINVARIANCEBYSECTORS == 1)) ? 8 : 1;
                uprintf("->bits:%d for level:%d -> [%d*%d=%d]%s\n", ArraySize(descriptor.bits[r]), r, v.length() / angleCount, angleCount, v.length(), v.c_str());
            }
            upop;
        }

        void vcprintf(const std::string &configName, const FeatureConfig &config) {
            upush;
            uprintf("featureConfig\n" );
            vprintf("descriptorPatterns",config.descriptorPatterns);
            vprintf("deepestLevel",config.deepestLevel);
            vprintf("lastLevel",config.lastLevel);
            vprintf("sampleLevel",config.sampleLevel);
            vprintf("scaleInvariance",config.scaleInvariance);
            vprintf("coordReferenceLevel",config.coordReferenceLevel);
            vprintf("refinementSteps",config.refinementSteps);
            vprintf("iterative",config.iterative);
            vprintf("stepFactor",config.stepFactor);
            vprintf("refinementStepPerLevelScaleFactor",config.refinementStepPerLevelScaleFactor);
            vprintf("descriptorSizeCoarse",config.descriptorSizeCoarse);
            vprintf("descriptorSizeFine",config.descriptorSizeFine);
            vprintf("descriptorLevelScale",config.descriptorLevelScale);
            vprintf("sourceToDestScaleRatio",config.sourceToDestScaleRatio);
            vprintf("scaleLevelChecksCount",config.scaleLevelChecksCount);
            vprintf("descriptorShapeSize",config.descriptorShapeSize);
            upop;
        }

        void unittest() {
            upush;
            printf("unittest start\n"); unitTick = 0;
            const std::string &configName = "c:/!mad/UnitTest/ESBIFT/esbift"; FILE *out = fopen((configName + ".txt").c_str(),"w"); fclose(out);
            Image image = imageFromFileNameExtension("c:/!mad/unittest.png", 0, 640, 0.0, 1.0, 1.0);
            preprocessImage(image);
            const int imagePyramidSize = setUpImagePyramid(image);
            vprintf(image);
            FeatureConfig config = createDefaultFeatureConfig(imagePyramidSize);
            vprintf(config);

            Array(Descriptor) descriptors; Array(Feature) features; Array(Feature) foundFeatures;
            ArrayResize(descriptors, 20);
            ArrayResize(features, ArraySize(descriptors));
            ArrayResize(foundFeatures, ArraySize(descriptors));

            uprintf("descriptorSize:%d\n", ArraySize(descriptors));
            uprintf("features:%d\n", ArraySize(features));
            uprintf("foundFeatures:%d\n", ArraySize(foundFeatures));

            int cn = 0;
            for (int i = 0; i < ArraySize(descriptors); ++i) {

                Descriptor &descriptor = descriptors[i];

                features[i].x = randomLike(i*99, 300);
                features[i].y = randomLike(i*777, 300);

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
                            int o = 0;
                            int nonValids = 0;
                            if (imageLevel < ArraySize(image.featureLayers)) {
                                const int descriptorChannels = channelElementCount();
                                Array(VECTORBOOL) &desc = descriptor.bits[descriptorLevel + (angle45Id * descriptorLevelCount) + (shapeId * descriptorLevelCount * angleCount)];
                                ArrayResize(desc, ArraySize(level.pairs) * descriptorChannels);
                                for (int j = 0; j < ArraySize(level.pairs); ++j) {
                                    //const DescriptorPair &pair = level.pairs[j];
                                    for (int n = 0; n < descriptorChannels; ++n) {
                                        desc[o++] = (int)randomLike(121+cn*10,4);
                                        cn++;
                                    }
                                }
                            }
                        }
                    }
                }
                vprintf("descriptor", descriptors[i]);
            }

            int currentImage = 0;
            for (int iterative = 0; iterative < 2; ++iterative) {
                for (int a = 0; a < Unit_Configuration::algorithmCount; ++a) {
                    Unit_Configuration::iterative = (iterative == 1) ? true : false;
                    Unit_Configuration::algorithm = a;
                    uprintf("tracking features for config:%d\n", currentImage);
                    for (int i = 0; i < ArraySize(features); ++i) {
                        threadPoolThread([&, i]() {
                            initThread();
                            const Feature f1 = Feature(image.featureLayers[0].width * 0.5, image.featureLayers[0].height * 0.5, features[i]);
                            const Feature f2 = features[i];
                            Feature feature = lookupFeature(image, descriptors[i], baseConfig.iterative ? f2 : f1, config);
                            foundFeatures[i] = feature;
                            });
                    }
                    threadPoolFinish();
                    char buffer[1000]; sprintf(buffer, "%s_png%04d.png", configName.c_str(), currentImage);
                    printf("\nsaving image to:%s\n", buffer);
                    savePixelLayerWithFeatures(image.pixelLayers[0], foundFeatures, buffer, 0);
                    currentImage++;
                }
            }
            uCompareWith("c:/!mad/UnitTest/ESBIFT/base1/esbift");
            upop;
        }

    }
    using namespace Unit;
    //-----------------
#endif