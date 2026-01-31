// The color space to use.
//#define YOZ
#define LINEAR
//
// Derivation/amplification of the standard output gamma correction of 1/2.4
#define OUTPUT_GAMMA 1.0
//
// Set this define for the opponend color space as in CIEDE2000
//#define OPPONENT_COLOR
//
#define KNEE_VALUE
// Use the PSI nonlinearity output function.
// Works only if LINEAR is set, then OUTPUT_GAMMA
// does not matter.
//#define USE_PSI
//
// The quantization matrices to use
// This is the DCTune setting.
#define AHUMADA
#define AHUMADA_EXPONENT 3.5f // usually 3.5
#define AHUMADA_FACTOR   1.0f // usually 1.0
//
// Shall the probability be weighted with the
// square error?
// Otherwise, probabilities are summed up.
//#define WEIGHT_MSE
//
// Shall the probability be weighted with the
// DeltaE color difference value? Otherwise,
// probabilities are summed up.
//#define WEIGHT_DELTA_E
//
// Daly instead of taubman masking
//#define DALY_MASKING
//
// Additional option to build another low-pass from the low-pass
#define EXTENDED_FILTER
//
// The input gamma exponent
#define INPUT_GAMMA 1.0
//
// The (inverse) base visibility.
#define BASE_VISIBILITY 0.08f
//
// Detection threshold 
#define DETECT_THREAS 1.0f //0.75 //1.0 // was: 1.5
