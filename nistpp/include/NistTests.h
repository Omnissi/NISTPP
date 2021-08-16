#pragma once
/* --------------------------------------------------------------------------
   Title       :  The NIST Statistical Test Suite

   Date        :  December 1999

   Programmer  :  Juan Soto

   Summary     :  For use in the evaluation of the randomness of bitstreams
                  produced by cryptographic random number generators.

   Package     :  Version 1.0

   Copyright   :  (c) 1999 by the National Institute Of Standards & Technology

   History     :  Version 1.0 by J. Soto, October 1999
                  Revised by J. Soto, November 1999
                  Revised by Larry Bassham, March 2008

   Keywords    :  Pseudorandom Number Generator (PRNG), Randomness, Statistical
                  Tests, Complementary Error functions, Incomplete Gamma
                  Function, Random Walks, Rank, Fast Fourier Transform,
                  Template, Cryptographically Secure PRNG (CSPRNG),
                  Approximate Entropy (ApEn), Secure Hash Algorithm (SHA-1),
                  Blum-Blum-Shub (BBS) CSPRNG, Micali-Schnorr (MS) CSPRNG,

   Source      :  David Banks, Elaine Barker, James Dray, Allen Heckert,
                  Stefan Leigh, Mark Levenson, James Nechvatal, Andrew Rukhin,
                  Miles Smid, Juan Soto, Mark Vangel, and San Vo.

   Technical
   Assistance  :  Larry Bassham, Ron Boisvert, James Filliben, Daniel Lozier,
                  and Bert Rust.

   Warning     :  Portability Issues.

   Limitation  :  Amount of memory allocated for workspace.

   Restrictions:  Permission to use, copy, and modify this software without
                  fee is hereby granted, provided that this entire notice is
                  included in all copies of any software which is or includes
                  a copy or modification of this software and in all copies
                  of the supporting documentation for such software.
   -------------------------------------------------------------------------- */

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
                              M A C R O S
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#include <cmath>
#include <string>
#include "cephes.h"
#include <UnderKit/BitData.h>

#include <array>

#define MAX(x,y)             ((x) <  (y)  ? (y)  : (x))
#define MIN(x,y)             ((x) >  (y)  ? (y)  : (x))
#define isNonPositive(x)     ((x) <= 0.e0 ?   1  : 0)
#define isPositive(x)        ((x) >  0.e0 ?   1 : 0)
#define isNegative(x)        ((x) <  0.e0 ?   1 : 0)
#define isGreaterThanOne(x)  ((x) >  1.e0 ?   1 : 0)
#define isZero(x)            ((x) == 0.e0 ?   1 : 0)
#define isOne(x)             ((x) == 1.e0 ?   1 : 0)

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
                         G L O B A L  C O N S T A N T S
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#define ALPHA							0.01	/* SIGNIFICANCE LEVEL */
#define MAXNUMOFTEMPLATES				148		/* APERIODIC TEMPLATES: 148=>temp_length=9 */
#define BOOL bool

//==============================================================================

#ifdef __dll__		//defined in stdafx.h
#define IMPEXP __declspec(dllexport)
#else 
//#define IMPEXP __declspec(dllimport)
#endif

struct TEST_PAR {
    const char  *Name;
    uint32_t	MinBitsNum;		//min number of bits in seq for test correct run
    uint32_t	DefBlockLen;	//default block length
    uint32_t	SubTestsNum;	//num of subtests in certain test
};
//const TEST_PAR TestPar[];

//const char Template9[MAXNUMOFTEMPLATES][9];

// ============= CNistTests class =======================
//works with CByteArray - one random bit in one byte
class CNistTests
{
public:
    enum TESTS {
        FREQUENCY	= 0
        , BLOCK_FREQUENCY
        , CUSUM
        , CUSUM_REVERSE
        , RUNS
        , LONGEST_RUN
        , RANK
        , FFT
        , NONPERIODIC
        , OVERLAPPING
        , UNIVERSAL
        , APP_ENTROPY
        , RND_EXCURSION
        , RND_EXCURSION_VAR
        , SERIAL
        , LINEARCOMPLEXITY
        , NUMOFTESTS = LINEARCOMPLEXITY+1
    };

public:
    CNistTests(uint32_t nBytes, const char *pData);
    /*virtual*/ ~CNistTests();

    int RunAllForShortKey();

    BOOL Frequency(double &p_value);
    BOOL BlockFrequency(double &p_value, int BlockLength = 128);
    BOOL CumulativeSums(double &p_value);
    BOOL CumulativeSumsReverse(double &p_value);
    BOOL Runs(double &p_value);
    BOOL LongestRunOfOnes(double &p_value);
    BOOL Rank(double &p_value);
    BOOL DiscreteFourierTransform(double &p_value);
    BOOL NonOverlappingTemplateMatchings(double &p_value, int BlockLength = 9);
    BOOL OverlappingTemplateMatchings(double &p_value, int BlockLength = 9);
    BOOL Universal(double &p_value);
    BOOL ApproximateEntropy(double &p_value, int BlockLength = 10);
    BOOL RandomExcursions(double &p_value);
    BOOL RandomExcursionsVariant(double &p_value);
    BOOL Serial(double &p_value, int BlockLength = 16);
    BOOL LinearComplexity(double &p_value1, int BlockLength = 500);

    std::string GetDllVersion();

private:
    double psi2(int m, int n);	//used in Serial

    const uint32_t	m_nBytes;
    const char	*m_pData;
};

// ============= CNistTests2 class =======================
//1) works with CBitData - one random bit in one bit of byte array
//2) calculates individual p_value & test result for all subtests

class CBitData;

class CNistTests2
{
public:
    enum TESTS {
        FREQUENCY	= 0
        , BLOCK_FREQUENCY
        , CUSUM
        , CUSUM_REVERSE
        , RUNS
        , LONGEST_RUN
        , RANK
        , FFT
        , NONPERIODIC
        , OVERLAPPING
        , UNIVERSAL
        , APP_ENTROPY
        , RND_EXCURSION
        , RND_EXCURSION_VAR
        , SERIAL
        , LINEARCOMPLEXITY
        , NUMOFTESTS = LINEARCOMPLEXITY+1
    };

    struct TEST_RESULT {
        BOOL  res;
        double  p_value;
    };

    struct ALL_TESTS_RESULT {
        BOOL   result;
        double  p_valueAv;
        std::array<double, NUMOFTESTS>    failed;
        double  time_s;
    };

public:
    CNistTests2(CBitData *pB);
    /*virtual*/ ~CNistTests2();

    void	GetConfidenceInterval(int sampleSize, double &threshold_min, double &threshold_max);
    ALL_TESTS_RESULT RunAll();

    //each test return num of (sub)tests passed without errors
    int Frequency(TEST_RESULT *res);
    int BlockFrequency(TEST_RESULT *res, int BlockLength = 128);
    int CumulativeSums(TEST_RESULT *res);
    int CumulativeSumsReverse(TEST_RESULT *res);
    int Runs(TEST_RESULT *res);
    int LongestRunOfOnes(TEST_RESULT *res);
    int Rank(TEST_RESULT *res);
    int DiscreteFourierTransform(TEST_RESULT *res);
    int NonOverlappingTemplateMatchings(TEST_RESULT *res, int BlockLength = 9);
    int OverlappingTemplateMatchings(TEST_RESULT *res, int BlockLength = 9);
    int Universal(TEST_RESULT *res);
    int ApproximateEntropy(TEST_RESULT *res, int BlockLength = 10);
    int RandomExcursions(TEST_RESULT *res);
    int RandomExcursionsVariant(TEST_RESULT *res);
    int Serial(TEST_RESULT *res, int BlockLength = 16);
    int LinearComplexity(TEST_RESULT *res, int BlockLength = 500);

    std::string GetDllVersion();

private:
    double psi2(int m, int n);	//used in Serial

    uint32_t	nBits;
    CBitData	&Bit;
};
