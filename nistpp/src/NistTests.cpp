
#include <cfloat>
#include <memory.h>
#include "nist.h"
#include "cephes.h"
#include "dfft.h"
#include "matrix.h"
#include "template9.h"	//for NonOverlappingTemplate test

#include <stdio.h>
#include <time.h>

#include <UnderKit/ctrlexception.h>

//TODO обработка исключений

// Num of Tests = CNistTests::NUMOFTESTS
const TEST_PAR TestPar[] = {
    //	Name				MinBitsNum	DefBlockLen  SubTestsNum
    {"Frequency",					100,	0,			1,}
    , {"BlockFrequency",			100,	128,		1,}
    , {"CumulativeSums",			100,	0,			1,}
    , {"CumulativeSumsReverse",		100,	0,			1,}
    , {"Runs",						100,	0,			1,}
    , {"LongestRun",				128,	0,			1,}
    , {"Rank",						40000,	0,			1,}
    , {"FFT",						1000,	0,			1,}
    , {"NonOverlappingTemplate",	100,	9,	MAXNUMOFTEMPLATES,}	//At the moment works with block length = 9 only
    , {"OverlappingTemplate",		100,	9,			1,}
    , {"Universal",					162000, 0,			1,}
    , {"ApproximateEntropy",		10,		10,			1,}
    , {"RandomExcursions",			100,	0,			8,}
    , {"RandomExcursionsVariant",	100,	0,			18,}
    , {"Serial",					10,		16,			2,}
    , {"LinearComplexity",			1000000, 500,		1}
};

CNistTests::CNistTests(const uint32_t nBytes, const char *pData) :
    m_nBytes(nBytes)
  , m_pData(pData)
{
}

CNistTests::~CNistTests()
= default;

int CNistTests::RunAllForShortKey()
{
    BOOL res = false;
    double p_value;
    int i;
    try {
        for( i = 0; i < NUMOFTESTS; i++ ) {
            if( m_nBytes <= TestPar[i].MinBitsNum )
                continue;
#ifdef	SELECTED
            if(i == NONPERIODIC || i == OVERLAPPING || i == RND_EXCURSION || i == RND_EXCURSION_VAR || i == APP_ENTROPY) //not good for short key
                continue;
#endif
            switch(i) {
            case FREQUENCY :
                res = Frequency(p_value);
                break;
            case BLOCK_FREQUENCY :
                res = BlockFrequency(p_value, TestPar[i].DefBlockLen);
                break;
            case CUSUM :
                res = CumulativeSums(p_value);
                break;
            case CUSUM_REVERSE :
                res = CumulativeSumsReverse(p_value);
                break;
            case RUNS :
                res = Runs(p_value);
                break;
            case LONGEST_RUN :
                res = LongestRunOfOnes(p_value);
                break;
            case RANK :
                res = Rank(p_value);
                break;
            case FFT :
                res = DiscreteFourierTransform(p_value);
                break;
            case NONPERIODIC :
                res = NonOverlappingTemplateMatchings(p_value, TestPar[i].DefBlockLen);
                break;
            case OVERLAPPING :
                res = OverlappingTemplateMatchings(p_value, TestPar[i].DefBlockLen);
                break;
            case UNIVERSAL :
                res = Universal(p_value);
                break;
            case APP_ENTROPY :
                res = ApproximateEntropy(p_value, TestPar[i].DefBlockLen);
                break;
            case RND_EXCURSION :
                res = RandomExcursions(p_value);
                break;
            case RND_EXCURSION_VAR :
                res = RandomExcursionsVariant(p_value);
                break;
            case SERIAL :
                res = Serial(p_value, TestPar[i].DefBlockLen);
                break;
            case LINEARCOMPLEXITY :
                res = LinearComplexity(p_value, TestPar[i].DefBlockLen);
                break;
            default: break;
            }
            if(!res)
                return i;
        }
    } catch(CCtrlException *Cx) {
        res = false;
        return i;
    }

    return -1;
}

BOOL CNistTests::Frequency(double &p_value)	//m_nBytes should be > 100
{
    if( m_nBytes < TestPar[FREQUENCY].MinBitsNum )
    {
        p_value = 0.0;
        char buff[100];
        std::snprintf(buff, sizeof (buff), " test! Nbits < %d!", TestPar[FREQUENCY].MinBitsNum);
        throw  CCtrlException( TestPar[FREQUENCY].Name, buff );
    }

    double	 f;
    double	 s_obs;
    double	 sum;
    double	 sqrt2 = 1.41421356237309504880;

    sum = 0.0;
    for( uint32_t i = 0; i < m_nBytes; i++ )
        sum += 2.0 * (double)m_pData[i] - 1;
    s_obs = fabs(sum) / sqrt( (double)m_nBytes );
    f = s_obs / sqrt2;
    p_value = cephes_erfc(f);
    //COMPUTATIONAL INFORMATION:
    //The m_nBytes'th partial sum = (int)sum
    //Sum/n
    return p_value >= ALPHA;
}

BOOL CNistTests::BlockFrequency(double &p_value, int BlockLength)	//m_nBytes should be > 100
{
    if( m_nBytes < TestPar[BLOCK_FREQUENCY].MinBitsNum )
    {
        p_value = 0.0;
        char buff[100];
        std::snprintf(buff, sizeof (buff), " test! Nbits < %d!", TestPar[BLOCK_FREQUENCY].MinBitsNum);
        throw  CCtrlException(   TestPar[BLOCK_FREQUENCY].Name, buff );
    }

    int		 i;
    int		 j;
    int		 N;
    int		 blockSum;
    double	 sum;
    double	 pi;
    double	 v;
    double	 chi_squared;

    N = m_nBytes / BlockLength; 		/* # OF SUBSTRING BLOCKS  */
    sum = 0.0;

    for( i = 0; i < N; i++ ) {
        blockSum = 0;
        for( j = 0; j < BlockLength; j++ )
            blockSum += m_pData[j+i*BlockLength];
        pi = (double)blockSum/(double)BlockLength;
        v = pi - 0.5;
        sum += v*v;
    }
    chi_squared = 4.0 * BlockLength * sum;
    p_value = cephes_igamc(N/2.0, chi_squared/2.0);
    //COMPUTATIONAL INFORMATION:
    //Chi^2 = chi_squared
    //# of substrings = N
    //block length = BlockLength
    //Note: (m_nBytes % BlockLength) bits were discarded
    return p_value >= ALPHA;
}

BOOL CNistTests::CumulativeSums(double &p_value)	//m_nBytes should be > 100
{
    if( m_nBytes < TestPar[CUSUM].MinBitsNum ) {
        p_value = 0.0;
        char buff[100];
        std::snprintf(buff, sizeof (buff), " test! Nbits < %d!", TestPar[CUSUM].MinBitsNum);
        throw  CCtrlException(   TestPar[CUSUM].Name, buff );
    }

    int		 S;
    int		 sup;
    int		 inf;
    int		 z;
    int		 /*zrev,*/ k;
    double	 sum1;
    double	 sum2;

    S = 0;
    sup = 0;
    inf = 0;
    for( uint32_t i = 0; i < m_nBytes; i++ ) {
        m_pData[i] ? S++ : S--;
        if ( S > sup )
            sup++;
        if ( S < inf )
            inf--;
        z = (sup > -inf) ? sup : -inf;
        //zrev = (sup-S > S-inf) ? sup-S : S-inf;
    }

    // forward
    sum1 = 0.0;
    for( k = (-(int)m_nBytes/z+1)/4; k <= ((int)m_nBytes/z-1)/4; k++ ) {
        sum1 += cephes_normal(((4*k+1)*z)/sqrt((double)m_nBytes));
        sum1 -= cephes_normal(((4*k-1)*z)/sqrt((double)m_nBytes));
    }
    sum2 = 0.0;
    for( k = (-(int)m_nBytes/z-3)/4; k <= ((int)m_nBytes/z-1)/4; k++ ) {
        sum2 += cephes_normal(((4*k+3)*z)/sqrt((double)m_nBytes));
        sum2 -= cephes_normal(((4*k+1)*z)/sqrt((double)m_nBytes));
    }
    p_value = 1.0 - sum1 + sum2;
    //COMPUTATIONAL INFORMATION:
    //The maximum partial sum = z
    if( isNegative(p_value) || isGreaterThanOne(p_value) )
    {
        //fprintf(stats[TEST_CUSUM], "\t\tWARNING:  P_VALUE IS OUT OF RANGE\n");
        char buff[100];
        std::snprintf(buff, sizeof (buff), " test! p_value = %f!", p_value);
        throw  CCtrlException(   TestPar[CUSUM].Name, buff );
    }

    return p_value >= ALPHA;
}

BOOL CNistTests::CumulativeSumsReverse(double &p_value)	//m_nBytes should be > 100
{
    if( m_nBytes < TestPar[CUSUM_REVERSE].MinBitsNum )
    {
        p_value = 0.0;
        char buff[100];
        std::snprintf(buff, sizeof (buff), " test! Nbits < %d!", TestPar[CUSUM_REVERSE].MinBitsNum);
        throw  CCtrlException(   TestPar[CUSUM_REVERSE].Name, buff );
    }

    int		 S;
    int		 sup;
    int		 inf;
    int		 /*z,*/ zrev;
    int		 k;
    double	 sum1;
    double	 sum2;

    S = 0;
    sup = 0;
    inf = 0;
    for( uint32_t i = 0; i < m_nBytes; i++ ) {
        m_pData[i] ? S++ : S--;
        if ( S > sup )
            sup++;
        if ( S < inf )
            inf--;
        //z = (sup > -inf) ? sup : -inf;
        zrev = (sup-S > S-inf) ? sup-S : S-inf;
    }

    // backwards
    sum1 = 0.0;
    for( k = (-(int)m_nBytes/zrev+1)/4; k <= ((int)m_nBytes/zrev-1)/4; k++ ) {
        sum1 += cephes_normal(((4*k+1)*zrev)/sqrt((double)m_nBytes));
        sum1 -= cephes_normal(((4*k-1)*zrev)/sqrt((double)m_nBytes));
    }
    sum2 = 0.0;
    for( k = (-(int)m_nBytes/zrev-3)/4; k <= ((int)m_nBytes/zrev-1)/4; k++ ) {
        sum2 += cephes_normal(((4*k+3)*zrev)/sqrt((double)m_nBytes));
        sum2 -= cephes_normal(((4*k+1)*zrev)/sqrt((double)m_nBytes));
    }
    p_value = 1.0 - sum1 + sum2;
    //COMPUTATIONAL INFORMATION:
    //The maximum partial sum = zrev
    if( isNegative(p_value) || isGreaterThanOne(p_value) )
    {
        //	fprintf(stats[TEST_CUSUM], "\t\tWARNING:  P_VALUE IS OUT OF RANGE\n");
        char buff[100];
        std::snprintf(buff, sizeof (buff), " test! p_value = %f!", p_value);
        throw  CCtrlException(   TestPar[CUSUM_REVERSE].Name, buff );
    }
    return p_value >= ALPHA;
}

BOOL CNistTests::Runs(double &p_value)	//m_nBytes should be > 100
{
    if( m_nBytes < TestPar[RUNS].MinBitsNum )
    {
        p_value = 0.0;
        char buff[100];
        std::snprintf(buff, sizeof (buff), " test! Nbits < %d!", TestPar[RUNS].MinBitsNum);
        throw  CCtrlException(   TestPar[RUNS].Name, buff );
    }

    int		S;
    uint32_t	k;
    double	 pi;
    double	 V;
    double	 erfc_arg;

    S = 0;
    for( k = 0; k < m_nBytes; k++ )
        if ( m_pData[k] )
            S++;
    pi = (double)S / (double)m_nBytes;

    if( fabs(pi - 0.5) > (2.0 / sqrt((double)m_nBytes)) )
    {
        //fprintf(stats[TEST_RUNS], "\t\tPI ESTIMATOR CRITERIA NOT MET! PI = %f\n", pi);
        p_value = 0.0;
        std::string str1;
        std::string str2;
        char buff1[100];
        char buff2[100];
        std::snprintf(buff1, sizeof (buff1), "in %s test!", TestPar[RUNS].Name);
        std::snprintf(buff2, sizeof (buff2), "Pi Estimator Criteria Not Met! Pi = %f", pi);
        throw  CCtrlException( buff1, buff2 );
    } 
        V = 1;
        for( k = 1; k < m_nBytes; k++ )
            if ( m_pData[k] != m_pData[k-1] )
                V++;

        erfc_arg = fabs(V - 2.0 * m_nBytes * pi * (1-pi)) / (2.0 * pi * (1-pi) * sqrt(2.0*m_nBytes));
        p_value = cephes_erfc(erfc_arg);
        //COMPUTATIONAL INFORMATION:
        //Pi = pi
        //V_n_obs (Total # of runs) = (int)V
        // fabs(V - 2.0 * m_nBytes * pi * (1-pi)) / (2.0 * pi * (1-pi) * sqrt(2.0*m_nBytes)) = erfc_arg
        if( isNegative(p_value) || isGreaterThanOne(p_value) ) {
            //	fprintf(stats[TEST_CUSUM], "\t\tWARNING:  P_VALUE IS OUT OF RANGE\n");
            char buff[100];
            std::snprintf(buff, sizeof (buff), " test! p_value = %f!", p_value);
            throw  CCtrlException(   TestPar[RUNS].Name, buff );
        }
    
    return p_value >= ALPHA;
}

BOOL CNistTests::LongestRunOfOnes(double &p_value)	//m_nBytes should be > 128
{
    if( m_nBytes < TestPar[LONGEST_RUN].MinBitsNum )
    {
        p_value = 0.0;
        char buff[100];
        std::snprintf(buff, sizeof (buff), " test! Nbits < %d!", TestPar[LONGEST_RUN].MinBitsNum);
        throw  CCtrlException(   TestPar[LONGEST_RUN].Name, buff );
    }

    double			 chi2;
    double			 pi[7];
    int				 run;
    int				 v_n_obs;
    int				 N;
    int				 i;
    int				 j;
    int				 K;
    int				 M;
    int				 V[7];
    unsigned int	nu[7] = { 0, 0, 0, 0, 0, 0, 0 };

    if ( m_nBytes < 6272 ) {
        K = 3;
        M = 8;
        V[0] = 1; V[1] = 2; V[2] = 3; V[3] = 4;
        pi[0] = 0.21484375;
        pi[1] = 0.3671875;
        pi[2] = 0.23046875;
        pi[3] = 0.1875;
    }
    else if ( m_nBytes < 750000 ) {
        K = 5;
        M = 128;
        V[0] = 4; V[1] = 5; V[2] = 6; V[3] = 7; V[4] = 8; V[5] = 9;
        pi[0] = 0.1174035788;
        pi[1] = 0.242955959;
        pi[2] = 0.249363483;
        pi[3] = 0.17517706;
        pi[4] = 0.102701071;
        pi[5] = 0.112398847;
    }
    else {
        K = 6;
        M = 10000;
        V[0] = 10; V[1] = 11; V[2] = 12; V[3] = 13; V[4] = 14; V[5] = 15; V[6] = 16;
        pi[0] = 0.0882;
        pi[1] = 0.2092;
        pi[2] = 0.2483;
        pi[3] = 0.1933;
        pi[4] = 0.1208;
        pi[5] = 0.0675;
        pi[6] = 0.0727;
    }

    N = m_nBytes / M;
    for( i = 0; i < N; i++ ) {
        v_n_obs = 0;
        run = 0;
        for( j = 0; j < M; j++ ) {
            if( m_pData[i*M+j] == 1 ) {
                run++;
                if( run > v_n_obs )
                    v_n_obs = run;
            }
            else
                run = 0;
        }
        if ( v_n_obs < V[0] )
            nu[0]++;
        for( j = 0; j <= K; j++ ) {
            if( v_n_obs == V[j] )
                nu[j]++;
        }
        if( v_n_obs > V[K] )
            nu[K]++;
    }

    chi2 = 0.0;
    for( i = 0; i <= K; i++ )
        chi2 += ((nu[i] - N * pi[i]) * (nu[i] - N * pi[i])) / (N * pi[i]);

    p_value = cephes_igamc((double)(K/2.0), chi2/2.0);
    //COMPUTATIONAL INFORMATION:
    //# of substrings = N
    //Substring Length = M
    //Chi^2 = chi2
    //fprintf(stats[TEST_LONGEST_RUN], "\t\t---------------------------------------------\n");
    //fprintf(stats[TEST_LONGEST_RUN], "\t\t      F R E Q U E N C Y\n");
    //fprintf(stats[TEST_LONGEST_RUN], "\t\t---------------------------------------------\n");

    //if ( K == 3 ) {
    //	fprintf(stats[TEST_LONGEST_RUN], "\t\t  <=1     2     3    >=4   P-value  Assignment");
    //	fprintf(stats[TEST_LONGEST_RUN], "\n\t\t %3d %3d %3d  %3d ", nu[0], nu[1], nu[2], nu[3]);
    //}
    //else if ( K == 5 ) {
    //	fprintf(stats[TEST_LONGEST_RUN], "\t\t<=4  5  6  7  8  >=9 P-value  Assignment");
    //	fprintf(stats[TEST_LONGEST_RUN], "\n\t\t %3d %3d %3d %3d %3d  %3d ", nu[0], nu[1], nu[2],
    //			nu[3], nu[4], nu[5]);
    //}
    //else {
    //	fprintf(stats[TEST_LONGEST_RUN],"\t\t<=10  11  12  13  14  15 >=16 P-value  Assignment");
    //	fprintf(stats[TEST_LONGEST_RUN],"\n\t\t %3d %3d %3d %3d %3d %3d  %3d ", nu[0], nu[1], nu[2],
    //			nu[3], nu[4], nu[5], nu[6]);
    //}
    if ( isNegative(p_value) || isGreaterThanOne(p_value) )
    {
        char buff[100];
        std::snprintf(buff, sizeof (buff), " test! p_value = %f!",  p_value);
        throw  CCtrlException(   TestPar[LONGEST_RUN].Name, buff );
    }

    return p_value >= ALPHA;
}

BOOL CNistTests::Rank(double &p_value)	//m_nBytes should be > 40000
{
    if( m_nBytes < TestPar[RANK].MinBitsNum )
    {
        p_value = 0.0;
        char buff[100];
        std::snprintf(buff, sizeof (buff), " test! Nbits < %d!", TestPar[RANK].MinBitsNum);
        throw  CCtrlException(   TestPar[RANK].Name, buff );
    }

    int		 N;
    int		 i;
    int		 j;
    int		 k;
    int		 r;
    double	 product;
    double	 chi_squared;
    double	 arg1;
    double	 p_32;
    double	 p_31;
    double	 p_30;
    double	 R;
    double	 F_32;
    double	 F_31;
    double	 F_30;
    char	**matrix = new char* [32];	//matrix 32x32
    for( i = 0; i < 32; i++ )
        matrix[i] = new char[32];

    N = m_nBytes/(32*32);
    if( isZero(N) ) {
        //fprintf(stats[TEST_RANK], "\t\tError: Insuffucient # Of Bits To Define An 32x32 (%dx%d) Matrix\n", 32, 32);
        p_value = 0.00;
        return false;
    }

    r = 32;					/* COMPUTE PROBABILITIES */
    product = 1;
    for( i = 0; i <= r-1; i++ )
        product *= ((1.e0-pow(2.0, i-32))*(1.e0-pow(2.0, i-32)))/(1.e0-pow(2.0, i-r));
    p_32 = pow(2.0, r*(32+32-r)-32*32) * product;

    r = 31;
    product = 1;
    for( i = 0; i <= r-1; i++ )
        product *= ((1.e0-pow(2.0, i-32))*(1.e0-pow(2.0, i-32)))/(1.e0-pow(2.0, i-r));
    p_31 = pow(2.0, r*(32+32-r)-32*32) * product;

    p_30 = 1 - (p_32+p_31);

    F_32 = 0;
    F_31 = 0;
    for( k = 0; k < N; k++ ) {			/* FOR EACH 32x32 MATRIX   */
        for( i = 0; i < 32; i++ )
            for( j = 0; j < 32; j++ )
                matrix[i][j] = m_pData[k*(32*32)+j+i*32];
        //#if (DISPLAY_MATRICES == 1)
        //		display_matrix(32, 32, matrix);
        //#endif
        R = computeRank(32, 32, matrix);
        if ( R == 32 )
            F_32++;			/* DETERMINE FREQUENCIES */
        if ( R == 31 )
            F_31++;
    }
    F_30 = (double)N - (F_32+F_31);

    chi_squared =(pow(F_32 - N*p_32, 2)/(double)(N*p_32) +
                  pow(F_31 - N*p_31, 2)/(double)(N*p_31) +
                  pow(F_30 - N*p_30, 2)/(double)(N*p_30));

    arg1 = -chi_squared/2.e0;
    p_value = exp(arg1);

    //fprintf(stats[TEST_RANK], "\t\t\t\tRANK TEST\n");
    //fprintf(stats[TEST_RANK], "\t\t---------------------------------------------\n");
    //fprintf(stats[TEST_RANK], "\t\tCOMPUTATIONAL INFORMATION:\n");
    //fprintf(stats[TEST_RANK], "\t\t---------------------------------------------\n");
    //fprintf(stats[TEST_RANK], "\t\t(a) Probability P_%d = %f\n", 32,p_32);
    //fprintf(stats[TEST_RANK], "\t\t(b)             P_%d = %f\n", 31,p_31);
    //fprintf(stats[TEST_RANK], "\t\t(c)             P_%d = %f\n", 30,p_30);
    //fprintf(stats[TEST_RANK], "\t\t(d) Frequency   F_%d = %d\n", 32,(int)F_32);
    //fprintf(stats[TEST_RANK], "\t\t(e)             F_%d = %d\n", 31,(int)F_31);
    //fprintf(stats[TEST_RANK], "\t\t(f)             F_%d = %d\n", 30,(int)F_30);
    //fprintf(stats[TEST_RANK], "\t\t(g) # of matrices    = %d\n", N);
    //fprintf(stats[TEST_RANK], "\t\t(h) Chi^2            = %f\n", chi_squared);
    //fprintf(stats[TEST_RANK], "\t\t(i) NOTE: %d BITS WERE DISCARDED.\n", n%(32*32));
    //fprintf(stats[TEST_RANK], "\t\t---------------------------------------------\n");

    for ( i = 0; i < 32; i++ )	/* DEALLOCATE MATRIX  */
        delete [] matrix[i];
    delete [] matrix;

    if ( isNegative(p_value) || isGreaterThanOne(p_value) )
    {
        //	fprintf(stats[TEST_RANK], "WARNING:  P_VALUE IS OUT OF RANGE.\n");
        char buff[100];
        std::snprintf(buff, sizeof (buff), " test! p_value = %f!", p_value);
        throw  CCtrlException(   TestPar[RANK].Name, buff );
    }

    return p_value >= ALPHA;
}

BOOL CNistTests::DiscreteFourierTransform(double &p_value)	//m_nBytes should be > 1000
{
    if( m_nBytes < TestPar[FFT].MinBitsNum )
    {
        p_value = 0.0;
        char buff[100];
        std::snprintf(buff, sizeof (buff), " test! Nbits < %d!", TestPar[FFT].MinBitsNum);
        throw  CCtrlException(   TestPar[FFT].Name, buff );
    }

    double	 upperBound;
    double	 percentile;
    double	 N_l;
    double	 N_o;
    double	 d;
    int		 count;
    int		 ifac[15];
    uint32_t	i;

    auto *X = (double*)calloc(m_nBytes, sizeof(double));
    auto *wsave = (double*)calloc(2*m_nBytes, sizeof(double));
    auto *m = (double*)calloc(m_nBytes/2+1, sizeof(double));
    if( X == nullptr || wsave == nullptr || m == nullptr )
    {
        free(X);
        free(wsave);
        free(m);
        p_value = 0.0;
        char buff[100];
        std::snprintf(buff, sizeof (buff), "%s test!", TestPar[FFT].Name);
        throw  CCtrlException(   nullptr, buff );
    }

    for( i = 0; i < m_nBytes; i++ )
        X[i] = 2*(int)m_pData[i] - 1;

    __ogg_fdrffti(m_nBytes, wsave, ifac);		/* INITIALIZE WORK ARRAYS */
    __ogg_fdrfftf(m_nBytes, X, wsave, ifac);	/* APPLY FORWARD FFT */

    m[0] = sqrt(X[0]*X[0]);	    /* COMPUTE MAGNITUDE */

    for( i = 0; i < m_nBytes/2; i++ )
        m[i+1] = sqrt(pow(X[2*i+1],2)+pow(X[2*i+2],2));
    count = 0;				       /* CONFIDENCE INTERVAL */
    upperBound = sqrt(2.995732274*m_nBytes);
    for( i = 0; i < m_nBytes/2; i++ )
        if ( m[i] < upperBound )
            count++;
    percentile = (double)count/(m_nBytes/2)*100;
    N_l = (double) count;       /* number of peaks less than h = sqrt(3*n) */
    N_o = (double) 0.95*m_nBytes/2.0;
    d = (N_l - N_o)/sqrt(m_nBytes/4.0*0.95*0.05);
    p_value = cephes_erfc(fabs(d)/sqrt(2.0));	//cephes_

    //fprintf(stats[TEST_FFT], "\t\t\t\tFFT TEST\n");
    //fprintf(stats[TEST_FFT], "\t\t-------------------------------------------\n");
    //fprintf(stats[TEST_FFT], "\t\tCOMPUTATIONAL INFORMATION:\n");
    //fprintf(stats[TEST_FFT], "\t\t-------------------------------------------\n");
    //fprintf(stats[TEST_FFT], "\t\t(a) Percentile = %f\n", percentile);
    //fprintf(stats[TEST_FFT], "\t\t(b) N_l        = %f\n", N_l);
    //fprintf(stats[TEST_FFT], "\t\t(c) N_o        = %f\n", N_o);
    //fprintf(stats[TEST_FFT], "\t\t(d) d          = %f\n", d);
    //fprintf(stats[TEST_FFT], "\t\t-------------------------------------------\n");

    free(X);
    free(wsave);
    free(m);

    return p_value >= ALPHA;
}

BOOL CNistTests::NonOverlappingTemplateMatchings(double &p_value, int BlockLength)
{
    if( m_nBytes < TestPar[NONPERIODIC].MinBitsNum )
    {
        p_value = 0.0;
        char buff[100];
        std::snprintf(buff, sizeof (buff), " test! Nbits < %d!", TestPar[NONPERIODIC].MinBitsNum);
        throw  CCtrlException(   TestPar[NONPERIODIC].Name, buff );
    }

    int		numOfTemplates[100] = {0, 0, 2, 4, 6, 12, 20, 40, 74, 148, 284, 568, 1116,
                                       2232, 4424, 8848, 17622, 35244, 70340, 140680, 281076, 562152};
    /*----------------------------------------------------------------------------
        NOTE:  Should additional templates lengths beyond 21 be desired, they must
        first be constructed, saved into files and then the corresponding
        number of nonperiodic templates for that file be stored in the BlockLength-th
        position in the numOfTemplates variable.
        ----------------------------------------------------------------------------*/
    //templates 17 - 21 are in binary format (".Z" extension): other reading needed...
    //if( BlockLength > 16 ) {
    //	AfxMessageBox("NONOVERLAPPING TEMPLATES TEST ABORTED!\nSelected BlockLength > 16. No correspondant template.");
    //	p_value = 0.0;
    //	return false;
    //}

    //std::string filepath;
    //GetAppDataFolder(filepath);
    //filepath.AppendFormat("template%d", BlockLength);
    //CStdioFile file;
    //if( !file.Open( filepath, CFile::modeRead | CFile::typeText ) ) {
    //	AfxMessageBox("NONOVERLAPPING TEMPLATES TEST ABORTED!\nRequired template file not existing.");
    //	p_value = 0.0;
    //	return false;
    //}

    BlockLength = 9;	//current version works with default template9 only!

    unsigned int	/*bit,*/ W_obs;
    unsigned int	/*bit,*/ nu[6];
    //FILE			*fp;
    //char			directory[100];
    double			 sum;
    double			 chi2;
    double			 lambda;
    double			 pi[6];
    double			 varWj;
    int				 i;
    int				 j;
    int				 jj;
    int				 k;
    int				 match;
    int				 SKIP;
    int				 M;
    int				 N;
    int				 K = 5;
    double			p_value_min = DBL_MAX;

    N = 8;
    M = m_nBytes/N;

    lambda = (M-BlockLength+1)/pow(2.0, BlockLength);
    if( isNegative(lambda) || isZero(lambda) ) {
        //AfxMessageBox("NONOVERLAPPING TEMPLATES TESTS ABORTED DUE TO: Lambda not being positive!");
        p_value = 0.0;
        throw  CCtrlException(   TestPar[NONPERIODIC].Name, " test! Lambda <= 0!" );
    }

    //if ( (Wj = (unsigned int*)calloc(N, sizeof(unsigned int))) == nullptr ) {
    //	fprintf(stats[TEST_NONPERIODIC], "\tNONOVERLAPPING TEMPLATES TESTS ABORTED DUE TO ONE OF THE FOLLOWING : \n");
    //	fprintf(stats[TEST_NONPERIODIC], "\tInsufficient memory for required work space.\n");
    //	return;
    //}

    varWj = M*(1.0/pow(2.0, BlockLength) - (2.0*BlockLength-1.0)/pow(2.0, 2.0*BlockLength));
    auto *Wj = new uint32_t[N];
    memset(Wj, 0, sizeof(uint32_t)*N);

    //char *sequence = new char[BlockLength];

    //sprintf(directory, "templates/template%d", BlockLength);
    //if ( ((isNegative(lambda)) || (isZero(lambda))) ||
    //	 ((fp = fopen(directory, "r")) == nullptr) ||
    //	 ((sequence = (BitSequence *) calloc(m, sizeof(BitSequence))) == nullptr) ) {
    //	fprintf(stats[TEST_NONPERIODIC], "\tNONOVERLAPPING TEMPLATES TESTS ABORTED DUE TO ONE OF THE FOLLOWING : \n");
    //	fprintf(stats[TEST_NONPERIODIC], "\tLambda (%f) not being positive!\n", lambda);
    //	fprintf(stats[TEST_NONPERIODIC], "\tTemplate file <%s> not existing\n", directory);
    //	fprintf(stats[TEST_NONPERIODIC], "\tInsufficient memory for required work space.\n");
    //	if ( sequence != nullptr )
    //		free(sequence);
    //}
    //else {
    //	fprintf(stats[TEST_NONPERIODIC], "\t\t  NONPERIODIC TEMPLATES TEST\n");
    //	fprintf(stats[TEST_NONPERIODIC], "-------------------------------------------------------------------------------------\n");
    //	fprintf(stats[TEST_NONPERIODIC], "\t\t  COMPUTATIONAL INFORMATION\n");
    //	fprintf(stats[TEST_NONPERIODIC], "-------------------------------------------------------------------------------------\n");
    //	fprintf(stats[TEST_NONPERIODIC], "\tLAMBDA = %f\tM = %d\tN = %d\tm = %d\tn = %d\n", lambda, M, N, m, n);
    //	fprintf(stats[TEST_NONPERIODIC], "-------------------------------------------------------------------------------------\n");
    //	fprintf(stats[TEST_NONPERIODIC], "\t\tF R E Q U E N C Y\n");
    //	fprintf(stats[TEST_NONPERIODIC], "Template   W_1  W_2  W_3  W_4  W_5  W_6  W_7  W_8    Chi^2   P_value Assignment Index\n");
    //	fprintf(stats[TEST_NONPERIODIC], "-------------------------------------------------------------------------------------\n");

    if( numOfTemplates[BlockLength] < MAXNUMOFTEMPLATES )
        SKIP = 1;
    else
        SKIP = (int)(numOfTemplates[BlockLength]/MAXNUMOFTEMPLATES);
    numOfTemplates[BlockLength] = (int)numOfTemplates[BlockLength]/SKIP;

    sum = 0.0;
    for( i = 0; i < 2; i++ ) {                      /* Compute Probabilities */
        pi[i] = exp(-lambda+i*log(lambda)-cephes_lgam(i+1));
        sum += pi[i];
    }
    pi[0] = sum;
    for( i = 2; i <= K; i++ ) {                      /* Compute Probabilities */
        pi[i-1] = exp(-lambda+i*log(lambda)-cephes_lgam(i+1));
        sum += pi[i-1];
    }
    pi[K] = 1 - sum;

    for( jj = 0; jj < MIN(MAXNUMOFTEMPLATES, numOfTemplates[BlockLength]); jj++ ) {
        sum = 0;

        ////for( k = 0; k < BlockLength; k++ ) {
        ////	fscanf(fp, "%d", &bit);
        ////	sequence[k] = bit;
        ////	fprintf(stats[TEST_NONPERIODIC], "%d", sequence[k]);
        ////}
        //std::string line;
        //if( !file.ReadString( line ) ) {
        //	AfxMessageBox("NONOVERLAPPING TEMPLATES TESTS ABORTED DUE TO: Unexpected end of template file!");
        //	delete [] sequence;
        //	delete [] Wj;
        //	p_value = 0.0;
        //	return false;
        //}
        //int kk = 0;
        //for( k = 0; k < line.GetLength(); k++ ) {
        //	char ch = line.GetAt(k);
        //	if(ch == '0' || ch == '1') {
        //		sequence[kk++] = ch - '0';
        //		if(kk > BlockLength) {
        //			AfxMessageBox("NONOVERLAPPING TEMPLATES TESTS ABORTED DUE TO: Line in template file > BlockLength!");
        //			delete [] sequence;
        //			delete [] Wj;
        //			p_value = 0.0;
        //			return false;
        //		}
        //		//fprintf(stats[TEST_NONPERIODIC], "%d", sequence[k]);
        //	}
        //}

        //fprintf(stats[TEST_NONPERIODIC], " ");
        for( k = 0; k <= K; k++ )
            nu[k] = 0;
        for( i = 0; i < N; i++ ) {
            W_obs = 0;
            for( j = 0; j < M-BlockLength+1; j++ ) {
                match = 1;
                for( k = 0; k < BlockLength; k++ ) {
                    //if( (int)sequence[k] != (int)m_pData[i*M+j+k] ) {
                    if( (int)Template9[jj][k] != (int)m_pData[i*M+j+k] ) {
                        match = 0;
                        break;
                    }
                }
                if( match == 1 )
                    W_obs++;
            }
            Wj[i] = W_obs;
        }
        sum = 0;
        chi2 = 0.0;                                   /* Compute Chi Square */
        for( i = 0; i < N; i++ ) {
            //if( m == 10 )
            //	fprintf(stats[TEST_NONPERIODIC], "%3d  ", Wj[i]);
            //else
            //	fprintf(stats[TEST_NONPERIODIC], "%4d ", Wj[i]);
            chi2 += pow(((double)Wj[i] - lambda)/pow(varWj, 0.5), 2);
        }
        p_value = cephes_igamc(N/2.0, chi2/2.0);

        if ( isNegative(p_value) || isGreaterThanOne(p_value) )
        {
            //			fprintf(stats[TEST_NONPERIODIC], "\t\tWARNING:  P_VALUE IS OUT OF RANGE.\n");
            char buff[100];
            std::snprintf(buff, sizeof (buff), " test! p_value[%d] = %f!", jj, p_value);
            throw  CCtrlException(   TestPar[NONPERIODIC].Name, buff );
        }

        //		fprintf(stats[TEST_NONPERIODIC], "%9.6f %f %s %3d\n", chi2, p_value, p_value < ALPHA ? "FAILURE" : "SUCCESS", jj);
        //		if ( SKIP > 1 )
        //			fseek(fp, (long)(SKIP-1)*2*m, SEEK_CUR);
        //		fprintf(results[TEST_NONPERIODIC], "%f\n", p_value); fflush(results[TEST_NONPERIODIC]);

        p_value_min = std::min(p_value, p_value_min);
    }
    //
    //fprintf(stats[TEST_NONPERIODIC], "\n"); fflush(stats[TEST_NONPERIODIC]);

    //fclose(fp);
    //file.Close();

    //if ( sequence != nullptr )
    //	free(sequence);
    //free(Wj);
    //delete [] sequence;
    delete [] Wj;

    p_value = p_value_min;
    return p_value >= ALPHA;
}

static double Pr(int u, double eta)	//used in OverlappingTemplateMatchings
{
    int		l;
    double	 sum;
    double	 p;

    if ( u == 0 )
        p = exp(-eta);
    else {
        sum = 0.0;
        for ( l=1; l<=u; l++ )
            sum += exp(-eta-u*log(2.0)+l*log(eta)-cephes_lgam(l+1)+cephes_lgam(u)-cephes_lgam(l)-cephes_lgam(u-l+1));
        p = sum;
    }
    return p;
}

BOOL CNistTests::OverlappingTemplateMatchings(double &p_value, int BlockLength)	//m = BlockLength, n = m_nBytes
{
    if( m_nBytes < TestPar[OVERLAPPING].MinBitsNum )
    {
        p_value = 0.0;
        char buff[100];
        std::snprintf(buff, sizeof (buff), " test! Nbits < %d!", TestPar[OVERLAPPING].MinBitsNum);
        throw  CCtrlException(   TestPar[OVERLAPPING].Name, buff );
    }

    int				 i;
    int				 k;
    int				 match;
    double			 W_obs;
    double			 eta;
    double			 sum;
    double			 chi2;
    double			 lambda;
    int				 M;
    int				 N;
    int				 j;
    int				 K = 5;
    unsigned int	nu[6] = { 0, 0, 0, 0, 0, 0 };
    double			pi[6] = { 0.143783, 0.139430, 0.137319, 0.124314, 0.106209, 0.348945 };

    M = 1032;
    N = m_nBytes/M;

    //if ( (sequence = (BitSequence *) calloc(m, sizeof(BitSequence))) == nullptr ) {
    //	fprintf(stats[TEST_OVERLAPPING], "\t\t    OVERLAPPING TEMPLATE OF ALL ONES TEST\n");
    //	fprintf(stats[TEST_OVERLAPPING], "\t\t---------------------------------------------\n");
    //	fprintf(stats[TEST_OVERLAPPING], "\t\tTEMPLATE DEFINITION:  Insufficient memory, Overlapping Template Matchings test aborted!\n");
    //}
    //else
    //	for ( i=0; i<m; i++ )
    //		sequence[i] = 1;
    char *sequence = new char[BlockLength];
    memset(sequence, 1, BlockLength);

    lambda = (double)(M-BlockLength+1)/pow(2.0,BlockLength);
    eta = lambda/2.0;
    sum = 0.0;
    for( i = 0; i < K; i++ ) {			/* Compute Probabilities */
        pi[i] = Pr(i, eta);
        sum += pi[i];
    }
    pi[K] = 1 - sum;

    for( i = 0; i < N; i++ ) {
        W_obs = 0;
        for( j = 0; j < M-BlockLength+1; j++ ) {
            match = 1;
            for( k = 0; k < BlockLength; k++ ) {
                if( sequence[k] != m_pData[i*M+j+k] )
                    match = 0;
            }
            if( match == 1 )
                W_obs++;
        }
        if( W_obs <= 4 )
            nu[(int)W_obs]++;
        else
            nu[K]++;
    }
    sum = 0;
    chi2 = 0.0;                                   /* Compute Chi Square */
    for( i = 0; i < K+1; i++ ) {
        chi2 += pow((double)nu[i] - (double)N*pi[i], 2)/((double)N*pi[i]);
        sum += nu[i];
    }
    p_value = cephes_igamc(K/2.0, chi2/2.0);

    //fprintf(stats[TEST_OVERLAPPING], "\t\t    OVERLAPPING TEMPLATE OF ALL ONES TEST\n");
    //fprintf(stats[TEST_OVERLAPPING], "\t\t-----------------------------------------------\n");
    //fprintf(stats[TEST_OVERLAPPING], "\t\tCOMPUTATIONAL INFORMATION:\n");
    //fprintf(stats[TEST_OVERLAPPING], "\t\t-----------------------------------------------\n");
    //fprintf(stats[TEST_OVERLAPPING], "\t\t(a) n (sequence_length)      = %d\n", n);
    //fprintf(stats[TEST_OVERLAPPING], "\t\t(b) m (block length of 1s)   = %d\n", m);
    //fprintf(stats[TEST_OVERLAPPING], "\t\t(c) M (length of substring)  = %d\n", M);
    //fprintf(stats[TEST_OVERLAPPING], "\t\t(d) N (number of substrings) = %d\n", N);
    //fprintf(stats[TEST_OVERLAPPING], "\t\t(e) lambda [(M-m+1)/2^m]     = %f\n", lambda);
    //fprintf(stats[TEST_OVERLAPPING], "\t\t(f) eta                      = %f\n", eta);
    //fprintf(stats[TEST_OVERLAPPING], "\t\t-----------------------------------------------\n");
    //fprintf(stats[TEST_OVERLAPPING], "\t\t   F R E Q U E N C Y\n");
    //fprintf(stats[TEST_OVERLAPPING], "\t\t  0   1   2   3   4 >=5   Chi^2   P-value  Assignment\n");
    //fprintf(stats[TEST_OVERLAPPING], "\t\t-----------------------------------------------\n");
    //fprintf(stats[TEST_OVERLAPPING], "\t\t%3d %3d %3d %3d %3d %3d  %f ",
    //	nu[0], nu[1], nu[2], nu[3], nu[4], nu[5], chi2);

    //free(sequence);
    delete [] sequence;

    if ( isNegative(p_value) || isGreaterThanOne(p_value) )
    {
        //	fprintf(stats[TEST_OVERLAPPING], "WARNING:  P_VALUE IS OUT OF RANGE.\n");
        char buff[100];
        std::snprintf(buff, sizeof (buff), " test! p_value = %f!", p_value);
        throw  CCtrlException(   TestPar[OVERLAPPING].Name, buff );
    }

    return p_value >= ALPHA;
}

BOOL CNistTests::Universal(double &p_value)	//m_nBytes should be > 162000
{
    if( m_nBytes < TestPar[UNIVERSAL].MinBitsNum )
    {
        p_value = 0.0;
        char buff[100];
        std::snprintf(buff, sizeof (buff), " test! Nbits < %d!", TestPar[UNIVERSAL].MinBitsNum);
        throw  CCtrlException(   TestPar[UNIVERSAL].Name, buff );
    }

    int		 i;
    int		 j;
    int		 p;
    int		 L;
    int		 Q;
    int		 K;
    double	 arg;
    double	 sqrt2;
    double	 sigma;
    double	 phi;
    double	 sum;
    double	 c;
    int64_t	decRep;
    double	expected_value[17] = { 0, 0, 0, 0, 0, 0, 5.2177052, 6.1962507, 7.1836656,
                                       8.1764248, 9.1723243, 10.170032, 11.168765,
                                       12.168070, 13.167693, 14.167488, 15.167379 };
    double   variance[17] = { 0, 0, 0, 0, 0, 0, 2.954, 3.125, 3.238, 3.311, 3.356, 3.384,
                              3.401, 3.410, 3.416, 3.419, 3.421 };

    /* * * * * * * * * ** * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
         * THE FOLLOWING REDEFINES L, SHOULD THE CONDITION:     n >= 1010*2^L*L       *
         * NOT BE MET, FOR THE BLOCK LENGTH L.                                        *
         * * * * * * * * * * ** * * * * * * * * * * * * * * * * * * * * * * * * * * * */
    L = 5;
    if ( m_nBytes >= 387840 )     L = 6;
    if ( m_nBytes >= 904960 )     L = 7;
    if ( m_nBytes >= 2068480 )    L = 8;
    if ( m_nBytes >= 4654080 )    L = 9;
    if ( m_nBytes >= 10342400 )   L = 10;
    if ( m_nBytes >= 22753280 )   L = 11;
    if ( m_nBytes >= 49643520 )   L = 12;
    if ( m_nBytes >= 107560960 )  L = 13;
    if ( m_nBytes >= 231669760 )  L = 14;
    if ( m_nBytes >= 496435200 )  L = 15;
    if ( m_nBytes >= 1059061760 ) L = 16;

    Q = 10*(int)pow(2.0, L);
    K = (int)(floor((double)(m_nBytes/L)) - (double)Q);	 		    /* BLOCKS TO TEST */

    p = (int)pow(2.0, L);
    //if( (L < 6) || (L > 16) || ((double)Q < 10*pow(2.0, L)) ||
    //	 ((T = (long *)calloc(p, sizeof(long))) == nullptr) ) {
    //	fprintf(stats[TEST_UNIVERSAL], "\t\tUNIVERSAL STATISTICAL TEST\n");
    //	fprintf(stats[TEST_UNIVERSAL], "\t\t---------------------------------------------\n");
    //	fprintf(stats[TEST_UNIVERSAL], "\t\tERROR:  L IS OUT OF RANGE.\n");
    //	fprintf(stats[TEST_UNIVERSAL], "\t\t-OR- :  Q IS LESS THAN %f.\n", 10*pow(2.0, L));
    //	fprintf(stats[TEST_UNIVERSAL], "\t\t-OR- :  Unable to allocate T.\n");
    //	return;
    //}
    //for( i = 0; i < p; i++ )
    //	T[i] = 0;
    auto *T = (int64_t*)calloc(p, sizeof(int64_t));
    if( T == nullptr )
    {
        free(T);
        p_value = 0.0;
        char buff[100];
        std::snprintf(buff, sizeof (buff), "%s test!", TestPar[UNIVERSAL].Name);
        throw  CCtrlException(   nullptr, buff );
    }
    //memset(T, 0, sizeof(int64_t)*p);

    /* COMPUTE THE EXPECTED:  Formula 16, in Marsaglia's Paper */
    c = 0.7 - 0.8/(double)L + (4 + 32/(double)L)*pow(K, -3/(double)L)/15;
    sigma = c * sqrt(variance[L]/(double)K);
    sqrt2 = sqrt(2.0);
    sum = 0.0;
    for( i = 1; i <= Q; i++ ) {		/* INITIALIZE TABLE */
        decRep = 0;
        for( j = 0; j < L; j++ )
            decRep += m_pData[(i-1)*L+j] * (int64_t)pow(2.0, L-1-j);
        T[decRep] = i;
    }
    for( i = Q+1; i <= Q+K; i++ ) { 	/* PROCESS BLOCKS */
        decRep = 0;
        for( j = 0; j < L; j++ )
            decRep += m_pData[(i-1)*L+j] * (int64_t)pow(2.0, L-1-j);
        sum += log((double)(i - T[decRep]))/log(2.0);
        T[decRep] = i;
    }
    phi = (double)(sum/(double)K);
    arg = fabs(phi-expected_value[L])/(sqrt2 * sigma);
    p_value = cephes_erfc(arg);	//cephes_

    //fprintf(stats[TEST_UNIVERSAL], "\t\tUNIVERSAL STATISTICAL TEST\n");
    //fprintf(stats[TEST_UNIVERSAL], "\t\t--------------------------------------------\n");
    //fprintf(stats[TEST_UNIVERSAL], "\t\tCOMPUTATIONAL INFORMATION:\n");
    //fprintf(stats[TEST_UNIVERSAL], "\t\t--------------------------------------------\n");
    //fprintf(stats[TEST_UNIVERSAL], "\t\t(a) L         = %d\n", L);
    //fprintf(stats[TEST_UNIVERSAL], "\t\t(b) Q         = %d\n", Q);
    //fprintf(stats[TEST_UNIVERSAL], "\t\t(c) K         = %d\n", K);
    //fprintf(stats[TEST_UNIVERSAL], "\t\t(d) sum       = %f\n", sum);
    //fprintf(stats[TEST_UNIVERSAL], "\t\t(e) sigma     = %f\n", sigma);
    //fprintf(stats[TEST_UNIVERSAL], "\t\t(f) variance  = %f\n", variance[L]);
    //fprintf(stats[TEST_UNIVERSAL], "\t\t(g) exp_value = %f\n", expected_value[L]);
    //fprintf(stats[TEST_UNIVERSAL], "\t\t(h) phi       = %f\n", phi);
    //fprintf(stats[TEST_UNIVERSAL], "\t\t(i) WARNING:  %d bits were discarded.\n", n-(Q+K)*L);
    //fprintf(stats[TEST_UNIVERSAL], "\t\t-----------------------------------------\n");

    free(T);

    if ( isNegative(p_value) || isGreaterThanOne(p_value) )
    {
        //	fprintf(stats[TEST_UNIVERSAL], "\t\tWARNING:  P_VALUE IS OUT OF RANGE\n");
        char buff[100];
        std::snprintf(buff, sizeof (buff), "test! p_value = %f!", p_value);
        throw  CCtrlException(   TestPar[UNIVERSAL].Name, buff );
    }

    return p_value >= ALPHA;
}

BOOL CNistTests::ApproximateEntropy(double &p_value, int BlockLength)	//m_nBytes should be > 10; BlockLength < (int)log2(m_nBytes) - 2
{
    if( m_nBytes < TestPar[APP_ENTROPY].MinBitsNum )
    {
        p_value = 0.0;
        char buff[100];
        std::snprintf(buff, sizeof (buff), " test! Nbits < %d!", TestPar[APP_ENTROPY].MinBitsNum);
        throw  CCtrlException(   TestPar[APP_ENTROPY].Name, buff );
    }
    //if ( BlockLength > (int)(log((double)seqLength)/log(2.0)-5) ) {
    //	fprintf(stats[TEST_APEN], "\t\tNote: The blockSize = %d exceeds recommended value of %d\n", m,
    //		MAX(1, (int)(log((double)seqLength)/log(2.0)-5)));
    //	fprintf(stats[TEST_APEN], "\t\tResults are inaccurate!\n");

    int		 i;
    int		 j;
    int		 k;
    int		 r;
    int		 blockSize;
    int		 seqLength;
    int		 powLen;
    int		 index;
    double	 sum;
    double	 numOfBlocks;
    double	 ApEn[2];
    double	 apen;
    double	 chi_squared;
    uint32_t	*P;

    seqLength = m_nBytes;
    r = 0;

    for( blockSize = BlockLength; blockSize <= BlockLength+1; blockSize++ ) {
        if( blockSize == 0 ) {
            ApEn[0] = 0.00;
            r++;
        } else {
            numOfBlocks = (double)seqLength;
            powLen = (int)pow(2.0, blockSize+1)-1;
            //if ( (P = (unsigned int*)calloc(powLen,sizeof(unsigned int)))== nullptr ) {
            //	fprintf(stats[TEST_APEN], "ApEn:  Insufficient memory available.\n");
            //	return;
            //}
            //for ( i=1; i<powLen-1; i++ )
            //	P[i] = 0;
            P = new uint32_t[powLen];
            memset(P, 0, sizeof(uint32_t)*powLen);
            for( i = 0; i < numOfBlocks; i++ ) { /* COMPUTE FREQUENCY */
                k = 1;
                for( j = 0; j < blockSize; j++ ) {
                    k <<= 1;
                    if( (int)m_pData[(i+j) % seqLength] == 1 )
                        k++;
                }
                P[k-1]++;
            }
            /* DISPLAY FREQUENCY */
            sum = 0.0;
            index = (int)pow(2.0, blockSize)-1;
            for( i = 0; i < (int)pow(2.0, blockSize); i++ ) {
                if ( P[index] > 0 )
                    sum += P[index]*log(P[index]/numOfBlocks);
                index++;
            }
            sum /= numOfBlocks;
            ApEn[r] = sum;
            r++;
            //free(P);
            delete [] P;
        }
    }
    apen = ApEn[0] - ApEn[1];
    chi_squared = 2.0*seqLength*(log(2.0) - apen);
    p_value = cephes_igamc(pow(2.0, BlockLength-1), chi_squared/2.0);

    //fprintf(stats[TEST_APEN], "\t\t(b) n (sequence length) = %d\n", seqLength);
    //fprintf(stats[TEST_APEN], "\t\t(c) Chi^2               = %f\n", chi_squared);
    //fprintf(stats[TEST_APEN], "\t\t(d) Phi(m)	       = %f\n", ApEn[0]);
    //fprintf(stats[TEST_APEN], "\t\t(e) Phi(m+1)	       = %f\n", ApEn[1]);
    //fprintf(stats[TEST_APEN], "\t\t(f) ApEn                = %f\n", apen);
    //fprintf(stats[TEST_APEN], "\t\t(g) Log(2)              = %f\n", log(2.0));
    //fprintf(stats[TEST_APEN], "\t\t--------------------------------------------\n");

    return p_value >= ALPHA;
}

BOOL CNistTests::RandomExcursions(double &p_value)
{
    if( m_nBytes < TestPar[RND_EXCURSION].MinBitsNum )
    {
        p_value = 0.0;
        char buff[100];
        std::snprintf(buff, sizeof (buff), " test! Nbits < %d!", TestPar[RND_EXCURSION].MinBitsNum);
        throw  CCtrlException(   TestPar[RND_EXCURSION].Name, buff );
    }

    int		 b;
    int		 i;
    int		 j;
    int		 k;
    int		 J;
    int		 x;
    int		 cycleStart;
    int		 cycleStop;
    int		stateX[8] = { -4, -3, -2, -1, 1, 2, 3, 4 };
    int		counter[8] = { 0, 0, 0, 0, 0, 0, 0, 0 };
    double	 sum;
    double	 constraint;
    double	 nu[6][8];
    double	pi[5][6] = { {0.0000000000, 0.00000000000, 0.00000000000, 0.00000000000, 0.00000000000, 0.0000000000},
                             {0.5000000000, 0.25000000000, 0.12500000000, 0.06250000000, 0.03125000000, 0.0312500000},
                             {0.7500000000, 0.06250000000, 0.04687500000, 0.03515625000, 0.02636718750, 0.0791015625},
                             {0.8333333333, 0.02777777778, 0.02314814815, 0.01929012346, 0.01607510288, 0.0803755143},
                             {0.8750000000, 0.01562500000, 0.01367187500, 0.01196289063, 0.01046752930, 0.0732727051} };
    double	p_value_min = DBL_MAX;

    int *S_k = (int*)calloc(m_nBytes, sizeof(int));
    int *cycle = (int*)calloc(MAX(1000, m_nBytes/100), sizeof(int));
    if( S_k == nullptr || cycle == nullptr )
    {
        free(S_k);
        free(cycle);
        p_value = 0.0;
        char buff[100];
        std::snprintf(buff, sizeof (buff), "%s test!", TestPar[RND_EXCURSION].Name);
        throw  CCtrlException(   nullptr, buff );
    }

    J = 0; 					/* DETERMINE CYCLES */
    S_k[0] = 2*(int)m_pData[0] - 1;
    for( i = 1; i < (int)m_nBytes; i++ ) {
        S_k[i] = S_k[i-1] + 2*m_pData[i] - 1;
        if( S_k[i] == 0 ) {
            J++;
            if ( J > (int)MAX(1000, m_nBytes/100) ) {
                free(S_k);
                free(cycle);
                p_value = 0.0;
                throw  CCtrlException(   TestPar[RND_EXCURSION].Name,
                                          " test! Exceeding the max number of cycles expected." );
            }
            cycle[J] = i;
        }
    }
    if( S_k[m_nBytes-1] != 0 )
        J++;
    cycle[J] = m_nBytes;

    //fprintf(stats[TEST_RND_EXCURSION], "\t\t\t  RANDOM EXCURSIONS TEST\n");
    //fprintf(stats[TEST_RND_EXCURSION], "\t\t--------------------------------------------\n");
    //fprintf(stats[TEST_RND_EXCURSION], "\t\tCOMPUTATIONAL INFORMATION:\n");
    //fprintf(stats[TEST_RND_EXCURSION], "\t\t--------------------------------------------\n");
    //fprintf(stats[TEST_RND_EXCURSION], "\t\t(a) Number Of Cycles (J) = %04d\n", J);
    //fprintf(stats[TEST_RND_EXCURSION], "\t\t(b) Sequence Length (n)  = %d\n", n);

    constraint = MAX(0.005*pow((double)m_nBytes, 0.5), 500);
    if (J < constraint) {
        //fprintf(stats[TEST_RND_EXCURSION], "\t\t---------------------------------------------\n");
        //fprintf(stats[TEST_RND_EXCURSION], "\t\tWARNING:  TEST NOT APPLICABLE.  THERE ARE AN\n");
        //fprintf(stats[TEST_RND_EXCURSION], "\t\t\t  INSUFFICIENT NUMBER OF CYCLES.\n");
        //fprintf(stats[TEST_RND_EXCURSION], "\t\t---------------------------------------------\n");
        //for(i = 0; i < 8; i++)
        //	fprintf(results[TEST_RND_EXCURSION], "%f\n", 0.0);
        free(S_k);
        free(cycle);
        p_value = 0.0;
        throw  CCtrlException(   TestPar[RND_EXCURSION].Name,
                                  " test! There are an insufficient number of cycles." );
    } 
        //fprintf(stats[TEST_RND_EXCURSION], "\t\t(c) Rejection Constraint = %f\n", constraint);
        //fprintf(stats[TEST_RND_EXCURSION], "\t\t-------------------------------------------\n");

        cycleStart = 0;
        cycleStop  = cycle[1];
        for( k = 0; k < 6; k++ )
            for( i = 0; i < 8; i++ )
                nu[k][i] = 0.;
        for( j = 1; j <= J; j++ ) {                           /* FOR EACH CYCLE */
            for( i = 0; i < 8; i++ )
                counter[i] = 0;
            for( i = cycleStart; i < cycleStop; i++ ) {
                if( (S_k[i] >= 1 && S_k[i] <= 4) || (S_k[i] >= -4 && S_k[i] <= -1) ) {
                    if( S_k[i] < 0 )
                        b = 4;
                    else
                        b = 3;
                    counter[S_k[i]+b]++;
                }
            }
            cycleStart = cycle[j]+1;
            if( j < J )
                cycleStop = cycle[j+1];

            for( i = 0; i < 8; i++ ) {
                if ( (counter[i] >= 0) && (counter[i] <= 4) )
                    nu[counter[i]][i]++;
                else if ( counter[i] >= 5 )
                    nu[5][i]++;
            }
        }

        free(S_k);
        free(cycle);

        for( i = 0; i < 8; i++ ) {
            x = stateX[i];
            sum = 0.;
            for( k = 0; k < 6; k++ )
                sum += pow(nu[k][i] - J*pi[(int)fabs((double)x)][k], 2) / (J*pi[(int)fabs((double)x)][k]);
            p_value = cephes_igamc(2.5, sum/2.0);

            if ( isNegative(p_value) || isGreaterThanOne(p_value) )
            {
                //	fprintf(stats[TEST_RND_EXCURSION], "WARNING:  P_VALUE IS OUT OF RANGE.\n");
                char buff[100];
                std::snprintf(buff, sizeof (buff), " test! p_value[%d] = %f!", i, p_value);
                throw  CCtrlException(   TestPar[NONPERIODIC].Name, buff );
            }

            //fprintf(stats[TEST_RND_EXCURSION], "%s\t\tx = %2d chi^2 = %9.6f p_value = %f\n",
            //		p_value < ALPHA ? "FAILURE" : "SUCCESS", x, sum, p_value);

            p_value_min = std::min(p_value, p_value_min);
        }
    

    p_value = p_value_min;
    return p_value >= ALPHA;
}

BOOL CNistTests::RandomExcursionsVariant(double &p_value)
{
    if( m_nBytes < TestPar[RND_EXCURSION_VAR].MinBitsNum )
    {
        p_value = 0.0;
        char buff[100];
        std::snprintf(buff, sizeof (buff), " test! Nbits < %d!", TestPar[RND_EXCURSION_VAR].MinBitsNum);
        throw  CCtrlException(   TestPar[RND_EXCURSION_VAR].Name, buff );
    }

    int		 i;
    int		 p;
    int		 J;
    int		 x;
    int		 constraint;
    int		 count;
    int		stateX[18] = { -9, -8, -7, -6, -5, -4, -3, -2, -1, 1, 2, 3, 4, 5, 6, 7, 8, 9 };
    double	p_value_min = DBL_MAX;

    int *S_k = (int*)calloc(m_nBytes, sizeof(int));
    if( S_k == nullptr )
    {
        free(S_k);
        p_value = 0.0;
        char buff[100];
        std::snprintf(buff, sizeof (buff), "%s test!", TestPar[RND_EXCURSION_VAR].Name);
        throw  CCtrlException(   nullptr, buff );
    }

    J = 0;
    S_k[0] = 2*(int)m_pData[0] - 1;
    for( i = 1; i < (int)m_nBytes; i++ ) {
        S_k[i] = S_k[i-1] + 2*m_pData[i] - 1;
        if ( S_k[i] == 0 )
            J++;
    }
    if ( S_k[m_nBytes-1] != 0 )
        J++;

    //fprintf(stats[TEST_RND_EXCURSION_VAR], "\t\t\tRANDOM EXCURSIONS VARIANT TEST\n");
    //fprintf(stats[TEST_RND_EXCURSION_VAR], "\t\t--------------------------------------------\n");
    //fprintf(stats[TEST_RND_EXCURSION_VAR], "\t\tCOMPUTATIONAL INFORMATION:\n");
    //fprintf(stats[TEST_RND_EXCURSION_VAR], "\t\t--------------------------------------------\n");
    //fprintf(stats[TEST_RND_EXCURSION_VAR], "\t\t(a) Number Of Cycles (J) = %d\n", J);
    //fprintf(stats[TEST_RND_EXCURSION_VAR], "\t\t(b) Sequence Length (n)  = %d\n", n);
    //fprintf(stats[TEST_RND_EXCURSION_VAR], "\t\t--------------------------------------------\n");

    constraint = (int)MAX(0.005*pow(m_nBytes, 0.5), 500);
    if (J < constraint) {
        //fprintf(stats[TEST_RND_EXCURSION_VAR], "\n\t\tWARNING:  TEST NOT APPLICABLE.  THERE ARE AN\n");
        //fprintf(stats[TEST_RND_EXCURSION_VAR], "\t\t\t  INSUFFICIENT NUMBER OF CYCLES.\n");
        //fprintf(stats[TEST_RND_EXCURSION_VAR], "\t\t---------------------------------------------\n");
        //for ( i=0; i<18; i++ )
        //	fprintf(results[TEST_RND_EXCURSION_VAR], "%f\n", 0.0);
        free(S_k);
        p_value = 0.0;
        throw  CCtrlException(   TestPar[RND_EXCURSION_VAR].Name,
                                  " test! There are an insufficient number of cycles." );
    } 
        for( p = 0; p <= 17; p++ ) {
            x = stateX[p];
            count = 0;
            for( i = 0;  i < (int)m_nBytes; i++ )
                if( S_k[i] == x )
                    count++;
            p_value = cephes_erfc(fabs((double)(count-J))/(sqrt(2.0*J*(4.0*fabs((double)x)-2))));	//cephes_

            if ( isNegative(p_value) || isGreaterThanOne(p_value) )
            {
                //	fprintf(stats[TEST_RND_EXCURSION_VAR], "\t\t(b) WARNING: P_VALUE IS OUT OF RANGE.\n");
                free(S_k);
                char buff[100];
                std::snprintf(buff, sizeof (buff), " test! p_value[%d] = %f!", p, p_value);
                throw  CCtrlException(   TestPar[RND_EXCURSION_VAR].Name, buff );
            }
            //fprintf(stats[TEST_RND_EXCURSION_VAR], "%s\t\t", p_value < ALPHA ? "FAILURE" : "SUCCESS");
            //fprintf(stats[TEST_RND_EXCURSION_VAR], "(x = %2d) Total visits = %4d; p-value = %f\n", x, count, p_value);

            p_value_min = std::min(p_value, p_value_min);
        }
    
    free(S_k);

    p_value = p_value_min;
    return p_value >= ALPHA;
}

double CNistTests::psi2(int m, int n)	//used in Serial
{
    int				 i;
    int				 j;
    int				 k;
    int				 powLen;
    double			 sum;
    double			 numOfBlocks;

    if ( (m == 0) || (m == -1) )
        return 0.0;
    numOfBlocks = n;
    powLen = (int)pow(2.0, m+1)-1;

    auto *P = (uint32_t*)calloc(powLen, sizeof(uint32_t));
    if( P == nullptr )
    {
        free(P);
        char buff[100];
        std::snprintf(buff, sizeof (buff), "%s test in psi2 func!", TestPar[SERIAL].Name);
        throw  CCtrlException(   nullptr, buff );
    }
    //for ( i=1; i<powLen-1; i++ )
    //	P[i] = 0;	  /* INITIALIZE NODES */
    memset(P, 0, sizeof(uint32_t)*powLen);

    for( i = 0; i < numOfBlocks; i++ ) {		 /* COMPUTE FREQUENCY */
        k = 1;
        for( j = 0; j < m; j++ ) {
            if( m_pData[(i+j)%n] == 0 )
                k *= 2;
            else if ( m_pData[(i+j)%n] == 1 )
                k = 2*k+1;
        }
        P[k-1]++;
    }
    sum = 0.0;
    for( i = (int)pow(2.0, m)-1; i < (int)pow(2.0, m+1)-1; i++ )
        sum += pow(P[i], 2.0);
    sum = (sum * pow(2.0, m)/(double)n) - (double)n;

    free(P);

    return sum;
}

BOOL CNistTests::Serial(double &p_value, int BlockLength)	//m_nBytes should be > 10; BlockLength < (int)log2(m_nBytes) - 2
{
    if( m_nBytes < TestPar[SERIAL].MinBitsNum )
    {
        p_value = 0.0;
        char buff[100];
        std::snprintf(buff, sizeof (buff), " test! Nbits < %d!", TestPar[SERIAL].MinBitsNum);
        throw  CCtrlException(   TestPar[SERIAL].Name, buff );
    }

    double	 p_value1;
    double	 p_value2;
    double	 psim0;
    double	 psim1;
    double	 psim2;
    double	 del1;
    double	 del2;

    psim0 = psi2(BlockLength, m_nBytes);
    psim1 = psi2(BlockLength-1, m_nBytes);
    psim2 = psi2(BlockLength-2, m_nBytes);
    del1 = psim0 - psim1;
    del2 = psim0 - 2.0*psim1 + psim2;
    p_value1 = cephes_igamc(pow(2.0, BlockLength-1)/2, del1/2.0);
    p_value2 = cephes_igamc(pow(2.0, BlockLength-2)/2, del2/2.0);

    //fprintf(stats[TEST_SERIAL], "\t\t\t       SERIAL TEST\n");
    //fprintf(stats[TEST_SERIAL], "\t\t---------------------------------------------\n");
    //fprintf(stats[TEST_SERIAL], "\t\t COMPUTATIONAL INFORMATION:		  \n");
    //fprintf(stats[TEST_SERIAL], "\t\t---------------------------------------------\n");
    //fprintf(stats[TEST_SERIAL], "\t\t(a) Block length    (m) = %d\n", m);
    //fprintf(stats[TEST_SERIAL], "\t\t(b) Sequence length (n) = %d\n", n);
    //fprintf(stats[TEST_SERIAL], "\t\t(c) Psi_m               = %f\n", psim0);
    //fprintf(stats[TEST_SERIAL], "\t\t(d) Psi_m-1             = %f\n", psim1);
    //fprintf(stats[TEST_SERIAL], "\t\t(e) Psi_m-2             = %f\n", psim2);
    //fprintf(stats[TEST_SERIAL], "\t\t(f) Del_1               = %f\n", del1);
    //fprintf(stats[TEST_SERIAL], "\t\t(g) Del_2               = %f\n", del2);
    //fprintf(stats[TEST_SERIAL], "\t\t---------------------------------------------\n");

    //fprintf(stats[TEST_SERIAL], "%s\t\tp_value1 = %f\n", p_value1 < ALPHA ? "FAILURE" : "SUCCESS", p_value1);
    //fprintf(results[TEST_SERIAL], "%f\n", p_value1);

    //fprintf(stats[TEST_SERIAL], "%s\t\tp_value2 = %f\n\n", p_value2 < ALPHA ? "FAILURE" : "SUCCESS", p_value2);
    //fprintf(results[TEST_SERIAL], "%f\n", p_value2);

    p_value = std::min(p_value1, p_value2);

    return p_value >= ALPHA;
}

BOOL CNistTests::LinearComplexity(double &p_value, int BlockLength)		//m_nBytes should be > 1000000
{
    if( m_nBytes < TestPar[LINEARCOMPLEXITY].MinBitsNum )
    {
        p_value = 0.0;
        char buff[100];
        std::snprintf(buff, sizeof (buff), " test! Nbits < %d!", TestPar[LINEARCOMPLEXITY].MinBitsNum);
        throw  CCtrlException(   TestPar[LINEARCOMPLEXITY].Name, buff );
    }

    int       i;
    int       ii;
    int       j;
    int       d;
    int       N;
    int       L;
    int       m;
    int       N_;
    int       parity;
    int       sign;
    int       K = 6;
    double    T_;
    double    mean;
    double    nu[7];
    double    chi2;
    double    pi[7] = { 0.01047, 0.03125, 0.12500, 0.50000, 0.25000, 0.06250, 0.020833 };

    N = (int)floor((double)(m_nBytes/BlockLength));

    char *B_ = (char*)calloc(BlockLength, sizeof(char));
    char *C  = (char*)calloc(BlockLength, sizeof(char));
    char *P  = (char*)calloc(BlockLength, sizeof(char));
    char *T  = (char*)calloc(BlockLength, sizeof(char));
    if( B_ == nullptr || C ==nullptr || P == nullptr || T == nullptr )
    {
        free(B_);
        free(P);
        free(C);
        free(T);
        char buff[100];
        std::snprintf(buff, sizeof (buff), "%s test!", TestPar[LINEARCOMPLEXITY].Name);
        throw  CCtrlException(   nullptr, buff );
    }

    //fprintf(stats[TEST_LINEARCOMPLEXITY], "-----------------------------------------------------\n");
    //fprintf(stats[TEST_LINEARCOMPLEXITY], "\tL I N E A R  C O M P L E X I T Y\n");
    //fprintf(stats[TEST_LINEARCOMPLEXITY], "-----------------------------------------------------\n");
    //fprintf(stats[TEST_LINEARCOMPLEXITY], "\tM (substring length)     = %d\n", M);
    //fprintf(stats[TEST_LINEARCOMPLEXITY], "\tN (number of substrings) = %d\n", N);
    //fprintf(stats[TEST_LINEARCOMPLEXITY], "-----------------------------------------------------\n");
    //fprintf(stats[TEST_LINEARCOMPLEXITY], "        F R E Q U E N C Y                            \n");
    //fprintf(stats[TEST_LINEARCOMPLEXITY], "-----------------------------------------------------\n");
    //fprintf(stats[TEST_LINEARCOMPLEXITY], "  C0   C1   C2   C3   C4   C5   C6    CHI2    P-value\n");
    //fprintf(stats[TEST_LINEARCOMPLEXITY], "-----------------------------------------------------\n");
    //fprintf(stats[TEST_LINEARCOMPLEXITY], "\tNote: %d bits were discarded!\n", n%M);

    for( i = 0; i < K+1; i++ )
        nu[i] = 0.00;
    for( ii = 0; ii < N; ii++ ) {
        for( i = 0; i < BlockLength; i++ ) {
            B_[i] = 0;
            C[i] = 0;
            T[i] = 0;
            P[i] = 0;
        }
        L = 0;
        m = -1;
        d = 0;
        C[0] = 1;
        B_[0] = 1;

        /* DETERMINE LINEAR COMPLEXITY */
        N_ = 0;
        while( N_ < BlockLength ) {
            d = (int)m_pData[ii*BlockLength+N_];
            for( i = 1; i <= L; i++ )
                d += C[i] * m_pData[ii*BlockLength+N_-i];
            d = d%2;
            if( d == 1 ) {
                for( i = 0; i < BlockLength; i++ ) {
                    T[i] = C[i];
                    P[i] = 0;
                }
                for( j = 0; j < BlockLength; j++ )
                    if( B_[j] == 1 )
                        P[j+N_-m] = 1;
                for( i = 0; i < BlockLength; i++ )
                    C[i] = (C[i] + P[i])%2;
                if( L <= N_/2 ) {
                    L = N_ + 1 - L;
                    m = N_;
                    for( i = 0; i < BlockLength; i++ )
                        B_[i] = T[i];
                }
            }
            N_++;
        }
        if( (parity = (BlockLength+1)%2) == 0 )
            sign = -1;
        else
            sign = 1;
        mean = BlockLength/2.0 + (9.0+sign)/36.0 - 1.0/pow(2.0, BlockLength) * (BlockLength/3.0 + 2.0/9.0);
        if( (parity = BlockLength%2) == 0 )
            sign = 1;
        else
            sign = -1;
        T_ = sign * (L - mean) + 2.0/9.0;

        if( T_ <= -2.5 )
            nu[0]++;
        else if ( T_ > -2.5 && T_ <= -1.5 )
            nu[1]++;
        else if ( T_ > -1.5 && T_ <= -0.5 )
            nu[2]++;
        else if ( T_ > -0.5 && T_ <= 0.5 )
            nu[3]++;
        else if ( T_ > 0.5 && T_ <= 1.5 )
            nu[4]++;
        else if ( T_ > 1.5 && T_ <= 2.5 )
            nu[5]++;
        else
            nu[6]++;
    }
    chi2 = 0.00;
    //for( i = 0; i < K+1; i++ )
    //	fprintf(stats[TEST_LINEARCOMPLEXITY], "%4d ", (int)nu[i]);
    for( i = 0; i < K+1; i++ )
        chi2 += pow(nu[i]-N*pi[i], 2) / (N*pi[i]);
    p_value = cephes_igamc(K/2.0, chi2/2.0);

    //fprintf(stats[TEST_LINEARCOMPLEXITY], "%9.6f%9.6f\n", chi2, p_value);
    //fprintf(results[TEST_LINEARCOMPLEXITY], "%f\n", p_value);

    free(B_);
    free(P);
    free(C);
    free(T);

    return p_value >= ALPHA;
}

std::string CNistTests::GetDllVersion()
{
    std::string	str = "Unavailable";
    //	char filepath[_MAX_PATH];
    //	HMODULE hModule = GetModuleHandle("Nist.dll");
    //	GetModuleFileName( hModule, filepath, _MAX_PATH );
    //	DWORD	size, buf;
    //	if( size = GetFileVersionInfoSize( filepath, &buf ) ) {
    //		void*	pData = malloc(size);
    //		if( GetFileVersionInfo( filepath, buf, size, pData ) ) {
    //			VS_FIXEDFILEINFO vsf;
    //			void*	pbuf;
    //            uint32_t	len;
    //			if( VerQueryValue( pData, "\\", &pbuf, &len ) ) {
    //				memcpy( &vsf, pbuf, sizeof(VS_FIXEDFILEINFO) );
    //				str.Format("%d.%d.%d.%d",
    //					(WORD)(vsf.dwProductVersionMS >> 16), (WORD)(vsf.dwProductVersionMS),
    //					(WORD)(vsf.dwProductVersionLS >> 16), (WORD)(vsf.dwProductVersionLS) );
    //			}
    // 		}
    //		free(pData);
    //	}
    return str;
}

//============================================================================================= CNistTest2 ===========================
CNistTests2::CNistTests2(CBitData *pB) :
    Bit(*pB)
  //, nBits(Bit.Nbits)
{
    //  Bit = *pB;
    nBits = Bit.Nbits;
}

CNistTests2::~CNistTests2()
= default;

void CNistTests2::GetConfidenceInterval(int sampleSize, double &threshold_min, double &threshold_max)
{
    double p_hat = 1.0 - ALPHA;
    threshold_max = (p_hat + 3.0 * sqrt((p_hat*ALPHA)/sampleSize));
    threshold_min = (p_hat - 3.0 * sqrt((p_hat*ALPHA)/sampleSize));
}

CNistTests2::ALL_TESTS_RESULT CNistTests2::RunAll()
{
    clock_t t0 = clock();
    clock_t t1;
    ALL_TESTS_RESULT result{};
    memset(&result, 0, sizeof(ALL_TESTS_RESULT));
    TEST_RESULT *res = nullptr;
    int		 nTests = 0;
    int		 nPassed = 0;
    int		 nFailed = 0;
    double	p_valAv = 0;
    int		 i;
    int		 j;
    try {
        for( i = 0; i < NUMOFTESTS; i++ ) {
            if( nBits <= TestPar[i].MinBitsNum )
                continue;
#ifdef	SELECTED
            if(i == NONPERIODIC || i == OVERLAPPING || i == RND_EXCURSION || i == RND_EXCURSION_VAR || i == APP_ENTROPY) //not good for short key
                continue;
#endif
            int nDone = 0;
            delete [] res;
            res = new TEST_RESULT[TestPar[i].SubTestsNum];
            memset(res, 0, sizeof(TEST_RESULT)*TestPar[i].SubTestsNum);
            switch(i) {
            case FREQUENCY :
                nDone = Frequency(res);
                break;
            case BLOCK_FREQUENCY :
                nDone = BlockFrequency(res, TestPar[i].DefBlockLen);
                break;
            case CUSUM :
                nDone = CumulativeSums(res);
                break;
            case CUSUM_REVERSE :
                nDone = CumulativeSumsReverse(res);
                break;
            case RUNS :
                nDone = Runs(res);
                break;
            case LONGEST_RUN :
                nDone = LongestRunOfOnes(res);
                break;
            case RANK :
                nDone = Rank(res);
                break;
            case FFT :
                nDone = DiscreteFourierTransform(res);
                break;
            case NONPERIODIC :
                nDone = NonOverlappingTemplateMatchings(res, TestPar[i].DefBlockLen);
                break;
            case OVERLAPPING :
                nDone = OverlappingTemplateMatchings(res, TestPar[i].DefBlockLen);
                break;
            case UNIVERSAL :
                nDone = Universal(res);
                break;
            case APP_ENTROPY :
                nDone = ApproximateEntropy(res, TestPar[i].DefBlockLen);
                break;
            case RND_EXCURSION :
                nDone = RandomExcursions(res);
                break;
            case RND_EXCURSION_VAR :
                nDone = RandomExcursionsVariant(res);
                break;
            case SERIAL :
                nDone = Serial(res, TestPar[i].DefBlockLen);
                break;
            case LINEARCOMPLEXITY :
                nDone = LinearComplexity(res, TestPar[i].DefBlockLen);
                break;
            default: break;
            }
            if(nDone == 0)
                result.failed[i] = -1;  //test was not applicable
            nTests += nDone;
            for( j = 0; j < nDone; j++ ) {
                if( res[j].res )
                    nPassed++;
                else {
                    nFailed++;
                    result.failed[i]++;
                }
                p_valAv += res[j].p_value;
            }
        }
        delete [] res;
    } catch(CCtrlException *Cx) {
        delete [] res;
    }
    p_valAv /= nTests;
    result.p_valueAv = p_valAv;

    double max;
    double min;
    GetConfidenceInterval(nTests, min, max);
    result.result = nPassed / (double)nTests >= min;
    t1 = clock();
    result.time_s = (double)(t1-t0) / CLOCKS_PER_SEC;

    return result;
}

int CNistTests2::Frequency(TEST_RESULT *res)	//nBites should be > 100
{
    if( nBits < TestPar[FREQUENCY].MinBitsNum || TestPar[FREQUENCY].SubTestsNum != 1 )
    {
        char buff[100];
        std::snprintf(buff, sizeof (buff), " test! Nbits < %d!", TestPar[FREQUENCY].MinBitsNum);
        throw  CCtrlException(   TestPar[FREQUENCY].Name, buff );
    }

    double	 f;
    double	 s_obs;
    double	 sum;
    double	 sqrt2 = 1.4142135623730950488016887242097;

    sum = 0.0;
    uint32_t i;
    for( i = 0; i < nBits; i++ )
        sum += 2.0 * (double)Bit[i] - 1;
    s_obs = fabs(sum) / sqrt( (double)nBits );
    f = s_obs / sqrt2;
    res[0].p_value = cephes_erfc(f);
    //COMPUTATIONAL INFORMATION:
    //The m_nBytes'th partial sum = (int)sum
    //Sum/n
    res[0].res = res[0].p_value >= ALPHA;
    return 1;
}

int CNistTests2::BlockFrequency(TEST_RESULT *res, int BlockLength)	//nBits should be > 100
{
    if( nBits < TestPar[BLOCK_FREQUENCY].MinBitsNum || TestPar[BLOCK_FREQUENCY].SubTestsNum != 1 )
    {
        char buff[100];
        std::snprintf(buff, sizeof (buff), " test! Nbits < %d!", TestPar[BLOCK_FREQUENCY].MinBitsNum);
        throw  CCtrlException(   TestPar[BLOCK_FREQUENCY].Name, buff );
    }

    uint32_t	 i;
    uint32_t	 j;
    uint32_t	 N;
    uint32_t	 blockSum;
    double	 sum;
    double	 pi;
    double	 v;
    double	 chi_squared;

    N = (uint32_t)(nBits / BlockLength); 		/* # OF SUBSTRING BLOCKS  */
    sum = 0.0;

    for( i = 0; i < N; i++ ) {
        blockSum = 0;
        for( j = 0; j < (uint32_t)BlockLength; j++ )
            blockSum += Bit[j+i*BlockLength];
        pi = (double)blockSum/(double)BlockLength;
        v = pi - 0.5;
        sum += v*v;
    }
    chi_squared = 4.0 * BlockLength * sum;
    res[0].p_value = cephes_igamc(N/2.0, chi_squared/2.0);
    //COMPUTATIONAL INFORMATION:
    //Chi^2 = chi_squared
    //# of substrings = N
    //block length = BlockLength
    //Note: (m_nBytes % BlockLength) bits were discarded
    res[0].res = res[0].p_value >= ALPHA;
    return 1;
}

int CNistTests2::CumulativeSums(TEST_RESULT *res)	//m_nBytes should be > 100
{
    if( nBits < TestPar[CUSUM].MinBitsNum || TestPar[CUSUM].SubTestsNum != 1 )
    {
        char buff[100];
        std::snprintf(buff, sizeof (buff), " test! Nbits < %d!", TestPar[CUSUM].MinBitsNum);
        throw  CCtrlException(   TestPar[CUSUM].Name, buff );
    }

    int		 S;
    int		 sup;
    int		 inf;
    int		 z;
    int		 /*zrev,*/ k;
    double	 sum1;
    double	 sum2;

    S = 0;
    sup = 0;
    inf = 0;
    for( uint32_t i = 0; i < nBits; i++ ) {
        Bit[i] ? S++ : S--;
        if ( S > sup )
            sup++;
        if ( S < inf )
            inf--;
        z = (sup > -inf) ? sup : -inf;
        //zrev = (sup-S > S-inf) ? sup-S : S-inf;
    }

    // forward
    sum1 = 0.0;
    for( k = (-(int)nBits/z+1)/4; k <= ((int)nBits/z-1)/4; k++ ) {
        sum1 += cephes_normal(((4*k+1)*z)/sqrt((double)nBits));
        sum1 -= cephes_normal(((4*k-1)*z)/sqrt((double)nBits));
    }
    sum2 = 0.0;
    for( k = (-(int)nBits/z-3)/4; k <= ((int)nBits/z-1)/4; k++ ) {
        sum2 += cephes_normal(((4*k+3)*z)/sqrt((double)nBits));
        sum2 -= cephes_normal(((4*k+1)*z)/sqrt((double)nBits));
    }
    res[0].p_value = 1.0 - sum1 + sum2;
    //COMPUTATIONAL INFORMATION:
    //The maximum partial sum = z
    if( isNegative(res[0].p_value) || isGreaterThanOne(res[0].p_value) )
    {
        char buff[100];
        std::snprintf(buff, sizeof (buff), " test! p_value = %f!", res[0].p_value);
        throw  CCtrlException(   TestPar[CUSUM].Name, buff );
    }

    res[0].res = res[0].p_value >= ALPHA;
    return 1;
}

int CNistTests2::CumulativeSumsReverse(TEST_RESULT *res)	//nBits should be > 100
{
    if( nBits < TestPar[CUSUM_REVERSE].MinBitsNum || TestPar[CUSUM_REVERSE].SubTestsNum != 1 )
    {
        char buff[100];
        std::snprintf(buff, sizeof (buff), " test! Nbits < %d!", TestPar[CUSUM_REVERSE].MinBitsNum);
        throw  CCtrlException(   TestPar[CUSUM_REVERSE].Name, buff );
    }

    int		 S;
    int		 sup;
    int		 inf;
    int		 /*z,*/ zrev;
    int		 k;
    double	 sum1;
    double	 sum2;

    S = 0;
    sup = 0;
    inf = 0;
    for( uint32_t i = 0; i < nBits; i++ ) {
        Bit[i] ? S++ : S--;
        if ( S > sup )
            sup++;
        if ( S < inf )
            inf--;
        //z = (sup > -inf) ? sup : -inf;
        zrev = (sup-S > S-inf) ? sup-S : S-inf;
    }

    // backwards
    sum1 = 0.0;
    for( k = (-(int)nBits/zrev+1)/4; k <= ((int)nBits/zrev-1)/4; k++ ) {
        sum1 += cephes_normal(((4*k+1)*zrev)/sqrt((double)nBits));
        sum1 -= cephes_normal(((4*k-1)*zrev)/sqrt((double)nBits));
    }
    sum2 = 0.0;
    for( k = (-(int)nBits/zrev-3)/4; k <= ((int)nBits/zrev-1)/4; k++ ) {
        sum2 += cephes_normal(((4*k+3)*zrev)/sqrt((double)nBits));
        sum2 -= cephes_normal(((4*k+1)*zrev)/sqrt((double)nBits));
    }
    res[0].p_value = 1.0 - sum1 + sum2;
    //COMPUTATIONAL INFORMATION:
    //The maximum partial sum = zrev
    if( isNegative(res[0].p_value) || isGreaterThanOne(res[0].p_value) )
    {
        char buff[100];
        std::snprintf(buff, sizeof (buff), " test! p_value = %f!", res[0].p_value);
        throw  CCtrlException(   TestPar[CUSUM_REVERSE].Name, buff );
    }
    res[0].res = res[0].p_value >= ALPHA;
    return 1;
}

int CNistTests2::Runs(TEST_RESULT *res)	//nBits should be > 100
{
    if( nBits < TestPar[RUNS].MinBitsNum || TestPar[RUNS].SubTestsNum != 1 )
    {
        char buff[100];
        std::snprintf(buff, sizeof (buff), " test! Nbits < %d!", TestPar[RUNS].MinBitsNum);
        throw  CCtrlException(   TestPar[RUNS].Name, buff );
    }

    int		S;
    uint32_t	k;
    double	 pi;
    double	 V;
    double	 erfc_arg;

    S = 0;
    for( k = 0; k < nBits; k++ )
        if ( Bit[k] )
            S++;
    pi = (double)S / (double)nBits;

    if( fabs(pi - 0.5) > (2.0 / sqrt((double)nBits)) )
    {
        std::string str1;
        std::string str2;
        char buff1[100];
        char buff2[100];
        std::snprintf(buff1, sizeof (buff1), "in %s test!", TestPar[RUNS].Name);
        std::snprintf(buff2, sizeof (buff2), "Pi Estimator Criteria Not Met! Pi = %F", pi);
        throw  CCtrlException( buff1, buff2 );
    } 
        V = 1;
        for( k = 1; k < nBits; k++ )
            if ( Bit[k] != Bit[k-1] )
                V++;

        erfc_arg = fabs(V - 2.0 * nBits * pi * (1-pi)) / (2.0 * pi * (1-pi) * sqrt(2.0*nBits));
        res[0].p_value = cephes_erfc(erfc_arg);
        //COMPUTATIONAL INFORMATION:
        //Pi = pi
        //V_n_obs (Total # of runs) = (int)V
        // fabs(V - 2.0 * nBits * pi * (1-pi)) / (2.0 * pi * (1-pi) * sqrt(2.0*nBits)) = erfc_arg
        if( isNegative(res[0].p_value) || isGreaterThanOne(res[0].p_value) )
        {
            //	fprintf(stats[TEST_CUSUM], "\t\tWARNING:  P_VALUE IS OUT OF RANGE\n");
            char buff[100];
            std::snprintf(buff, sizeof (buff), " test! p_value = %f!", res[0].p_value);
            throw  CCtrlException(   TestPar[RUNS].Name, buff );
        }
    
    res[0].res = res[0].p_value >= ALPHA;
    return 1;
}

int CNistTests2::LongestRunOfOnes(TEST_RESULT *res)	//nBits should be > 128
{
    if( nBits < TestPar[LONGEST_RUN].MinBitsNum || TestPar[LONGEST_RUN].SubTestsNum != 1 )
    {
        char buff[100];
        std::snprintf(buff, sizeof (buff), " test! Nbits < %d!", TestPar[LONGEST_RUN].MinBitsNum);
        throw  CCtrlException(   TestPar[LONGEST_RUN].Name, buff );
    }

    double			 chi2;
    double			 pi[7];
    int				 run;
    int				 v_n_obs;
    int				 N;
    int				 i;
    int				 j;
    int				 K;
    int				 M;
    int				 V[7];
    unsigned int	nu[7] = { 0, 0, 0, 0, 0, 0, 0 };

    if ( nBits < 6272 ) {
        K = 3;
        M = 8;
        V[0] = 1; V[1] = 2; V[2] = 3; V[3] = 4;
        pi[0] = 0.21484375;
        pi[1] = 0.3671875;
        pi[2] = 0.23046875;
        pi[3] = 0.1875;
    }
    else if ( nBits < 750000 ) {
        K = 5;
        M = 128;
        V[0] = 4; V[1] = 5; V[2] = 6; V[3] = 7; V[4] = 8; V[5] = 9;
        pi[0] = 0.1174035788;
        pi[1] = 0.242955959;
        pi[2] = 0.249363483;
        pi[3] = 0.17517706;
        pi[4] = 0.102701071;
        pi[5] = 0.112398847;
    }
    else {
        K = 6;
        M = 10000;
        V[0] = 10; V[1] = 11; V[2] = 12; V[3] = 13; V[4] = 14; V[5] = 15; V[6] = 16;
        pi[0] = 0.0882;
        pi[1] = 0.2092;
        pi[2] = 0.2483;
        pi[3] = 0.1933;
        pi[4] = 0.1208;
        pi[5] = 0.0675;
        pi[6] = 0.0727;
    }

    N = nBits / M;
    for( i = 0; i < N; i++ ) {
        v_n_obs = 0;
        run = 0;
        for( j = 0; j < M; j++ ) {
            if( Bit[i*M+j] == 1 ) {
                run++;
                if( run > v_n_obs )
                    v_n_obs = run;
            }
            else
                run = 0;
        }
        if ( v_n_obs < V[0] )
            nu[0]++;
        for( j = 0; j <= K; j++ ) {
            if( v_n_obs == V[j] )
                nu[j]++;
        }
        if( v_n_obs > V[K] )
            nu[K]++;
    }

    chi2 = 0.0;
    for( i = 0; i <= K; i++ )
        chi2 += ((nu[i] - N * pi[i]) * (nu[i] - N * pi[i])) / (N * pi[i]);

    res[0].p_value = cephes_igamc((double)(K/2.0), chi2/2.0);
    //COMPUTATIONAL INFORMATION:
    //# of substrings = N
    //Substring Length = M
    //Chi^2 = chi2
    //fprintf(stats[TEST_LONGEST_RUN], "\t\t---------------------------------------------\n");
    //fprintf(stats[TEST_LONGEST_RUN], "\t\t      F R E Q U E N C Y\n");
    //fprintf(stats[TEST_LONGEST_RUN], "\t\t---------------------------------------------\n");

    //if ( K == 3 ) {
    //	fprintf(stats[TEST_LONGEST_RUN], "\t\t  <=1     2     3    >=4   P-value  Assignment");
    //	fprintf(stats[TEST_LONGEST_RUN], "\n\t\t %3d %3d %3d  %3d ", nu[0], nu[1], nu[2], nu[3]);
    //}
    //else if ( K == 5 ) {
    //	fprintf(stats[TEST_LONGEST_RUN], "\t\t<=4  5  6  7  8  >=9 P-value  Assignment");
    //	fprintf(stats[TEST_LONGEST_RUN], "\n\t\t %3d %3d %3d %3d %3d  %3d ", nu[0], nu[1], nu[2],
    //			nu[3], nu[4], nu[5]);
    //}
    //else {
    //	fprintf(stats[TEST_LONGEST_RUN],"\t\t<=10  11  12  13  14  15 >=16 P-value  Assignment");
    //	fprintf(stats[TEST_LONGEST_RUN],"\n\t\t %3d %3d %3d %3d %3d %3d  %3d ", nu[0], nu[1], nu[2],
    //			nu[3], nu[4], nu[5], nu[6]);
    //}
    if ( isNegative(res[0].p_value) || isGreaterThanOne(res[0].p_value) )
    {
        //	fprintf(stats[TEST_LONGEST_RUN], "WARNING:  P_VALUE IS OUT OF RANGE.\n");
        char buff[100];
        std::snprintf(buff, sizeof (buff), " test! p_value = %f!", res[0].p_value);
        throw  CCtrlException(   TestPar[LONGEST_RUN].Name, buff );
    }

    res[0].res = res[0].p_value >= ALPHA;
    return 1;
}

int CNistTests2::Rank(TEST_RESULT *res)	//nBits should be > 40000
{
    if( nBits < TestPar[RANK].MinBitsNum || TestPar[RANK].SubTestsNum != 1 )
    {
        char buff[100];
        std::snprintf(buff, sizeof (buff), " test! Nbits < %d!", TestPar[RANK].MinBitsNum);
        throw  CCtrlException(   TestPar[RANK].Name, buff );
    }

    int		 N;
    int		 i;
    int		 j;
    int		 k;
    int		 r;
    double	 product;
    double	 chi_squared;
    double	 arg1;
    double	 p_32;
    double	 p_31;
    double	 p_30;
    double	 R;
    double	 F_32;
    double	 F_31;
    double	 F_30;
    char	**matrix = new char* [32];	//matrix 32x32
    for( i = 0; i < 32; i++ )
        matrix[i] = new char[32];

    N = nBits/(32*32);
    if( isZero(N) ) {
        //fprintf(stats[TEST_RANK], "\t\tError: Insuffucient # Of Bits To Define An 32x32 (%dx%d) Matrix\n", 32, 32);
        return 0;
    }

    r = 32;					/* COMPUTE PROBABILITIES */
    product = 1;
    for( i = 0; i <= r-1; i++ )
        product *= ((1.e0-pow(2.0, (double)(i-32)))*(1.e0-pow(2.0, (double)(i-32))))/(1.e0-pow(2.0, (double)(i-r)));
    p_32 = pow(2.0, (double)(r*(32+32-r)-32*32)) * product;

    r = 31;
    product = 1;
    for( i = 0; i <= r-1; i++ )
        product *= ((1.e0-pow(2.0, (double)(i-32)))*(1.e0-pow(2.0, (double)(i-32))))/(1.e0-pow(2.0, (double)(i-r)));
    p_31 = pow(2.0, (double)(r*(32+32-r)-32*32)) * product;

    p_30 = 1 - (p_32+p_31);

    F_32 = 0;
    F_31 = 0;
    for( k = 0; k < N; k++ ) {			/* FOR EACH 32x32 MATRIX   */
        for( i = 0; i < 32; i++ )
            for( j = 0; j < 32; j++ )
                matrix[i][j] = Bit[k*(32*32)+j+i*32];
        //#if (DISPLAY_MATRICES == 1)
        //		display_matrix(32, 32, matrix);
        //#endif
        R = computeRank(32, 32, matrix);
        if ( R == 32 )
            F_32++;			/* DETERMINE FREQUENCIES */
        if ( R == 31 )
            F_31++;
    }
    F_30 = (double)N - (F_32+F_31);

    chi_squared =(pow(F_32 - N*p_32, 2)/(double)(N*p_32) +
                  pow(F_31 - N*p_31, 2)/(double)(N*p_31) +
                  pow(F_30 - N*p_30, 2)/(double)(N*p_30));

    arg1 = -chi_squared/2.e0;
    res[0].p_value = exp(arg1);

    //fprintf(stats[TEST_RANK], "\t\t\t\tRANK TEST\n");
    //fprintf(stats[TEST_RANK], "\t\t---------------------------------------------\n");
    //fprintf(stats[TEST_RANK], "\t\tCOMPUTATIONAL INFORMATION:\n");
    //fprintf(stats[TEST_RANK], "\t\t---------------------------------------------\n");
    //fprintf(stats[TEST_RANK], "\t\t(a) Probability P_%d = %f\n", 32,p_32);
    //fprintf(stats[TEST_RANK], "\t\t(b)             P_%d = %f\n", 31,p_31);
    //fprintf(stats[TEST_RANK], "\t\t(c)             P_%d = %f\n", 30,p_30);
    //fprintf(stats[TEST_RANK], "\t\t(d) Frequency   F_%d = %d\n", 32,(int)F_32);
    //fprintf(stats[TEST_RANK], "\t\t(e)             F_%d = %d\n", 31,(int)F_31);
    //fprintf(stats[TEST_RANK], "\t\t(f)             F_%d = %d\n", 30,(int)F_30);
    //fprintf(stats[TEST_RANK], "\t\t(g) # of matrices    = %d\n", N);
    //fprintf(stats[TEST_RANK], "\t\t(h) Chi^2            = %f\n", chi_squared);
    //fprintf(stats[TEST_RANK], "\t\t(i) NOTE: %d BITS WERE DISCARDED.\n", n%(32*32));
    //fprintf(stats[TEST_RANK], "\t\t---------------------------------------------\n");

    for ( i = 0; i < 32; i++ )	/* DEALLOCATE MATRIX  */
        delete [] matrix[i];
    delete [] matrix;

    if ( isNegative(res[0].p_value) || isGreaterThanOne(res[0].p_value) )
    {
        char buff[100];
        std::snprintf(buff, sizeof (buff), " test! p_value = %f!", res[0].p_value);
        throw  CCtrlException(   TestPar[RANK].Name, buff );
    }

    res[0].res = res[0].p_value >= ALPHA;
    return 1;
}

int CNistTests2::DiscreteFourierTransform(TEST_RESULT *res)	//nBits should be > 1000
{
    if( nBits < TestPar[FFT].MinBitsNum || TestPar[FFT].SubTestsNum != 1 )
    {
        char buff[100];
        std::snprintf(buff, sizeof (buff), " test! Nbits < %d!", TestPar[FFT].MinBitsNum);
        throw  CCtrlException(   TestPar[FFT].Name, buff );
    }

    double	 upperBound;
    double	 percentile;
    double	 N_l;
    double	 N_o;
    double	 d;
    int		 count;
    int		 ifac[15];
    uint32_t	i;

    double *X = nullptr;
    double *wsave = nullptr;
    double *m = nullptr;
    X = (double*)calloc(nBits, sizeof(double));
    if( X != nullptr ) {
        wsave = (double*)calloc(2*nBits, sizeof(double));
        if( wsave != nullptr )
            m = (double*)calloc(nBits/2+1, sizeof(double));
    }
    if( X == nullptr || wsave == nullptr || m == nullptr )
    {
        free(X);
        free(wsave);
        free(m);
        char buff[100];
        std::snprintf(buff, sizeof (buff), "%s test!", TestPar[FFT].Name);
        throw  CCtrlException(   nullptr, buff );
    }

    for( i = 0; i < nBits; i++ )
        X[i] = 2*(int)Bit[i] - 1;

    __ogg_fdrffti(nBits, wsave, ifac);		/* INITIALIZE WORK ARRAYS */
    __ogg_fdrfftf(nBits, X, wsave, ifac);	/* APPLY FORWARD FFT */

    m[0] = sqrt(X[0]*X[0]);	    /* COMPUTE MAGNITUDE */

    for( i = 0; i < nBits/2; i++ )
        m[i+1] = sqrt(pow(X[2*i+1],2)+pow(X[2*i+2],2));
    count = 0;				       /* CONFIDENCE INTERVAL */
    upperBound = sqrt(2.995732274*nBits);
    for( i = 0; i < nBits/2; i++ )
        if ( m[i] < upperBound )
            count++;
    percentile = (double)count/(nBits/2)*100;
    N_l = (double) count;       /* number of peaks less than h = sqrt(3*n) */
    N_o = (double) 0.95*nBits/2.0;
    d = (N_l - N_o)/sqrt(nBits/4.0*0.95*0.05);
    res[0].p_value = cephes_erfc(fabs(d)/sqrt(2.0));	//cephes_

    //fprintf(stats[TEST_FFT], "\t\t\t\tFFT TEST\n");
    //fprintf(stats[TEST_FFT], "\t\t-------------------------------------------\n");
    //fprintf(stats[TEST_FFT], "\t\tCOMPUTATIONAL INFORMATION:\n");
    //fprintf(stats[TEST_FFT], "\t\t-------------------------------------------\n");
    //fprintf(stats[TEST_FFT], "\t\t(a) Percentile = %f\n", percentile);
    //fprintf(stats[TEST_FFT], "\t\t(b) N_l        = %f\n", N_l);
    //fprintf(stats[TEST_FFT], "\t\t(c) N_o        = %f\n", N_o);
    //fprintf(stats[TEST_FFT], "\t\t(d) d          = %f\n", d);
    //fprintf(stats[TEST_FFT], "\t\t-------------------------------------------\n");

    free(X);
    free(wsave);
    free(m);

    res[0].res = res[0].p_value >= ALPHA;
    return 1;
}

int CNistTests2::NonOverlappingTemplateMatchings(TEST_RESULT *res, int BlockLength)
{
    if( nBits < TestPar[NONPERIODIC].MinBitsNum || TestPar[NONPERIODIC].SubTestsNum != MAXNUMOFTEMPLATES )
    {
        char buff[100];
        std::snprintf(buff, sizeof (buff), " test! Nbits < %d!", TestPar[NONPERIODIC].MinBitsNum);
        throw  CCtrlException(   TestPar[NONPERIODIC].Name, buff );
    }

    int		numOfTemplates[100] = {0, 0, 2, 4, 6, 12, 20, 40, 74, 148, 284, 568, 1116,
                                       2232, 4424, 8848, 17622, 35244, 70340, 140680, 281076, 562152};
    /*----------------------------------------------------------------------------
        NOTE:  Should additional templates lengths beyond 21 be desired, they must
        first be constructed, saved into files and then the corresponding
        number of nonperiodic templates for that file be stored in the BlockLength-th
        position in the numOfTemplates variable.
        ----------------------------------------------------------------------------*/
    //templates 17 - 21 are in binary format (".Z" extension): other reading needed...
    //if( BlockLength > 16 ) {
    //	AfxMessageBox("NONOVERLAPPING TEMPLATES TEST ABORTED!\nSelected BlockLength > 16. No correspondant template.");
    //	p_value = 0.0;
    //	return false;
    //}

    //std::string filepath;
    //GetAppDataFolder(filepath);
    //filepath.AppendFormat("template%d", BlockLength);
    //CStdioFile file;
    //if( !file.Open( filepath, CFile::modeRead | CFile::typeText ) ) {
    //	AfxMessageBox("NONOVERLAPPING TEMPLATES TEST ABORTED!\nRequired template file not existing.");
    //	p_value = 0.0;
    //	return false;
    //}

    BlockLength = 9;	//current version works with default template9 only!

    unsigned int	/*bit,*/ W_obs;
    unsigned int	/*bit,*/ nu[6];
    //FILE			*fp;
    //char			directory[100];
    double			 sum;
    double			 chi2;
    double			 lambda;
    double			 pi[6];
    double			 varWj;
    int				 i;
    int				 j;
    int				 jj;
    int				 k;
    int				 match;
    int				 SKIP;
    int				 N;
    int				 M;
    int				 K = 5;
    double			p_value_min = DBL_MAX;

    N = 8;
    M = nBits/N;

    lambda = (M-BlockLength+1)/pow(2.0, BlockLength);
    if( isNegative(lambda) || isZero(lambda) ) {
        //AfxMessageBox("NONOVERLAPPING TEMPLATES TESTS ABORTED DUE TO: Lambda not being positive!");
        throw  CCtrlException(   TestPar[NONPERIODIC].Name, " test! Lambda <= 0!" );
    }

    //if ( (Wj = (unsigned int*)calloc(N, sizeof(unsigned int))) == nullptr ) {
    //	fprintf(stats[TEST_NONPERIODIC], "\tNONOVERLAPPING TEMPLATES TESTS ABORTED DUE TO ONE OF THE FOLLOWING : \n");
    //	fprintf(stats[TEST_NONPERIODIC], "\tInsufficient memory for required work space.\n");
    //	return;
    //}

    varWj = M*(1.0/pow(2.0, BlockLength) - (2.0*BlockLength-1.0)/pow(2.0, 2.0*BlockLength));
    auto *Wj = new uint32_t[N];
    memset(Wj, 0, sizeof(uint32_t)*N);

    //char *sequence = new char[BlockLength];

    //sprintf(directory, "templates/template%d", BlockLength);
    //if ( ((isNegative(lambda)) || (isZero(lambda))) ||
    //	 ((fp = fopen(directory, "r")) == nullptr) ||
    //	 ((sequence = (BitSequence *) calloc(m, sizeof(BitSequence))) == nullptr) ) {
    //	fprintf(stats[TEST_NONPERIODIC], "\tNONOVERLAPPING TEMPLATES TESTS ABORTED DUE TO ONE OF THE FOLLOWING : \n");
    //	fprintf(stats[TEST_NONPERIODIC], "\tLambda (%f) not being positive!\n", lambda);
    //	fprintf(stats[TEST_NONPERIODIC], "\tTemplate file <%s> not existing\n", directory);
    //	fprintf(stats[TEST_NONPERIODIC], "\tInsufficient memory for required work space.\n");
    //	if ( sequence != nullptr )
    //		free(sequence);
    //}
    //else {
    //	fprintf(stats[TEST_NONPERIODIC], "\t\t  NONPERIODIC TEMPLATES TEST\n");
    //	fprintf(stats[TEST_NONPERIODIC], "-------------------------------------------------------------------------------------\n");
    //	fprintf(stats[TEST_NONPERIODIC], "\t\t  COMPUTATIONAL INFORMATION\n");
    //	fprintf(stats[TEST_NONPERIODIC], "-------------------------------------------------------------------------------------\n");
    //	fprintf(stats[TEST_NONPERIODIC], "\tLAMBDA = %f\tM = %d\tN = %d\tm = %d\tn = %d\n", lambda, M, N, m, n);
    //	fprintf(stats[TEST_NONPERIODIC], "-------------------------------------------------------------------------------------\n");
    //	fprintf(stats[TEST_NONPERIODIC], "\t\tF R E Q U E N C Y\n");
    //	fprintf(stats[TEST_NONPERIODIC], "Template   W_1  W_2  W_3  W_4  W_5  W_6  W_7  W_8    Chi^2   P_value Assignment Index\n");
    //	fprintf(stats[TEST_NONPERIODIC], "-------------------------------------------------------------------------------------\n");

    if( numOfTemplates[BlockLength] < MAXNUMOFTEMPLATES )
        SKIP = 1;
    else
        SKIP = numOfTemplates[BlockLength]/MAXNUMOFTEMPLATES;
    numOfTemplates[BlockLength] = numOfTemplates[BlockLength]/SKIP;

    sum = 0.0;
    for( i = 0; i < 2; i++ ) {                      /* Compute Probabilities */
        pi[i] = exp(-lambda+(double)i*log(lambda)-cephes_lgam((double)i+1));
        sum += pi[i];
    }
    pi[0] = sum;
    for( i = 2; i <= K; i++ ) {                      /* Compute Probabilities */
        pi[i-1] = exp(-lambda+(double)i*log(lambda)-cephes_lgam((double)i+1));
        sum += pi[i-1];
    }
    pi[K] = 1 - sum;

    int nSubtests = 0;
    for( jj = 0; jj < MIN(MAXNUMOFTEMPLATES, numOfTemplates[BlockLength]); jj++ ) {
        sum = 0;

        ////for( k = 0; k < BlockLength; k++ ) {
        ////	fscanf(fp, "%d", &bit);
        ////	sequence[k] = bit;
        ////	fprintf(stats[TEST_NONPERIODIC], "%d", sequence[k]);
        ////}
        //std::string line;
        //if( !file.ReadString( line ) ) {
        //	AfxMessageBox("NONOVERLAPPING TEMPLATES TESTS ABORTED DUE TO: Unexpected end of template file!");
        //	delete [] sequence;
        //	delete [] Wj;
        //	p_value = 0.0;
        //	return false;
        //}
        //int kk = 0;
        //for( k = 0; k < line.GetLength(); k++ ) {
        //	char ch = line.GetAt(k);
        //	if(ch == '0' || ch == '1') {
        //		sequence[kk++] = ch - '0';
        //		if(kk > BlockLength) {
        //			AfxMessageBox("NONOVERLAPPING TEMPLATES TESTS ABORTED DUE TO: Line in template file > BlockLength!");
        //			delete [] sequence;
        //			delete [] Wj;
        //			p_value = 0.0;
        //			return false;
        //		}
        //		//fprintf(stats[TEST_NONPERIODIC], "%d", sequence[k]);
        //	}
        //}

        //fprintf(stats[TEST_NONPERIODIC], " ");
        for( k = 0; k <= K; k++ )
            nu[k] = 0;
        for( i = 0; i < N; i++ ) {
            W_obs = 0;
            for( j = 0; j < M-BlockLength+1; j++ ) {
                match = 1;
                for( k = 0; k < BlockLength; k++ ) {
                    //if( (int)sequence[k] != (int)Bit[i*M+j+k] ) {
                    if( (int)Template9[jj][k] != (int)Bit[i*M+j+k] ) {
                        match = 0;
                        break;
                    }
                }
                if( match == 1 )
                    W_obs++;
            }
            Wj[i] = W_obs;
        }
        sum = 0;
        chi2 = 0.0;                                   /* Compute Chi Square */
        for( i = 0; i < N; i++ ) {
            //if( m == 10 )
            //	fprintf(stats[TEST_NONPERIODIC], "%3d  ", Wj[i]);
            //else
            //	fprintf(stats[TEST_NONPERIODIC], "%4d ", Wj[i]);
            chi2 += pow(((double)Wj[i] - lambda)/pow(varWj, 0.5), 2);
        }
        res[jj].p_value = cephes_igamc(N/2.0, chi2/2.0);

        if ( isNegative(res[jj].p_value) || isGreaterThanOne(res[jj].p_value) )
        {
            //			fprintf(stats[TEST_NONPERIODIC], "\t\tWARNING:  P_VALUE IS OUT OF RANGE.\n");
            char buff[100];
            std::snprintf(buff, sizeof (buff), "  test! p_value[%d] = %f!", jj, res[jj].p_value);
            throw  CCtrlException(   TestPar[NONPERIODIC].Name, buff );
        }

        //		fprintf(stats[TEST_NONPERIODIC], "%9.6f %f %s %3d\n", chi2, p_value, p_value < ALPHA ? "FAILURE" : "SUCCESS", jj);
        //		if ( SKIP > 1 )
        //			fseek(fp, (long)(SKIP-1)*2*m, SEEK_CUR);
        //		fprintf(results[TEST_NONPERIODIC], "%f\n", p_value); fflush(results[TEST_NONPERIODIC]);

        res[jj].res = res[jj].p_value >= ALPHA;
        nSubtests++;
    }
    //
    //fprintf(stats[TEST_NONPERIODIC], "\n"); fflush(stats[TEST_NONPERIODIC]);

    //fclose(fp);
    //file.Close();

    //if ( sequence != nullptr )
    //	free(sequence);
    //free(Wj);
    //delete [] sequence;
    delete [] Wj;

    return nSubtests;
}
//
//static double Pr(int u, double eta)	//used in OverlappingTemplateMatchings
//{
//	int		l;
//	double	sum, p;
//
//	if ( u == 0 )
//		p = exp(-eta);
//	else {
//		sum = 0.0;
//		for ( l=1; l<=u; l++ )
//			sum += exp(-eta-u*log(2.0)+l*log(eta)-cephes_lgam(l+1)+cephes_lgam(u)-cephes_lgam(l)-cephes_lgam(u-l+1));
//		p = sum;
//	}
//	return p;
//}

int CNistTests2::OverlappingTemplateMatchings(TEST_RESULT *res, int BlockLength)	//m = BlockLength, n = nBits
{
    if( nBits < TestPar[OVERLAPPING].MinBitsNum || TestPar[OVERLAPPING].SubTestsNum != 1 )
    {
        char buff[100];
        std::snprintf(buff, sizeof (buff), " test! Nbits < %d!", TestPar[OVERLAPPING].MinBitsNum);
        throw  CCtrlException(   TestPar[OVERLAPPING].Name, buff );
    }

    int				 i;
    int				 k;
    int				 match;
    double			 W_obs;
    double			 eta;
    double			 sum;
    double			 chi2;
    double			 lambda;
    int				 M;
    int				 N;
    int				 j;
    int				 K = 5;
    unsigned int	nu[6] = { 0, 0, 0, 0, 0, 0 };
    double			pi[6] = { 0.143783, 0.139430, 0.137319, 0.124314, 0.106209, 0.348945 };

    M = 1032;
    N = nBits/M;

    //if ( (sequence = (BitSequence *) calloc(m, sizeof(BitSequence))) == nullptr ) {
    //	fprintf(stats[TEST_OVERLAPPING], "\t\t    OVERLAPPING TEMPLATE OF ALL ONES TEST\n");
    //	fprintf(stats[TEST_OVERLAPPING], "\t\t---------------------------------------------\n");
    //	fprintf(stats[TEST_OVERLAPPING], "\t\tTEMPLATE DEFINITION:  Insufficient memory, Overlapping Template Matchings test aborted!\n");
    //}
    //else
    //	for ( i=0; i<m; i++ )
    //		sequence[i] = 1;
    char *sequence = new char[BlockLength];
    memset(sequence, 1, BlockLength);

    lambda = (double)(M-BlockLength+1)/pow(2.0,BlockLength);
    eta = lambda/2.0;
    sum = 0.0;
    for( i = 0; i < K; i++ ) {			/* Compute Probabilities */
        pi[i] = Pr(i, eta);
        sum += pi[i];
    }
    pi[K] = 1 - sum;

    for( i = 0; i < N; i++ ) {
        W_obs = 0;
        for( j = 0; j < M-BlockLength+1; j++ ) {
            match = 1;
            for( k = 0; k < BlockLength; k++ ) {
                if( sequence[k] != Bit[i*M+j+k] )
                    match = 0;
            }
            if( match == 1 )
                W_obs++;
        }
        if( W_obs <= 4 )
            nu[(int)W_obs]++;
        else
            nu[K]++;
    }
    sum = 0;
    chi2 = 0.0;                                   /* Compute Chi Square */
    for( i = 0; i < K+1; i++ ) {
        chi2 += pow((double)nu[i] - (double)N*pi[i], 2)/((double)N*pi[i]);
        sum += nu[i];
    }
    res[0].p_value = cephes_igamc(K/2.0, chi2/2.0);

    //fprintf(stats[TEST_OVERLAPPING], "\t\t    OVERLAPPING TEMPLATE OF ALL ONES TEST\n");
    //fprintf(stats[TEST_OVERLAPPING], "\t\t-----------------------------------------------\n");
    //fprintf(stats[TEST_OVERLAPPING], "\t\tCOMPUTATIONAL INFORMATION:\n");
    //fprintf(stats[TEST_OVERLAPPING], "\t\t-----------------------------------------------\n");
    //fprintf(stats[TEST_OVERLAPPING], "\t\t(a) n (sequence_length)      = %d\n", n);
    //fprintf(stats[TEST_OVERLAPPING], "\t\t(b) m (block length of 1s)   = %d\n", m);
    //fprintf(stats[TEST_OVERLAPPING], "\t\t(c) M (length of substring)  = %d\n", M);
    //fprintf(stats[TEST_OVERLAPPING], "\t\t(d) N (number of substrings) = %d\n", N);
    //fprintf(stats[TEST_OVERLAPPING], "\t\t(e) lambda [(M-m+1)/2^m]     = %f\n", lambda);
    //fprintf(stats[TEST_OVERLAPPING], "\t\t(f) eta                      = %f\n", eta);
    //fprintf(stats[TEST_OVERLAPPING], "\t\t-----------------------------------------------\n");
    //fprintf(stats[TEST_OVERLAPPING], "\t\t   F R E Q U E N C Y\n");
    //fprintf(stats[TEST_OVERLAPPING], "\t\t  0   1   2   3   4 >=5   Chi^2   P-value  Assignment\n");
    //fprintf(stats[TEST_OVERLAPPING], "\t\t-----------------------------------------------\n");
    //fprintf(stats[TEST_OVERLAPPING], "\t\t%3d %3d %3d %3d %3d %3d  %f ",
    //	nu[0], nu[1], nu[2], nu[3], nu[4], nu[5], chi2);

    //free(sequence);
    delete [] sequence;

    if ( isNegative(res[0].p_value) || isGreaterThanOne(res[0].p_value) )
    {
        //	fprintf(stats[TEST_OVERLAPPING], "WARNING:  P_VALUE IS OUT OF RANGE.\n");
        char buff[100];
        std::snprintf(buff, sizeof (buff), " test! p_value = %f!", res[0].p_value);
        throw  CCtrlException(   TestPar[OVERLAPPING].Name, buff );
    }

    res[0].res = res[0].p_value >= ALPHA;
    return 1;
}

int CNistTests2::Universal(TEST_RESULT *res)	//nBits should be > 162000
{
    if( nBits < TestPar[UNIVERSAL].MinBitsNum || TestPar[UNIVERSAL].SubTestsNum != 1 )
    {
        char buff[100];
        std::snprintf(buff, sizeof (buff), " test! Nbits < %d!", TestPar[UNIVERSAL].MinBitsNum);
        throw  CCtrlException(   TestPar[UNIVERSAL].Name, buff );
    }

    int		 i;
    int		 j;
    int		 p;
    int		 L;
    int		 Q;
    int		 K;
    double	 arg;
    double	 sqrt2;
    double	 sigma;
    double	 phi;
    double	 sum;
    double	 c;
    int64_t	decRep;
    double	expected_value[17] = { 0, 0, 0, 0, 0, 0, 5.2177052, 6.1962507, 7.1836656,
                                       8.1764248, 9.1723243, 10.170032, 11.168765,
                                       12.168070, 13.167693, 14.167488, 15.167379 };
    double   variance[17] = { 0, 0, 0, 0, 0, 0, 2.954, 3.125, 3.238, 3.311, 3.356, 3.384,
                              3.401, 3.410, 3.416, 3.419, 3.421 };

    /* * * * * * * * * ** * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
         * THE FOLLOWING REDEFINES L, SHOULD THE CONDITION:     n >= 1010*2^L*L       *
         * NOT BE MET, FOR THE BLOCK LENGTH L.                                        *
         * * * * * * * * * * ** * * * * * * * * * * * * * * * * * * * * * * * * * * * */
    L = 5;
    if ( nBits >= 387840 )     L = 6;
    if ( nBits >= 904960 )     L = 7;
    if ( nBits >= 2068480 )    L = 8;
    if ( nBits >= 4654080 )    L = 9;
    if ( nBits >= 10342400 )   L = 10;
    if ( nBits >= 22753280 )   L = 11;
    if ( nBits >= 49643520 )   L = 12;
    if ( nBits >= 107560960 )  L = 13;
    if ( nBits >= 231669760 )  L = 14;
    if ( nBits >= 496435200 )  L = 15;
    if ( nBits >= 1059061760 ) L = 16;

    Q = 10*(int)pow(2.0, L);
    K = (int)(floor((double)(nBits/L)) - (double)Q);	 		    /* BLOCKS TO TEST */

    p = (int)pow(2.0, L);
    //if( (L < 6) || (L > 16) || ((double)Q < 10*pow(2.0, L)) ||
    //	 ((T = (long *)calloc(p, sizeof(long))) == nullptr) ) {
    //	fprintf(stats[TEST_UNIVERSAL], "\t\tUNIVERSAL STATISTICAL TEST\n");
    //	fprintf(stats[TEST_UNIVERSAL], "\t\t---------------------------------------------\n");
    //	fprintf(stats[TEST_UNIVERSAL], "\t\tERROR:  L IS OUT OF RANGE.\n");
    //	fprintf(stats[TEST_UNIVERSAL], "\t\t-OR- :  Q IS LESS THAN %f.\n", 10*pow(2.0, L));
    //	fprintf(stats[TEST_UNIVERSAL], "\t\t-OR- :  Unable to allocate T.\n");
    //	return;
    //}
    //for( i = 0; i < p; i++ )
    //	T[i] = 0;
    auto *T = (int64_t*)calloc(p, sizeof(int64_t));
    if( T == nullptr )
    {
        free(T);
        char buff[100];
        std::snprintf(buff, sizeof (buff), "%s test!", TestPar[UNIVERSAL].Name);
        throw  CCtrlException(   nullptr, buff );
    }
    //memset(T, 0, sizeof(int64_t)*p);

    /* COMPUTE THE EXPECTED:  Formula 16, in Marsaglia's Paper */
    c = 0.7 - 0.8/(double)L + (4 + 32/(double)L)*pow(K, -3/(double)L)/15;
    sigma = c * sqrt(variance[L]/(double)K);
    sqrt2 = sqrt(2.0);
    sum = 0.0;
    for( i = 1; i <= Q; i++ ) {		/* INITIALIZE TABLE */
        decRep = 0;
        for( j = 0; j < L; j++ )
            decRep += Bit[(i-1)*L+j] * (int64_t)pow(2.0, L-1-j);
        T[decRep] = i;
    }
    for( i = Q+1; i <= Q+K; i++ ) { 	/* PROCESS BLOCKS */
        decRep = 0;
        for( j = 0; j < L; j++ )
            decRep += Bit[(i-1)*L+j] * (int64_t)pow(2.0, L-1-j);
        sum += log((double)(i - T[decRep]))/log(2.0);
        T[decRep] = i;
    }
    phi = (double)(sum/(double)K);
    arg = fabs(phi-expected_value[L])/(sqrt2 * sigma);
    res[0].p_value = cephes_erfc(arg);	//cephes_

    //fprintf(stats[TEST_UNIVERSAL], "\t\tUNIVERSAL STATISTICAL TEST\n");
    //fprintf(stats[TEST_UNIVERSAL], "\t\t--------------------------------------------\n");
    //fprintf(stats[TEST_UNIVERSAL], "\t\tCOMPUTATIONAL INFORMATION:\n");
    //fprintf(stats[TEST_UNIVERSAL], "\t\t--------------------------------------------\n");
    //fprintf(stats[TEST_UNIVERSAL], "\t\t(a) L         = %d\n", L);
    //fprintf(stats[TEST_UNIVERSAL], "\t\t(b) Q         = %d\n", Q);
    //fprintf(stats[TEST_UNIVERSAL], "\t\t(c) K         = %d\n", K);
    //fprintf(stats[TEST_UNIVERSAL], "\t\t(d) sum       = %f\n", sum);
    //fprintf(stats[TEST_UNIVERSAL], "\t\t(e) sigma     = %f\n", sigma);
    //fprintf(stats[TEST_UNIVERSAL], "\t\t(f) variance  = %f\n", variance[L]);
    //fprintf(stats[TEST_UNIVERSAL], "\t\t(g) exp_value = %f\n", expected_value[L]);
    //fprintf(stats[TEST_UNIVERSAL], "\t\t(h) phi       = %f\n", phi);
    //fprintf(stats[TEST_UNIVERSAL], "\t\t(i) WARNING:  %d bits were discarded.\n", n-(Q+K)*L);
    //fprintf(stats[TEST_UNIVERSAL], "\t\t-----------------------------------------\n");

    free(T);

    if ( isNegative(res[0].p_value) || isGreaterThanOne(res[0].p_value) )
    {
        //	fprintf(stats[TEST_UNIVERSAL], "\t\tWARNING:  P_VALUE IS OUT OF RANGE\n");
        char buff[100];
        std::snprintf(buff, sizeof (buff), " test! p_value = %f!", res[0].p_value);
        throw  CCtrlException(   TestPar[UNIVERSAL].Name, buff );
    }

    res[0].res = res[0].p_value >= ALPHA;
    return 1;
}

int CNistTests2::ApproximateEntropy(TEST_RESULT *res, int BlockLength)	//nBits should be > 10; BlockLength < (int)log2(nBits) - 2
{
    if( nBits < TestPar[APP_ENTROPY].MinBitsNum || TestPar[APP_ENTROPY].SubTestsNum != 1 )
    {
        char buff[100];
        std::snprintf(buff, sizeof (buff), " test! Nbits < %d!", TestPar[APP_ENTROPY].MinBitsNum);
        throw  CCtrlException(   TestPar[APP_ENTROPY].Name, buff );
    }
    //if ( BlockLength > (int)(log((double)seqLength)/log(2.0)-5) ) {
    //	fprintf(stats[TEST_APEN], "\t\tNote: The blockSize = %d exceeds recommended value of %d\n", m,
    //		MAX(1, (int)(log((double)seqLength)/log(2.0)-5)));
    //	fprintf(stats[TEST_APEN], "\t\tResults are inaccurate!\n");

    int	 i;
    int	 j;
    int	 k;
    int	 r;
    int	 blockSize;
    int	 seqLength;
    int	 /*powLen,*/ index;
    uint32_t	powLen;
    double	 sum;
    double	 numOfBlocks;
    double	 ApEn[2];
    double	 apen;
    double	 chi_squared;
    uint32_t	*P;

    seqLength = nBits;
    r = 0;

    for( blockSize = BlockLength; blockSize <= BlockLength+1; blockSize++ ) {
        if( blockSize == 0 ) {
            ApEn[0] = 0.00;
            r++;
        } else {
            numOfBlocks = (double)seqLength;
            powLen = (uint32_t)pow(2.0, (double)(blockSize+1))-1;
            //if ( (P = (unsigned int*)calloc(powLen,sizeof(unsigned int)))== nullptr ) {
            //	fprintf(stats[TEST_APEN], "ApEn:  Insufficient memory available.\n");
            //	return;
            //}
            //for ( i=1; i<powLen-1; i++ )
            //	P[i] = 0;
            P = new uint32_t[powLen];
            memset(P, 0, sizeof(uint32_t)*powLen);
            for( i = 0; i < numOfBlocks; i++ ) { /* COMPUTE FREQUENCY */
                k = 1;
                for( j = 0; j < blockSize; j++ ) {
                    k <<= 1;
                    if( (int)Bit[(i+j) % seqLength] == 1 )
                        k++;
                }
                P[k-1]++;
            }
            /* DISPLAY FREQUENCY */
            sum = 0.0;
            index = (int)pow(2.0, (double)blockSize)-1;
            for( i = 0; i < (int)pow(2.0, (double)blockSize); i++ ) {
                if ( P[index] > 0 )
                    sum += P[index]*log(P[index]/numOfBlocks);
                index++;
            }
            sum /= numOfBlocks;
            ApEn[r] = sum;
            r++;
            //free(P);
            delete [] P;
        }
    }
    apen = ApEn[0] - ApEn[1];
    chi_squared = 2.0*seqLength*(log(2.0) - apen);
    res[0].p_value = cephes_igamc(pow(2.0, BlockLength-1), chi_squared/2.0);

    //fprintf(stats[TEST_APEN], "\t\t(b) n (sequence length) = %d\n", seqLength);
    //fprintf(stats[TEST_APEN], "\t\t(c) Chi^2               = %f\n", chi_squared);
    //fprintf(stats[TEST_APEN], "\t\t(d) Phi(m)	       = %f\n", ApEn[0]);
    //fprintf(stats[TEST_APEN], "\t\t(e) Phi(m+1)	       = %f\n", ApEn[1]);
    //fprintf(stats[TEST_APEN], "\t\t(f) ApEn                = %f\n", apen);
    //fprintf(stats[TEST_APEN], "\t\t(g) Log(2)              = %f\n", log(2.0));
    //fprintf(stats[TEST_APEN], "\t\t--------------------------------------------\n");

    res[0].res = res[0].p_value >= ALPHA;
    return 1;
}

int CNistTests2::RandomExcursions(TEST_RESULT *res)
{
    if( nBits < TestPar[RND_EXCURSION].MinBitsNum || TestPar[RND_EXCURSION].SubTestsNum != 8 )
    {
        char buff[100];
        std::snprintf(buff, sizeof (buff), " test! Nbits < %d!", TestPar[RND_EXCURSION].MinBitsNum);
        throw  CCtrlException(   TestPar[RND_EXCURSION].Name, buff );
    }

    int		nSubtests = 0;
    int		 b;
    int		 i;
    int		 j;
    int		 k;
    int		 J;
    int		 x;
    int		 cycleStart;
    int		 cycleStop;
    int		stateX[8] = { -4, -3, -2, -1, 1, 2, 3, 4 };
    int		counter[8] = { 0, 0, 0, 0, 0, 0, 0, 0 };
    double	 sum;
    double	 constraint;
    double	 nu[6][8];
    double	pi[5][6] = { {0.0000000000, 0.00000000000, 0.00000000000, 0.00000000000, 0.00000000000, 0.0000000000},
                             {0.5000000000, 0.25000000000, 0.12500000000, 0.06250000000, 0.03125000000, 0.0312500000},
                             {0.7500000000, 0.06250000000, 0.04687500000, 0.03515625000, 0.02636718750, 0.0791015625},
                             {0.8333333333, 0.02777777778, 0.02314814815, 0.01929012346, 0.01607510288, 0.0803755143},
                             {0.8750000000, 0.01562500000, 0.01367187500, 0.01196289063, 0.01046752930, 0.0732727051} };
    double	p_value_min = DBL_MAX;

    int *S_k = (int*)calloc(nBits, sizeof(int));
    int *cycle = (int*)calloc(MAX(1000, nBits/100), sizeof(int));
    if( S_k == nullptr || cycle == nullptr )
    {
        free(S_k);
        free(cycle);
        char buff[100];
        std::snprintf(buff, sizeof (buff), "%s test!", TestPar[RND_EXCURSION].Name);
        throw  CCtrlException(   nullptr, buff );
    }

    J = 0; 					/* DETERMINE CYCLES */
    S_k[0] = 2*(int)Bit[0] - 1;
    for( i = 1; i < (int)nBits; i++ ) {
        S_k[i] = S_k[i-1] + 2*Bit[i] - 1;
        if( S_k[i] == 0 ) {
            J++;
            if ( J > (int)MAX(1000, nBits/100) ) {
                free(S_k);
                free(cycle);
                throw  CCtrlException(   TestPar[RND_EXCURSION].Name,
                                          " test! Exceeding the max number of cycles expected." );
            }
            cycle[J] = i;
        }
    }
    if( S_k[nBits-1] != 0 )
        J++;
    cycle[J] = nBits;

    //fprintf(stats[TEST_RND_EXCURSION], "\t\t\t  RANDOM EXCURSIONS TEST\n");
    //fprintf(stats[TEST_RND_EXCURSION], "\t\t--------------------------------------------\n");
    //fprintf(stats[TEST_RND_EXCURSION], "\t\tCOMPUTATIONAL INFORMATION:\n");
    //fprintf(stats[TEST_RND_EXCURSION], "\t\t--------------------------------------------\n");
    //fprintf(stats[TEST_RND_EXCURSION], "\t\t(a) Number Of Cycles (J) = %04d\n", J);
    //fprintf(stats[TEST_RND_EXCURSION], "\t\t(b) Sequence Length (n)  = %d\n", n);

    constraint = MAX(0.005*pow((double)nBits, 0.5), 500);
    if (J < constraint) {
        //fprintf(stats[TEST_RND_EXCURSION], "\t\t---------------------------------------------\n");
        //fprintf(stats[TEST_RND_EXCURSION], "\t\tWARNING:  TEST NOT APPLICABLE.  THERE ARE AN\n");
        //fprintf(stats[TEST_RND_EXCURSION], "\t\t\t  INSUFFICIENT NUMBER OF CYCLES.\n");
        //fprintf(stats[TEST_RND_EXCURSION], "\t\t---------------------------------------------\n");
        //for(i = 0; i < 8; i++)
        //	fprintf(results[TEST_RND_EXCURSION], "%f\n", 0.0);
        free(S_k);
        free(cycle);
        //throw  CCtrlException(   0, TestPar[RND_EXCURSION].Name,
        //									" test! There are an insufficient number of cycles." );
        return 0;
    } 
        //fprintf(stats[TEST_RND_EXCURSION], "\t\t(c) Rejection Constraint = %f\n", constraint);
        //fprintf(stats[TEST_RND_EXCURSION], "\t\t-------------------------------------------\n");

        cycleStart = 0;
        cycleStop  = cycle[1];
        for( k = 0; k < 6; k++ )
            for( i = 0; i < 8; i++ )
                nu[k][i] = 0.;
        for( j = 1; j <= J; j++ ) {                           /* FOR EACH CYCLE */
            for( i = 0; i < 8; i++ )
                counter[i] = 0;
            for( i = cycleStart; i < cycleStop; i++ ) {
                if( (S_k[i] >= 1 && S_k[i] <= 4) || (S_k[i] >= -4 && S_k[i] <= -1) ) {
                    if( S_k[i] < 0 )
                        b = 4;
                    else
                        b = 3;
                    counter[S_k[i]+b]++;
                }
            }
            cycleStart = cycle[j]+1;
            if( j < J )
                cycleStop = cycle[j+1];

            for( i = 0; i < 8; i++ ) {
                if ( (counter[i] >= 0) && (counter[i] <= 4) )
                    nu[counter[i]][i]++;
                else if ( counter[i] >= 5 )
                    nu[5][i]++;
            }
        }

        free(S_k);
        free(cycle);

        for( i = 0; i < 8; i++ ) {
            x = stateX[i];
            sum = 0.;
            for( k = 0; k < 6; k++ )
                sum += pow(nu[k][i] - J*pi[(int)fabs((double)x)][k], 2) / (J*pi[(int)fabs((double)x)][k]);
            res[i].p_value = cephes_igamc(2.5, sum/2.0);

            if ( isNegative(res[i].p_value) || isGreaterThanOne(res[i].p_value) )
            {
                //	fprintf(stats[TEST_RND_EXCURSION], "WARNING:  P_VALUE IS OUT OF RANGE.\n");
                char buff[100];
                std::snprintf(buff, sizeof (buff), " test! p_value[%d] = %f!", i, res[i].p_value);
                throw  CCtrlException(   TestPar[NONPERIODIC].Name, buff );
            }

            //fprintf(stats[TEST_RND_EXCURSION], "%s\t\tx = %2d chi^2 = %9.6f p_value = %f\n",
            //		p_value < ALPHA ? "FAILURE" : "SUCCESS", x, sum, p_value);

            res[i].res = res[i].p_value >= ALPHA;
            nSubtests++;
        }
    

    return nSubtests;
}

int CNistTests2::RandomExcursionsVariant(TEST_RESULT *res)
{
    if( nBits < TestPar[RND_EXCURSION_VAR].MinBitsNum || TestPar[RND_EXCURSION_VAR].SubTestsNum != 18 )
    {
        char buff[100];
        std::snprintf(buff, sizeof (buff), " test! Nbits < %d!", TestPar[RND_EXCURSION_VAR].MinBitsNum);
        throw  CCtrlException(   TestPar[RND_EXCURSION_VAR].Name, buff );
    }

    int		 i;
    int		 p;
    int		 J;
    int		 x;
    int		 constraint;
    int		 count;
    int		stateX[18] = { -9, -8, -7, -6, -5, -4, -3, -2, -1, 1, 2, 3, 4, 5, 6, 7, 8, 9 };
    double	p_value_min = DBL_MAX;

    int *S_k = (int*)calloc(nBits, sizeof(int));
    if( S_k == nullptr )
    {
        free(S_k);
        char buff[100];
        std::snprintf(buff, sizeof (buff), "%s test!", TestPar[RND_EXCURSION_VAR].Name);
        throw  CCtrlException(   nullptr, buff );
    }

    J = 0;
    S_k[0] = 2*(int)Bit[0] - 1;
    for( i = 1; i < (int)nBits; i++ ) {
        S_k[i] = S_k[i-1] + 2*Bit[i] - 1;
        if ( S_k[i] == 0 )
            J++;
    }
    if ( S_k[nBits-1] != 0 )
        J++;

    //fprintf(stats[TEST_RND_EXCURSION_VAR], "\t\t\tRANDOM EXCURSIONS VARIANT TEST\n");
    //fprintf(stats[TEST_RND_EXCURSION_VAR], "\t\t--------------------------------------------\n");
    //fprintf(stats[TEST_RND_EXCURSION_VAR], "\t\tCOMPUTATIONAL INFORMATION:\n");
    //fprintf(stats[TEST_RND_EXCURSION_VAR], "\t\t--------------------------------------------\n");
    //fprintf(stats[TEST_RND_EXCURSION_VAR], "\t\t(a) Number Of Cycles (J) = %d\n", J);
    //fprintf(stats[TEST_RND_EXCURSION_VAR], "\t\t(b) Sequence Length (n)  = %d\n", n);
    //fprintf(stats[TEST_RND_EXCURSION_VAR], "\t\t--------------------------------------------\n");

    int nSubtests = 0;
    constraint = (int)MAX(0.005*pow((double)nBits, 0.5), 5e2);
    if (J < constraint) {
        //fprintf(stats[TEST_RND_EXCURSION_VAR], "\n\t\tWARNING:  TEST NOT APPLICABLE.  THERE ARE AN\n");
        //fprintf(stats[TEST_RND_EXCURSION_VAR], "\t\t\t  INSUFFICIENT NUMBER OF CYCLES.\n");
        //fprintf(stats[TEST_RND_EXCURSION_VAR], "\t\t---------------------------------------------\n");
        //for ( i=0; i<18; i++ )
        //	fprintf(results[TEST_RND_EXCURSION_VAR], "%f\n", 0.0);
        free(S_k);
        //throw  CCtrlException(   0, TestPar[RND_EXCURSION_VAR].Name,
        //									" test! There are an insufficient number of cycles." );
        return 0;
    } 
        for( p = 0; p <= 17; p++ ) {
            x = stateX[p];
            count = 0;
            for( i = 0;  i < (int)nBits; i++ )
                if( S_k[i] == x )
                    count++;
            res[p].p_value = cephes_erfc(fabs((double)(count-J))/(sqrt(2.0*J*(4.0*fabs((double)x)-2))));	//cephes_

            if ( isNegative(res[p].p_value) || isGreaterThanOne(res[p].p_value) )
            {
                free(S_k);
                char buff[100];
                std::snprintf(buff, sizeof (buff), " test! p_value[%d] = %f!", p, res[p].p_value);
                throw  CCtrlException(   TestPar[RND_EXCURSION_VAR].Name, buff );
            }
            //fprintf(stats[TEST_RND_EXCURSION_VAR], "%s\t\t", p_value < ALPHA ? "FAILURE" : "SUCCESS");
            //fprintf(stats[TEST_RND_EXCURSION_VAR], "(x = %2d) Total visits = %4d; p-value = %f\n", x, count, p_value);
            res[p].res = res[p].p_value >= ALPHA;
            nSubtests++;
        }
    
    free(S_k);

    return nSubtests;
}

double CNistTests2::psi2(int m, int n)	//used in Serial
{
    int				 i;
    int				 j;
    int				 k;
    int				 powLen;
    double			 sum;
    double			 numOfBlocks;

    if ( (m == 0) || (m == -1) )
        return 0.0;
    numOfBlocks = n;
    powLen = (int)pow(2.0, m+1)-1;

    auto *P = (uint32_t*)calloc(powLen, sizeof(uint32_t));
    if( P == nullptr )
    {
        free(P);
        char buff[100];
        std::snprintf(buff, sizeof (buff), "%s test in psi2 func!", TestPar[SERIAL].Name);
        throw  CCtrlException(   nullptr, buff );
    }
    //for ( i=1; i<powLen-1; i++ )
    //	P[i] = 0;	  /* INITIALIZE NODES */
    memset(P, 0, sizeof(uint32_t)*powLen);

    for( i = 0; i < numOfBlocks; i++ ) {		 /* COMPUTE FREQUENCY */
        k = 1;
        for( j = 0; j < m; j++ ) {
            if( Bit[(i+j)%n] == 0 )
                k *= 2;
            else if ( Bit[(i+j)%n] == 1 )
                k = 2*k+1;
        }
        P[k-1]++;
    }
    sum = 0.0;
    for( i = (int)pow(2.0, m)-1; i < (int)pow(2.0, m+1)-1; i++ )
        sum += pow(P[i], 2.0);
    sum = (sum * pow(2.0, m)/(double)n) - (double)n;

    free(P);

    return sum;
}

int CNistTests2::Serial(TEST_RESULT *res, int BlockLength)	//nBits should be > 10; BlockLength < (int)log2(nBits) - 2
{
    if( nBits < TestPar[SERIAL].MinBitsNum || TestPar[SERIAL].SubTestsNum != 2 )
    {
        char buff[100];
        std::snprintf(buff, sizeof (buff), " test! Nbits < %d!", TestPar[SERIAL].MinBitsNum);
        throw  CCtrlException(   TestPar[SERIAL].Name, buff );
    }

    double	 psim0;
    double	 psim1;
    double	 psim2;
    double	 del1;
    double	 del2;

    psim0 = psi2(BlockLength, nBits);
    psim1 = psi2(BlockLength-1, nBits);
    psim2 = psi2(BlockLength-2, nBits);
    del1 = psim0 - psim1;
    del2 = psim0 - 2.0*psim1 + psim2;
    res[0].p_value = cephes_igamc(pow(2.0, BlockLength-1)/2, del1/2.0);
    res[0].res = res[0].p_value >= ALPHA;
    res[1].p_value = cephes_igamc(pow(2.0, BlockLength-2)/2, del2/2.0);
    res[1].res = res[1].p_value >= ALPHA;

    //fprintf(stats[TEST_SERIAL], "\t\t\t       SERIAL TEST\n");
    //fprintf(stats[TEST_SERIAL], "\t\t---------------------------------------------\n");
    //fprintf(stats[TEST_SERIAL], "\t\t COMPUTATIONAL INFORMATION:		  \n");
    //fprintf(stats[TEST_SERIAL], "\t\t---------------------------------------------\n");
    //fprintf(stats[TEST_SERIAL], "\t\t(a) Block length    (m) = %d\n", m);
    //fprintf(stats[TEST_SERIAL], "\t\t(b) Sequence length (n) = %d\n", n);
    //fprintf(stats[TEST_SERIAL], "\t\t(c) Psi_m               = %f\n", psim0);
    //fprintf(stats[TEST_SERIAL], "\t\t(d) Psi_m-1             = %f\n", psim1);
    //fprintf(stats[TEST_SERIAL], "\t\t(e) Psi_m-2             = %f\n", psim2);
    //fprintf(stats[TEST_SERIAL], "\t\t(f) Del_1               = %f\n", del1);
    //fprintf(stats[TEST_SERIAL], "\t\t(g) Del_2               = %f\n", del2);
    //fprintf(stats[TEST_SERIAL], "\t\t---------------------------------------------\n");

    //fprintf(stats[TEST_SERIAL], "%s\t\tp_value1 = %f\n", p_value1 < ALPHA ? "FAILURE" : "SUCCESS", p_value1);
    //fprintf(results[TEST_SERIAL], "%f\n", p_value1);

    //fprintf(stats[TEST_SERIAL], "%s\t\tp_value2 = %f\n\n", p_value2 < ALPHA ? "FAILURE" : "SUCCESS", p_value2);
    //fprintf(results[TEST_SERIAL], "%f\n", p_value2);

    return 2;
}

int CNistTests2::LinearComplexity(TEST_RESULT *res, int BlockLength)		//nBits should be > 1000000
{
    if( nBits < TestPar[LINEARCOMPLEXITY].MinBitsNum || TestPar[LINEARCOMPLEXITY].SubTestsNum != 1 )
    {
        char buff[100];
        std::snprintf(buff, sizeof (buff), " test! Nbits < %d!", TestPar[LINEARCOMPLEXITY].MinBitsNum);
        throw  CCtrlException(   TestPar[LINEARCOMPLEXITY].Name, buff );
    }

    int       i;
    int       ii;
    int       j;
    int       d;
    int       N;
    int       L;
    int       m;
    int       N_;
    int       parity;
    int       sign;
    int       K = 6;
    double    T_;
    double    mean;
    double    nu[7];
    double    chi2;
    double    pi[7] = { 0.01047, 0.03125, 0.12500, 0.50000, 0.25000, 0.06250, 0.020833 };

    N = (int)floor((double)(nBits/BlockLength));

    char *B_ = (char*)calloc(BlockLength, sizeof(char));
    char *C  = (char*)calloc(BlockLength, sizeof(char));
    char *P  = (char*)calloc(BlockLength, sizeof(char));
    char *T  = (char*)calloc(BlockLength, sizeof(char));
    if( B_ == nullptr || C ==nullptr || P == nullptr || T == nullptr ) {
        free(B_);
        free(P);
        free(C);
        free(T);
        char buff[100];
        std::snprintf(buff, sizeof (buff), "%s test!", TestPar[LINEARCOMPLEXITY].Name);
        throw  CCtrlException(   nullptr, buff );
    }

    //fprintf(stats[TEST_LINEARCOMPLEXITY], "-----------------------------------------------------\n");
    //fprintf(stats[TEST_LINEARCOMPLEXITY], "\tL I N E A R  C O M P L E X I T Y\n");
    //fprintf(stats[TEST_LINEARCOMPLEXITY], "-----------------------------------------------------\n");
    //fprintf(stats[TEST_LINEARCOMPLEXITY], "\tM (substring length)     = %d\n", M);
    //fprintf(stats[TEST_LINEARCOMPLEXITY], "\tN (number of substrings) = %d\n", N);
    //fprintf(stats[TEST_LINEARCOMPLEXITY], "-----------------------------------------------------\n");
    //fprintf(stats[TEST_LINEARCOMPLEXITY], "        F R E Q U E N C Y                            \n");
    //fprintf(stats[TEST_LINEARCOMPLEXITY], "-----------------------------------------------------\n");
    //fprintf(stats[TEST_LINEARCOMPLEXITY], "  C0   C1   C2   C3   C4   C5   C6    CHI2    P-value\n");
    //fprintf(stats[TEST_LINEARCOMPLEXITY], "-----------------------------------------------------\n");
    //fprintf(stats[TEST_LINEARCOMPLEXITY], "\tNote: %d bits were discarded!\n", n%M);

    for( i = 0; i < K+1; i++ )
        nu[i] = 0.00;
    for( ii = 0; ii < N; ii++ ) {
        for( i = 0; i < BlockLength; i++ ) {
            B_[i] = 0;
            C[i] = 0;
            T[i] = 0;
            P[i] = 0;
        }
        L = 0;
        m = -1;
        d = 0;
        C[0] = 1;
        B_[0] = 1;

        /* DETERMINE LINEAR COMPLEXITY */
        N_ = 0;
        while( N_ < BlockLength ) {
            d = (int)Bit[ii*BlockLength+N_];
            for( i = 1; i <= L; i++ )
                d += C[i] * Bit[ii*BlockLength+N_-i];
            d = d%2;
            if( d == 1 ) {
                for( i = 0; i < BlockLength; i++ ) {
                    T[i] = C[i];
                    P[i] = 0;
                }
                for( j = 0; j < BlockLength; j++ )
                    if( B_[j] == 1 )
                        P[j+N_-m] = 1;
                for( i = 0; i < BlockLength; i++ )
                    C[i] = (C[i] + P[i])%2;
                if( L <= N_/2 ) {
                    L = N_ + 1 - L;
                    m = N_;
                    for( i = 0; i < BlockLength; i++ )
                        B_[i] = T[i];
                }
            }
            N_++;
        }
        if( (parity = (BlockLength+1)%2) == 0 )
            sign = -1;
        else
            sign = 1;
        mean = BlockLength/2.0 + (9.0+sign)/36.0 - 1.0/pow(2.0, BlockLength) * (BlockLength/3.0 + 2.0/9.0);
        if( (parity = BlockLength%2) == 0 )
            sign = 1;
        else
            sign = -1;
        T_ = sign * (L - mean) + 2.0/9.0;

        if( T_ <= -2.5 )
            nu[0]++;
        else if ( T_ > -2.5 && T_ <= -1.5 )
            nu[1]++;
        else if ( T_ > -1.5 && T_ <= -0.5 )
            nu[2]++;
        else if ( T_ > -0.5 && T_ <= 0.5 )
            nu[3]++;
        else if ( T_ > 0.5 && T_ <= 1.5 )
            nu[4]++;
        else if ( T_ > 1.5 && T_ <= 2.5 )
            nu[5]++;
        else
            nu[6]++;
    }
    chi2 = 0.00;
    //for( i = 0; i < K+1; i++ )
    //	fprintf(stats[TEST_LINEARCOMPLEXITY], "%4d ", (int)nu[i]);
    for( i = 0; i < K+1; i++ )
        chi2 += pow(nu[i]-N*pi[i], 2) / (N*pi[i]);
    res[0].p_value = cephes_igamc(K/2.0, chi2/2.0);
    res[0].res = res[0].p_value >= ALPHA;
    //fprintf(stats[TEST_LINEARCOMPLEXITY], "%9.6f%9.6f\n", chi2, p_value);
    //fprintf(results[TEST_LINEARCOMPLEXITY], "%f\n", p_value);

    free(B_);
    free(P);
    free(C);
    free(T);

    return 1;
}

std::string CNistTests2::GetDllVersion()
{
    std::string	str = "Unavailable";
    //	char filepath[_MAX_PATH];
    //	HMODULE hModule = GetModuleHandle("Nist.dll");
    //	GetModuleFileName( hModule, filepath, _MAX_PATH );
    //	DWORD	size, buf;
    //	if( size = GetFileVersionInfoSize( filepath, &buf ) ) {
    //		void*	pData = malloc(size);
    //		if( GetFileVersionInfo( filepath, buf, size, pData ) ) {
    //			VS_FIXEDFILEINFO vsf;
    //			void*	pbuf;
    //            uint32_t	len;
    //			if( VerQueryValue( pData, "\\", &pbuf, &len ) ) {
    //				memcpy( &vsf, pbuf, sizeof(VS_FIXEDFILEINFO) );
    //				str.Format("%d.%d.%d.%d",
    //					(WORD)(vsf.dwProductVersionMS >> 16), (WORD)(vsf.dwProductVersionMS),
    //					(WORD)(vsf.dwProductVersionLS >> 16), (WORD)(vsf.dwProductVersionLS) );
    //			}
    // 		}
    //		free(pData);
    //	}
    return str;
}
