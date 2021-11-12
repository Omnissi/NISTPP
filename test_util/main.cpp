#include <nistpp/tests.h>
#include <nistpp/bits_storage.h>

#include <boost/filesystem.hpp>

#include <iostream>
#include <chrono>

class time
{
public:
    time() = default;
    ~time() noexcept = default;

    /**
     * @brief Запуск секундомера.
     */
    void start()
    {
        m_begin = std::chrono::high_resolution_clock::now();
    }

    /**
     * @brief Перезапуск секундомера и получение времени, прошедшего с последнего старта или рестарта.
     * @return Время, прошедшее с последнего стрта/рестара, в миллисекундах.
     *
     * @todo Заменить на любое время (секунды, часы и т.д.)
     */
    uint32_t restart()
    {
        auto tmp = elapsed();
        start();

        return tmp;
    }

    /**
     * @brief Получение времени, прошедшего с последнего старта или рестарта.
     * @return Время, прошедшее с последнего стрта/рестара, в миллисекундах.
     */
    uint32_t elapsed()
    {
        typedef std::chrono::duration<uint32_t, std::ratio<1, 1000>> milli;
        return std::chrono::duration_cast<milli>(std::chrono::high_resolution_clock::now() - m_begin).count();
    }

private:
    std::chrono::high_resolution_clock::time_point m_begin{};
};

int main(int argc, char* argv[])
{
    if(argc < 4)
    {
        std::cerr << "File not set!" << std::endl;
        return -1;
    }

    boost::filesystem::path path(argv[1]);
    if( !boost::filesystem::exists( path ) )
    {
        std::cerr << "File [" << path.c_str() << "] doesn't exist!" << std::endl;
    }
    std::ifstream stream;
    stream.open(path.string(), std::ios::binary);

    if(stream.fail() || stream.bad())
    {
        std::cerr <<  "Can't open file: " << stream.failbit  << std::endl;
    }

    auto fileSize = boost::filesystem::file_size(path);

    const size_t blockSize = static_cast<std::size_t>(atoi(argv[3]));
    std::vector<uint8_t> tmp(static_cast<std::size_t>(atoi(argv[2]) / 8));

    if(blockSize * tmp.size() > fileSize)
    {
        std::cerr << "Not enough data!" << std::endl;
        return -1;
    }

    class time res_t;
    res_t.start();
    for(size_t i = 0; i < blockSize; ++i)
    {
//        std::cout << "Reading file..." << std::endl;
        std::size_t size = 0;
        class time t;
        t.start();
        while(size < tmp.size())
        {
            auto err = stream.readsome(reinterpret_cast<char*>(tmp.data()) + size, static_cast<ssize_t>(tmp.size() - size));
            if(err < 0)
            {
                std::cerr << "Error read sequence" << std::endl;
            }

            size += static_cast<std::size_t>(err);
        }
//        std::cout << "Complete read file: " << t.restart() << std::endl;

//        std::cout << "Start parsing..." << std::endl;
        nistpp::BitsStorage bitsStorage(tmp);
//        std::cout << "Complete parsing: " << t.restart() << std::endl;

        std::vector<std::pair<std::string, std::function<nistpp::return_t()>>> tests =
        {
                    {"freq",        [&](){ return nistpp::FrequencyTest(bitsStorage); }},
                    {"BlkFreq",     [&](){ return nistpp::BlockFrequencyTest(bitsStorage, 128); }},
                    {"Runs",        [&](){ return nistpp::RunsTest(bitsStorage); }},
                    {"LongRuns",    [&](){ return nistpp::LongestRunOfOnesTest(bitsStorage); }},
                    {"Runk",        [&](){ return nistpp::RankTest(bitsStorage); }},
                    {"FFT",         [&](){ return nistpp::FftTest(bitsStorage); }},
                    {"NonOverlap",  [&](){
                        std::vector<double> P;
                        auto res = nistpp::NonOverlappingTemplateTest(bitsStorage, 9, P);
//                        for(auto& el : P)
//                        {
//                            std::cout << "    " << el << std::endl;
//                        }
                        return res; }},
                    {"Overlap",     [&](){ return nistpp::OverlappingTemplateTest(bitsStorage, 9); }},
                    {"Universal",   [&](){ return nistpp::UniversalTest(bitsStorage); }},
                    {"Linear",      [&](){ return nistpp::LinearComplexityTest(bitsStorage, 500); }},
                    {"Serial",      [&](){ return nistpp::SerialTest(bitsStorage, 16); }},
                    {"Approximate", [&](){ return nistpp::ApproximateEntropyTest(bitsStorage, 10); }},
                    {"Cusum",       [&](){ return nistpp::CumulativeSumsTest(bitsStorage); }},
                    {"Rand",        [&](){
                        std::array<double, 8> P{};
                        nistpp::return_t res;
                        try
                        {
                            res = nistpp::RandomExcursionsTest(bitsStorage, P);
                        }
                        catch(const std::exception& ex)
                        {
                            std::cerr << "Rand failed: " << ex.what() << std::endl;
                            std::fill(P.begin(), P.end(), 0);
                        }

//                        for(auto& el : P)
//                        {
//                            std::cout << "    " << el << std::endl;
//                        }
                        return res;
                    }},
                    {"RandVar",     [&](){
                        std::array<double, 18> P{};
                        nistpp::return_t res;
                        try
                        {
                            res = nistpp::RandomExcursionsVariantTest(bitsStorage, P);
                        }
                        catch(const std::exception& ex)
                        {
                            std::cerr << "RandVar failed: " << ex.what() << std::endl;
                            std::fill(P.begin(), P.end(), 0);
                        }
//                        for(auto& el : P)
//                        {
//                            std::cout << "    " << el << std::endl;
//                        }
                        return res;
                    }},
    };

        t.restart();
#pragma omp parallel for
        for(auto& el : tests)
        {
            el.second();
        }

        std::cout << "Complete iteration: " << i << ". Time: " << t.restart() << std::endl;
    }
    std::cout << "Testing complete! Result time: " << res_t.elapsed() << std::endl;

    return 0;
}
