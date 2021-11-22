#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/search/views/minimiser_hash.hpp>
#include <vector>

int convert2base4(auto value)
{

    auto converted_val = 0;
    int i = 0;
    while (value > 0)
    {
        converted_val += (pow(10, i) * (value % 4));
        value = value / 4;
        i++;
    }
    return converted_val;
}

void testMinimizer()
{
    using namespace seqan3::literals;
    seqan3::dna4_vector text = "CCACGTCGACGGTT"_dna4;
    std::vector<seqan3::dna4> text2{"CCACGTCGACGGTT"_dna4};

    // Here a consecutive shape with size 4 (so the k-mer size is 4) and a window size of 8 is used.
    auto minimisers = text | seqan3::views::minimiser_hash(seqan3::shape{seqan3::ungapped{4}},
                                                           seqan3::window_size{8});
    // results in: [10322096095657499224, 10322096095657499142, 10322096095657499224]
    // representing the k-mers [GACG, TCGA, GACG]
    seqan3::debug_stream << minimisers << '\n';
    seqan3::debug_stream << *minimisers.begin();

    // getting the hash for the first element
    uint64_t seed = 0x8F3F73B5CF1C9ADE;
    uint64_t first_elem = *minimisers.begin();
    auto hash_first_elem = first_elem ^ seed;
    seqan3::debug_stream << "recover kmer " << hash_first_elem << std::endl;

    // getting the hash for all elements
    auto hash_values = minimisers | std::views::transform([seed](uint64_t i)
                                                          { return i ^ seed; });
    seqan3::debug_stream << hash_values << '\n'; // results in: [182, 216, 134]

    //auto convert_values = hash_values | std::views::transform(convert2base4);
    //seqan3::debug_stream << convert_values << '\n'; // results in: [182, 216, 134]

    // getting the hash using normal for
    for (auto &&h : minimisers)
    {
        auto hash = (h ^ seed);
        auto converted = convert2base4(hash);
        seqan3::debug_stream << hash << std::endl;
        seqan3::debug_stream << "to base 4 " << converted << std::endl;
    }

    // speficify seed when computing the minimiser
    auto example_a = text | seqan3::views::minimiser_hash(seqan3::shape{seqan3::ungapped{4}},
                                                          seqan3::window_size{4},
                                                          seqan3::seed{0}); // sets lexicographic
    auto convert_values2 = example_a | std::views::transform([](auto i){ return convert2base4(i); });

    seqan3::debug_stream << example_a << std::endl;
    seqan3::debug_stream << "to base 4 " << convert_values2 << std::endl;
}