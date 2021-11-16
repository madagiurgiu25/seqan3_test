#include <vector>
#include <iostream>
#include <seqan3/std/ranges>  
 
// use the seqan3::view functionality
#include <seqan3/alphabet/views/all.hpp>
#include <seqan3/alphabet/views/complement.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/argument_parser/all.hpp> 

// used for the bitpacked sequence functionality
#include <seqan3/alphabet/container/bitpacked_sequence.hpp>     // include bitpacked sequence
#include <seqan3/alphabet/nucleotide/dna4.hpp>         

void testRanges(){

    std::vector vec{1, 2, 3, 4, 5, 6};
    auto v = std::views::reverse(vec);
    
    // why cant I print std::cout << v.begin()
    std::cout << *v.begin() << '\n';

    // filter out uneven numbers
    auto even = [](int i) { return 0 == i % 2; };
    auto square = [](int i) { return i * i; };

    for (int i : vec 
               | std::views::filter(even) 
               | std::views::transform(square)) {
        std::cout << i << ' ';
    }

    auto v2 = vec | std::views::filter(even) | std::views::transform(square);

    // v = vec | std::views::filter([] (auto const i) { return i % 2 == 0; }) | std::views::transform([] (auto const i) { return i*i; });
    
    for (int t : v2){
        std::cout << t << " ";
    }
    std::cout << *v2.begin() << '\n';
}

int complementSequence(int argc, char * argv[]){

    // We use the seqan3::argument_parser which was introduced in the second chapter
    // of the tutorial: "Parsing command line arguments with SeqAn".
    seqan3::argument_parser myparser{"Assignment-3", argc, argv}; // initialize
    std::string s{};
 
    myparser.add_positional_option(s, "Please specify the DNA string.");
 
    try
    {
       myparser.parse();
    }
    catch (seqan3::argument_parser_error const & ext) // the user did something wrong
    {
       std::cerr << "[PARSER ERROR]" << ext.what() << '\n'; // you can customize your error message
       return 0;
    }
 
    // read the string into an seqan3::dna5 object
    auto s_as_dna = s | seqan3::views::char_to<seqan3::dna5>;

    seqan3::debug_stream << "Original: " << s_as_dna << '\n';
    
    // this checks for characters not in the alphabet
    auto s_as_dna2 = s | std::views::transform([] (char const c)
    {
       return seqan3::assign_char_strictly_to(c, seqan3::dna5{});
    });

    // seqan3::debug_stream << "Cleaned: " << s_as_dna << '\n';
    seqan3::debug_stream << "Reverse:  " << (s_as_dna  | std::views::reverse ) << '\n';
    seqan3::debug_stream << "Complement:  " << (s_as_dna  | seqan3::views::complement) << '\n';
    seqan3::debug_stream << "RevComp:  " << (s_as_dna | std::views::reverse | seqan3::views::complement) << '\n';
    seqan3::debug_stream << "Frames:   " << (s_as_dna | seqan3::views::translate(seqan3::translation_frames::six_frames)) << '\n';

    return 1;
}

int bitpackedSequence(int argc, char * argv[]){

    using namespace seqan3::literals;
 
    seqan3::argument_parser myparser("Vector-implementations-comparison", argc, argv);
    size_t size{};
    bool use_bitvector{};
    myparser.add_positional_option(size, "Size of vector");
    myparser.add_flag(use_bitvector, 'b', "bitvector", "Use bitvector instead of vector");
 
    try
    {
         myparser.parse();
    }
    catch (seqan3::argument_parser_error const & ext)                     // catch user errors
    {
        seqan3::debug_stream << "[Error] " << ext.what() << "\n";
        return -1;
    }
 
    if (use_bitvector)
    {
        seqan3::bitpacked_sequence<seqan3::dna4> vector{size, 'A'_dna4};
        // vector.resize(size, 'A'_dna4);
        seqan3::debug_stream << "Allocated seqan3::bitpacked_sequence<seqan3::dna4> of size "
                             << vector.size() << '\n';
    }
    else
    {
        std::vector<seqan3::dna4> vector{size};
        // vector.resize(size, 'A'_dna4);
        seqan3::debug_stream << "Allocated std::vector<seqan3::dna4> of size " << vector.size() << '\n';
    }

    return 0;
}