#include <array> // std::array
#include <string> // std::string
#include <vector> // std::vector
#include <iostream>
#include <seqan3/alphabet/all.hpp>
#include <seqan3/argument_parser/all.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp> 


using namespace seqan3::literals;

int readSequence(int argc, char * argv[]){

    std::string input{};
    seqan3::argument_parser parser("GC-Content", argc, argv);
    parser.add_positional_option(input, "Specify an input sequence.");

    try
    {
        parser.parse();
    }
    catch (seqan3::argument_parser_error const & ext) // the input is invalid
    {
        seqan3::debug_stream << "[PARSER ERROR] " << ext.what() << '\n';
        return -1;
    }

    std::vector<seqan3::dna5> sequence{};
    for (char c : input){
        sequence.push_back( seqan3::assign_char_to(c, seqan3::dna5{}) );
    }

    std::array<size_t, seqan3::dna5::alphabet_size> count{}; // default initialised with zeroes
 
    // Increase the symbol count according to the sequence.
    for (seqan3::dna5 symbol : sequence){
        ++count[symbol.to_rank()];
        seqan3::debug_stream << symbol.to_char() << " " << symbol.to_rank() << std::endl;
    }

    seqan3::debug_stream << std::endl << "A C G T" << std::endl;
    for (size_t v: count){
        seqan3::debug_stream << v << " ";
    }
 
    // Calculate the GC content: (#G + #C) / (#A + #T + #G + #C).
    size_t gc = count['C'_dna5.to_rank()] + count['G'_dna5.to_rank()];
    size_t atgc = input.size() - count['N'_dna5.to_rank()];
    float gc_content = 1.0f * gc / atgc;
 
    seqan3::debug_stream << "The GC content of " << sequence << " is " << 100 * gc_content << "%.\n";
    return 0;
}