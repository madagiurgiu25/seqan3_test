#include <string>                         // for std::string
#include <sstream>
#include <vector>
#include <cstring>  
#include <numeric>
#include <seqan3/argument_parser/all.hpp> // for argument_parser
#include <seqan3/core/debug_stream.hpp>   // for debug_stream

#include <seqan3/io/all.hpp>
#include <vector>                                         // for std::vector
#include <tuple>                                          // for std::tie
#include <seqan3/alignment/aligned_sequence/all.hpp>      // for alignment stream operator <<
#include <seqan3/alignment/configuration/all.hpp>         // for all configs in the seqan3::align_cfg namespace
#include <seqan3/alignment/pairwise/all.hpp>              // for seqan3::align_pairwise
#include <seqan3/alphabet/nucleotide/dna5.hpp>            // for dna5 datastrucutre
#include <seqan3/io/sequence_file/all.hpp>                // for sequence_file_input and sequence_file_output
#include <seqan3/std/filesystem>   
#include <seqan3/std/charconv>            // includes std::from_chars


int read_from_commandline(int argc, char * argv[]){
    // Create a buffer for the input

    std::string input{};
    std::string input2{};
    int age{30};
    bool verborse;

    // Initialise the Argument Parser and add an option.
    seqan3::argument_parser parser("My-SeqAn3-Tutorial", argc, argv);

    // parser.add_option(age, 'a', "user-age", "Please specify your age.");
    parser.add_positional_option(input, "Input");
    parser.add_positional_option(input2, "Output");
    parser.add_flag(verborse, 'v', "verborsse", "verborse description", seqan3::option_spec::required);
 
    try
    {
        // Parse the given arguments and catch possible errors.
        parser.parse();
        seqan3::debug_stream << age << std::endl << input << std::endl << input2 << std::endl << verborse << std::endl;

    }
    catch (seqan3::argument_parser_error const & ext)
    {
        seqan3::debug_stream << "[PARSER ERROR] " << ext.what() << '\n';
        return -1;
    }

    return 0;

}

void read_fasta(int argc, char * argv[]){

    std::string input{};

    seqan3::argument_parser parser("My-Fasta-Tutorial", argc, argv);
    parser.add_positional_option(input, "Input fasta file");

    try
    {
        // Parse the given arguments and catch possible errors.
        parser.parse();
        // seqan3::debug_stream << age << std::endl << input << std::endl << input2 << std::endl << verborse << std::endl;
        seqan3::debug_stream << input << std::endl;

        std::filesystem::path tmp_dir = std::filesystem::temp_directory_path(); // get the tmp directory
        // Initialise a file input object with a FastA file.
        seqan3::sequence_file_input file_in1{tmp_dir/"seq.fasta"};

        // read file using the arguments
        // Initialise a file input object with a FastA file.
        seqan3::sequence_file_input file_in2{input};
    
        // Retrieve the sequences and ids.
        seqan3::debug_stream << tmp_dir << std::endl;
        for (auto & record : file_in1)
        {
            seqan3::debug_stream << "ID:     " << record.id() << '\n';
            seqan3::debug_stream << "SEQ:    " << record.sequence() << '\n';
            // seqan3::debug_stream << "Empty Qual." << qual << '\n';  // qual is empty for FastAfiles
        }

        seqan3::debug_stream << std::endl << input << std::endl;
        std::vector<seqan3::dna5_vector> sequences;
        for (auto const & [seq, id, qual] : file_in2)
        {
            seqan3::debug_stream << "ID:     " << id << '\n';
            seqan3::debug_stream << "SEQ:    " << seq << '\n';
            seqan3::debug_stream << "Empty Qual." << qual << '\n';  // qual is empty for FastAfiles

            sequences.push_back(seq);

        }

        seqan3::debug_stream << sequences << '\n';


    }
    catch (seqan3::argument_parser_error const & ext)
    {
        seqan3::debug_stream << "[PARSER ERROR] " << ext.what() << '\n';
    }

}

void align_two_sequences(){

    using namespace seqan3::literals;
 
    auto tmp_dir = std::filesystem::temp_directory_path();
    std::string filename{tmp_dir/"seq_2.fasta"};
    {
        // Create a /tmp/seq.fasta file.
        seqan3::sequence_file_output file_out{filename};
 
        file_out.emplace_back("ACGTGATG"_dna5, std::string{"seq1"});
        file_out.emplace_back("AGTGATACT"_dna5, std::string{"seq2"});
    }
 
    // Initialise a file input object and a vector
    seqan3::sequence_file_input file_in{filename};
    std::vector<seqan3::dna5_vector> sequences;
 
    for (auto & record : file_in)
    {
        sequences.push_back(record.sequence());
    }
 
    // Call a global pairwise alignment with edit distance and traceback.
    for (auto && res : align_pairwise(std::tie(sequences[0], sequences[1]),
                                      seqan3::align_cfg::method_global{} |
                                      seqan3::align_cfg::edit_scheme |
                                      seqan3::align_cfg::output_alignment{} |
                                      seqan3::align_cfg::output_score{}))
    {
        // Print the resulting score and the alignment.
        seqan3::debug_stream << res.score() << '\n';      // => -4
        seqan3::debug_stream << res.alignment() << '\n';  // =>       0     .    :
                                                          //            ACGTGATG--
                                                          //            | |||||
                                                          //            A-GTGATACT
    }
    std::filesystem::remove(filename);

}

