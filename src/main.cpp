#include <string>                         // for std::string
#include <sstream>
#include <vector>
#include <cstring>  
#include <numeric>
#include <seqan3/argument_parser/all.hpp> // for argument_parser
#include <seqan3/core/debug_stream.hpp>   // for debug_stream
// #include <seqan3/io/sequence_file/input.hpp> // for sequence_file_input
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

#include "starter.hpp"
#include "gccontent.hpp"
#include "ranges.hpp"

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

// std::vector<std::string> split(const std::string s, char delimiter)
// {
//    std::vector<std::string> tokens;
//    std::string token;
//    std::istringstream tokenStream(s);

// //    std::cout << "enter split function" << std::endl;

//    while (std::getline(tokenStream, token, delimiter))
//    {
//       tokens.push_back(token);
//     //   std::cout << token << std::endl;
//    }
//    return tokens;
// }

// void print_data(std::vector<std::vector<std::string>> &data, std::vector<std::string> headers){

//     std::cout << "Size data " << data.size() << " x " << data[0].size() << std::endl;
//     if (!headers.empty()){
//         // display header
//         for(auto const& col: headers){
//             std::cout << col << " ";
//         }
//     }

//     for(auto const& row: data) {
//         for(auto const& col: row){
//             std::cout << col << " ";
//         }
//         std::cout << std::endl;
//     }

// }

// template <typename number_type, typename range_type>
// number_type to_number(range_type && range)
// {
//     std::string str;
//     number_type num;
//     std::ranges::copy(range, std::cpp20::back_inserter(str));
//     auto res = std::from_chars(&str[0], &str[0] + str.size(), num);
 
//     if (res.ec != std::errc{})
//     {
//         seqan3::debug_stream << "Could not cast '" << range << "' to a valid number\n";
//         throw std::invalid_argument{"CAST ERROR"};
//     }
//     return num;
// }

// template <typename T>
// void aggregation_operator(const std::vector<T> &v, std::string aggr){
//     if (aggr == "median")
//             seqan3::debug_stream << ([&v] () { std::sort(v.begin(), v.end()); return v[v.size()/2]; })() << '\n';
//         else if (aggr == "mean")
//             seqan3::debug_stream << ([&v] () { double sum{}; for (auto i : v) sum += i; return sum / v.size(); })()
//                                  << '\n';
//         else
//             seqan3::debug_stream << "I do not know the aggregation method " << aggr << '\n';
// }

// void aggregate_data_by_seasons(std::vector<std::vector<std::string>> &data, std::vector<uint8_t> sn, std::string aggr){
//     // assume we know seasons are at column 0
//     std::vector<uint8_t> sv{};

//     for(auto const& row: data) {
//         auto it = std::next(row.begin()); 
//         if (std::find(sn.begin(), sn.end(), to_number<uint8_t>(*it)) != sn.end()){
//                 sv.push_back(to_number<double>(*std::next(it)));
//         }
//     }

//     seqan3::debug_stream << aggr << " views for seasons is: " << '\n';
//     aggregation_operator(sv, aggr);
  
// }

// void aggregate_data_by_year(std::vector<std::vector<std::string>> &data, uint32_t yr, std::string aggr){
//     // assume we know year are at column 3
//     std::vector<uint32_t> yv{};

//     for(auto const& row: data) {
//         auto it = std::next(row.begin(), 3); 
//         if (to_number<uint32_t>(*it) >= yr)
//                 yv.push_back(to_number<double>(*std::next(it)));
//     }

//     seqan3::debug_stream << aggr << " views for >=year is: " << '\n';
//     aggregation_operator(yv, aggr);

// }

// void initialise_argument_parser(seqan3::argument_parser & parser){
//     parser.info.app_name = "Cersei";
//     parser.info.short_description = "Aggregate average US. Game of Thrones viewers by season.";
//     parser.info.version = "1.0.0";
// }
 
// void run_program(std::filesystem::path & path, uint32_t yr, std::string & aggr_by, bool hd_is_set)
// {
//     std::ifstream file{path.string()};
 
//     if (file.is_open())
//     {
//         std::vector<double> v;
//         std::string line;
 
//         if (hd_is_set)
//             std::getline(file, line); // ignore first line
 
//         while (std::getline(file, line))
//         {
//             auto splitted_line = line | std::views::split('\t');
//             auto it = std::next(splitted_line.begin(), 3); // move to 4th column
 
//             if (to_number<uint32_t>(*it) >= yr)
//                 v.push_back(to_number<double>(*std::next(it)));
//         }
 
//         if (aggr_by == "median")
//             seqan3::debug_stream << ([&v] () { std::sort(v.begin(), v.end()); return v[v.size()/2]; })() << '\n';
//         else if (aggr_by == "mean")
//             seqan3::debug_stream << ([&v] () { double sum{}; for (auto i : v) sum += i; return sum / v.size(); })()
//                                  << '\n';
//         else
//             seqan3::debug_stream << "I do not know the aggregation method " << aggr_by << '\n';
//     }
//     else
//     {
//         seqan3::debug_stream << "Error: Cannot open file for reading.\n";
//     }
// }
// // -----------------------------------------------------------------------------

// struct cmd_arguments
// {
//     std::filesystem::path file_path{};
//     // std::string file_path{};
//     std::vector<uint8_t> seasons{};
//     std::uint32_t year{};
//     std::string aggregate_by{"mean"};
//     bool header{};
// };

// void parser_homework(int argc, char * argv[]){
    
//     cmd_arguments args;
//     seqan3::argument_parser myparser{"Game-of-Parsing", argc, argv};        // initialise myparser

//     initialise_argument_parser(myparser);
//     myparser.add_positional_option(args.file_path, "Input");
//     myparser.add_flag(args.header, 'H', "header", "File includes header");
//     myparser.add_option(args.seasons, 's', "season", "Choose the seasons to aggregate.");
//     myparser.add_option(args.year, 'y', "year", "Choose the year to aggregate.");
//     myparser.add_option(args.aggregate_by, 'a', "aggregate-by", "Choose your method of aggregation: mean or median.");
 
//     try
//     {
//         myparser.parse();                                                  // trigger command line parsing
//         seqan3::debug_stream << args.header << " " << args.file_path << std::endl;

//         // method from the course
//         run_program(args.file_path, args.year, args.aggregate_by, args.header);

//         // my method file
//         std::ifstream myfile(args.file_path.string());
//         std::vector<std::string> headers{};
//         std::vector<std::vector<std::string>> data{};
//         char delimiter = '\t';

//         if ( myfile.is_open() ) { 
//             std::string myline;
//             while ( myfile ) {

//                 std::getline (myfile, myline);
//                 std::vector<std::string> row = split(myline, delimiter);
                
//                 if (headers.empty()){
//                     if (args.header == 1){
//                         headers = row;
//                     }else{
//                         data.push_back(row);
//                     }
//                 }else{
//                     if (row.size() > 0){ // avoid adding last line
//                         data.push_back(row);
//                     }
                    
//                 }
//             }
//         }

//         // display data
//         print_data(data, headers);

//         // compute seasons and year
//         // aggregate_data_by_seasons(data, args.seasons, args.aggregate_by);
//         // aggregate_data_by_year(data, args.year, args.aggregate_by);


//     }
//     catch (seqan3::argument_parser_error const & ext)                     // catch user errors
//     {
//         seqan3::debug_stream << "[Winter has come] " << ext.what() << "\n"; // customise your error message
//     }
// }



int main(int argc, char * argv[])
{
    // 1. Read from command line
    // read_from_commandline(argc, argv);

    // 2. Read fasta seq
    // read_fasta(argc, argv);

    // 3. Align sequences
    // align_two_sequences();

    // 3. Parse US data
    // parser_homework(argc, argv);

    // 4. GC content
    // readSequence(argc, argv);

    // 5. Ranges
    // testRanges();
    // complementSequence(argc, argv);
    bitpackedSequence(argc, argv);


    
    return 0;
}