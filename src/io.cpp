#include <filesystem>
#include <string>
#include <sstream>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/io/sequence_file/all.hpp>  // for sequence_file_input and sequence_file_output

#include "io.hpp"

#include <string>

enum bioinformatics_format {
    fasta,
    embl,
    sam,
    vienna,
    fastq,
    noformat
};

bioinformatics_format hashit (std::string const& inString) {
    if (inString == "fasta") return fasta;
    if (inString == "embl") return embl;
    if (inString == "sam") return sam;
    if (inString == "vienna") return vienna;
    if (inString == "fastq") return fastq;
    return noformat;
}

void readFromFile(std::string filename, std::string format){

    // auto filename = std::filesystem::current_path() / "my.fasta";

    std::cout << format << " " << hashit(format) << std::endl;

    seqan3::sequence_file_input fin_from_filename{filename};
 
    for (auto & record : fin_from_filename)
    {
        seqan3::debug_stream << "ID:  " << record.id() << '\n';
        seqan3::debug_stream << "SEQ: " << record.sequence() << '\n';
        // a quality field also exists, but is not printed, because we know it's empty for FASTA files.
    }
    
    // switch (hashit(format)) {
    //     case bioinformatics_format::fasta:
    //         seqan3::sequence_file_input fin_from_filename{filename};
 
    //         for (auto & record : fin_from_filename)
    //         {
    //             seqan3::debug_stream << "ID:  " << record.id() << '\n';
    //             seqan3::debug_stream << "SEQ: " << record.sequence() << '\n';
    //             // a quality field also exists, but is not printed, because we know it's empty for FASTA files.
    //         }
    //         break;
    //     case bioinformatics_format::embl:
    //         std::cout << "embl" << std::endl;
    //         break;
    //     case bioinformatics_format::sam:
    //         std::cout << "sam" << std::endl;
    //         break;
    //     case bioinformatics_format::vienna:
    //         std::cout << "vienna" << std::endl;
    //         break;
    //     case bioinformatics_format::fastq:
    //         std::cout << "fastq" << std::endl;
    //         break;
    //     default:
    //         break;
    // }
}

void testFasta(){
    
    auto input = R"(> TEST1
ACGT
> Test2
AGGCTGA
> Test3
GGAGTATAATATATATATATATAT)";

    seqan3::sequence_file_input fin{std::istringstream{input}, seqan3::format_fasta{}};
 
    for (auto & record : fin)
    {
        seqan3::debug_stream << "ID:  " << record.id() << '\n';
        seqan3::debug_stream << "SEQ: " << record.sequence() << '\n';
        // a quality field also exists, but is not printed, because we know it's empty for FASTA files.
    }
}

void testFastaVector(){
    
    auto input = R"(> TEST1
ACGT
> Test2
AGGCTGA
> Test3
GGAGTATAATATATATATATATAT)";

    seqan3::sequence_file_input fin{std::istringstream{input}, seqan3::format_fasta{}};
    using record_type = typename decltype(fin)::record_type;
    std::vector<record_type> entries{};
    // record_type rec = std::move(*fin.begin()); // avoid copying

    for (auto & record : fin)
    {
        entries.push_back(record);
        seqan3::debug_stream << "ID:  " << record.id() << '\n';
        seqan3::debug_stream << "SEQ: " << record.sequence() << '\n';
        // a quality field also exists, but is not printed, because we know it's empty for FASTA files.
    }

    for (auto & record: entries){
        seqan3::debug_stream <<  record << "\n";
    }

    seqan3::debug_stream << "backwards\n";
    // But you can also do this:
    std::ranges::copy(fin, std::cpp20::back_inserter(entries));
    seqan3::debug_stream << entries << '\n';
}

void testFastaVectorEfficient(){

       auto input = R"(> TEST1
ACGT
> Test2
AGGCTGA
> Test3
GGAGTATAATATATATATATATAT)";

    seqan3::sequence_file_input fin{std::istringstream{input}, seqan3::format_fasta{}};
    using record_type = typename decltype(fin)::record_type;
    std::vector<record_type> entries{};


    for (auto & record : fin)
    {
        entries.push_back(std::move(record));
    }
}

void testFastaVectorEfficient2(){

       auto input = R"(> TEST1
ACGT
> Test2
AGGCTGA
> Test3
GGAGTATAATATATATATATATAT)";

    seqan3::sequence_file_input fin{std::istringstream{input}, seqan3::format_fasta{}};
    using record_type = typename decltype(fin)::record_type;
    std::vector<record_type> entries{};


    for (auto & record : fin)
    {
        entries.push_back(record);
    }
}
