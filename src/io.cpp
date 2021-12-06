#include <filesystem>
#include <string>
#include <sstream>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/io/sequence_file/all.hpp>  // for sequence_file_input and sequence_file_output

#include <numeric> // std::accumulate
#include <seqan3/std/ranges>
 
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/all.hpp>

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

void testFastq(){

    std::cout << "Current path: " << std::filesystem::current_path() << '\n';
    seqan3::sequence_file_input fin{std::filesystem::current_path() / "my.fastq"};
 
    // `&&` is important because seqan3::views::chunk returns temporaries!
    for (auto && records : fin | seqan3::views::chunk(10))
    {
        // `records` contains 10 elements (or less at the end)
        seqan3::debug_stream << "Taking the next 10 sequences:\n";
        seqan3::debug_stream << "ID:  " << (*records.begin()).id() << '\n'; // prints first ID in batch
    }
}

void filterFastq(std::string filename){

    seqan3::sequence_file_input fin{filename};

    // std::views::filter takes a function object (a lambda in this case) as input that returns a boolean
    auto minimum_quality_filter = std::views::filter([] (auto const & rec)
    {
        auto qualities = rec.base_qualities()
                       | std::views::transform([] (auto quality) { return seqan3::to_phred(quality); });
 
        auto sum = std::accumulate(qualities.begin(), qualities.end(), 0);
        return sum / std::ranges::size(qualities) >= 2; // minimum average quality >= 40
    });
 
    // // `&&` is important because seqan3::views::chunk returns temporaries!
    // for (auto && records : fin | seqan3::views::chunk(10) | minimum_quality_filter)
    // {
    //     // `records` contains 10 elements (or less at the end)
    //     seqan3::debug_stream << "Taking the next 10 sequences:\n";
    //     seqan3::debug_stream << "ID:  " << (records.id() << '\n'; // prints first ID in batch
    // }

    for (auto && rec : fin | minimum_quality_filter | seqan3::views::chunk(10))
    {
        seqan3::debug_stream << "ID: " << *(rec.begin()) << '\n';
    }
}

void filterSelect(std::string filename){
 
    seqan3::sequence_file_input fin{filename};
 
    auto length_filter = std::views::filter([] (auto const & rec)
    {
        return std::ranges::size(rec.sequence()) >= 5;
    });
 
    // Store all IDs into a vector:
    std::vector<std::string> ids{};
    for (auto & record : fin | length_filter | std::views::take(2))
    {
        ids.push_back(std::move(record.id()));
    }
 
    seqan3::debug_stream << ids << '\n';
}