/**
 * identify_barcode.cpp
 *
 * C++ port of scripts/python/identify_barcodes.py for the SPIDR pipeline.
 * Identifies bead IDs (DPM) from read 1 and round barcodes (Y/ODD/EVEN)
 * from read 2, appending the barcode string to each FASTQ header.
 *
 * Reads/writes gzip-compressed FASTQ files using zlib.
 *
 * Usage (mirrors Python CLI exactly):
 *   identify_barcode \
 *     --input_read1  <r1.fq.gz>  --input_read2  <r2.fq.gz> \
 *     --output_read1 <out1.fq.gz> --output_read2 <out2.fq.gz> \
 *     --read1_format 'DPM' \
 *     --read2_format 'Y|SPACER|ODD|SPACER|EVEN|SPACER|ODD|SPACER|EVEN|SPACER|ODD' \
 *     --read1_start_offset 0  --read2_start_offset 0 \
 *     --config config_6_rounds_mTOR.txt
 */

#include <algorithm>
#include <cstring>
#include <fstream>
#include <iostream>
#include <queue>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <zlib.h>

// ---------------------------------------------------------------------------
// Constants – must match the Python script
// ---------------------------------------------------------------------------
static constexpr int SPACER         = 7;
static constexpr int LAXITY         = 6;
static constexpr int BARCODE_LENGTH = 17;
static const std::string NOT_FOUND  = "NOT_FOUND";

// ---------------------------------------------------------------------------
// GzipReader – line-oriented reader for .fq.gz files
// ---------------------------------------------------------------------------
class GzipReader {
public:
    explicit GzipReader(const std::string& path) {
        file_ = gzopen(path.c_str(), "rb");
        if (!file_)
            throw std::runtime_error("Cannot open for reading: " + path);
        gzbuffer(file_, 1 << 20); // 1 MiB read buffer
    }
    ~GzipReader() { if (file_) gzclose(file_); }
    GzipReader(const GzipReader&)            = delete;
    GzipReader& operator=(const GzipReader&) = delete;

    // Returns false on EOF. Strips trailing \r\n.
    bool readline(std::string& line) {
        if (!gzgets(file_, buf_, sizeof(buf_)))
            return false;
        line = buf_;
        while (!line.empty() && (line.back() == '\n' || line.back() == '\r'))
            line.pop_back();
        return true;
    }

private:
    gzFile file_ = nullptr;
    char   buf_[1 << 17]; // 128 KiB – large enough for any FASTQ line
};

// ---------------------------------------------------------------------------
// GzipWriter – buffered gzip writer for .fq.gz files
// ---------------------------------------------------------------------------
class GzipWriter {
public:
    explicit GzipWriter(const std::string& path) {
        file_ = gzopen(path.c_str(), "wb");
        if (!file_)
            throw std::runtime_error("Cannot open for writing: " + path);
        gzbuffer(file_, 1 << 20);
    }
    ~GzipWriter() { if (file_) gzclose(file_); }
    GzipWriter(const GzipWriter&)            = delete;
    GzipWriter& operator=(const GzipWriter&) = delete;

    void write(const std::string& s) {
        gzwrite(file_, s.data(), static_cast<unsigned>(s.size()));
    }

private:
    gzFile file_ = nullptr;
};

// ---------------------------------------------------------------------------
// Hamming-distance neighbourhood
//
// Replicates the Python BFS in identify_barcodes.py:
//   get_all_seqs_within_k(seq, k)
//
// Important subtlety: the BFS explores paths of exactly k mutation *steps*,
// not Hamming distance k.  When k≥2 the same position may be mutated more
// than once, so the returned set includes sequences at actual Hamming
// distances 0…k.  The caller also adds the original sequence separately
// (matching: hamming_hashmap[name].add(seq)).
// ---------------------------------------------------------------------------
static std::unordered_set<std::string>
get_all_seqs_within_k(const std::string& seq, int k) {
    std::unordered_set<std::string> result;

    struct State { std::string seq; int dist; };
    std::queue<State> q;
    q.push({seq, 0});

    while (!q.empty()) {
        auto [cur, dist] = std::move(q.front());
        q.pop();

        if (dist == k) {
            result.insert(cur);
            continue;
        }
        for (size_t pos = 0; pos < cur.size(); ++pos) {
            for (char nuc : {'A', 'C', 'G', 'T'}) {
                if (nuc != cur[pos]) {
                    std::string mut = cur;
                    mut[pos] = nuc;
                    q.push({std::move(mut), dist + 1});
                }
            }
        }
    }
    return result;
}

// Build an inverted lookup: mutated_sequence → barcode_name
// entries = [(name, sequence, tolerance), ...]
using LookupMap = std::unordered_map<std::string, std::string>;

static LookupMap build_lookup(
        const std::vector<std::tuple<std::string, std::string, int>>& entries) {
    LookupMap lookup;
    for (const auto& [name, seq, tol] : entries) {
        auto neighbours = get_all_seqs_within_k(seq, tol);
        neighbours.insert(seq); // Python also does: hamming_hashmap[name].add(seq)
        for (const auto& nb : neighbours)
            lookup[nb] = name;
    }
    return lookup;
}

// ---------------------------------------------------------------------------
// Config file parsing
//
// The file has 3 header lines (READ1=…, READ2=…, blank) then TSV rows:
//   type  name  sequence  tolerance
// This mirrors pandas read_csv(..., skiprows=3, usecols=[0,1,2,3]).
// ---------------------------------------------------------------------------
struct ConfigEntry {
    std::string type, name, sequence;
    int tolerance;
};

static std::vector<ConfigEntry> parse_config(const std::string& path) {
    std::ifstream f(path);
    if (!f) throw std::runtime_error("Cannot open config: " + path);

    std::vector<ConfigEntry> entries;
    std::string line;
    int line_num = 0;
    while (std::getline(f, line)) {
        if (++line_num <= 3) continue; // skip READ1, READ2, blank
        if (line.empty())   continue;

        std::istringstream ss(line);
        ConfigEntry e;
        std::string tol_str;
        if (std::getline(ss, e.type, '\t') &&
            std::getline(ss, e.name, '\t') &&
            std::getline(ss, e.sequence, '\t') &&
            std::getline(ss, tol_str))
        {
            e.tolerance = std::stoi(tol_str);
            entries.push_back(std::move(e));
        }
    }
    return entries;
}

// ---------------------------------------------------------------------------
// Split a '|'-delimited format string, trimming whitespace from each token
// ---------------------------------------------------------------------------
static std::vector<std::string> split_format(const std::string& s) {
    std::vector<std::string> result;
    std::istringstream ss(s);
    std::string tok;
    while (std::getline(ss, tok, '|')) {
        size_t a = tok.find_first_not_of(" \t");
        size_t b = tok.find_last_not_of(" \t");
        if (a != std::string::npos)
            result.push_back(tok.substr(a, b - a + 1));
    }
    return result;
}

// ---------------------------------------------------------------------------
// find_bead_id  (read 1 / DPM barcodes – exact match only, tolerance=0)
// ---------------------------------------------------------------------------
static std::string find_bead_id(
        const std::string& read,
        int start_offset,
        const std::vector<int>& bead_lengths,   // sorted unique lengths
        const LookupMap&        bead_hashmap)    // seq → name
{
    for (int len : bead_lengths) {
        if (start_offset + len > (int)read.size()) continue;
        std::string candidate = read.substr(start_offset, len);
        auto it = bead_hashmap.find(candidate);
        if (it != bead_hashmap.end())
            return it->second;
    }
    return NOT_FOUND;
}

// ---------------------------------------------------------------------------
// find_terminal_barcode  (Y barcodes at the 5' end of read 2)
//
// Mirrors Python find_terminal_barcode().  Returns {name, new_start}.
// new_start = start + term_size of the match (or last iterated length on
// NOT_FOUND, matching Python's post-loop value of `term_size`).
// ---------------------------------------------------------------------------
static std::pair<std::string, int> find_terminal_barcode(
        const std::string& read,
        const LookupMap&   term_lookup,
        int                start,
        const std::vector<int>& possible_lengths) // sorted
{
    int last_size = possible_lengths.empty() ? 0 : possible_lengths.back();
    for (int term_size : possible_lengths) {
        last_size = term_size;
        if (start + term_size > (int)read.size()) continue;
        auto it = term_lookup.find(read.substr(start, term_size));
        if (it != term_lookup.end())
            return {it->second, start + term_size};
    }
    return {NOT_FOUND, start + last_size};
}

// ---------------------------------------------------------------------------
// find_nonterminal_barcode  (ODD / EVEN barcodes in read 2)
//
// Mirrors Python find_nonterminal_barcode().
// Searches offsets [0, LAXITY] and all possible barcode sizes.
// Returns {name, start + BARCODE_LENGTH} on success or NOT_FOUND.
// ---------------------------------------------------------------------------
static std::pair<std::string, int> find_nonterminal_barcode(
        const std::string& read,
        const LookupMap&   lookup,
        int                start,
        const std::vector<int>& possible_lengths) // sorted ascending
{
    const int max_bc_size = possible_lengths.back();

    for (int offset = 0; offset <= LAXITY; ++offset) {
        const int window_start = start + offset;
        for (int bc_size : possible_lengths) {
            if (window_start + bc_size > (int)read.size())
                return {NOT_FOUND, start + BARCODE_LENGTH};

            auto it = lookup.find(read.substr(window_start, bc_size));
            if (it != lookup.end())
                return {it->second, start + BARCODE_LENGTH};

            if (bc_size == max_bc_size && offset == LAXITY)
                return {NOT_FOUND, start + BARCODE_LENGTH};
        }
    }
    return {NOT_FOUND, start + BARCODE_LENGTH}; // unreachable in practice
}

// ---------------------------------------------------------------------------
// find_barcodes  (orchestrate all barcode finding on read 2)
//
// Mirrors Python find_barcodes() including the layout-iteration logic:
//   layout = read2_format + ["END"]
//   barcode_type = layout.pop(0)
//   while len(layout) > 0:
//       ...
//       barcode_type = layout.pop(0)
// ---------------------------------------------------------------------------
static std::vector<std::string> find_barcodes(
        const std::string& read,
        const LookupMap&   even_lookup,
        const LookupMap&   odd_lookup,
        const LookupMap&   term_lookup,
        const std::vector<std::string>& read2_format,
        int read2_start_offset,
        const std::vector<int>& term_lengths,
        const std::vector<int>& odd_lengths,
        const std::vector<int>& even_lengths)
{
    int start = read2_start_offset;
    std::vector<std::string> barcodes;

    std::vector<std::string> layout = read2_format;
    layout.push_back("END");

    // Pop the first element (idx advances to 1); loop while idx < layout.size()
    // i.e. while layout still has remaining elements after the current one.
    size_t idx = 1;
    const std::string* barcode_type = &layout[0];

    while (idx < layout.size()) {
        if (*barcode_type == "Y") {
            auto [bc, ns] = find_terminal_barcode(read, term_lookup, start, term_lengths);
            barcodes.push_back(std::move(bc));
            start = ns;
        } else if (*barcode_type == "SPACER") {
            start += SPACER;
        } else if (*barcode_type == "END") {
            break; // safety – normally reached by loop condition
        } else if (*barcode_type == "ODD") {
            auto [bc, ns] = find_nonterminal_barcode(read, odd_lookup, start, odd_lengths);
            barcodes.push_back(std::move(bc));
            start = ns;
        } else if (*barcode_type == "EVEN") {
            auto [bc, ns] = find_nonterminal_barcode(read, even_lookup, start, even_lengths);
            barcodes.push_back(std::move(bc));
            start = ns;
        } else {
            throw std::runtime_error("Invalid barcode type: " + *barcode_type);
        }

        // If fewer than BARCODE_LENGTH bases remain, no more barcodes are possible
        if ((int)read.size() - start < BARCODE_LENGTH)
            break;

        barcode_type = &layout[idx++];
    }
    return barcodes;
}

// ---------------------------------------------------------------------------
// pad_barcodes  – pad with NOT_FOUND to reach expected_length
// ---------------------------------------------------------------------------
static std::vector<std::string> pad_barcodes(
        std::vector<std::string> barcodes, int expected_length)
{
    if ((int)barcodes.size() > expected_length)
        throw std::runtime_error("More barcodes found than expected");
    while ((int)barcodes.size() < expected_length)
        barcodes.push_back(NOT_FOUND);
    return barcodes;
}

// ---------------------------------------------------------------------------
// CLI argument parsing
// ---------------------------------------------------------------------------
struct Args {
    std::string input_read1, input_read2;
    std::string output_read1, output_read2;
    std::string read1_format, read2_format;
    int read1_start_offset = 0;
    int read2_start_offset = 0;
    std::string config;
};

static void usage(const char* prog) {
    std::cerr
        << "Usage: " << prog << "\n"
        << "  --input_read1 <r1.fq.gz>    --input_read2 <r2.fq.gz>\n"
        << "  --output_read1 <out1.fq.gz> --output_read2 <out2.fq.gz>\n"
        << "  --read1_format 'DPM'\n"
        << "  --read2_format 'Y|SPACER|ODD|SPACER|EVEN|...'\n"
        << "  --read1_start_offset <int>   (default 0)\n"
        << "  --read2_start_offset <int>   (default 0)\n"
        << "  --config <config.txt>\n";
}

static Args parse_args(int argc, char* argv[]) {
    Args a;
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        auto next = [&]() -> std::string {
            if (i + 1 >= argc) {
                std::cerr << "Missing value for " << arg << "\n";
                std::exit(1);
            }
            return argv[++i];
        };
        if      (arg == "--input_read1")        a.input_read1        = next();
        else if (arg == "--input_read2")        a.input_read2        = next();
        else if (arg == "--output_read1")       a.output_read1       = next();
        else if (arg == "--output_read2")       a.output_read2       = next();
        else if (arg == "--read1_format")       a.read1_format       = next();
        else if (arg == "--read2_format")       a.read2_format       = next();
        else if (arg == "--read1_start_offset") a.read1_start_offset = std::stoi(next());
        else if (arg == "--read2_start_offset") a.read2_start_offset = std::stoi(next());
        else if (arg == "--config")             a.config             = next();
        else if (arg == "--show_progress_bar")  { next(); /* ignored */ }
        else {
            std::cerr << "Unknown argument: " << arg << "\n";
            usage(argv[0]);
            std::exit(1);
        }
    }
    if (a.input_read1.empty() || a.input_read2.empty() ||
        a.output_read1.empty() || a.output_read2.empty() ||
        a.read1_format.empty() || a.read2_format.empty() ||
        a.config.empty()) {
        std::cerr << "Missing required argument(s).\n";
        usage(argv[0]);
        std::exit(1);
    }
    return a;
}

// ---------------------------------------------------------------------------
// main
// ---------------------------------------------------------------------------
int main(int argc, char* argv[]) {
    Args args = parse_args(argc, argv);

    // Parse format strings (e.g. "Y|SPACER|ODD" → ["Y","SPACER","ODD"])
    auto read1_format = split_format(args.read1_format);
    auto read2_format = split_format(args.read2_format);

    int num_r1_barcodes = 0;
    for (const auto& b : read1_format) if (b != "SPACER") ++num_r1_barcodes;
    int num_r2_barcodes = 0;
    for (const auto& b : read2_format) if (b != "SPACER") ++num_r2_barcodes;
    const int total_barcode_length = num_r1_barcodes + num_r2_barcodes;

    // ---- Parse config ----
    std::cerr << "Parsing config..." << std::endl;
    auto config_entries = parse_config(args.config);

    // Separate entries by type
    LookupMap bead_hashmap; // exact: seq → name (no Hamming expansion for DPM)
    std::vector<std::tuple<std::string,std::string,int>> odd_entries, even_entries, term_entries;
    std::unordered_set<int> bead_len_set;
    std::vector<int>        bead_lengths;

    for (const auto& e : config_entries) {
        if (e.type == "DPM") {
            bead_hashmap[e.sequence] = e.name;
            if (bead_len_set.insert((int)e.sequence.size()).second)
                bead_lengths.push_back((int)e.sequence.size());
        } else if (e.type == "ODD") {
            odd_entries.emplace_back(e.name, e.sequence, e.tolerance);
        } else if (e.type == "EVEN") {
            even_entries.emplace_back(e.name, e.sequence, e.tolerance);
        } else if (e.type == "Y") {
            term_entries.emplace_back(e.name, e.sequence, e.tolerance);
        }
    }
    std::sort(bead_lengths.begin(), bead_lengths.end());

    // ---- Build Hamming lookup maps ----
    std::cerr << "Building Hamming hashmaps (this may take a few seconds)..." << std::endl;
    auto odd_lookup  = build_lookup(odd_entries);
    auto even_lookup = build_lookup(even_entries);
    auto term_lookup = build_lookup(term_entries);

    // Collect unique barcode lengths for each type (sorted ascending)
    auto unique_lengths = [](const std::vector<std::tuple<std::string,std::string,int>>& v) {
        std::unordered_set<int> seen;
        std::vector<int> lens;
        for (const auto& [n, s, t] : v) {
            if (seen.insert((int)s.size()).second)
                lens.push_back((int)s.size());
        }
        std::sort(lens.begin(), lens.end());
        return lens;
    };

    auto odd_lengths  = unique_lengths(odd_entries);
    // NOTE: mirrors the Python code which uses odd_df for both odd and even lengths
    auto even_lengths = unique_lengths(odd_entries);
    auto term_lengths = unique_lengths(term_entries);

    // ---- Open I/O ----
    std::cerr << "Processing reads..." << std::endl;
    GzipReader r1_in(args.input_read1);
    GzipReader r2_in(args.input_read2);
    GzipWriter r1_out(args.output_read1);
    GzipWriter r2_out(args.output_read2);

    std::string h1, seq1, plus1, qual1;
    std::string h2, seq2, plus2, qual2;
    long long counter = 0;

    while (r1_in.readline(h1)    && r1_in.readline(seq1)  &&
           r1_in.readline(plus1) && r1_in.readline(qual1) &&
           r2_in.readline(h2)    && r2_in.readline(seq2)  &&
           r2_in.readline(plus2) && r2_in.readline(qual2))
    {
        ++counter;
        if (counter % 100000 == 0)
            std::cerr << "Processed " << counter << " reads" << std::endl;

        // Strip space and everything after (matching Python header_r1.split(" ")[0])
        {
            auto p = h1.find(' ');
            if (p != std::string::npos) h1.resize(p);
            p = h2.find(' ');
            if (p != std::string::npos) h2.resize(p);
        }

        // Bead ID from read 1
        std::string bead_id = find_bead_id(
            seq1, args.read1_start_offset, bead_lengths, bead_hashmap);

        // Round barcodes from read 2
        auto barcodes = find_barcodes(
            seq2, even_lookup, odd_lookup, term_lookup,
            read2_format, args.read2_start_offset,
            term_lengths, odd_lengths, even_lengths);

        // Combine bead id + round barcodes, pad to total_barcode_length
        std::vector<std::string> all_barcodes;
        all_barcodes.reserve(total_barcode_length);
        all_barcodes.push_back(std::move(bead_id));
        for (auto& bc : barcodes) all_barcodes.push_back(std::move(bc));
        all_barcodes = pad_barcodes(std::move(all_barcodes), total_barcode_length);

        // Build barcode string: [BC1][BC2]...[BCN]
        std::string barcode_str;
        barcode_str.reserve(total_barcode_length * 20);
        for (const auto& bc : all_barcodes) {
            barcode_str += '[';
            barcode_str += bc;
            barcode_str += ']';
        }

        // Write output FASTQ records
        // Format: header::barcode_str \n seq \n + \n qual \n
        r1_out.write(h1 + "::" + barcode_str + "\n" + seq1 + "\n" + plus1 + "\n" + qual1 + "\n");
        r2_out.write(h2 + "::" + barcode_str + "\n" + seq2 + "\n" + plus2 + "\n" + qual2 + "\n");
    }

    std::cerr << "Done. Processed " << counter << " reads total." << std::endl;
    return 0;
}
