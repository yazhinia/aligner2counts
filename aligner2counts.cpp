#include <iostream>
#include <stdio.h>
#include <string>
#include <cmath>
#include <fstream>
#include <cstring>
#include <sstream>
#include <vector>
#include <chrono>
#include <algorithm>
#include <functional>
#include <filesystem>
#include <unordered_map>


const unsigned int FLAG_0 = 0;      // 0 primary alignment
const unsigned int FLAG_1 = 1;      // 0x1 template having multiple segments in sequencing
const unsigned int PAIRED_FLAG = 2;      // 0x2 each segment properly aligned according to the aligner (paired)
const unsigned int FLAG_3 = 4;      // 0x4 segment unmapped
const unsigned int FLAG_4 = 8;      // 0x8 next segment in the template unmapped
const unsigned int FLAG_5 = 16;     // 0x10 SEQ being reverse complemented
const unsigned int FLAG_6 = 32;     // 0x20 SEQ of the next segment in the template being reverse complemented
const unsigned int FORWARD_FLAG = 64;     // 0x40 the first segment in the template (forward read R1)
const unsigned int REVERSE_FLAG = 128;    // 0x80 the last segment in the template (reverse read R2)
const unsigned int FLAG_9 = 256;    // 0x100 secondary alignment
const unsigned int FLAG_10 = 512;   // 0x200 not passing filters, such as platform/vendor quality controls
const unsigned int FLAG_11 = 1024;  // 0x400 PCR or optical duplicate
const unsigned int FLAG_12 = 2048;  // 0x800 supplementary alignment

struct al_data{
    std::string read;
    std::string contig;
    unsigned int pair_dir; 
    bool paired;
    float sequence_identity;
};

struct PairHash {
    template <class T1, class T2>
    std::size_t operator () (const std::pair<T1, T2>& p) const {
        auto h1 = std::hash<T1>{}(p.first);
        auto h2 = std::hash<T2>{}(p.second);

        // A simple combination function
        return h1 ^ h2;
    }
};

using AlnParser = std::vector<al_data>;
using AlnStat = std::pair<float, float>;
using ContigsMap = std::unordered_map<std::string, float>;
using ContigLinks = std::unordered_map<std::pair<std::string, std::string>, float, PairHash>;

void addpairlinks(
    ContigLinks& contiglinks,
    const std::string& c1,
    const std::string& c2,
    float frac_contigs_mapped) {
    
    std::string first = (c1 < c2) ? c1 : c2;
    std::string second = (c1 < c2) ? c2 : c1;
<<<<<<< HEAD
    std::pair<std::string, std::string> linkPair = std::make_pair(first, second);
    contiglinks[linkPair] += frac_contigs_mapped;
=======
    contiglinks[std::make_pair(first, second)] += frac_contigs_mapped;
>>>>>>> 2debd9a (update)
}

void construct_map(
    ContigsMap& contigs_map,
    std::string line,
    unsigned int minlength) {

    float init_count = 0.0f;
    if (line.compare(0, 3, "@SQ") == 0) {
        std::string header, id, length;
        std::istringstream iss(line);
        iss >> header >> id >> length;
        
        // Extract contig id and its length
        std::string contig_id = id.substr(3);
        unsigned int contig_length = std::stoi(length.substr(3));

        // Store only if contig length is at least minlength
        if (contig_length >= minlength) {
            contigs_map[contig_id] =init_count;
        }
    }
}

std::pair<int, int> get_alnpos(std::string &cigar_str) {
    int start, end, pos1, pos2 = 0;
    char op1, op2;
    // bool softclip_flag = false;

    if (sscanf(cigar_str.c_str(), "%d%c", &pos1, &op1) == 2) {
        // eg: 150M
        start = 0;
        end = pos1;
        // softclip_flag = true;
    } else {
        if (sscanf(cigar_str.c_str(), "%d%c%d%c", &pos1, &op1, &pos2, &op2) == 4) {
            // eg: 140M10S
            if(op1 == 'M') {
                start = 0;
                end = pos1;
                // if(op2 == 'S') {
                //     softclip_flag = true;
                // }
            } else {
                // eg: 10S140M
                start = pos1;
                end = start + pos2;
                // if(op1 == 'S') {
                //     softclip_flag = true;
                // }
            }
        }
    }
    // } else {
    //     std::cout << "should not enter this part where cigar does not have perfect or soft clip\n";
    //     softclip_flag = true;
    // }

    return std::make_pair(start, end);
}


std::pair<float, float> get_seqid_alncov(std::pair<int, int> &alnpos, std::string &qual_str, std::string &md_str) {
    
    unsigned int matches = 0;
    unsigned int mismatches = 0;
    std::string mismatch_string;
    std::istringstream ss(md_str);
    unsigned int future_matches = 0;
    ss >> future_matches;

    int start = alnpos.first;
    int end = alnpos.second;

    // int start = std::get<0>(alnpos);
    // int end = std::get<1>(alnpos);
    // if(!std::get<2>(alnpos)) {
    //     // Alignment starts from the beginning
    //     // It means hard clip
    //     // In the current version we ignore hard clip alignment
    //     // Hence, it this condition is not necessary
    //     start = 0;
    //     end = end - start;
    // }

    unsigned int alignment_length = end - start;
    std::cout << start << " " << end << " start and end position \n";
    //  Sequence identity only considers aligned region (soft clip region is ignored)
    //  Alignment coverage is calculated w.r.t full read length
    for (unsigned int i = start; i < end; i++) {
        bool is_match;
        if (future_matches > 0) {
            future_matches--;
            is_match = true;
        } else {
            is_match = false;
            char tt;
            ss >> tt;
            int t_match;
            if(ss >> t_match) {
                future_matches = t_match;
            } else {
                //do nothing
            }
        }
        // Quality check per nt position
        // 33 is the ASCII offset
        if (qual_str[i] >= 20 + 33) {
            if(is_match){
                matches++;
            } else {
                mismatches++;
            }
        }
    }
    float seq_id, alignment_coverage;
    std::cout << matches << " " << mismatches << " " << matches+mismatches << " matches and mismatches \n";
    if (! (matches+mismatches) > 0) {
        std::cerr << "Zero matches and mismatched position. Something wrong with the alignment\n";
        exit(1);
    }
    std::cout << qual_str.length() << " quality string length \n";
    if (!qual_str.length() > 0) {
        std::cerr << "Zero quality string length. Something wrong with the alignment \n";
        exit(1);
    }

    seq_id = (static_cast<float>(matches) * 100) / (matches+mismatches);
    alignment_coverage = (static_cast<float>(alignment_length) * 100) / (qual_str.length() + start); // account for variable read length
    // std::cout << matches << " matches " << mismatches << " mismatches " <<  seq_id << " " << alignment_coverage << " " << matches+mismatches << " in alnstats \n";
    return std::make_pair(seq_id, alignment_coverage);
}

void storealnstats(AlnParser &parsedaln, AlnStat &alnstats, std::string &currentread_id, std::string &contig_id, unsigned int &pair_dir, bool & proper_pair) {
    auto ic = std::find_if(parsedaln.rbegin(), parsedaln.rend(),
                [&](const al_data& a) {return a.contig == contig_id;});
    if (ic != parsedaln.rend()){ 
        // if mate pair alignment is found
        if (ic->pair_dir != pair_dir){
            float new_identity = (ic->sequence_identity + alnstats.first) / 2.0f;
            if (ic->sequence_identity < new_identity) {
                ic->sequence_identity = new_identity;
            }
            ic->paired = true;
        } else { // update sequence identity if the same read (of the same direction) mapped to the same contig with higher sequence identity
            if (ic->sequence_identity < alnstats.first) {
                ic->sequence_identity = alnstats.first;
            }
        }
    } else { // read mapped to different contig 
        parsedaln.push_back({currentread_id, contig_id, pair_dir, proper_pair, alnstats.first});
    }
}

void counting(AlnParser &parsedaln, ContigsMap &contigs_map, ContigsMap &unique_map, ContigsMap &cross_map, ContigLinks &contig_links, float &sequenceidentity, std::ofstream &eachcount) {
    if (parsedaln.empty()) return;
    auto max_it = std::max_element(parsedaln.begin(), parsedaln.end(),
        [](const al_data& a,const al_data& b)
        { return a.sequence_identity < b.sequence_identity;});

    // for (size_t r = 0; r < parsedaln.size(); r++) {
    //     std::cout << parsedaln[r].read << " " << parsedaln[r].contig << " " << parsedaln[r].pair_dir << " " << parsedaln[r].paired << " " << parsedaln[r].sequence_identity << " complete set \n";
    //     // eachcount << parsedaln[r].read << " " << parsedaln[r].contig << " " << parsedaln[r].pair_dir << " " << parsedaln[r].paired << " " << parsedaln.size() << " " << parsedaln[r].sequence_identity << "\n";
    // }

<<<<<<< HEAD
    if (it->sequence_identity >= sequenceidentity) {
        parsedaln.erase(std::remove_if(parsedaln.begin(),parsedaln.end(), [&](const al_data& a) { return it->sequence_identity > a.sequence_identity;}), parsedaln.end());
=======
    if (max_it->sequence_identity >= sequenceidentity) {
        parsedaln.erase(
            std::remove_if(parsedaln.begin(),parsedaln.end(), [&](const al_data& a) { 
                return max_it->sequence_identity > a.sequence_identity;}),
            parsedaln.end()
            );
>>>>>>> 2debd9a (update)

        if (parsedaln.size() > 0) {
            // Check if there are any elements with paired == true
            bool has_paired = std::any_of(parsedaln.begin(), parsedaln.end(), [](const al_data& a) {
                return a.paired;
            });

            // Keep only paired end read mapping (paired == true) if present and remove single read mapping
            if (has_paired) {
                parsedaln.erase(
                    std::remove_if(parsedaln.begin(), parsedaln.end(), [](const al_data& a) {
                    return !a.paired;}),
                    parsedaln.end()
                );
            }
            // Keep only the first 10 elements
            if (parsedaln.size() > 10) {
                parsedaln.resize(10);
            }

            float frac_contigs_mapped = 1.0f / parsedaln.size();
<<<<<<< HEAD
            float value_unpaired = 0.0f, value_added_to_paired = 0.0f, value_added_per_paired = 0.0f;
            if (paired_count > 0) {
                value_unpaired = frac_contigs_mapped / 2.0f;
                value_added_to_paired = value_unpaired * non_paired_count;
                value_added_per_paired = value_added_to_paired / paired_count;
            } else {
                value_unpaired = frac_contigs_mapped; // likely discordant alignment
            }

            std::vector<std::string> contiglist;
            for (size_t r = 0; r < parsedaln.size(); r++) {
                contiglist.push_back(parsedaln[r].contig);
                auto it = contigs_map.find(parsedaln[r].contig);
                if (parsedaln[r].paired) {
                    it->second = it->second + frac_contigs_mapped + value_added_per_paired;
                    // std::cout << parsedaln[r].read << " " << parsedaln[r].contig << " " << parsedaln[r].pair_dir << " " << parsedaln[r].paired << " " << parsedaln[r].sequence_identity << " " << frac_contigs_mapped + value_added_per_paired << " after filter set \n";
                    eachcount << parsedaln[r].read << " " << parsedaln[r].contig << " " << parsedaln[r].pair_dir << " " << parsedaln[r].paired << " " << parsedaln.size() << " " << frac_contigs_mapped + value_added_per_paired << "\n";
                } else {
                    it->second = it->second + value_unpaired;
                    // std::cout << parsedaln[r].read << " " << parsedaln[r].contig << " " << parsedaln[r].pair_dir << " " << parsedaln[r].paired << " " << parsedaln[r].sequence_identity << " " << value_unpaired << " after filter set \n";
                    eachcount << parsedaln[r].read << " " << parsedaln[r].contig << " " << parsedaln[r].pair_dir << " " << parsedaln[r].paired << " " << parsedaln.size() << " " << value_unpaired << "\n";
                }
            }

            // add shared reads links to contigs
            for (std::size_t i = 0; i < contiglist.size(); ++i) {
                for (std::size_t j = i + 1; j < contiglist.size(); ++j) {
                    auto it1 = contigs_map.find(contiglist[i]);
                    auto it2 = contigs_map.find(contiglist[j]);
                    addpairlinks(contig_links, contiglist[i], contiglist[j], frac_contigs_mapped);
=======

            // Update counts in contigs_map and unique_map or cross_map
            if (parsedaln.size() == 1) {
                contigs_map[parsedaln[0].contig] += 1.0f;
                unique_map[parsedaln[0].contig] += 1.0f;
            } else {    
                std::vector<std::string> contiglist;
                for (const auto &aln: parsedaln) {
                    contiglist.push_back(aln.contig);
                    contigs_map[aln.contig] += frac_contigs_mapped;
                    cross_map[aln.contig] += frac_contigs_mapped;
                    eachcount << aln.read << " " << aln.contig << " " << aln.pair_dir << " " << aln.paired << " " << parsedaln.size() << " " << aln.sequence_identity << "\n";
                }
                // add shared reads links to contigs
                for (std::size_t i = 0; i < contiglist.size(); ++i) {
                    for (std::size_t j = i + 1; j < contiglist.size(); ++j) {;
                        addpairlinks(contig_links, contiglist[i], contiglist[j], frac_contigs_mapped);
                    }
>>>>>>> 2debd9a (update)
                }
                // clear contiglist data of current read alignments
                contiglist.clear();
            }
        }
    } else {
        parsedaln.clear();
    }
}

void fractionate_countlinks(ContigLinks &contiglinks, ContigsMap &contigs_map, std::string &tmp_dir, std::string &outname) {
    std::ofstream linkfile;
    linkfile.open(tmp_dir + '/' + outname + "_countlinks");
<<<<<<< HEAD
    float count_1 = 1.0f, count_2 = 1.0f, normalize_factor = 1.0f;
    for (auto const & k: contiglinks) {
        float frac_links = k.second;
        if (contigs_map.find(k.first.first) != contigs_map.end()) {
            count_1 = contigs_map[k.first.first];
        } else {
            std::cerr << k.first.first << " is not present in contigs_map constructed. Revisit the code \n";
        }
        if (contigs_map.find(k.first.second) != contigs_map.end()) {
            count_2 = contigs_map[k.first.second];
        } else {
            std::cerr << k.first.second << " is not present in contigs_map constructed. Revisit the code \n";
        }
        normalize_factor = (count_1 < count_2) ? count_1 : count_2;
        if (frac_links > normalize_factor) {
            frac_links = normalize_factor;
        }
        frac_links /= normalize_factor;

        // output file format: (i) pair_contig 1, (ii) pair_contig 2, 
        // (iii) fraction of value obtained from the split-mapped read, (iv) iii normalized by the total count of shortest contig in the pair
        linkfile << k.first.first << " " << k.first.second << " " << k.second << " " << frac_links << "\n";
=======
    if (!linkfile.is_open()) {
        std::cerr << "Error opening file: " << tmp_dir + '/' + outname + "_countlinks" << "\n";
        return;
>>>>>>> 2debd9a (update)
    }

    for (auto const & k: contiglinks) {
        const std::string &contig1 = k.first.first;
        const std::string &contig2 = k.first.second;
        float frac_links = k.second;

        auto it1 = contigs_map.find(contig1);
        auto it2 = contigs_map.find(contig2);

        if (it1 == contigs_map.end() || it2 == contigs_map.end()) {
            std::cerr << "One of the contigs (" << contig1 << ", " << contig2 << ") is not present in contigs_map. Revisit the code.\n";
            continue;
        }

        float count_1 = it1->second;
        float count_2 = it2->second;
        float normalize_factor = std::min(count_1, count_2);
        if (frac_links > normalize_factor) {
            frac_links = normalize_factor;
            std::cout << "fraction link count is greater than the minimum total counts of a pair. Revisit the code \n";
        }
        frac_links /= normalize_factor;

        // output file format: (i) pair_contig 1, (ii) pair_contig 2, 
        // (iii) fraction of value obtained from the split-mapped read, (iv) iii normalized by the total count of shortest contig in the pair
        linkfile << contig1 << " " << contig2 << " " << k.second << " " << frac_links << "\n";
    }

    linkfile.close();
}

void process_header_line(
    const std::string& line, 
    ContigsMap& contigs_map,
    ContigsMap& unique_map,
    ContigsMap& cross_map,
    unsigned int minlength) {
    if ((line.rfind("@HD",0) == 0)) { // check if input alignment is unsorted
        if ((line.find("SO:unsorted") == std::string::npos) && line.find("SS") != std::string::npos) {
            std::cerr << "Input alignment file is sorted. Please provide unsorted file or ordered by read ids\n";
            exit(1);
        }
    } else { construct_map(contigs_map, line, minlength); 
        unique_map = contigs_map;
        cross_map = contigs_map;
    };    
}

void process_alignment_line(
    const std::string& line,
    ContigsMap& contigs_map,
    ContigsMap& unique_map,
    ContigsMap& cross_map,
    ContigLinks& contiglinks,
    AlnParser& parsedaln,
    float sequenceidentity,
    float read_coverage,
    std::ofstream& eachcount) {
    if (contigs_map.size() == 0) {
        std::cerr << "Input sam/bam file doesn't has header. Please provide input file with header \n";
        exit(1);
    }
    std::string currentread_id, contig_id, cigar_str, qual_str, md_str, field;
    unsigned int bitflag;

    std::istringstream iss(line);
    iss >> currentread_id >> bitflag >> contig_id;

    // No reference or not of minimum length, don't process the alignment;
    if ((contig_id == "*") || (contigs_map.find(contig_id) == contigs_map.end())) return;

    iss >> field >> field >> cigar_str;
    if (cigar_str.find_first_of("IDNPH*") != std::string::npos) {
        // indel, hard clipping, don't process the alignment;
        if (bitflag & PAIRED_FLAG && !parsedaln.empty()) {
            parsedaln.erase(std::remove_if(parsedaln.begin(), parsedaln.end(),
            [contig_id](const al_data& element) {
                return element.contig == contig_id;
            }),
            parsedaln.end());
        }
        return; 
    }

    iss >> field >> field >> field >> field >> qual_str >> field >> field;
    if (field.find("XS") != std::string::npos) {
        // suboptimal alignment score is given;
        iss >> field >> field >> field >> field >> field >> md_str;
    } else {
        iss >> field >> field >> field >> field >> md_str;
    }

    md_str = md_str.substr(5, md_str.length());

    auto alnpos = get_alnpos(cigar_str);

    // clip at both ends, don't process the alignment;
    if (std::get<1>(alnpos) == 0) return;

    auto alnstats = get_seqid_alncov(alnpos, qual_str, md_str);

    if ((!std::isnan(alnstats.first)) && (alnstats.second >= read_coverage)) {
        // get read direction
        unsigned int pair_dir = (bitflag & FORWARD_FLAG) ? 1: (bitflag & REVERSE_FLAG) ? 2: 0;
        if (pair_dir == 0){
            std::cerr << "input reads are not paired\n";
            exit(1);
        }
        
        // get if alignment is paired
        bool proper_pair = bitflag & PAIRED_FLAG;

        if (parsedaln.empty() || parsedaln.rbegin()->read != currentread_id) {
            // new read alignment
            if (!parsedaln.empty()) {
                counting(parsedaln, contigs_map, unique_map, cross_map, contiglinks, sequenceidentity, eachcount);
                parsedaln.clear();
            }
            // store new read alignment
            parsedaln.push_back({currentread_id, contig_id, pair_dir, proper_pair, alnstats.first});
        } else { // continue store alignment of current read
            storealnstats(parsedaln, alnstats, currentread_id, contig_id, pair_dir, proper_pair);
        }
    }
}

void write_counts(
    const std::string& tmp_dir,
    const std::string& outname,
    const ContigsMap& map,
    const std::string& suffix) {
    std::ofstream outfile(tmp_dir + '/' + outname + suffix);
    for (const auto& k : map) {
        outfile << k.first << " " << outname << " " << k.second << "\n";
    }
    outfile.close();
}

// TO DO LIST
// Flags in sam files are considered only as per bowtie2 output. Need to generalize or check compatibility for bwa-mem and minimap2 output
// 

int main(int argc, char *argv[]) {
    auto start = std::chrono::high_resolution_clock::now();

    // never use scanf or any library which would use scanf, but it makes input much faster
    std::ios::sync_with_stdio(false);

    if (argc == 1 || strcmp(argv[1], "--help") == 0 || strcmp(argv[1], "-h") == 0) {
        std::cerr << "aligner command | ./aligner2counts working_dir outputname minlength[optional, default=1000]" << "\n";
        return 1;
    }
    std::string tmp_dir = argv[1];
    std::string outname = (argc > 2) ? argv[2] : "";
    unsigned int minlength = (argc > 3) ? std::stoul(argv[3]) : 1000;

    ContigsMap contigs_map;
    ContigsMap unique_map;
    ContigsMap cross_map;
    ContigLinks contiglinks;
    AlnParser parsedaln;
    float sequenceidentity = 97.0f;
    float read_coverage = 99.0f;

    std::ofstream samfile(tmp_dir + '/' + outname + ".sam");
    std::ofstream eachcount(tmp_dir + '/' + outname + "_eachcount");

    std::string line;
    while(std::getline(std::cin, line)) {
        samfile << line << "\n"; // write all line before applying any conditions and proceed with continue
        if (line.empty()) continue;

<<<<<<< HEAD
                std::string field;
                iss >> bitflag;
                iss >> contig_id;

                if ((contig_id == "*") || (contigs_map.find(contig_id) == contigs_map.end())) {
                    // No reference or not of minimum length, don't process the alignment;
                    continue;
                }

                iss >> field >> field;
                iss >> cigar_str;
                if (cigar_str.find_first_of("IDNP*") != std::string::npos) {
                    // indel, don't process the alignment;
                    if (bitflag & PAIRED_FLAG) {
                        if (!parsedaln.empty()) {
                            parsedaln.erase(std::remove_if(parsedaln.begin(), parsedaln.end(),
                            [contig_id](const al_data& element) {
                                return element.contig == contig_id;
                            }),
                            parsedaln.end());
                        }  
                    }
                    continue;
                }

                iss >> field >> field >> field >> field;
                iss >> qual_str;
                iss >> field >> field;
                if (field.find("XS") != std::string::npos) {
                    // suboptimal alignment score is given;
                    iss >> field >> field >> field >> field >> field >> md_str;
                } else {
                    iss >> field >> field >> field >> field >> md_str;
                }

                md_str = md_str.substr(5, md_str.length());

                auto alnpos = get_alnpos(cigar_str);

                if (std::get<1>(alnpos) == 0) {
                    // clip at both ends, don't process the alignment;
                    continue;
                }

                auto alnstats = get_seqid_alncov(alnpos, qual_str, md_str);

                if ((!std::isnan(alnstats.first)) && (alnstats.second >= 70.0f)) {
                    // get read direction
                    unsigned int pair_flag = 0;
                    if (bitflag & FORWARD_FLAG) {
                        pair_flag = 1;
                    } else if (bitflag & REVERSE_FLAG) {
                        pair_flag = 2;
                    } else {
                        std::cerr << "input reads are not paired\n";
                        return 1;
                    }
                    
                    // get if alignment is paired
                    bool proper_pair = false;
                    if (bitflag & PAIRED_FLAG) {
                        proper_pair = true;
                    } else {;}

                    if (parsedaln.empty()) {
                        parsedaln.push_back({currentread_id, contig_id, pair_flag, proper_pair, false, alnstats.first});
                    } else { // store alignments of current read
                        if (parsedaln.rbegin()->read == currentread_id) {
                            storealnstats(parsedaln,alnstats, currentread_id, contig_id, pair_flag, proper_pair);
                        } else { // new read alignment
                            counting(parsedaln, contigs_map, contiglinks, sequenceidentity, eachcount);
                            parsedaln.clear(); // clean alignment data for current read
                            parsedaln.push_back({currentread_id, contig_id, pair_flag, proper_pair, false, alnstats.first});
                        }
                    }
                }
            }
=======
        if (line[0] == '@') {
            process_header_line(line, contigs_map, unique_map, cross_map, minlength);
        } else {
            process_alignment_line(line, contigs_map, unique_map, cross_map, contiglinks, parsedaln, sequenceidentity, read_coverage, eachcount);
>>>>>>> 2debd9a (update)
        }
    }
    // while end
    counting(parsedaln, contigs_map, unique_map, cross_map, contiglinks, sequenceidentity, eachcount); // process last read alignment(s)
    parsedaln.clear();

    fractionate_countlinks(contiglinks, contigs_map, tmp_dir, outname);

    write_counts(tmp_dir, outname, contigs_map, "_count");
    write_counts(tmp_dir, outname, unique_map, "_uniqcount");
    write_counts(tmp_dir, outname, cross_map, "_crosscount");

    samfile.close();
    eachcount.close();

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    std::cout << "Alignment file processed in " << duration.count() << " seconds\n";

    return EXIT_SUCCESS;
}
