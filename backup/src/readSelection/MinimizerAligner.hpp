
#ifndef MDBG_METAG_MINIMIZERALIGNER
#define MDBG_METAG_MINIMIZERALIGNER


#include "../Commons.hpp"
//#include <cstdint>
//#include <memory>
//#include <string>
//#include <vector>



class MinimizerAligner {

    public:
    

    typedef vector<pair<int64_t, int64_t>> Alignment;
    int64_t kNegativeInfinity;

    std::int8_t m_;
    std::int8_t n_;
    std::int8_t g_;
    std::int8_t e_;
    std::int8_t q_;
    std::int8_t c_;
    std::vector<int8_t> sequence_profile;
    std::vector<int32_t> M;
    //vector<MinimizerType> rank_to_node;

    //MinimizerAligner(){

    //}
    MinimizerAligner(int8_t m, int8_t n, int8_t g, int8_t e){
        m_ = m;
        n_ = n;
        g_ = g;
        e_ = e;
        q_ = g;
        c_ = e;

        kNegativeInfinity = std::numeric_limits<std::int64_t>::min() + 1024;
    }

    /*
    void setup(int8_t m, int8_t n, int8_t g, int8_t e){
        m_ = m;
        n_ = n;
        g_ = g;
        e_ = e;
        q_ = g;
        c_ = e;

        kNegativeInfinity = std::numeric_limits<std::int64_t>::min() + 1024;
    }
    */

    ~MinimizerAligner(){
    }


    Alignment performAlignment(const vector<MinimizerType>& referenceMinimizers, const vector<MinimizerType>& queryMinimizers){

        //sequence_profile.clear();
        //M.clear();

        size_t referenceSize = referenceMinimizers.size();
        size_t querySize = queryMinimizers.size();

        //size_t sequence_len = queryMinimizers.size();
        if(referenceMinimizers.size() == 0 || queryMinimizers.size() == 0) return {};

        realloc(querySize + 1, referenceSize + 1, referenceMinimizers.size());
        initialize(referenceMinimizers, queryMinimizers);
        return align(referenceMinimizers, queryMinimizers);
    }

    private: 
    
    void realloc(uint64_t matrix_width, uint64_t matrix_height, uint64_t num_codes) {
        //if (node_id_to_rank.size() < matrix_height - 1) {
        //    node_id_to_rank.resize(matrix_height - 1, 0);
        //}
        if (sequence_profile.size() < num_codes * matrix_width) {
            sequence_profile.resize(num_codes * matrix_width, 0);
        }
        if (M.size() < matrix_height * matrix_width) {
            M.resize(matrix_width * matrix_height, 0);
            //H = M.data();
            //F = nullptr;
            //E = nullptr;
        }
    }

    //void clear(){
      //  cout << "Can't manage to clear without seg fault" << endl;
        //sequence_profile.clear();
        //M.clear();
        //H.clear();

    //}
    

    void initialize(const vector<MinimizerType>& referenceMinimizers, const vector<MinimizerType>& queryMinimizers) {
            
        std::uint64_t matrix_width = queryMinimizers.size() + 1;
        std::uint64_t matrix_height = referenceMinimizers.size() + 1;

        //int64_t num_codes = referenceMinimizers.size();

        for (std::uint64_t i = 0; i < referenceMinimizers.size(); ++i) {
            std::int64_t c = referenceMinimizers[i];//graph.decoder(i);
            sequence_profile[i * matrix_width] = 0;
            for (std::uint64_t j = 0; j < queryMinimizers.size(); ++j) {
            sequence_profile[i * matrix_width + (j + 1)] =
                (c == queryMinimizers[j] ? m_ : n_);
            }
        }

        //const auto& rank_to_node = graph.rank_to_node();
        //for (std::uint64_t i = 0; i < rank_to_node.size(); ++i) {
        //    node_id_to_rank[rank_to_node[i]->id] = i;
        //}

        // initialize secondary matrices
        M[0] = 0;


        //for (std::uint64_t i = 1; i < matrix_height; ++i) {
        //  H[i * matrix_width] = i * g_;
        //}
        for (std::uint64_t j = 1; j < matrix_width; ++j) {
            M[j] = 0 ;//j * g_;
        }
        
        for (std::uint64_t i = 1; i < matrix_height; ++i) {
            M[i * matrix_width] = 0;
        }

        //for (std::uint64_t j = 1; j < matrix_width; ++j) {
            //H[j] = 0;
        //}
    }

    
    Alignment align(const vector<MinimizerType>& referenceMinimizers, const vector<MinimizerType>& queryMinimizers){
        
        //size_t referenceSize = referenceMinimizers.size();
        size_t querySize = queryMinimizers.size();

        std::uint64_t matrix_width = querySize + 1;
        //const auto& rank_to_node = graph.rank_to_node();

        std::int64_t max_score = kNegativeInfinity;

        std::uint64_t max_i = 0;
        std::uint64_t max_j = 0;

        /*
        auto update_max_score = [&max_score, &max_i, &max_j] (
            std::int64_t* H_row,
            std::uint64_t i,
            std::uint64_t j) -> void {
                if (max_score < H_row[j]) {
                max_score = H_row[j];
                max_i = i;
                max_j = j;
                }
                return;
        };
        */

        auto update_max_score_ov = [&max_score, &max_i, &max_j] (
            auto* H_row,
            auto i,
            auto j) -> void {
                if (max_score <= H_row[j]) {
                    max_score = H_row[j];
                    max_i = i;
                    max_j = j;
                }
                return;
        };
        
        // alignment
        for (size_t ii=0; ii<referenceMinimizers.size(); ii++) {

            const auto& char_profile =
                &(sequence_profile[ii * matrix_width]);

            std::uint64_t i = ii +1;//node_id_to_rank[it->id] + 1;
            std::uint64_t pred_i = ii;//it->inedges.empty() ? 0 : node_id_to_rank[it->inedges[0]->tail->id] + 1;

            auto* H_row = &(M[i * matrix_width]);
            auto* H_pred_row = &(M[pred_i * matrix_width]);

            // update H
            for (std::uint64_t j = 1; j < matrix_width; ++j) {
            H_row[j] = std::max(
                H_pred_row[j - 1] + char_profile[j],
                H_pred_row[j] + g_);
            }
            /*
            // check other predeccessors
            for (std::uint64_t p = 1; p < it->inedges.size(); ++p) {
            pred_i = node_id_to_rank[it->inedges[p]->tail->id] + 1;

            H_pred_row = &(H[pred_i * matrix_width]);

            for (std::uint64_t j = 1; j < matrix_width; ++j) {
                H_row[j] = std::max(
                    H_pred_row[j - 1] + char_profile[j],
                    std::max(
                        H_row[j],
                        H_pred_row[j] + g_));
            }
            }
            */


            for (std::uint64_t j = 1; j < matrix_width; ++j) {
                H_row[j] = std::max(H_row[j - 1] + g_, H_row[j]);

                //if (type_ == AlignmentType::kSW) {
                //    H_row[j] = std::max(H_row[j], (std::int64_t)0);
                //    update_max_score(H_row, i, j);
                //} else if (type_ == AlignmentType::kNW &&
                //    it->outedges.empty() && j == matrix_width - 1) {
                //    update_max_score(H_row, i, j);
                //} else if (type_ == AlignmentType::kOV && it->outedges.empty()) {
                //    update_max_score_ov(H_row, i, j);
                //}
                if(ii == referenceMinimizers.size()-1){
                    update_max_score_ov(H_row, i, j);
                }
            }

            //cout << "On doit prendre le max score dasn derniere ligne OU derniere colonne non ?" << endl;

            //if (type_ == AlignmentType::kOV && it->outedges.empty()) {
            //  for (std::uint64_t i = 1; i < matrix_width; ++i) {
            //    update_max_score_ov(H_row, i, j);
            //  }
            //}

        }


        std::uint64_t maxI = 0;
        std::uint64_t maxJ = 0;
        
        std::int64_t maxScore = kNegativeInfinity;
        //std::cout << maxScore << std::endl;

        for(size_t i=1; i<matrix_width; i++){

        int lastColScore = M[matrix_width*referenceMinimizers.size()+i];

        if(lastColScore >= maxScore){
            maxI = referenceMinimizers.size();
            maxJ = i;
            //std::cout << "Best col score: " << matrix_width*rank_to_node.size()+i << " " << H[matrix_width*rank_to_node.size()+i] << std::endl;
            maxScore = lastColScore;
            //std::cout << maxI << " " << maxJ << " " << maxScore << std::endl;
        }
        //std::cout << "Last col score: " << matrix_width*rank_to_node.size()+i << " " << H[matrix_width*rank_to_node.size()+i] << std::endl;
        }

        for(size_t i=1; i<referenceMinimizers.size()+1; i++){

        int lastRowScore = M[(matrix_width-1)+(matrix_width*i)];

        if(lastRowScore >= maxScore){

            maxI = i;
            maxJ = matrix_width-1;

            //std::cout << "Best row score: " << (matrix_width-1)+(matrix_width*i) << " " << H[(matrix_width-1)+(matrix_width*i)] << std::endl;
            maxScore = lastRowScore;
            //std::cout << maxI << " " << maxJ << " " << maxScore << std::endl;
        }

        //std::cout << "Last row score: " << (matrix_width-1)+(matrix_width*i) << " " << H[(matrix_width-1)+(matrix_width*i)] << std::endl;
        }
        
        max_i = maxI;
        max_j = maxJ;
        max_score = maxScore;


        if (max_i == 0 && max_j == 0) {
            return {};
        }
        //if (score) {
        //    *score = max_score;
        //}



        // backtrack
        Alignment alignment;
        std::uint64_t i = max_i;
        std::uint64_t j = max_j;
        //std::cout << "Max: " << i << " " << j << " " << max_score << std::endl;
        
        //auto sw_condition = [this, &i, &j, &matrix_width] () -> bool {
        //    return (H[i * matrix_width + j] == 0) ? false : true;
        //};
        //auto nw_condition = [&i, &j] () -> bool {
        //    return (i == 0 && j == 0) ? false : true;
        //};
        auto ov_condition = [&i, &j] () -> bool {
            return (i == 0 || j == 0) ? false : true;
        };

        std::uint64_t prev_i = 0;
        std::uint64_t prev_j = 0;
        while (ov_condition()) {

            auto H_ij = M[i * matrix_width + j];
            //std::cout << "Miam: " << i << " " << j << " " << H_ij << std::endl;

            bool predecessor_found = false;

            if (i != 0 && j != 0) {
                int64_t ii = i - 1;
                //const auto& it = rank_to_node[i - 1];
                std::int64_t match_cost =
                    sequence_profile[ii * matrix_width + j];

                uint64_t pred_i = ii;
                //std::uint64_t pred_i = it->inedges.empty() ?
                //    0 : node_id_to_rank[it->inedges[0]->tail->id] + 1;

                //std::uint64_t pred_i = ii==0 ?
                //    0 : node_id_to_rank[it->inedges[0]->tail->id] + 1;
                if (H_ij == M[pred_i * matrix_width + (j - 1)] + match_cost) {
                    prev_i = pred_i;
                    prev_j = j - 1;
                    predecessor_found = true;
                    //std::cout << "\tRA1 " << prev_i << " " << prev_j << std::endl;
                } else {
                    /*
                    for (std::uint64_t p = 1; p < it->inedges.size(); ++p) {
                    std::uint64_t pred_i =
                        node_id_to_rank[it->inedges[p]->tail->id] + 1;

                    if (H_ij == H[pred_i * matrix_width + (j - 1)] + match_cost) {
                        prev_i = pred_i;
                        prev_j = j - 1;
                        predecessor_found = true;
                        break;
                    }
                    }
                    */
                    //std::cout << "\tRA2 " << prev_i << " " << prev_j << std::endl;
                }
            }

            //getchar();
            
            if (!predecessor_found && i != 0) {
                int64_t ii = i - 1;
                //const auto& it = rank_to_node[i - 1];

                std::uint64_t pred_i = ii; //ii == 0 ? 0 : i +1;
                //node_id_to_rank[it->inedges[0]->tail->id] + 1;

                if (H_ij == M[pred_i * matrix_width + j] + g_) {
                    prev_i = pred_i;
                    prev_j = j;
                    predecessor_found = true;
                } else {
                    /*
                    for (std::uint64_t p = 1; p < it->inedges.size(); ++p) {
                        std::uint64_t pred_i =
                            node_id_to_rank[it->inedges[p]->tail->id] + 1;

                        if (H_ij == H[pred_i * matrix_width + j] + g_) {
                            prev_i = pred_i;
                            prev_j = j;
                            predecessor_found = true;
                            break;
                        }
                    }
                    */
                }
            }

            if (!predecessor_found && H_ij == M[i * matrix_width + j - 1] + g_) {  // NOLINT
                prev_i = i;
                prev_j = j - 1;
                predecessor_found = true;
            }
            

            alignment.emplace_back(
                i == prev_i ? -1 : i - 1,
                j == prev_j ? -1 : j - 1);

            i = prev_i;
            j = prev_j;
        }

        std::reverse(alignment.begin(), alignment.end());
        return alignment;
    }
    


};


#endif
