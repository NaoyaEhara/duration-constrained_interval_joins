#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include <utility>
#include <chrono>
#include <random>
#include <array>
#include <unistd.h>
#include <sys/resource.h>
#include <numeric>
using namespace std;

struct Interval
{
    int id, start, end, duration;
};

struct Group
{
    vector<Interval> Queries;
    int start;
    int end;
    int start_margin;
    int end_margin;
    int max_start;
};

using IntervalPair = std::pair<std::vector<Interval>, std::vector<Interval>>;

bool check_overlap(const Interval& r, const Interval& s, float param, unsigned long long* result)
{
    int threshold = param;
    int overlap_start = std::max(r.start, s.start);
    int overlap_end = std::min(r.end, s.end);
    int overlap_length = overlap_end - overlap_start;
    if (overlap_length >= threshold) {
        *result += (r.start ^ s.start);
        return true;
    }
    return false;
}

double process_mem_usage() {
    long rss = 0L;
    FILE* fp = nullptr;
    if ((fp = fopen("/proc/self/statm", "r")) == nullptr) return 0.0;
    if (fscanf(fp, "%*s%ld", &rss) != 1) { fclose(fp); return 0.0; }
    fclose(fp);
    long page_size_kb = sysconf(_SC_PAGESIZE) / 1024; // page size in KB
    double mem_mb = rss * page_size_kb / 1024.0;      // convert to MB
    return mem_mb;
}

vector<Interval> NextSubseq(const vector<Interval> &r, int h, int b){
    if ((h + b) >= r.size()){
        auto first = r.cbegin() + h;
        auto last = r.cend();  //h_ - 1 + 1 = h_
 
        vector<Interval> vec_last(first, last);
        return vec_last;
    }
    int h_ = h + b;
    if(r[h].start == r[h_].start){
        while((h_ < r.size()) && (r[h].start == r[h_].start)){
            h_ = h_ + 1;
        }
    }
    else{
        while(r[h_-1].start == r[h_].start){
            h_--;
        }
    }

    auto first = r.cbegin() + h;
    auto last = r.cbegin() + h_; 
 
    vector<Interval> vec(first, last);
    return vec;
}

inline unsigned long long overlapping_query(const vector<vector<IntervalPair>> &grid, float overlap_duration, const vector<int> &col_maxstart,
                                            const vector<vector<int>> &cell_minend, const vector<vector<int>> &cell_maxend, const vector<vector<int>> &cell_minstart, int start_index_common, const Group Group){
    unsigned long long result = 0;
    int index_i = 0;
    int grouped_queries_index;
    int max_grouped_queries_index = Group.Queries.size();
    Interval q;
    int threshold = overlap_duration;
    int start_index;
    double elapsed;
    std::chrono::system_clock::time_point start, end;
    
    while(index_i < start_index_common){
        int last_checked_cell_index = cell_maxend[index_i].size();
        int maxcell_index = cell_maxend[index_i].size();
        grouped_queries_index = max_grouped_queries_index - 1;
        while(grouped_queries_index >= 0){
            q = Group.Queries[grouped_queries_index];
            int threshold_end = q.start + threshold;
            auto cell_end_it = std::lower_bound(cell_maxend[index_i].begin(), cell_maxend[index_i].end(), threshold_end);
            int index_j = std::distance(cell_maxend[index_i].begin(), cell_end_it);
            while((index_j < maxcell_index)){
                if(last_checked_cell_index == index_j){
                    break;
                } 
                else if(cell_minend[index_i][index_j] >= threshold_end){
                    int tmp_cell_index = index_j;
                    while(index_j < last_checked_cell_index){
                        for(Interval r : grid[index_i][index_j].second){
                            int tmp_index = grouped_queries_index;
                            while(tmp_index >= 0){
                                q = Group.Queries[tmp_index];
                                result += (r.start ^ q.start);
                                tmp_index--;
                            }
                        }
                        index_j++;
                    }
                    last_checked_cell_index = tmp_cell_index;
                    break;
                }else{
                    auto cell_end_it = grid[index_i][index_j].second.end();
                    auto cell_it = std::lower_bound(grid[index_i][index_j].second.begin(), cell_end_it, threshold_end,
                        [](const Interval &a, int b) {     
                        return a.end < b;
                    });
                    while((cell_it < cell_end_it)){
                        result += ((*cell_it).start ^ q.start);
                        cell_it++;
                    }
                }
                index_j++;
            }
            grouped_queries_index--;
        }
        index_i++;
    }

    // normal comparison
    grouped_queries_index = 0;
    start_index = start_index_common;
    q = Group.Queries[grouped_queries_index];
    int cell_threshold_start;
    int threshold_end = q.start + threshold;
    int threshold_start = q.end - threshold;
    auto col_thst_it = std::lower_bound(col_maxstart.begin(), col_maxstart.end(), threshold_start);
    int threshold_start_index = std::distance(col_maxstart.begin(), col_thst_it);
    
    while((start_index < col_maxstart.size()) && (start_index <= threshold_start_index)){
        if(start_index < threshold_start_index){
            auto cell_thend_it = std::lower_bound(cell_maxend[start_index].begin(), cell_maxend[start_index].end(), threshold_end);
            int index_k = std::distance(cell_maxend[start_index].begin(), cell_thend_it);
            int max_column_size = cell_maxend[start_index].size();
            while((index_k < max_column_size)){
                if(cell_minend[start_index][index_k] >= Group.end_margin){
                    if (max_grouped_queries_index > 1){
                        for(Interval r : grid[start_index][index_k].first){
                            for(Interval g : Group.Queries){
                                result += (r.start ^ g.start);
                            }
                        }
                    }else{
                        for(Interval r : grid[start_index][index_k].first){
                            result += (r.start ^ q.start);
                        }
                    }
                    index_k++;
                    continue;
                }else if(cell_minstart[start_index][index_k] > (Group.end_margin - threshold)){
                    index_k++;
                    continue;
                }
                else{
                    auto cell_it = grid[start_index][index_k].first.begin();
                    auto cell_end_it = grid[start_index][index_k].first.end();
                    int cell_min_end = cell_minend[start_index][index_k];
                    cell_threshold_start = cell_min_end - threshold;
                    cell_threshold_start = std::max(q.start, cell_threshold_start);
                    
                    if(cell_min_end < threshold_end){
                        while((cell_it->start <= cell_threshold_start) && (cell_it < cell_end_it)){
                            
                            if(cell_it->duration < overlap_duration){
                                cell_it++;
                                continue;
                            }
                            if((Group.max_start + threshold) <= cell_it->end){
                                for (Interval r : Group.Queries){
                                    result += (cell_it->start ^ r.start);
                                }
                                cell_it++;
                                continue;
                            }
                            
                            int s_upper_threshold = cell_it->end - threshold;
                            grouped_queries_index = 0;
                            while((Group.Queries[grouped_queries_index].start <= s_upper_threshold) && (grouped_queries_index < max_grouped_queries_index)){
                                result += (cell_it->start ^ Group.Queries[grouped_queries_index].start);
                                grouped_queries_index++;
                            }
                            cell_it++;
                        }
                    }else{
                        while((cell_it->start <= cell_threshold_start) && (cell_it < cell_end_it)){
                            if((Group.max_start + threshold) <= cell_it->end){
                                for (Interval r : Group.Queries){
                                    result += (cell_it->start ^ r.start);
                                }
                                cell_it++;
                                continue;
                            }
                            
                            result += ((*cell_it).start ^ q.start);
                            if (max_grouped_queries_index > 1){
                                int s_upper_threshold = cell_it->end - threshold;
                                grouped_queries_index = 1;
                                while((Group.Queries[grouped_queries_index].start <= s_upper_threshold) && (grouped_queries_index < max_grouped_queries_index)){
                                    result += (cell_it->start ^ Group.Queries[grouped_queries_index].start);
                                    grouped_queries_index++;
                                }
                            }
                            cell_it++;
                        }
                    }

                    cell_threshold_start = std::min(cell_maxend[start_index][index_k], Group.end_margin) - threshold;
                    while((cell_it < cell_end_it) && (cell_it->start <=  cell_threshold_start)){
                        if(cell_it->duration < overlap_duration){
                            cell_it++;
                            continue;
                        }
                        
                        if(Group.end_margin <= cell_it->end){
                            for (Interval r : Group.Queries){
                                result += (cell_it->start ^ r.start);
                            }
                            cell_it++;
                            continue;
                        }
                        
                        int s_lower_threshold = cell_it->start + threshold;
                        int s_upper_threshold = cell_it->end - threshold;
                        grouped_queries_index = 0;
                        Interval tmp = Group.Queries[grouped_queries_index];
                        while((tmp.start <= cell_it->start) && (grouped_queries_index < max_grouped_queries_index)){
                            if(tmp.end >= s_lower_threshold){
                                result += (cell_it->start ^ tmp.start);
                            }
                            tmp = Group.Queries[++grouped_queries_index];
                        }
                        while((tmp.start <= s_upper_threshold) && (grouped_queries_index < max_grouped_queries_index)){
                            result += (cell_it->start ^ tmp.start);
                            tmp = Group.Queries[++grouped_queries_index];
                        }
                        cell_it++;
                    }
                }
                index_k++;
            }
        }else{
            auto cell_thend_it = std::lower_bound(cell_maxend[start_index].begin(), cell_maxend[start_index].end(), threshold_end);
            int index_l = std::distance(cell_maxend[start_index].begin(), cell_thend_it);
            int max_column_size = cell_minstart[start_index].size(); 
            while((index_l < max_column_size)){
                if (cell_minstart[start_index][index_l] > (Group.end_margin - threshold)) {
                    index_l++;
                    continue;
                }

                auto cell_it = grid[start_index][index_l].first.begin();
                auto cell_end_it = grid[start_index][index_l].first.end();
                int cell_min_end = cell_minend[start_index][index_l];
                cell_threshold_start = cell_min_end - threshold;
                cell_threshold_start = std::max(q.start, cell_threshold_start);
                cell_threshold_start = std::min(cell_threshold_start, threshold_start);

                while((cell_it->start <= cell_threshold_start) && (cell_it < cell_end_it)){                 
                    if((Group.max_start + threshold) <= cell_it->end){
                        for (Interval r : Group.Queries){
                            result += (cell_it->start ^ r.start);
                        }
                        cell_it++;
                        continue;
                    }
                    int s_upper_threshold = cell_it->end - threshold;
                    grouped_queries_index = 0;
                    while((Group.Queries[grouped_queries_index].start <= s_upper_threshold) && (grouped_queries_index < max_grouped_queries_index)){
                        result += (cell_it->start ^ Group.Queries[grouped_queries_index].start);
                        grouped_queries_index++;
                    }
                    cell_it++;
                }

                cell_threshold_start = std::min(cell_maxend[start_index][index_l], Group.end_margin) - threshold;

                while((cell_it < cell_end_it) && (cell_it->start <= cell_threshold_start)){
                    if(cell_it->duration < overlap_duration){
                        cell_it++;
                        continue;
                    }
                    grouped_queries_index = 0;
                    int s_lower_threshold = cell_it->start + threshold;
                    int s_upper_threshold = cell_it->end - threshold;
                    Interval tmp = Group.Queries[grouped_queries_index];
                    while((tmp.start < cell_it->start) && (grouped_queries_index < max_grouped_queries_index)){
                        if(tmp.end >= s_lower_threshold){
                            result += (cell_it->start ^ tmp.start);
                        }
                        tmp = Group.Queries[++grouped_queries_index];
                    }
                    while((tmp.start <= s_upper_threshold) && (grouped_queries_index < max_grouped_queries_index)){
                        result += (cell_it->start ^ tmp.start);
                        tmp = Group.Queries[++grouped_queries_index];
                    }
                    cell_it++;
                }
                index_l++;
            }
        }
        start_index++;
    }
    return result;
}

int main(int argc, char* argv[])
{
    if (argc < 7) {
        cerr << "Usage: " << argv[0] << " overlap_duration page_size relation_R relation_S max_batch_size margin_param " << endl;
        return 1;
    }
    float overlap_duration = std::stof(argv[1]);
    int page_size = std::stoi(argv[2]);
    string file1 = argv[3];
    string file2 = argv[4];
    int max_batch_size = std::stoi(argv[5]);
    float margin_param = std::stof(argv[6]);

    //load queries
    vector<Interval> queries;
    ifstream ifs(file1);
    int qstart, qend;
    int qduration;
    int numQueries = 0;
    while (ifs >> qstart >> qend)
    {
        qduration = qend - qstart;
        Interval q = {numQueries, qstart, qend, qduration};
        queries.emplace_back(q);
        numQueries++;
    }
    ifs.close();

    int margin = overlap_duration * margin_param;
    cout << "ours : calculated margin is " << margin << endl;
    if (margin <= 0) margin = 1; // at least 1

    //load relations
    vector<Interval> relations;
    int rstart, rend;
    int rduration;
    ifstream inp(file2);
    int numRecords = 0;

    if (!inp)
    {
        cerr << "Error - cannot open data file " << endl;
        exit(1);
    }

    while (inp >> rstart >> rend)
    {
        if (rstart > rend)
        {
            cerr << endl
                 << "Error - start is after end for interval [" << rstart << ".." << rend << "]" << endl
                 << endl;
            exit(1);
        }
        rduration = rend - rstart;
        Interval r = {numRecords, rstart, rend, rduration};
        relations.emplace_back(r);
        numRecords++;
    }
    inp.close();

    //query & indexing
    vector<int> data_size;
    vector<double> query_time;
    vector<double> total_query_time;
    vector<double> total_result_size;
    unsigned long long result = 0;
    std::chrono::system_clock::time_point start, end;
    double elapsed;
    double g_elapsed;

    double mem = process_mem_usage();
    cout << "ours : making index, size is " << relations.size() << endl;

    // aray parameters
    int max_col_height = relations.size() / (page_size * page_size) + 10; // 10 is margin
    using IntervalPair = std::pair<std::vector<Interval>, std::vector<Interval>>;
    std::vector<std::vector<IntervalPair>> grid(
        max_col_height,
        std::vector<IntervalPair>(page_size * 10, IntervalPair())
    );
    vector<int> col_maxstart;
    vector<vector<int>> cell_minend(max_col_height,vector<int>());
    vector<vector<int>> cell_maxend(max_col_height,vector<int>());
    vector<vector<int>> cell_minstart(max_col_height,vector<int>());
    vector<Interval> column;
    vector<Interval> cell;
    
    start = std::chrono::system_clock::now();

    // sort by start, ascending
    std::sort(relations.begin(), relations.end(), [](const Interval &a, const Interval &b) {
        return a.start < b.start;
    });

    // structuring
    int h = 0; // position in r
    int i = 0; // column index
    while(h < relations.size()){
        column.clear();
        column = NextSubseq(relations,h,page_size*page_size);
        
        auto max_ = std::max_element(column.begin(), column.end(), [](const Interval &a, const Interval &b) {
            return a.start < b.start;
        });
        col_maxstart.emplace_back(max_->start);

        std::sort(column.begin(), column.end(), [](const Interval &a, const Interval &b) {
            return a.end < b.end;
        });

        int k = 0;
        int j = 0;

        while(k < column.size()){
            cell.clear();
            cell = NextSubseq(column,k,page_size);

            std::sort(cell.begin(), cell.end(), [](const Interval &a, const Interval &b) {
                return a.start < b.start;
            });

            grid[i][j].first = cell;

            std::sort(cell.begin(), cell.end(), [](const Interval &a, const Interval &b) {
                return a.end < b.end;
            });

            grid[i][j].second = cell;

            auto min_end = std::min_element(cell.begin(), cell.end(), [](const Interval &a, const Interval &b) {
                return a.end < b.end;
            });
            cell_minend[i].emplace_back(min_end->end);

            auto max_end = std::max_element(cell.begin(), cell.end(), [](const Interval &a, const Interval &b) {
                return a.end < b.end;
            });
            cell_maxend[i].emplace_back(max_end->end);

            auto min_start = std::min_element(cell.begin(), cell.end(), [](const Interval &a, const Interval &b) {
                return a.start < b.start;
            });
            cell_minstart[i].emplace_back(min_start->start);

            j++;
            k = k + cell.size();
        }
        h = h + column.size();
        i++;
    }

    mem = process_mem_usage() - mem; // after build
    //std::cerr << "memory_usage_MB=" << mem << "\n";
    
    std::sort(queries.begin(), queries.end(), [](const Interval &a, const Interval &b) {
        return a.start < b.start;
    });

    end = std::chrono::system_clock::now();  
    elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end-start).count();

    cout << "ours : making index time is " << elapsed  / 1000000 << " secs" << endl;
    cout << "ours : make ok" << endl;
    
    // query processing
    vector<Group> Groups;
    int query_count = 0;
    int tmp_query_count = 0;
    result = 0;
    elapsed = 0;
    g_elapsed = 0;

    while(query_count < queries.size()){
        Groups.clear();
        Interval q;
        Interval current_query; 

        start = std::chrono::system_clock::now();

        q = queries[query_count];
        if (q.duration < overlap_duration){
            query_count++;
            tmp_query_count++;
            continue;
        }
        
        auto st_it = std::lower_bound(col_maxstart.begin(), col_maxstart.end(), q.start);
        int st_index = std::distance(col_maxstart.begin(), st_it);

        while((query_count < queries.size()) && (std::lower_bound(col_maxstart.begin(), col_maxstart.end(), queries[query_count].start) == st_it)){
            current_query = queries[query_count];
            if (current_query.duration < overlap_duration){
                tmp_query_count++;
                query_count++;
                continue;
            }
            tmp_query_count++;
            if (tmp_query_count >= max_batch_size){
                break;
            }
            bool found_group = false;
            for (auto& group : Groups) {
                if (group.start_margin >= current_query.start && group.end <= current_query.end 
                    && group.end_margin >= current_query.end) {
                    group.Queries.emplace_back(current_query);
                    if (current_query.start > group.max_start) {
                        group.max_start = current_query.start;
                    }
                    found_group = true;
                    break;
                }
            }
            if (!found_group) {
                Group new_group;


                auto col_thst_it = std::lower_bound(col_maxstart.begin(), col_maxstart.end(), current_query.end - overlap_duration);
                int threshold_start_index = std::distance(col_maxstart.begin(), col_thst_it);

                auto col_thmst_it = std::lower_bound(col_maxstart.begin(), col_maxstart.end(), (current_query.end + margin - overlap_duration));
                int threshold_m_start_index = std::distance(col_maxstart.begin(), col_thmst_it);

                int tmp = current_query.end + margin;
                if (threshold_start_index != threshold_m_start_index){
                    tmp = col_maxstart[threshold_start_index] + overlap_duration;
                }
                
                tmp = std::min(tmp, current_query.end + margin);
                
                new_group.start = current_query.start;
                new_group.end = current_query.end;
                new_group.start_margin = current_query.start + margin;
                new_group.end_margin = tmp;
                new_group.max_start = current_query.start;
                new_group.Queries.emplace_back(current_query);

                Groups.emplace_back(new_group);
            }
            query_count++;
        }

        end = std::chrono::system_clock::now();
        g_elapsed += std::chrono::duration_cast<std::chrono::microseconds>(end-start).count();

        start = std::chrono::system_clock::now();
        for (int groups_count = 0; groups_count < Groups.size(); groups_count++){
            result += overlapping_query(grid, overlap_duration, col_maxstart, cell_minend, cell_maxend, cell_minstart, st_index, Groups[groups_count]);
        }
        end = std::chrono::system_clock::now();  
        elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end-start).count();
        query_time.emplace_back(elapsed);

        tmp_query_count = 0;
    }
    

    cout << query_count << " queries are processed" << endl;
    cout << "ours : query processing is end" << endl;
    cout << "ours : grouping time is " << g_elapsed / 1000000 << " secs" << endl;

    total_query_time.emplace_back(accumulate(query_time.begin(), query_time.end(), 0.0) / 1000000);
    total_result_size.emplace_back(result);
    data_size.emplace_back(relations.size());

    cout << "total query time: " << total_query_time[0] << " secs" << endl;

    query_time.clear();
    grid.clear();
    column.clear();
    cell.clear();
    cell_minend.clear();
    cell_maxend.clear();
    cell_minstart.clear();
    col_maxstart.clear();

    // make output file
    std::string base = file2;
    size_t slash_pos = base.find_last_of('/');
    if (slash_pos != std::string::npos) base = base.substr(slash_pos + 1);
    size_t dot_pos = base.find_last_of('.');
    if (dot_pos != std::string::npos) base = base.substr(0, dot_pos);

    int overlap_int = static_cast<int>(overlap_duration);
    std::ostringstream oss;
    oss.width(2);
    oss.fill('0');
    oss << overlap_int;
    std::string out_filename = "batch_" + base + "_" + oss.str() + ".csv";
    std::ofstream ofs(out_filename);
    for (int i = 0; i < data_size.size(); i++){
        ofs << i << "," << data_size[i] << "," << total_query_time[i] << "," << result << endl;
    }
    ofs.close();

}

