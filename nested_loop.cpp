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
using namespace std;

struct Interval
{
    int id, start, end, duration;
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

inline unsigned long long overlapping_query(const vector<vector<IntervalPair>> &grid, const Interval q, float overlap_duration, const vector<int> &col_maxstart,
                       const vector<vector<int>> &cell_minend, const vector<vector<int>> &cell_maxend, const vector<vector<int>> &cell_minstart){
    
    // output pairs whose start is less than q.start and end is higher than threshold_end
    unsigned long long result = 0;
    int index_i = 0;
    int threshold = overlap_duration;
    int threshold_start = q.end - threshold;
    int threshold_end = q.start + threshold;
    auto start_it = std::lower_bound(col_maxstart.begin(), col_maxstart.end(), q.start);
    int start_index = std::distance(col_maxstart.begin(), start_it);
    while(index_i < start_index){
        int maxcell_index = cell_maxend[index_i].size();
        auto cell_end_it = std::lower_bound(cell_maxend[index_i].begin(), cell_maxend[index_i].end(), threshold_end);
        int index_j = std::distance(cell_maxend[index_i].begin(), cell_end_it);
        while((index_j < maxcell_index)){
            if(cell_minend[index_i][index_j] >= threshold_end){
                while(index_j < maxcell_index){
                    for(Interval r : grid[index_i][index_j].second){
                        result += (r.start ^ q.start);
                    }
                    index_j++;
                }
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
        index_i++;
    }

    // normal comparison
    int cell_threshold_start;
    auto col_st_it = std::lower_bound(col_maxstart.begin(), col_maxstart.end(), threshold_start);
    int threshold_start_index = std::distance(col_maxstart.begin(), col_st_it);
    while((start_index < col_maxstart.size()) && (start_index <= threshold_start_index)){
        if(start_index < threshold_start_index){
            auto cell_end_it = std::lower_bound(cell_maxend[start_index].begin(), cell_maxend[start_index].end(), threshold_end);
            int index_k = std::distance(cell_maxend[start_index].begin(), cell_end_it);
            int max_column_size = cell_maxend[start_index].size();
            while((index_k < max_column_size)){
                if(cell_minend[start_index][index_k] >= q.end){
                    while(index_k < max_column_size){
                        for(Interval r : grid[start_index][index_k].first){
                            result += (r.start ^ q.start);
                        }
                        index_k++;
                    }
                    break;
                }else if(cell_minstart[start_index][index_k] > threshold_start){
                    index_k++;
                    continue;
                }
                else{
                    auto it_temp = grid[start_index][index_k].first.begin();
                    auto end_it = grid[start_index][index_k].first.end();                    
                    int min_end = cell_minend[start_index][index_k];
                    cell_threshold_start = min_end - threshold;
                    cell_threshold_start = std::max(q.start, cell_threshold_start);
                    
                    if(min_end < threshold_end){
                        while((it_temp->start <= cell_threshold_start) && (it_temp < end_it)){
                            check_overlap(*it_temp, q, overlap_duration, &result);                                
                            it_temp++;
                        }
                    }else{
                        while((it_temp->start <= cell_threshold_start) && (it_temp < end_it)){
                            result += ((*it_temp).start ^ q.start);
                            it_temp++;
                        }
                    }
                    
                    cell_threshold_start = std::min(cell_maxend[start_index][index_k], q.end) - threshold;
                    while((it_temp < end_it) && (it_temp->start <=  cell_threshold_start)){
                        check_overlap(*it_temp, q, overlap_duration, &result);
                        it_temp++;
                    }
                }
                index_k++;
            }
        }else{
            auto cell_end_it = std::lower_bound(cell_maxend[start_index].begin(), cell_maxend[start_index].end(), threshold_end);
            int index_l = std::distance(cell_maxend[start_index].begin(), cell_end_it);
            int max_column_size = cell_minstart[start_index].size(); 
            while((index_l < max_column_size)){ 
                if (cell_minstart[start_index][index_l] > threshold_start) {
                    index_l++;
                    continue;
                }
                
                auto it_temp = grid[start_index][index_l].first.begin();
                auto end_it = grid[start_index][index_l].first.end();
                int min_end = cell_minend[start_index][index_l];

                if(min_end < threshold_end){                        
                    auto second_end_it = grid[start_index][index_l].second.end();
                    auto cell_end_it = std::lower_bound(grid[start_index][index_l].second.begin(), second_end_it, threshold_end,
                        [](const Interval &a, int b) {     
                        return a.end < b;
                    });
                    while(cell_end_it < second_end_it){
                        check_overlap(*cell_end_it, q, overlap_duration, &result);
                        cell_end_it++;
                    }
                    index_l++;
                    continue;                        
                }else{
                    cell_threshold_start = min_end - threshold;
                    cell_threshold_start = std::max(q.start, cell_threshold_start);
                    cell_threshold_start = std::min(cell_threshold_start, threshold_start); 
                    while((it_temp->start <= cell_threshold_start) && (it_temp < end_it)){
                        result += ((*it_temp).start ^ q.start);
                        it_temp++;
                    }
                }
                
                cell_threshold_start = std::min(cell_maxend[start_index][index_l], q.end) - threshold;
                while((it_temp < end_it) && (it_temp->start <=  cell_threshold_start)){
                    check_overlap(*it_temp, q, overlap_duration, &result);
                    it_temp++;
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
    if (argc < 5) {
        cerr << "Usage: " << argv[0] << " overlap_duration page_size relation_R relation_S" << endl;
        return 1;
    }
    float overlap_duration = std::stof(argv[1]);
    int page_size = std::stoi(argv[2]);
    string file1 = argv[3];
    string file2 = argv[4];
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

    //query & indexing
    vector<int> data_size;
    vector<double> query_time;
    vector<double> total_query_time;
    vector<double> total_result_size;
    double elapsed;
    unsigned long long result = 0;

    double mem = process_mem_usage();
    cout << "ours : making index: slice size is " << relations.size() << endl;

    // aray parameters
    int max_col_height = relations.size() / (page_size * page_size) + 10; //10 is margin
    using IntervalPair = std::pair<std::vector<Interval>, std::vector<Interval>>;
    std::vector<std::vector<IntervalPair>> grid(
        max_col_height,
        std::vector<IntervalPair>(page_size * 10, IntervalPair())
    );
    vector<int> col_minstart;
    vector<int> col_maxstart;
    vector<vector<int>> cell_minend(max_col_height,vector<int>());
    vector<vector<int>> cell_maxend(max_col_height,vector<int>());
    vector<vector<int>> cell_minstart(max_col_height,vector<int>());
    vector<Interval> column;
    vector<Interval> cell;
    
    std::chrono::system_clock::time_point  start, end;
    start = std::chrono::system_clock::now();

    // sort by start, ascending
    std::sort(relations.begin(), relations.end(), [](const Interval &a, const Interval &b) {
        return a.start < b.start;
    });

    // index structuring
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

    mem = process_mem_usage() - mem;
    //std::cerr << "memory_usage_MB=" << mem << "\n";  

    std::sort(queries.begin(), queries.end(), [](const Interval &a, const Interval &b) {
        return a.start < b.start;
    });

    end = std::chrono::system_clock::now();  
    elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end-start).count();

    cout << "ours : making index time is " << elapsed  / 1000000 << " secs" << endl;
    cout << "ours : make ok" << endl;
    
    // query processing
    int query_count = 0;
    result = 0;
    cout << "ours : queries size: " << queries.size() << endl;
    while(query_count < queries.size()){
        Interval q;            
        q = queries[query_count];
        if(q.duration < overlap_duration){
            query_count++;
            continue;
        }

        std::chrono::system_clock::time_point  start, end;
        start = std::chrono::system_clock::now();

        result += overlapping_query(grid, q, overlap_duration, col_maxstart, cell_minend, cell_maxend, cell_minstart);

        end = std::chrono::system_clock::now();  
        elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end-start).count();

        query_time.emplace_back(elapsed);            
        query_count++;
    }

    cout << query_count << " queries are processed" << endl;
    cout << "ours : query processing is end" << endl;

    total_query_time.emplace_back(accumulate(query_time.begin(), query_time.end(), 0.0) / 1000000);
    total_result_size.emplace_back(result);
    data_size.emplace_back(relations.size());

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
    std::string out_filename = "nested_loop_" + base + "_" + oss.str() + ".csv";
    std::ofstream ofs(out_filename);
    for (int i = 0; i < data_size.size(); i++){
        ofs << i << "," << data_size[i] << "," << total_query_time[i] << "," << result << endl;
    }
    ofs.close();

}

