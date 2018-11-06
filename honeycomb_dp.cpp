/**
 * For compiling -->
 *      g++ -std=c++11  honeycomb_dp.cpp -o honeycomb_dp
 *
 * For running -->
 *      ./honeycomb_dp ./input.csv 15 ./output.csv
 *      Here
 *          @arg1 is input file name
 *          @arg2 is max_radius (i.e number of bounding circles/rings)
 *          @arg2 is output file name
 *
 * Python execution:
 *      import os
 *      os.system('./honeycomb_dp ./input.csv 15 output.csv')
 */

#include <iostream>
#include <vector>
#include <string>
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <sstream>
#include <iterator>
#include <fstream>

using namespace std;

typedef vector<int>               integer_vector;
typedef vector<vector<int> >      integer_vector_vector;
typedef vector<double>              vd;
typedef vector<vector<double> >     vvd;

#define PB          push_back
#define SZ(x)       (int)x.size()


template<typename Out>
void split(const std::string &s, char delim, Out result) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        *(result++) = stoi(item);
    }
}

std::vector<int> split(const std::string &s, char delim) {
    std::vector<int> elems;
    split(s, delim, std::back_inserter(elems));
    return elems;
}

void read_csv(string filename,integer_vector_vector &v)
{
    string str;
    ifstream fin;
    fin.open(filename);
    while(getline(fin, str)){
        v.push_back(split(str, ','));
    }
    fin.close();
}

void convert_to_2D(integer_vector_vector &s,integer_vector_vector &t)
{
//    v[0] = # (Green)
//    v[1] = Y (Red)
//    v[2] = X (Blue)
//    v[3] = value

    int xMin = -579;
    int xMax = 1399;
    int yMin = -2370;
    int yMax = -433;
//    xMin = xMax = s[0][2];
//    yMin = yMax = s[0][1];
//
//    for (int i = 0; i < (int)s.size(); i++){
//        xMin = min(xMin, s[i][2]);
//        xMax = max(xMax, s[i][2]);
//        yMin = min(yMin, s[i][1]);
//        yMax = max(yMax, s[i][1]);
//    }

//    cout << "xMin" << xMin << endl;
//    cout << "xMax" << xMax << endl;
//    cout << "yMin" << yMin << endl;
//    cout << "yMax" << yMax << endl;
  integer_vector_vector v(xMax - xMin + 1, vector<int>(yMax - yMin + 1, 0));
    for (int i = 0; i < (int)s.size(); i++){
        int x = abs(s[i][2] - xMax);
        int y = s[i][1] - yMin;
//        v[x][y] = s[i][3];
        v[x][y] = 1;
    }
    t = v;
}


void aggr_rect(integer_vector_vector &arr)
{
    for (int i = 1; i < SZ(arr); i++){
        int cur = 0;
        for (int j = 1; j < SZ(arr[0]); j++){
            cur += arr[i][j];
            arr[i][j] = arr[i-1][j] + cur;
        }
    }
}

void aggr_upper(integer_vector_vector &arr)
{
    for (int j = 1; j < SZ(arr[0]); j++){
        int cur = 0;
        for (int i = 1; (i < SZ(arr)) &&  (i <= j); i++) {
            cur += arr[i][j];
            arr[i][j] = arr[i-1][j-1] + cur;
        }
    }
}

void aggr_lower(integer_vector_vector &arr)
{
    for (int i = 1; i < SZ(arr); i++){
        int cur = 0;
        for (int j = 1; (j < SZ(arr[0])) && (j <= i); j++){
            cur += arr[i][j];
            arr[i][j] = arr[i-1][j-1] + cur;
        }
    }
}

// computes the sum of a square of side length r
int calc_sqr(int lx, int ly, int r,integer_vector_vector &vc_rect)
{
    return  vc_rect[lx][ly] - vc_rect[lx - r][ly] - vc_rect[lx][ly - r] + vc_rect[lx - r][ly - r];
}

// computes lower triangle of side length r just below the diagonal
int calc_lower_tri(int lx, int ly, int r,integer_vector_vector &vc_rect,integer_vector_vector &vc_lower)
{
    return vc_lower[lx][ly] - vc_lower[lx - r][ly -r] - (vc_rect[lx][ly - r] - vc_rect[lx - r][ly - r]); // remove +1
}

// computes upper triagle of (r - 1) in any location of grid
int calc_upper_tri(int lx, int ly, int r,integer_vector_vector &vc_rect,integer_vector_vector &vc_upper,integer_vector_vector &vc_lower)
{
    if (lx <= ly){   // triangle is above the diagonal
        return vc_upper[lx][ly] - vc_upper[lx - r][ly - r] - (vc_rect[lx - r][ly] - vc_rect[lx - r][ly - r]);
    }
    // upper triangle is either below the diagonal or intersecting the diagonal
    return calc_sqr(lx, ly, r, vc_rect) - calc_lower_tri(lx, ly - 1, r - 1, vc_rect, vc_lower);
}

/*
 * computes the sum surrounded by r rows, r = 1 will look at the cells below (inlcuding it's own cell):
 * 00000
 * 01110
 * 01x10
 * 01110
 * 00000
 */
int calc_point_sum(int x, int y, int r,integer_vector_vector &vc_rect,integer_vector_vector &vc_upper,integer_vector_vector &vc_lower)
{
    int ret = 0;

    // get the value of bounding square
    ret = calc_sqr(x + r, y + r, 2 * r + 1, vc_rect);

    // Subtract the upper triangle upper right corner
    ret -= calc_upper_tri(x - 1, y + r, r, vc_rect, vc_upper, vc_lower);

    // Subtract lower left lower_triangle (equivalent to subtracting rectangle which
    //      contains tower_tri then adding the upper triangle of that rectangle)
    ret -= calc_sqr(x + r, y - 1, r, vc_rect);
    ret += calc_upper_tri(x + r - 1, y - 1, r - 1, vc_rect, vc_upper, vc_lower);

    return ret;
}

void convolution_hcdp(integer_vector_vector &a, int MX_RAD, int rad,integer_vector_vector &vc_rect,integer_vector_vector &vc_upper,integer_vector_vector &vc_lower,integer_vector_vector &res)
{
    // calculating the counts of each point
    res.resize(SZ(a));
    for(int i = 0; i < SZ(a); i++)
        for (int j = 0; j < SZ(a[0]); j++)
            res[i].PB(calc_point_sum(i + MX_RAD + 1, j + MX_RAD + 1, rad, vc_rect, vc_upper, vc_lower));
}

void init_mat(integer_vector_vector &a, int MX_RAD,integer_vector_vector &vc_rect,integer_vector_vector &vc_upper,integer_vector_vector &vc_lower)
{
    int sz_x, sz_y;
    sz_x= SZ(a) + 2 * MX_RAD;
    sz_y = SZ(a[0]) + 2 * MX_RAD;

    // Initialize arrays
  integer_vector_vector vT (sz_x + 2, vector<int>(sz_y + 2, 0));
    for (int i = 0; i < SZ(a); i++)
        for (int j = 0; j < SZ(a[0]); j++)
            vT[i + MX_RAD + 1][j + MX_RAD + 1] = a[i][j];

    vc_upper = vc_lower = vc_rect = vT;

    // Calculate cumulative sums
    aggr_rect(vc_rect);
    aggr_upper(vc_upper);
    aggr_lower(vc_lower);
}

double compute_avg(integer_vector_vector &v,integer_vector_vector &a,integer_vector_vector &c, int pen_cnt)
{
    // Sum of cells which has penguin
    int cCnt_occupied = 0;
    for (int i = 0; i < SZ(a); i++)
        for (int j = 0; j < SZ(a[0]); j++) {
            if(a[i][j] && c[i][j])
                cCnt_occupied += (c[i][j] - 1);
        }

    // Sum of all input cells
    int xMin, xMax, yMin, yMax, cCnt_existing = 0;
    xMin = xMax = v[0][2];
    yMin = yMax = v[0][1];
    for (int i = 0; i < (int)v.size(); i++){
        xMin = min(xMin, v[i][2]);
        xMax = max(xMax, v[i][2]);
        yMin = min(yMin, v[i][1]);
        yMax = max(yMax, v[i][1]);
    }
    for (int i = 0; i < (int)v.size(); i++){
        int x = abs(v[i][2] - xMax);
        int y = v[i][1] - yMin;
        if (c[x][y])
            cCnt_existing += (c[x][y] - 1);
    }

    // Sum of entire grid cells
    long long cCnt_nbr = 0;
    for (int i = 0; i < SZ(c); i++)
        for (int j = 0; j < SZ(c[0]); j++)
            if(c[i][j])
                cCnt_nbr += (c[i][j] - 1);

//    cout << "occ: " << cCnt_occupied << " inp: " << cCnt_existing << " nbr: " << cCnt_nbr << endl;
//
//    cout << "occ/tot: " << cCnt_occupied/(double)SZ(v) << " occ/pen: " << cCnt_occupied/(double)pen_cnt << endl;
//    cout << "inp/tot: " << cCnt_existing/(double)SZ(v) << " inp/pen: " << cCnt_existing/(double)pen_cnt << endl;
//    cout << "nbr/tot: " << cCnt_nbr/(double)SZ(v) << " nbr/pen: " << cCnt_nbr/(double)pen_cnt << endl;
    return cCnt_occupied/(double)pen_cnt;
}

void compute(integer_vector_vector &v,integer_vector_vector &a, int MX_RAD, string filename)
{
  integer_vector_vector vc_rect, vc_upper, vc_lower;
    init_mat(a, MX_RAD, vc_rect, vc_upper, vc_lower);

    double pen_cnt = 0;
    for(int i = 0; i < SZ(v); i ++)
        if(v[i][3] == 1)
            pen_cnt ++;

//    cout << "Notion: \n\tNumber of input cells (tot): " << SZ(v) << endl;
//    cout << "\tNumber of Penguins (pen): " << pen_cnt << endl;
//    cout << "\tSum of cells which has penguin: \"occ\" \n\tSum of all input cells: \"inp\" \n\tSum of entire grid cells: \"nbr\""<< endl;

    string str;
    ofstream ofs;
    ofs.open (filename, std::ofstream::out | std::ofstream::app);

    double pre_avg = 0.0, cur_avg = 0.0;
    ofs << "0.0";
    for (int rad = 1; rad <= MX_RAD; ++rad) {
      integer_vector_vector c;
        convolution_hcdp(a, MX_RAD, rad, vc_rect, vc_upper, vc_lower, c);
//        cout << "\n****** STATS R" << rad << " ******" << endl;
        cur_avg = compute_avg(v, a, c, pen_cnt);
        ofs << ", " << cur_avg - pre_avg;

        pre_avg = cur_avg;
    }
    ofs << endl;
}

int main(int argc, char *argv[])
{
    if(argc < 2){
        cout << "Please Provide input file name!\nAborting..." << endl;
        return 0;
    }

    int MX_RAD;
    if(argc < 3){
        cout << "Max Radius not set! Setting max radius size to 1" << endl;
        MX_RAD = 1;
    }
    else MX_RAD = atoi(argv[2]);

    if (argc < 4){
        cout << "Output file Name name not available!\nAborting..." << endl;
        return 0;
    }

//    freopen(argv[1], "r", stdin);
//    double x;
//  integer_vector_vector a;
//
//    int rows, cols;
//    cin >> rows >> cols;
//    a.resize(rows);
//    for (int i = 0; i < rows; i++)
//        for (int j = 0; j < cols; j++){
//            cin >> x;
//            a[i].PB(x);
//        }

//    cout << "Reading Inputs from file... ";


    string next_file;
    ifstream file_list;
    file_list.open(string(argv[1]));
    while (getline(file_list, next_file)){
    cout << "processing " << next_file << endl;

  integer_vector_vector v, a;
    read_csv(next_file, v);
//    for ( const std::vector<int> &v1 : v )
//    {
//        for ( int x : v1 ){ std::cout << x << '|';}
//        std::cout << std::endl;
//    }
//    cout << "Finished." << endl;

//       for ( const std::vector<int> &v1 : v )
//    {
//        for ( int x : v1 ){ std::cout << x << '|';}
//        std::cout << std::endl;
//    }

    convert_to_2D(v, a);

    compute(v, a, MX_RAD, argv[3]);

//    cout << "Output Array:" << endl;
//    for (int i = 0; i < SZ(c); i++){
//        for (int j = 0; j < SZ(c[0]); j++)
//            cout << c[i][j] << "\t";
//        cout << endl;
//    }
       };
    return 0;
}
