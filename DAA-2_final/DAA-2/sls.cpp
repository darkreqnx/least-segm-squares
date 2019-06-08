/*!
 * \file
  Class containing all functions required for the Segmented Least Squares Algorithm
*/
#include "sls.h"
#include <bits/stdc++.h>

#define MAXN 1000
#define INF numeric_limits<double>::infinity()

using namespace std;

static vector<Graph> gvector;
static Graph gfirst;
static int Count = 0;
static Graph x;
static stack<Point> hull;

SLS::SLS()
{

}

/*!
 * \brief A vector of type Point(int x, int y) that contains the input points the way they have been entered
 */
static vector <Point> points;

/*!
 * \brief A simple function that compares two points wrt their x-coordinate values
 * \param a is a Point
 * \param b is a Point
 * \return bool, 1 if <, 0 otherwise
 */
bool compareX(Point a, Point b)
{
    if (a.x == b.x)
        return (a.y <= b.y);

    return (a.x <= b.x);
}

/*!
 * \brief Main function that precomputes the error table and finds the optimum cost and number of segments
 * \param cost is a double, denotes the cost of creating a new line segment
 * \return vector<Line> containing all the segments denoted by their starting and ending points
 */
vector<Line> leastSegmSquares(double cost)
{

    int n, i, j, k, interval;
    long long x_sum, y_sum, xy_sum, x2_sum, num, denom;
    double buffer, minimum;

    vector<Point> S;
    Point buff;
    buff.x = points[0].x, buff.y = points[0].y;
    S.push_back(buff);
    for (auto iter = points.begin(); iter != points.end(); iter++)
    {
        S.push_back(*iter);
    }
    n = points.size();
    S.push_back(points[n - 1]);

    // sums_all_x[i] is sum(S[j].x) for 1 <= j <= i
    // sums_all_y[i] is sum(S[j].y) for 1 <= j <= i
    // sums_all_xy[i] is sum(S[j].x * S[j].y) for 1 <= j <= i
    // sums_all_x2[i] is sum(S[j].x * S[j].x) for 1 <= j <= i
    long long sums_all_x[n + 1], sums_all_y[n + 1], sums_all_xy[n + 1], sums_all_x2[n + 1];

    // slopes[i][j] is the slope of the segment that is best fit to
    // the points {S[i], S[i+1], ..., S[j]}

    // intercept[i][j] is the y-intercept of the segment that is best fit to
    // the points {S[i], S[i+1], ..., S[j]}

    // err[i][j] is the square error of the segment that is best fit to
    // the points {S[i], S[i+1], ..., S[j]}
    double slopes[n + 1][n + 1], intercept[n + 1][n + 1], err[n + 1][n + 1];

    // OPT[i] is the optimal solution (minimum cost) for the points {S[1], S[2], ..., S[i]}
    double OPT[n + 1];

    // [opt_segment[i], i] is the last segment in the optimal solution
    // for the points {S[1], S[2], ..., S[i]}
    int opt_segment[n + 1];

    // precompute the error terms
    sums_all_x[0] = sums_all_y[0] = sums_all_xy[0] = sums_all_x2[0] = 0;

    j = 1;
    for (auto iter = S.begin() + 1; iter != S.end() - 1; iter++)
    {
        Point buff = (*iter);
        sums_all_x[j] = sums_all_x[j - 1] + buff.x;
        sums_all_y[j] = sums_all_y[j - 1] + buff.y;
        sums_all_xy[j] = sums_all_xy[j - 1] + buff.x * buff.y;
        sums_all_x2[j] = sums_all_x2[j - 1] + buff.x * buff.x;

        for (i = 1; i <= j; i++)
        {
            interval = j - i + 1;
            x_sum = sums_all_x[j] - sums_all_x[i - 1];
            y_sum = sums_all_y[j] - sums_all_y[i - 1];
            xy_sum = sums_all_xy[j] - sums_all_xy[i - 1];
            x2_sum = sums_all_x2[j] - sums_all_x2[i - 1];

            num = (interval * xy_sum) - (x_sum * y_sum);
            denom = (interval * x2_sum) - (x_sum * x_sum);

            if (num == 0)
                slopes[i][j] = 0.0;
            else
            {
                slopes[i][j] = (denom == 0) ? INF : (num / double(denom));
            }
            intercept[i][j] = (y_sum - slopes[i][j] * x_sum) / double(interval);

            for (k = i, err[i][j] = 0.0; k <= j; k++)
            {
                if (denom != 0)
                    buffer = S[k].y - slopes[i][j] * S[k].x - intercept[i][j];
                else if (denom == 0)
                    buffer = S[k].x + (x_sum / double(interval));

                err[i][j] += buffer * buffer;
            }
        }
        j++;
    }

    // find the cost of the optimal solution
    OPT[0] = opt_segment[0] = 0;
    for (j = 1; j <= n; j++)
    {
        for (i = 1, minimum = INF, k = 0; i <= j; i++)
        {
            buffer = err[i][j] + OPT[i - 1];
            if (buffer < minimum)
            {
                minimum = buffer;
                k = i;
            }
        }
        OPT[j] = minimum + cost;
        opt_segment[j] = k;
    }

    printf("\nCost of the optimal solution : %lf\n", OPT[n]);

    // find the optimal solution
    vector<int> segments;
    for (i = n, j = opt_segment[n]; i > 0; i = j - 1, j = opt_segment[i])
    {
        segments.push_back(i);
        segments.push_back(j);
    }

    vector<Line> edges;

    cout << "Optimal segments: " << endl;
    while (!segments.empty())
        {

            i = segments.back();
            segments.pop_back();
            j = segments.back();
            segments.pop_back();

            struct Line l;

            l.p1= S[i];
            l.p2 = S[j + 1];
            edges.push_back(l);

//            edges.push_back(make_pair(S[i], S[j + 1]));

            cout << "y = (" << slopes[i][j] << ")x + " << intercept[i][j]
                 << " from (" << S[i].x << "," << S[i].y << ")"
                 << " to (" << S[j + 1].x << "," << S[j + 1].y << ")"
                 << "; squared error: " << err[i][j] << endl;
        }

    return edges;
}

/*!
 * \brief Utility function for tokenising and storing cost and points when reading from a txt, not required otherwise when using the GUI
 * \param str_line is a string
 */
void saveTokenize(string str_line)
{
    Point newPoint;

    char line[1024];
    strcpy(line, str_line.c_str());

    char *token = strtok(line, " ");
    char *x = token;
    token = strtok(NULL, " ");
    char *y = token;

    string str_x(x);
    string str_y(y);

    newPoint.x = stoi(str_x);
    newPoint.y = stoi(str_y);

    points.push_back(newPoint);
}

int SLS::SLS_algo(double cost){
    gvector.clear();
    Count =0;

    Graph& g= getGraph();

    gfirst = g;
    Graph gfinal = gfirst;
    points = gfirst.getNormalPoints();

    sort(points.begin(), points.end(), compareX);


    vector<Line> slsLines = leastSegmSquares(cost);

    for(auto line: slsLines){
        gfinal.addBoldLine(line.p1.x, line.p1.y,line.p2.x, line.p2.y);
        cout<<line.p1.x<<" " <<line.p1.y<<" "<<line.p2.x<<" "<< line.p2.y<< endl;
    }
    gfinal.Message = "End of SLS-Algo!!";
    gvector.push_back(gfinal);Count++;
    return Count;
}

int SLS::Iteration(int x){
    if(x<Count){
        Graph& g = getGraph();
        g = gvector[x];
        return 1;
    }else{
        return 0;
    }
}
