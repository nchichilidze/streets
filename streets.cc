#include <iostream>
#include <string>
#include <list>
#include <vector>
#include <queue>
#include <utility>
#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <math.h>
#include <algorithm>
#include <iomanip>
#include <limits>


using namespace std;

struct Vertex {
public:
    Vertex() : x(-1), y(-1) {}
    Vertex(float newX, float newY) : x(newX), y(newY) {}
    float x; float y;
};

int order, size;

/* graph structure built with these 2 vectors */
vector <list <int> > network;
vector <Vertex> vertices;

/* queries */
void query1();
void query2();
void query3();
void query4(int v);
void query5(int v, int h);
void query6(int s, int e);
void query7(int s, int e);

/* other functions*/
void readGraphFile(string filename);
void readXYZFile(string filename);
void execute(int query);
void maxAndAverage(string ma);
string printList(list<int> l, bool sort);
list<int> retracePath(int v1, int v2, vector<int> p);
float distance (Vertex v1, Vertex v2);
float pathLength(list<int> path);


/* output order and size */
void query1() {
    cout << "n= " << order << "; m= " << size << "." << endl;
}

/* output max degree */
void query2() {
    maxAndAverage("max");
}

/* output average degree */
void query3() {
    maxAndAverage("avg");
}

/* direct neighbours */
void query4(int v) {
    list<int> neighbours = network[v];
    cout <<"N(" << v << ")=" << printList(neighbours, true) << "." << endl;
}

/* neighbours h hops away */
void query5(int v, int h) {
    list<int> ngb;
    list<int>::const_iterator it;
    list<int>::const_iterator nit;
    ngb = network[v];
    for (int i = 1; i < h; i++) {
        list<int> temp;
        for(it = ngb.begin(); it != ngb.end(); it++) { /* it.over current vertices */
            list<int> cngb = network[*it];
            for (nit = cngb.begin(); nit != cngb.end(); nit++) { /* it.over neighbours */
                if(!(find(temp.begin(), temp.end(), *nit) != temp.end()))
                    temp.push_back(*nit);
            }
        }
        ngb = temp;
    }
    cout << "N(" << v << "," << h << ")= " + printList(ngb,true) << "."<<endl;
}

/* wighted shortest path: Dijkstra algorithm implementation */
typedef pair<float,int> ifpair;
void query6(int v1, int v2) {
    priority_queue< ifpair, vector<ifpair> , greater<ifpair> > pq;
    vector<float> dist(order + 1 , std::numeric_limits<float>::max());
    vector<int> prev(order + 1 , 0);
    pq.push(make_pair (0.0,v1));
    dist[v1] = 0.0;
    list <int>::iterator i;
    /* looping until pq is empty */
    while (!pq.empty()) {
        int u = pq.top().second;
        pq.pop();
        list<int> ngb = network[u];
        for (i = ngb.begin(); i != ngb.end(); i++) {
            int v = *i;
            float weight = distance(vertices[u], vertices[v]);
            if (dist[v] > dist[u] + weight) {
                dist[v] = dist[u] + weight;
                pq.push(make_pair(dist[v], v));
                prev[v] = u;
            }
        }
    }
    list<int> path = retracePath (v1, v2, prev);
    if (path.empty()) cout << "Could not find path. " << endl;
    else {
        // cout << printList(path) << endl;
        float length = pathLength(path);
        cout << "d(" << v1 << "," << v2 <<")= "<< length << setprecision(6) << fixed << "."<<endl;
    }
}

/* unweighted shortest path: BFS algorithm implementation */
void query7(int v1, int v2) {
    list <int> q; q.push_back(v1);
    list<int>::const_iterator it;
    vector <int> prev (order + 1, 0);
    vector <bool> visited (order + 1, false);
    visited[v1] = true;
    while (!q.empty()) {
        int node = q.front(); q.pop_front();
        list <int> ngb = network[node];
        for (it = ngb.begin(); it != ngb.end(); it++) {
            int next = *it;
            if (!visited[next]) {
                q.push_back(next);
                prev[next] = node; visited[next] = true;
            }
            if (next == v2) break;
        }
    }
    list <int> path = retracePath(v1,v2,prev);
    if (path.empty()) cout << "Could not find path. " << endl;
    else {
       // cout << printList(path, false) << endl;
        cout << "ed(" << v1 << "," << v2 <<")= " << path.size() - 1 <<"." << endl;
    }
}


int main (int argc, char *argv[]) {
    string fileName; cin >> fileName;
    readGraphFile(fileName + ".osm.graph");
    readXYZFile(fileName + ".osm.xyz");
    int query; cin >> query; execute(query);
    exit(0);
}

void maxAndAverage(string ma) {
    int maxVertex = -1; int maxDegree = -1;
    float avgDegree = 0; float sum = 0;
    for (int i = 1; i < network.size(); i++) {
        int currentVertex = i;
        int currentDegree = network[currentVertex].size();
        if (currentDegree > maxDegree) {
            maxDegree = currentDegree; maxVertex = currentVertex;
        }
        sum += currentDegree;
    }
    
    avgDegree = float (sum / order);
    
    if (ma == "max")
        cout << "v=" << maxVertex <<"; |N(v)|=" << maxDegree << "." << endl;
    else if (ma == "avg")
        cout << "avg |N(v)|= " << fixed << setprecision(6) << avgDegree <<"." << endl;
}

string printList(list<int> lon, bool sort) {
    string result;
    if (sort) lon.sort();
    list<int>::const_iterator lit_1;
    for(lit_1 = lon.begin(); lit_1 != lon.end(); lit_1++) {
        result += " " + to_string(*lit_1)  ;
    }
    return result;
}

list<int> retracePath (int v1, int v2, vector<int> prev) {
    list <int> path;
    for (int i = v2; i != 0; i = prev[i]) {
        path.push_back(i);
    }
    path.reverse();
    if (path.front() != v1) {
        list<int> emptyList;
        return emptyList;
    } else { return path; }
}

void readGraphFile(string fileName) {
    int n, m;
    ifstream file;
    file.open(fileName);
    //if (file == NULL) cout << "Error: can not read graph file." << endl;
    string line;
    while (getline(file,line)) {
        if (isdigit(line[0])) {
            istringstream s(line);
            s >> n >> m; break;
        }
    }
    order = n; size = m;
    list<int> emptyList;
    vector <list <int> > graph (order + 1, emptyList);
    network = graph;
    int v = 1; /* initialize vertex & neighbour */
    int ngb;
    while(getline(file,line)) {
        istringstream s(line);
        list <int> *list_ptr = &network[v];
        while (s >> ngb) { /* add edge (v, ngb) */
            if (!(find(list_ptr->begin(), list_ptr->end(), ngb) != list_ptr->end()))
                list_ptr->push_back(ngb);
        }
        v++;
    }
    file.close();
}

void readXYZFile (string fileName) {
    vector <Vertex> result;
    ifstream xyzFile;
    xyzFile.open(fileName);
    //fstream xyzFile(fileName);
   // if (xyzFile == NULL) cout << "Error: can not read xyz file. " << endl;
    float x = 0; float y = 0; float z = 0;
    result.push_back(*new Vertex(-1,-1));
    while (xyzFile >> x >> y >> z) {
        result.push_back(*new Vertex(x,y));
    }
    xyzFile.close();
    vertices = result;
}

void execute(int query) {
    switch(query) {
        case 1: query1(); break;
        case 2: query2(); break;
        case 3: query3(); break;
        case 4: { int v; cin >> v; query4(v); break; }
        case 5: { int v; int h; cin >> v ; cin >> h; query5(v,h); break; }
        case 6: { int s; int e; cin >> s ; cin >> e; query6(s,e); break; }
        case 7: { int s; int e; cin >> s ; cin >> e; query7(s,e); break; }
    }
}


float distance (Vertex v1, Vertex v2) {
    float x1 = v1.x; float x2 = v2.x;
    float y1 = v1.y; float y2 = v2.y;
    return sqrt ( (x2 - x1) * (x2-x1) + (y2-y1) * (y2-y1) );
}

float pathLength(list<int> path) {
    vector<int> result;
    float totDist = 0.0;
    float curDist = 0.0;
    list<int>::const_iterator it;
    for (it = path.begin() ; it != path.end(); it++) {
        result.push_back(*it);
    }
    for (int i = 0; i < result.size() - 1; i++) {
        Vertex v1 = vertices[result[i]];
        Vertex v2 = vertices[result[i+1]];
        curDist = distance(v1,v2);
        totDist = totDist + curDist;
    }
    cout << setprecision(6) << fixed;
    return totDist;
}

/*
 float distance(int s, int e) {
 Vertex *v1 = &vertices[s];
 Vertex *v2 = &vertices[e];
 int x1 = v1->getX(); int y1 = v1->getY();
 int x2 = v2->getY(); int y2 = v2->getX();
 return sqrt((x2 - x1) * (x2-x1) + (y2-y1) * (y2-y1));
 }
 */
