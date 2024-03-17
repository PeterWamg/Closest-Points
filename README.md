# Closest Pair of Points

## Introduction
given an array of n points in the plane, and the problem is to find out the closest pair of points in the array.
## Description
Closest Pair of Points problem revolves around identifying the pair of points with the smallest Euclidean distance among a given set of points in a two-dimensional plane. We are provided with a set of n points, denoted as ｛p1,p2,...,pn｝where each point pi is defined by its(x,y) coordinates.
The task is to find the pair of points (p i,p j) within this set, minimizing their Euclidean distance ||pi-pj||.

![Image Title](01.png)

## Gold
Identify the pair of points with the smallest Euclidean distance within a given set, facilitating efficient spatial analysis and computational geometry applications.

## Soultion

### Brute Force

#### Algorithmic principle
Solve the problem by exhaustively listing all possible solutions, and then find the optimal solution that meets the conditions. Explore every possible combination of points and select the pair with the smallest distance, ensuring a thorough search for the closest pair within the dataset.

#### Pseudocode

```
ClosestPairBruteForce(P)
    minDistance = ∞
    closestPair = None
    
    for each point p1 in P:
        for each point p2 in P:
            if p1 ≠ p2:
                distance = EuclideanDistance(p1, p2)
                if distance < minDistance:
                    minDistance = distance
                    closestPair = (p1, p2)
    
    return closestPair

```
#### Operate
- The algorithm systematically examines all possible pairs of points within the given set.
- For each pair of points, it calculates the Euclidean distance between them.
- It keeps track of the pair with the smallest distance encountered so far.
- After comparing distances between all pairs, it selects the pair with the minimum distance as the optimal solution.

  
#### Algorithm Implement

```
// Function to calculate the Euclidean distance between two points
double distance(const Point& p1, const Point& p2) {
    return sqrt(pow(p1.x - p2.x, 2) + pow(p1.y - p2.y, 2));
}

// Brute Force method to find the closest pair of points
pair<Point, Point> closestPairBruteForce(const vector<Point>& points) {
    double minDistance = numeric_limits<double>::max();
    pair<Point, Point> closestPair;

    for (size_t i = 0; i < points.size(); ++i) {
        for (size_t j = i + 1; j < points.size(); ++j) {
            double dist = distance(points[i], points[j]);
            if (dist < minDistance) {
                minDistance = dist;
                closestPair = make_pair(points[i], points[j]);
            }
        }
    }

    return closestPair;
}

```

#### Time complexity

In the Brute Force method, we need to compare n given points pairwise to calculate the distance between them. Therefore, the total time complexity is O (n<sup>2</sup>), where n is the number of points.

### 1-D version

#### Algorithmic principle
exploit the sorted nature of the points along the one-dimensional space. By sorting the points initially, the algorithm ensures that the closest pair of points will be adjacent to each other in the sorted list, thereby simplifying the process of finding the closest pair. 

#### Pseudocode

```
ClosestPair1D(P)
    Sort points in array P
    minDistance = +∞
    closestPair = None
    for i = 1 to n-1 do
        distance = P[i+1] - P[i]
        if distance < minDistance then
            minDistance = distance
            closestPair = (P[i], P[i+1])
    return closestPair
```
#### Operate
- Sort the given points according to their coordinates in one-dimensional space.
- Initialize the minimum distance to positive infinity.
- Initialize variables to store the closest pair of points.
- Traverse the sorted point set.
     - Calculate the distance between adjacent points.
     - If the distance is less than the minimum distance, update the minimum distance and the closest point pair.
- Return the closest point pair.

  
#### Algorithm Implement

```
// Function to calculate the distance between two points
double distance(const Point& p1, const Point& p2) {
    return abs(p1.x - p2.x);
}

// Function to find the closest pair of points in 1D space
pair<Point, Point> closestPair1D(vector<Point>& points) {
    // Sort points based on their coordinates
    sort(points.begin(), points.end(), [](const Point& p1, const Point& p2) {
        return p1.x < p2.x;
    });

    // Initialize minimum distance to positive infinity
    double minDistance = numeric_limits<double>::infinity();
    pair<Point, Point> closestPair;

    // Iterate through sorted points and find closest pair
    for (size_t i = 0; i < points.size() - 1; ++i) {
        double dist = distance(points[i], points[i + 1]);
        if (dist < minDistance) {
            minDistance = dist;
            closestPair = make_pair(points[i], points[i + 1]);
        }
    }

    return closestPair;
}

```

#### Time complexity

The time complexity of the 1-D version of the Closed Pair of Points algorithm is O (n log n), where n is the number of points.






