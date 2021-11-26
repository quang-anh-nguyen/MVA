// Imagine++ project
// Project:  Fundamental
// Author:   Pascal Monasse

#include "./Imagine/Features.h"
#include <Imagine/Graphics.h>
#include <Imagine/LinAlg.h>
#include <vector>
#include <cstdlib>
#include <ctime>
using namespace Imagine;
using namespace std;

static const float BETA = 0.01f; // Probability of failure

struct Match {
    float x1, y1, x2, y2;
};

// Display SIFT points and fill vector of point correspondences
void algoSIFT(Image<Color,2> I1, Image<Color,2> I2,
              vector<Match>& matches) {
    // Find interest points
    SIFTDetector D;
    D.setFirstOctave(-1);
    Array<SIFTDetector::Feature> feats1 = D.run(I1);
    drawFeatures(feats1, Coords<2>(0,0));
    cout << "Im1: " << feats1.size() << flush;
    Array<SIFTDetector::Feature> feats2 = D.run(I2);
    drawFeatures(feats2, Coords<2>(I1.width(),0));
    cout << " Im2: " << feats2.size() << flush;

    const double MAX_DISTANCE = 100.0*100.0;
    for(size_t i=0; i < feats1.size(); i++) {
        SIFTDetector::Feature f1=feats1[i];
        for(size_t j=0; j < feats2.size(); j++) {
            double d = squaredDist(f1.desc, feats2[j].desc);
            if(d < MAX_DISTANCE) {
                Match m;
                m.x1 = f1.pos.x();
                m.y1 = f1.pos.y();
                m.x2 = feats2[j].pos.x();
                m.y2 = feats2[j].pos.y();
                matches.push_back(m);
            }
        }
    }
}

int computeIterNb(int currentIterNb, int bestInlierNb, int sampleNb, int matchesNb) {
    if (bestInlierNb == 0)
        return currentIterNb;
    int bound = ceil(log(BETA) / log(1 - pow(((float)bestInlierNb) / matchesNb, sampleNb)));
    if (bound>currentIterNb)
        bound = currentIterNb;
    return bound;
}

vector<Match> getSample(int sampleNb, vector<Match> matches) {
    int matchesNb = matches.size();
    if (sampleNb>matchesNb)
        cerr << "Not enough matches to estimate models." << endl;
    int sampleId[sampleNb];
    vector<Match> sampleMatches;
    for (int i=0; i<sampleNb; i++) {
        int chosen = rand() % matchesNb;
        sampleId[i] = chosen;
        sampleMatches.push_back(matches[chosen]);
    }

    return sampleMatches;
}

FMatrix<float,3,3> sampleComputeF(vector<Match>& sampleMatches) {
    if (sampleMatches.size()<8) {
        cerr << "Not sufficient points to compute F." << endl;
    }

    // Normalizaton matrix
    FVector<float,3> N;
    N[0] = 1e-3; N[1] = 1e-3; N[2]=1.;

    // Linear system
    FMatrix<float,9,9> A;
    for (int i=0; i<8; i++) {

        float x1 = sampleMatches[i].x1 * N[0];
        float y1 = sampleMatches[i].y1 * N[1];
        float x2 = sampleMatches[i].x2 * N[0];
        float y2 = sampleMatches[i].y2 * N[1];

        A(i,0) = x1*x2;
        A(i,1) = x1*y2;
        A(i,2) = x1;
        A(i,3) = y1*x2;
        A(i,4) = y1*y2;
        A(i,5) = y1;
        A(i,6) = x2;
        A(i,7) = y2;
        A(i,8) = 1.;
    }

    for (int j=0; j<9; j++)
        A(8,j) = 0.;

    FVector<float,9> S;
    FMatrix<float,9,9> U, Vt;
    svd(A, U, S, Vt);
    FMatrix<float,3,3> F;
    for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++) {
            F(i,j) = Vt.getRow(8)[3*i+j];
        }
    }

    // Enforce rank(F)=2
    FVector<float,3> SF;
    FMatrix<float,3,3> UF, VtF;
    svd(F, UF, SF, VtF);
    SF[2] = 0;

    F = UF * Diagonal(SF) * VtF;
    F = Diagonal(N) * F * Diagonal(N);
    return F;
}

vector<int> getInliers(vector<Match> matches, FMatrix<float,3,3> currentF, float distMax) {

    vector<int> inliers;
    DoublePoint3 point1, point2;
    FVector<float,3> epiline;

    for (int i=0; i<(int)matches.size(); i++) {

        point1[0] = matches[i].x1; point1[1] = matches[i].y1; point1[2] = 1.;
        point2[0] = matches[i].x2; point2[1] = matches[i].y2; point2[2] = 1.;

        // Calculate x'Fx normalized
        epiline = transpose(currentF) * point1;
        float dist = abs(epiline[0]*point2[0] + epiline[1]*point2[1] + epiline[2]) / sqrt(pow(epiline[0], 2) + pow(epiline[1], 2));
        if (dist<distMax) {
            inliers.push_back(i);
        }
    }

    return inliers;
}

// RANSAC algorithm to compute F from point matches (8-point algorithm)
// Parameter matches is filtered to keep only inliers as output.
FMatrix<float,3,3> computeF(vector<Match>& matches) {
    const float distMax = 1.5f; // Pixel error for inlier/outlier discrimination
    int Niter=100000; // Adjusted dynamically
    FMatrix<float,3,3> bestF;
    vector<int> bestInliers;
    // --------------- TODO ------------
    // DO NOT FORGET NORMALIZATION OF POINTS

    const int sampleSize = 8;
    int currentIter = 0;
    FMatrix<float,3,3> currentF;
    vector<int> inliers;
    vector<Match> sampleMatches;
    while (currentIter<Niter) {

        currentIter++;
        sampleMatches = getSample(sampleSize, matches);
        currentF = sampleComputeF(sampleMatches);
        inliers = getInliers(matches, currentF, distMax);
        cout << currentIter << "/" << Niter << "\t" << inliers.size() << " inliers" << endl;

        if (inliers.size()>bestInliers.size()) {
            bestInliers = inliers;
            bestF = currentF;
            if (inliers.size()>50) /* avoid overflow recomputing Niter */ {
                Niter = computeIterNb(Niter, bestInliers.size(), 8, matches.size());
                cout << "Better model with " << inliers.size() << " inliers, reducing to " << Niter << " iterations" << endl;
            }
        }
    }

    cout << "Best model with " << bestInliers.size() << " inliers" << endl;

    // Refine F with all inliers
    int inliersNb = bestInliers.size();
    Matrix<float> A(inliersNb, 9);
    FVector<float,3> N;
    N[0] = 1e-3; N[1] = 1e-3; N[2] = 1.;
    for (int i=0; i<inliersNb; i++) {

        float x1 = matches[bestInliers[i]].x1 * N[0];
        float y1 = matches[bestInliers[i]].y1 * N[1];
        float x2 = matches[bestInliers[i]].x2 * N[0];
        float y2 = matches[bestInliers[i]].y2 * N[1];

        A(i,0) = x1*x2;
        A(i,1) = x1*y2;
        A(i,2) = x1;
        A(i,3) = y1*x2;
        A(i,4) = y1*y2;
        A(i,5) = y1;
        A(i,6) = x2;
        A(i,7) = y2;
        A(i,8) = 1.;
    }

    Vector<float> S(9);
    Matrix<float> U(inliersNb, 9), Vt(9, 9);
    svd(A, U, S, Vt);

    for (int i=0; i<3; i++) {
        for (int j=0; j<3; j++) {
            bestF(i,j) = Vt.getRow(8)[3*i+j];
        }
    }

    FVector<float,3> SF;
    FMatrix<float,3,3> UF, VtF;
    svd(bestF, UF, SF, VtF);
    SF[2] = 0;
    bestF = UF * Diagonal(SF) * VtF;
    bestF = Diagonal(N) * bestF * Diagonal(N);


    // Updating matches with inliers only
    vector<Match> all=matches;
    matches.clear();
    for(size_t i=0; i<bestInliers.size(); i++)
        matches.push_back(all[bestInliers[i]]);
    return bestF;
}

// Expects clicks in one image and show corresponding line in other image.
// Stop at right-click.
void displayEpipolar(Image<Color> I1, Image<Color> I2,
                     const FMatrix<float,3,3>& F) {
    while(true) {
        int x,y;
        if(getMouse(x,y) == 3)
            break;
        // --------------- TODO ------------
        int w = I1.width();
        Color c(rand()%256, rand()%256, rand()%256);
        if (x < w) {
            DoublePoint3 point1;
            point1[0] = x; point1[1] = y; point1[2] = 1.;
            drawCircle(x, y, 5, c);
            FVector<float,3> line = transpose(F) * point1;
            drawLine(w, -line[2]/line[1], 2*w, (-line[2]-line[0]*w)/line[1], c);
        } else {
            DoublePoint3 point2;
            point2[0] = x - w; point2[1] = y; point2[2] = 1.;
            drawCircle(x, y, 5, c);
            FVector<float,3> line = F * point2;
            drawLine(0, -line[2]/line[1], w, (-line[2]-line[0]*w)/line[1], c);
        }
    }
}

int main(int argc, char* argv[])
{
    srand((unsigned int)time(0));

    const char* s1 = argc>1? argv[1]: srcPath("im1.jpg");
    const char* s2 = argc>2? argv[2]: srcPath("im2.jpg");

    // Load and display images
    Image<Color,2> I1, I2;
    if( ! load(I1, s1) ||
        ! load(I2, s2) ) {
        cerr<< "Unable to load images" << endl;
        return 1;
    }
    int w = I1.width();
    openWindow(2*w, I1.height());
    display(I1,0,0);
    display(I2,w,0);

    vector<Match> matches;
    algoSIFT(I1, I2, matches);
    const int n = (int)matches.size();
    cout << " matches: " << n << endl;
    drawString(100,20,std::to_string(n)+ " matches",RED);
    click();
    
    FMatrix<float,3,3> F = computeF(matches);
    cout << "F="<< endl << F;

    // Redisplay with matches
    display(I1,0,0);
    display(I2,w,0);
    for(size_t i=0; i<matches.size(); i++) {
        Color c(rand()%256,rand()%256,rand()%256);
        fillCircle(matches[i].x1+0, matches[i].y1, 2, c);
        fillCircle(matches[i].x2+w, matches[i].y2, 2, c);        
    }
    drawString(100, 20, to_string(matches.size())+"/"+to_string(n)+" inliers", RED);
    click();


    // Redisplay without SIFT points
    display(I1,0,0);
    display(I2,w,0);
    displayEpipolar(I1, I2, F);

    endGraphics();
    return 0;
}
