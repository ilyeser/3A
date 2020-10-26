// Imagine++ project
// Project:  Fundamental
// Author:   Pascal Monasse
// Date:     2013/10/08

#include "./Imagine/Features.h"
#include <Imagine/Graphics.h>
#include <Imagine/LinAlg.h>
#include <vector>
#include <cstdlib>
#include <ctime>
using namespace Imagine;
using namespace std;

static const bool refine = true;  //CHOOSE TO REFINE OR NOT TO REFINE

static const float BETA = 0.01;
static const float LOG_BETA = logf(BETA);

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

// RANSAC algorithm to compute F from point matches (8-point algorithm)
// Parameter matches is filtered to keep only inliers as output.
FMatrix<float,3,3> computeF(vector<Match>& matches) {
    srand(time(0));
    float distMax = 1.5f; // Threshold for inliers/outliers discrimination
    int Niter=100000; // Adjusted dynamically
    FMatrix<float,3,3> bestF;
    vector<int> bestInliers;
    
    // -------------------------------------------------------
    // DO NOT FORGET NORMALIZATION OF POINTS
    
    vector<int> currentInliers;
    vector<int> usedIndices;
    vector<Match> currentMatches;

    vector<bool> markTab(matches.size() , false);

    FMatrix<float , 9 , 9> A;       //Matrix A
    FMatrix<float , 3 , 3> F;

    FMatrix<float , 9 , 9> U;        //Matrix U of the SVD
    FMatrix<float , 9 , 9> Vt;       //Matrix Vt
    FVector<float , 9> S;            //Matrix Sigma

    FMatrix<float , 3 , 3> Uf;       //SVD matrices for F 2-rank forcing
    FMatrix<float , 3 , 3> Vft;
    FVector<float , 3> Sf;
    FMatrix<float , 3 , 3> DSf;

    FMatrix<float , 3 , 3> N(0.);    //Normalization matric
    N(0,0) = 0.001;
    N(1,1) = 0.001;
    N(2,2) = 1;

    cout << "RANSAC algorithm..." << endl;
    int iter = 0;
    while (iter < Niter)   //Main loop of the RANSAC algorithm
    {
    currentMatches.clear();
    usedIndices.clear();


    while (currentMatches.size() < 8)  //choose 8 matches randomly
    {
        int r = intRandom(0,matches.size()-1);
        while (markTab[r])
        {
            r = intRandom(0,matches.size()-1);
        }
        usedIndices.push_back(r);
        markTab[r] = true;
        currentMatches.push_back(matches[r]);
    }
    for (int i = 0 ; i < 8 ; i++)
    {
        markTab[usedIndices[i]] = false;   //re-initialize markTab to false for the next loop
    }

    vector<Match> normalizedMatches;     //Normalization of the coordinates
    for (int i = 0 ; i < 8 ; i++)
    {
        Match m = currentMatches[i];
        m.x1 = 0.001 * m.x1;
        m.x2 = 0.001 * m.x2;
        m.y1 = 0.001 * m.y1;
        m.y2 = 0.001 * m.y2;
        normalizedMatches.push_back(m);
    }
    
    
    for (int i = 0 ; i < 8 ; i++)    //Filling each A_i
    {
        A(i,0) = normalizedMatches[i].x1 * normalizedMatches[i].x2;
        A(i,1) = normalizedMatches[i].x1 * normalizedMatches[i].y2;
        A(i,2) = normalizedMatches[i].x1;

        A(i,3) = normalizedMatches[i].y1 * normalizedMatches[i].x2;
        A(i,4) = normalizedMatches[i].y1 * normalizedMatches[i].y2;
        A(i,5) = normalizedMatches[i].y1;

        A(i,6) = normalizedMatches[i].x2;
        A(i,7) = normalizedMatches[i].y2;
        A(i,8) = 1;
    }
    for (int i = 0 ; i < 9 ; i++)    //Adding a zero line in A
    {
        A(8,i) = 0.;
    }

    svd(A , U , S , Vt);            //SVD decomposition of A

    for (int k = 0 ; k < 3 ; k++)
    {
        for (int l = 0 ; l < 3 ; l++)
        {
            F(k,l) = Vt(7 , 3*k+l);     //Computing F
        }
    }

    //We need now to force F to be 2-rank
    
    svd(F , Uf , Sf , Vft);                //SVD decomposition of F
    DSf = Diagonal(Sf);
    DSf(2,2) = 0;                          //We force sigma_3 to be 0 in the SVD

    F = Uf * (DSf * Vft);                   //then recompute F

    
    F = (N*F)*N;                            //Normalization

    int nbInliers = 0;
    currentInliers.clear();
    for (size_t i = 0 ; i < matches.size() ; i++)       //Inliers/outliers discrimination
    {
            FVector<float, 3> x1(matches[i].x1,matches[i].y1,1);
            FVector<float, 3> x2(matches[i].x2,matches[i].y2,1);
            FVector<float, 3> Fx1 = transpose(F)*x1;
            float Dist = abs((x2*Fx1)/(sqrt(pow(Fx1[0],2)+pow(Fx1[1],2))));

            if (Dist < distMax)
            {
                currentInliers.push_back(i);
                ++nbInliers;
            }
    }

    if (nbInliers > bestInliers.size())                 //Test : do we keep those inliers ?
    {
            bestInliers = currentInliers;
            cout << "Better inliers found, number = " << bestInliers.size() << endl;
            bestF = F;
            float q = 1.0*nbInliers/(1.0*matches.size());

            int nextNiter;
            nextNiter = (LOG_BETA/log(1- pow( q , 8))+1);
            Niter = nextNiter;
            if (Niter < 30){Niter = 30;}
            cout << "New max iterations : " << Niter << endl;
            
    }

    iter = iter + 1;
    } //end of the main loop


    //Refining with least square minimization (the step are basically the same as in the main loop)
    if (refine){
    Matrix<float> A_ls (bestInliers.size() , 9);
    Matrix<float> U_ls (bestInliers.size() , bestInliers.size());
    Matrix<float> Vt_ls (9,9);
    Vector<float> S_ls (9);

    for (int i = 0 ; i < bestInliers.size() ; i++)
    {
        Match normalizedMatch;
        normalizedMatch.x1 = 0.001*matches[bestInliers[i]].x1;
        normalizedMatch.x2 = 0.001*matches[bestInliers[i]].x2;
        normalizedMatch.y1 = 0.001*matches[bestInliers[i]].y1;
        normalizedMatch.y2 = 0.001*matches[bestInliers[i]].y2;
        
        A_ls(i,0) = normalizedMatch.x1 * normalizedMatch.x2;
        A_ls(i,1) = normalizedMatch.x1 * normalizedMatch.y2;
        A_ls(i,2) = normalizedMatch.x1;
        
        A_ls(i,3) = normalizedMatch.y1 * normalizedMatch.x2;
        A_ls(i,4) = normalizedMatch.y1 * normalizedMatch.y2;
        A_ls(i,5) = normalizedMatch.y1;

        A_ls(i,6) = normalizedMatch.x2;
        A_ls(i,7) = normalizedMatch.y2;
        A_ls(i,8) = 1;
    }
    svd(A_ls , U_ls , S_ls , Vt_ls);

    for (int k = 0 ; k < 3 ; k++)
    {
        for (int l = 0 ; l < 3 ; l++)
        {
            bestF(k,l) = Vt_ls(7,3*k+l);
        }
    }

    //Forcing 2-rank for bestF
    svd(bestF , Uf , Sf , Vft);
    DSf = Diagonal(Sf);
    DSf(2,2) = 0;
    bestF = Uf * (DSf * Vft);

    bestF = N * (bestF * N);
    }
    cout << "End of the RANSAC algorithm" << endl;
    cout << "Number of inliers : " << bestInliers.size() << endl;
    //-------------------------------------------------------
    
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
    float a;
    float b;
    while(true) {
        int x,y;
        if(getMouse(x,y) == 3)
            break;
        // ------------------------------------------------------
        int x1, y1, x2, y2;

        FVector<double,3> X;

        X[0] = x;
        X[1] = y;
        X[2] = 1;

        if (x<I1.width()){
            X = transpose(F)*X; //epipolar line of X
            a = -X[0]/X[1];
            b = -X[2]/X[1];
            x1 = 0;
            x2 = I2.width();

            y1 = a*x1 + b;
            y2 = a*x2 + b;
            x1 += I1.width();
            x2 += I1.width();

        }
        else if(x>I1.width()){
            X[0] = x - I1.width();
            X = F*X;
            a = -X[0]/X[1];
            b = -X[2]/X[1];
            x1 = 0;
            x2 = I1.width()-1;

            y1 = a*x1 + b;
            y2 = a*x2 + b;
        }
        Color c(rand()%32 + 128+64+32 , 0 , rand()%50 + 50);
        drawLine(x+5, y, x-5 , y , c , 2);
        drawLine(x, y+5, x , y-5 , c , 2);
        drawLine(x1,y1,x2,y2,c);
    //-------------------------------------------------------
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
    cout << " matches: " << matches.size() << endl;
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
    click();

    // Redisplay without SIFT points
    display(I1,0,0);
    display(I2,w,0);
    displayEpipolar(I1, I2, F);

    endGraphics();
    return 0;
}