// Imagine++ project
// Project:  Panorama
// Author:   Pascal Monasse
// Date:     2013/10/08

#include <Imagine/Graphics.h>
#include <Imagine/Images.h>
#include <Imagine/LinAlg.h>
#include <vector>
#include <sstream>
using namespace Imagine;
using namespace std;

// Record clicks in two images, until right button click
void getClicks(Window w1, Window w2,
               vector<IntPoint2>& pts1, vector<IntPoint2>& pts2) {

    //FIRST PERSONAL CODE
    IntPoint2 clickedPoint;
    Window clickedWindow;
    int pointNumber = 0;
    int subWin = 0;
    int button;

    cout << "Click of one of the two pictures." << endl << "A number will appear beside each clicked point." << endl << "Please make sure that each points labelled with the same number match with the same point of the physical world." << endl << "Also make sure to choose at least 4 points of each picture." << endl;
    button = anyGetMouse(clickedPoint, clickedWindow, subWin);   //Records the type of the click (3 for right click)

    while (button != 3)  //While the user does not right click, the clicks recording process continues
    {
        if (clickedWindow == w1)  //If the user clicks on the first window
        {
            setActiveWindow(w1);
            pts1.push_back(clickedPoint);   //Records the point in the first points vector
            pointNumber = pts1.size();
            drawLine(clickedPoint.x() + 5 , clickedPoint.y() , clickedPoint.x() - 5 , clickedPoint.y() , GREEN , 2);    //Draws a cross on the clicked point
            drawLine(clickedPoint.x() , clickedPoint.y() + 5 , clickedPoint.x() , clickedPoint.y() - 5 , GREEN , 2);
            drawString(clickedPoint.x() + 2 , clickedPoint.y() + 15 , to_string(pointNumber) , GREEN , 12 , 00 , false , true);
            cout << "Recorded point n. " << pointNumber << " located at coordinates x,y = " << clickedPoint << "on the first picture." << endl;
        }
        else if (clickedWindow == w2)  //The process if the same for the second window
        {
            setActiveWindow(w2);
            pts2.push_back(clickedPoint);
            pointNumber = pts2.size();
            drawLine(clickedPoint.x() + 5 , clickedPoint.y() , clickedPoint.x() - 5 , clickedPoint.y() , GREEN , 2);
            drawLine(clickedPoint.x() , clickedPoint.y() + 5 , clickedPoint.x() , clickedPoint.y() - 5 , GREEN , 2);
            drawString(clickedPoint.x() + 2 , clickedPoint.y() + 15 , to_string(pointNumber) , GREEN , 12 , 0 , false , true);
            cout << "Recorded point n. " << pointNumber << " located at coordinates x,y = " << clickedPoint << "on the second picture." << endl;
        }
        else
        {
            cout << "Please click on one of the two pictures !" << endl ;
        }
        button = anyGetMouse(clickedPoint, clickedWindow, subWin);
    }
    cout << "End of the point recording process." << endl << "NEXT STEP" << endl;
    //END
}

// Return homography compatible with point matches
Matrix<float> getHomography(const vector<IntPoint2>& pts1,
                            const vector<IntPoint2>& pts2) {
    size_t n = min(pts1.size(), pts2.size());
    if(n<4) {
        cout << "Not enough correspondences: " << n << endl;
        return Matrix<float>::Identity(3);
    }
    Matrix<double> A(2*n,8);
    Vector<double> B(2*n);
    
    //SECOND PERSONAL CODE
    for (int i = 0 ; i < n ; i++)
    {
        IntPoint2 point1 = pts1[i];
        IntPoint2 point2 = pts2[i];

        //Filling even rows of A according to the courses
        A(2*i,0) = point1.x();
        A(2*i,1) = point1.y();
        A(2*i,2) = 1;
        A(2*i,3) = 0;
        A(2*i,4) = 0;
        A(2*i,5) = 0;
        A(2*i,6) = - point2.x() * point1.x();
        A(2*i,7) = - point2.x() * point1.y();

        //Filling odd rows of A according to the courses
        A(2*i+1,0) = 0;
        A(2*i+1,1) = 0;
        A(2*i+1,2) = 0;
        A(2*i+1,3) = point1.x();
        A(2*i+1,4) = point1.y();
        A(2*i+1,5) = 1;
        A(2*i+1,6) = - point2.y() * point1.x();
        A(2*i+1,7) = - point2.y() * point1.y();

        //Idem for B
        B[2*i] = point2.x();
        B[2*i+1] = point2.y();
    }
    //END

    B = linSolve(A, B);
    Matrix<float> H(3, 3);
    H(0,0)=B[0]; H(0,1)=B[1]; H(0,2)=B[2];
    H(1,0)=B[3]; H(1,1)=B[4]; H(1,2)=B[5];
    H(2,0)=B[6]; H(2,1)=B[7]; H(2,2)=1;

    // Sanity check
    for(size_t i=0; i<n; i++) {
        float v1[]={(float)pts1[i].x(), (float)pts1[i].y(), 1.0f};
        float v2[]={(float)pts2[i].x(), (float)pts2[i].y(), 1.0f};
        Vector<float> x1(v1,3);
        Vector<float> x2(v2,3);
        x1 = H*x1;
        cout << x1[1]*x2[2]-x1[2]*x2[1] << ' '
             << x1[2]*x2[0]-x1[0]*x2[2] << ' '
             << x1[0]*x2[1]-x1[1]*x2[0] << endl;
    }
    return H;
}

// Grow rectangle of corners (x0,y0) and (x1,y1) to include (x,y)
void growTo(float& x0, float& y0, float& x1, float& y1, float x, float y) {
    if(x<x0) x0=x;
    if(x>x1) x1=x;
    if(y<y0) y0=y;
    if(y>y1) y1=y;    
}

// Panorama construction
void panorama(const Image<Color,2>& I1, const Image<Color,2>& I2,
              Matrix<float> H) {
    Vector<float> v(3);
    float x0=0, y0=0, x1=I2.width(), y1=I2.height();

    v[0]=0; v[1]=0; v[2]=1;
    v=H*v; v/=v[2];
    growTo(x0, y0, x1, y1, v[0], v[1]);

    v[0]=I1.width(); v[1]=0; v[2]=1;
    v=H*v; v/=v[2];
    growTo(x0, y0, x1, y1, v[0], v[1]);

    v[0]=I1.width(); v[1]=I1.height(); v[2]=1;
    v=H*v; v/=v[2];
    growTo(x0, y0, x1, y1, v[0], v[1]);

    v[0]=0; v[1]=I1.height(); v[2]=1;
    v=H*v; v/=v[2];
    growTo(x0, y0, x1, y1, v[0], v[1]);

    cout << "x0 x1 y0 y1=" << x0 << ' ' << x1 << ' ' << y0 << ' ' << y1<<endl;

    Image<Color> I(int(x1-x0), int(y1-y0));
    setActiveWindow( openWindow(I.width(), I.height()) );
    I.fill(WHITE);
    
    //THIRD PERSONAL CODE
    Vector<float> xyz(3);    //Vector that will contain the coordinates of the current point 
                            //point studied in the double-loop
    Matrix<float> inv = inverse(H);    //If I is based on I2, we will need the reverse of H
                                       //to compute the coordonates of the point of I1 after
                                       //the application of the homography

    for (int x = 0 ; x < I.width() ; x++)
    {
        for (int y = 0 ; y < I.height() ; y++)
        {
            xyz[0] = x+x0;
            xyz[1] = y+y0;
            xyz[2] = 1;
            if (xyz[0] >= 0 && xyz[0] < I2.width() && xyz[1] >= 0 && xyz[1] < I2.height())  //Careful to the range of I
            {
                I(x,y) = I2(xyz[0] , xyz[1]);
            }

            xyz = inv * xyz;
            xyz /= xyz[2];
            if (xyz[0] >= 0 && xyz[0] < I1.width() && xyz[1] >= 0 && xyz[1] < I1.height())
            {
                I(x,y) = I1.interpolate(xyz[0] , xyz[1]);
            }
        }
    }
    //END
    display(I,0,0);
}

// Main function
int main(int argc, char* argv[]) {
    const char* s1 = argc>1? argv[1]: srcPath("image0006.jpg");
    const char* s2 = argc>2? argv[2]: srcPath("image0007.jpg");

    // Load and display images
    Image<Color> I1, I2;
    if( ! load(I1, s1) ||
        ! load(I2, s2) ) {
        cerr<< "Unable to load the images" << endl;
        return 1;
    }
    Window w1 = openWindow(I1.width(), I1.height(), s1);
    display(I1,0,0);
    Window w2 = openWindow(I2.width(), I2.height(), s2);
    setActiveWindow(w2);
    display(I2,0,0);

    // Get user's clicks in images
    vector<IntPoint2> pts1, pts2;
    getClicks(w1, w2, pts1, pts2);

    vector<IntPoint2>::const_iterator it;
    cout << "pts1="<<endl;
    for(it=pts1.begin(); it != pts1.end(); it++)
        cout << *it << endl;
    cout << "pts2="<<endl;
    for(it=pts2.begin(); it != pts2.end(); it++)
        cout << *it << endl;

    // Compute homography
    Matrix<float> H = getHomography(pts1, pts2);
    cout << "H=" << H/H(2,2);

    // Apply homography
    panorama(I1, I2, H);

    endGraphics();
    return 0;
}
