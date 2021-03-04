#ifndef FUNCTIONS
#define FUNCTIONS

// std header
#include <iostream>
#include <string>

// pcl header
#include <pcl/register_point_struct.h>
#include <pcl/register_point_struct.h>
#include <pcl/io/pcd_io.h>
#include <pcl/io/vtk_lib_io.h>
#include <pcl/filters/voxel_grid.h>
#include <pcl/point_types.h>
#include <pcl/console/parse.h>
#include <pcl/visualization/vtk.h>

// opencv header
#include <opencv2/core/utility.hpp>
#include <opencv2/core/utility.hpp>
#include <opencv2/imgproc/imgproc_c.h> // cvFindContours
#include <opencv2/imgproc.hpp>
#include <opencv2/rgbd.hpp>
#include <opencv2/videoio.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/core.hpp>

// vtk header
#include <vtkGeneralTransform.h>
#include <vtkTransform.h>
#include <vtkCellLocator.h>

#include "boost.h"

using namespace std;
using namespace pcl;

#define EPS 0.00001
const double pi = 3.14159265358979f;

struct ScanParameters
{
    double fov_vert;          // vertical fov in radian
    double fov_hor;           // horizontal fov in radian
    int    res_vert;          // vertical resolution
    int    res_hor;           // horizontal  resolution
    double distance;          // distance in mm
    double screen_dist;       // screen distance in mm
    int    rotInterval;       // rotation interval in degree
};

class Timer
{
public:
    Timer() : start_(0), time_(0) {}

    void start()
    {
        start_ = cv::getTickCount();
    }

    void stop()
    {
        CV_Assert(start_ != 0);
        int64 end = cv::getTickCount();
        time_ += end - start_;
        start_ = 0;
    }

    double time()
    {
        double ret = time_ / cv::getTickFrequency();
        time_ = 0;
        return ret;
    }

private:
    int64 start_, time_;
};

// Copy of cv_mouse from cv_utilities
class Mouse
{
public:
    static void start(const std::string& a_img_name)
    {
        cv::setMouseCallback(a_img_name.c_str(), Mouse::cv_on_mouse, 0);
    }
    static int event(void)
    {
        int l_event = m_event;
        m_event = -1;
        return l_event;
    }
    static int x(void)
    {
        return m_x;
    }
    static int y(void)
    {
        return m_y;
    }

private:
    static void cv_on_mouse(int a_event, int a_x, int a_y, int, void*)
    {
        m_event = a_event;
        m_x = a_x;
        m_y = a_y;
    }

    static int m_event;
    static int m_x;
    static int m_y;
};
int Mouse::m_event;
int Mouse::m_x;
int Mouse::m_y;


void PrintUsage(){
    cout<<"Usage: ./scanner [options] model.ply viewpoints.ply"<<endl;
    cout<<"[options]"<<endl;
    cout<<"\t-dist <double> : scan distance in mm (default:3000)"<<endl;
    cout<<"                   or put three numbers for multiple distances [start] [interval] [#]"<<endl;
    cout<<"\t-rot  <int>    : rotation interval in degree (default:10)"<<endl;
    cout<<"\t-path <string> : directory path for result files (default: .)"<<endl;
    cout<<"\t-cv   [depth/color/both] : print openCV matrix in YML"<<endl;
    cout<<"\t-ply  [0/1]    : print point cloud w/ color in PLY [0:object fixed/1:camera fixed]"<<endl;
    cout<<"\t-pcd  [0/1]    : quick print of point cloud w/o color in PCD [0:object fixed/1:camera fixed]"<<endl;
    exit(1);
}
void PrintKeyUsage(){
    cout<<"[Key Usage]"<<endl;
    cout<<"\tw : write template"<<endl;
    cout<<"\tv : to next view"<<endl;
    cout<<"\tf : show the number of failed templates"<<endl;
}

void EstimateColor(vtkPolyData* poly, vtkIdType cellId, double* point, double* rgb){
    double a[3], b[3], c[3];
    poly->GetPoint(poly->GetCell(cellId)->GetPointId(0), a);
    poly->GetPoint(poly->GetCell(cellId)->GetPointId(1), b);
    poly->GetPoint(poly->GetCell(cellId)->GetPointId(2), c);
    double ab[3], ac[3], ba[3], bc[3];
    vtkMath::Subtract(b, a, ab);
    vtkMath::Subtract(c, a, ac);
    vtkMath::Subtract(a, b, ba);
    vtkMath::Subtract(c, b, bc);

    double n0[3], n[3], end[3], l[3], d, l0p0[3];
    //a (p0=b, l0=a, l=ap, n0=(baXbc), n=n0Xbc)
    vtkMath::Subtract(point, a, l); //l=ap
    vtkMath::Cross(ba, bc, n0);     //n0
    vtkMath::Cross(n0, bc, n);      //n
    vtkMath::Subtract(b, a, l0p0);  //l0p0=p0-l0
    d = vtkMath::Dot(l0p0, n)/vtkMath::Dot(l, n);
    vtkMath::MultiplyScalar(l, d);
    vtkMath::Add(a, l, end);
    double a_w = sqrt(vtkMath::Distance2BetweenPoints(point, end)/vtkMath::Distance2BetweenPoints(a, end));
    //b (p0=a, l0=b, l=bp, n0=(abXac), n=n0Xac)
    vtkMath::Subtract(point, b, l); //l=bp
    vtkMath::Cross(ab, ac, n0);     //n0
    vtkMath::Cross(n0, ac, n);      //n
    vtkMath::Subtract(a, b, l0p0);  //l0p0=p0-l0
    d = vtkMath::Dot(l0p0, n)/vtkMath::Dot(l, n);
    vtkMath::MultiplyScalar(l, d);
    vtkMath::Add(b, l, end);
    double b_w = sqrt(vtkMath::Distance2BetweenPoints(point, end)/vtkMath::Distance2BetweenPoints(b, end));
    //c (p0=a, l0=c, l=cp, n0=(baXbc), n=n0Xba)
    vtkMath::Subtract(point, c, l); //l=cp
    vtkMath::Cross(ba, bc, n0);     //n0
    vtkMath::Cross(n0, ba, n);      //n
    vtkMath::Subtract(a, c, l0p0);  //l0p0=p0-l0
    d = vtkMath::Dot(l0p0, n)/vtkMath::Dot(l, n);
    vtkMath::MultiplyScalar(l, d);
    vtkMath::Add(c, l, end);
    double c_w = sqrt(vtkMath::Distance2BetweenPoints(point, end)/vtkMath::Distance2BetweenPoints(c, end));

    double rgb_a[3] = {poly->GetAttributes(0)->GetArray(0)->GetComponent(poly->GetCell(cellId)->GetPointId(0),0),
                       poly->GetAttributes(0)->GetArray(0)->GetComponent(poly->GetCell(cellId)->GetPointId(0),1),
                       poly->GetAttributes(0)->GetArray(0)->GetComponent(poly->GetCell(cellId)->GetPointId(0),2)};
    double rgb_b[3] = {poly->GetAttributes(0)->GetArray(0)->GetComponent(poly->GetCell(cellId)->GetPointId(1),0),
                       poly->GetAttributes(0)->GetArray(0)->GetComponent(poly->GetCell(cellId)->GetPointId(1),1),
                       poly->GetAttributes(0)->GetArray(0)->GetComponent(poly->GetCell(cellId)->GetPointId(1),2)};
    double rgb_c[3] = {poly->GetAttributes(0)->GetArray(0)->GetComponent(poly->GetCell(cellId)->GetPointId(2),0),
                       poly->GetAttributes(0)->GetArray(0)->GetComponent(poly->GetCell(cellId)->GetPointId(2),1),
                       poly->GetAttributes(0)->GetArray(0)->GetComponent(poly->GetCell(cellId)->GetPointId(2),2)};

    vtkMath::MultiplyScalar(rgb_a, a_w);
    vtkMath::MultiplyScalar(rgb_b, b_w);
    vtkMath::MultiplyScalar(rgb_c, c_w);
    vtkMath::Add(rgb_a, rgb_b, rgb);
    vtkMath::Add(rgb, rgb_c, rgb);
}

void drawResponse(const std::vector<cv::linemod::Template>& templates,
    int num_modalities, cv::Mat& dst, cv::Point offset, int T)
{
    static const cv::Scalar COLORS[5] = { CV_RGB(0, 0, 255),
                                          CV_RGB(0, 255, 0),
                                          CV_RGB(255, 255, 0),
                                          CV_RGB(255, 140, 0),
                                          CV_RGB(255, 0, 0) };

    for (int m = 0; m < num_modalities; ++m)
    {
        // NOTE: Original demo recalculated max response for each feature in the TxT
        // box around it and chose the display color based on that response. Here
        // the display color just depends on the modality.
        cv::Scalar color = COLORS[m];

        for (int i = 0; i < (int)templates[m].features.size(); ++i)
        {
            cv::linemod::Feature f = templates[m].features[i];
            cv::Point pt(f.x + offset.x, f.y + offset.y);
            cv::circle(dst, pt, T / 2, color);
        }
    }
}

static void writeLinemod(const cv::Ptr<cv::linemod::Detector>& detector, const std::string& filename)
{
    cv::FileStorage fs(filename, cv::FileStorage::WRITE);
    detector->write(fs);

    std::vector<cv::String> ids = detector->classIds();
    fs << "classes" << "[";
    for (int i = 0; i < (int)ids.size(); ++i)
    {
        fs << "{";
        detector->writeClass(ids[i], fs);
        fs << "}"; // current class
    }
    fs << "]"; // classes
}

#endif
