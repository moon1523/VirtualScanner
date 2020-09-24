#ifndef FUNCTIONS
#define FUNCTIONS

#include <iostream>
#include <pcl/register_point_struct.h>
#include <pcl/io/pcd_io.h>
#include <pcl/io/vtk_lib_io.h>
#include <pcl/filters/voxel_grid.h>
#include <pcl/point_types.h>
#include <pcl/console/parse.h>
#include <pcl/visualization/vtk.h>
#include <opencv2/core/utility.hpp>
#include <opencv2/imgproc/imgproc_c.h> // cvFindContours
#include <opencv2/imgproc.hpp>
#include <opencv2/rgbd.hpp>
#include <opencv2/videoio.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/core.hpp>

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

#endif
