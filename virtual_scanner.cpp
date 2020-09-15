/*
 * Software License Agreement (BSD License)
 *
 *  Copyright (c) 2010, Willow Garage, Inc.
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions
 *  are met:
 *
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above
 *     copyright notice, this list of conditions and the following
 *     disclaimer in the documentation and/or other materials provided
 *     with the distribution.
 *   * Neither the name of the copyright holder(s) nor the names of its
 *     contributors may be used to endorse or promote products derived
 *     from this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 *  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 *  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 *  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 *  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 *  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 *  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 *  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *
 * $Id$
 *
 */

/**
  * \author Radu Bogdan Rusu
  *
  * @b virtual_scanner takes in a .ply or a .vtk file of an object model, and virtually scans it
  * in a raytracing fashion, saving the end results as PCD (Point Cloud Data) files. In addition,
  * it noisifies the PCD models, and downsamples them.
  * The viewpoint can be set to 1 or multiple views on a sphere.
  */
#include <string>
#include <pcl/register_point_struct.h>
#include <pcl/io/pcd_io.h>
#include <pcl/io/vtk_lib_io.h>
#include <pcl/filters/voxel_grid.h>
#include <pcl/point_types.h>
#include <pcl/console/parse.h>
#include <pcl/visualization/vtk.h>
#include "boost.h"
#include "vtkGeneralTransform.h"
#include "vtkTransform.h"
#include "vtkCellLocator.h"

using namespace pcl;

#define EPS 0.00001

//struct ScanParameters
//{
//    int nr_scans;             // number of steps for sweep movement
//    int nr_points_in_scans;   // number of laser beam measurements per scan
//    double vert_res;          // vertical resolution (res. of sweep) in degrees
//    double hor_res;           // horizontal  resolution (of laser beam) in degrees
//    double max_dist;          // maximum distance in units.
//};

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

using namespace std;
void PrintUsage(){
    cout<<"Usage: ./scanner [options] model.ply viewpoints.ply"<<endl;
    cout<<"[options]"<<endl;
    cout<<"\t-dist       <double>  : scan distance in mm (default:3000)"<<endl;
    cout<<"                          or put three numbers for multiple distances [start] [interval] [#]"<<endl;
    cout<<"\t-rot        <int>     : rotation interval in degree (default:10)"<<endl;
    cout<<"\t-path       <string>  : directory path for result files (default: .)"<<endl;
    exit(1);
}

void EstimateColor(vtkPolyData* data, vtkIdType cellId, double* point, double* rgb);

int main (int argc, char** argv)
{
    //arguements (file names)
    string modelFileName = argv[argc-2];
    string viewFileName  = argv[argc-1];

    //set scan parameters
    ScanParameters scan;
    scan.fov_hor  = 90.*pi/180.;
    scan.fov_vert = 59.*pi/180.;
    scan.res_hor  = 1280;
    scan.res_vert  = 720;
    scan.distance = 3000; //mm
    scan.rotInterval = 10;
    scan.screen_dist = 5000; //mm
    vector<double> scandists;

    //arguments2
    string path(".");
    for(int i=1;i<argc-2;i++){
        if(string(argv[i])=="-dist"){
            stringstream ss(argv[++i]);
            vector<double> inputs;
            double input;
            while(ss>>input) inputs.push_back(input);
            if(inputs.size()==1) scan.distance = inputs[0];
            else
                for(int n=0;n<inputs[2];n++)
                    scandists.push_back(inputs[0]+inputs[1]*n);
        }
        else if(string(argv[i])=="-rot")  scan.rotInterval = atoi(argv[++i]);
        else if(string(argv[i])=="-path") path = argv[++i];

        else PrintUsage();
    }
    if(scandists.size()==0) scandists.push_back(scan.distance);

    //read ply file
    vtkSmartPointer<vtkPolyData> data;
    vtkPLYReader* reader = vtkPLYReader::New ();
    reader->SetFileName (modelFileName.c_str ());
    reader->Update ();
    data = reader->GetOutput ();

    // Build a spatial locator for our dataset
    vtkSmartPointer<vtkCellLocator> tree = vtkSmartPointer<vtkCellLocator>::New ();
    tree->SetDataSet (data);
    tree->CacheCellBoundsOn ();
    tree->SetTolerance (0.0);
    tree->SetNumberOfCellsPerBucket (1);
    tree->AutomaticOn ();
    tree->BuildLocator ();
    tree->Update ();

    //read viewpoint ply file
    vtkSmartPointer<vtkPolyData> viewPointData;
    vtkPLYReader* viewPointReader = vtkPLYReader::New ();
    viewPointReader->SetFileName (viewFileName.c_str ());
    viewPointReader->Update ();
    viewPointData = viewPointReader->GetOutput ();

    //get viewpoints
    vector<double*> viewRays;
    for(vtkIdType i=0;i<viewPointData->GetNumberOfCells();i++){
        double a[3], b[3], c[3];
        viewPointData->GetPoint(viewPointData->GetCell(i)->GetPointId(0), a);
        viewPointData->GetPoint(viewPointData->GetCell(i)->GetPointId(1), b);
        viewPointData->GetPoint(viewPointData->GetCell(i)->GetPointId(2), c);

        double ab[3], ac[3];
        double* normal = new double[3];
        vtkMath::Subtract(a,b,ab);
        vtkMath::Subtract(a,c,ac);
        vtkMath::Cross(ab, ac, normal);
        vtkMath::Normalize(normal);
        viewRays.push_back(normal);
    }

    // Virtual camera parameters
    //double eye0[3]     = {0.0, 0.0, 0.0};
    double viewray0[3] = {0.0, 0.0, -1.0};
    double up0[3]      = {0.0, 1.0, 0.0}; //normalize

    // Screen parameters
    double hor_len  = (scan.distance+scan.screen_dist) * 2. * tan(scan.fov_hor*0.5);
    double vert_len = (scan.distance+scan.screen_dist) * 2. * tan(scan.fov_vert*0.5);
    double hor_interval  = hor_len  / (double) scan.res_hor;
    double vert_interval = vert_len / (double) scan.res_vert;
    cout<<hor_len<<"*"<<vert_len<<endl<<hor_interval<<"*"<<vert_interval<<endl;

    //main loop start
    if(path!=".") system(("mkdir "+path).c_str());
    ofstream log(path+"/labels.txt");
    vtkTransform* tr1 = vtkTransform::New ();
    vtkTransform* tr = vtkTransform::New ();
    vtkMatrix4x4* mat = vtkMatrix4x4::New();
    vtkMatrix4x4* mat1 = vtkMatrix4x4::New();
    for(int i=0;i<viewRays.size();i++){ //viewpoint loop

        log<<"viewpoint"<<i<<" "<<viewRays[i][0]<<" "<<viewRays[i][1]<<" "<<viewRays[i][2]<<endl;
        //set viewray and upVector
        double angle = vtkMath::AngleBetweenVectors(viewRays[i], viewray0); //radian
        double axis[3], up1[3];
        vtkMath::Cross(viewray0,viewRays[i], axis);
        tr1->Identity ();
        tr1->RotateWXYZ(angle*180/pi, axis); //requires angle in degree
        tr1->TransformPoint(up0, up1);
        tr1->Inverse();
        tr1->GetMatrix(mat1);
        vtkMath::Normalize(up1);

        for(int deg = 0; deg<360; deg+=scan.rotInterval){ //rotation loop

            //set up/right vectors
            double up[3], right[3];
            tr->Identity ();
            tr->RotateWXYZ((double)deg, viewRays[i]);
            tr->TransformPoint(up1, up);
            tr->Inverse();
            tr->GetMatrix(mat);
            vtkMatrix4x4::Multiply4x4(mat1,mat, mat);
            tr->Identity();
            tr->SetMatrix(mat);

            vtkMath::Normalize(up);
            vtkMath::Cross(viewRays[i], up, right);
            vtkMath::Normalize(right);

            log<<to_string(i)+" "+to_string(deg)<<" "<<"upVector: "<<up[0]<<" "<<up[1]<<" "<<up[2]<<endl;
            log<<*mat<<endl;

            for(double dist:scandists){ //distance loop
                cout<<"\rStart Scanning... -view: "<<i<<"/"<<viewRays.size()<<" -deg: "<<deg<<" -dist: "<<dist<<flush;
                scan.distance = dist;
                //Start!!
                double eye[3];
                eye[0] = -viewRays[i][0]*scan.distance;
                eye[1] = -viewRays[i][1]*scan.distance;
                eye[2] = -viewRays[i][2]*scan.distance;

                // Generate screen
                double screen_center[3];
                for(int n=0;n<3;n++) screen_center[n] = viewRays[i][n] * scan.screen_dist;
                double screen_leftup[3];
                for(int n=0;n<3;n++) screen_leftup[n] = screen_center[n]+up[n]*vert_len*0.5 - right[n]*hor_len*0.5;

                pcl::PointCloud<pcl::PointXYZ> screenCloud;
                for(int vert = 0 ; vert<scan.res_vert;vert++){
                    double left[3] = {screen_leftup[0] - up[0]*vert_interval*(vert+0.5),
                                      screen_leftup[1] - up[1]*vert_interval*(vert+0.5),
                                      screen_leftup[2] - up[2]*vert_interval*(vert+0.5)};
                    for(int hor = 0 ; hor<scan.res_hor;hor++){
                        double point[3] = {left[0]+right[0]*hor_interval*(hor+0.5),
                                           left[1]+right[1]*hor_interval*(hor+0.5),
                                           left[2]+right[2]*hor_interval*(hor+0.5)};
                        screenCloud.push_back(pcl::PointXYZ(point[0], point[1], point[2]));
                    }
                }

                // Scanning
                pcl::PointCloud<pcl::PointXYZ> cloud;
                vector<double> xx, yy, zz;
                vector<int> rr, gg, bb;
                double p_coords[3], x[3], t, rgb[3];
                int subId;
                for(int vert = 0 ; vert<scan.res_vert;vert++){
                    for(int hor = 0 ; hor<scan.res_hor;hor++){
                        vtkIdType cellId;
                        double point[3];
                        point[0] = screenCloud.at(vert*scan.res_hor+hor).x;
                        point[1] = screenCloud.at(vert*scan.res_hor+hor).y;
                        point[2] = screenCloud.at(vert*scan.res_hor+hor).z;
                        if (tree->IntersectWithLine (eye, point, 0, t, x, p_coords, subId, cellId))
                        {
                            EstimateColor(data, cellId, x, rgb);
                            rr.push_back(floor(rgb[0]+0.5));gg.push_back(floor(rgb[1]+0.5));bb.push_back(floor(rgb[2]+0.5));
                            xx.push_back(x[0]);yy.push_back(x[1]);zz.push_back(x[2]);

                            //vtkMath::Subtract(x,data->GetCenter(),x);
                            tr->TransformPoint(x,x);
                            cloud.push_back (PointXYZ(x[0], x[1], x[2]));
                        }
                        //else cloud.push_back (PointXYZ(0,0,0));
                        else{
                            rr.push_back(0);gg.push_back(0);bb.push_back(0);
                            xx.push_back(point[0]);yy.push_back(point[1]);zz.push_back(point[2]);
                            //vtkMath::Subtract(point,data->GetCenter(),point);
                            tr->TransformPoint(point,point);
                            cloud.push_back (PointXYZ(point[0],point[1],point[2]));
                        }
                    } // Horizontal
                } // Vertical

                pcl::PCDWriter writer;
                writer.writeBinaryCompressed(path+"/"+to_string(i)+"_r"+to_string(deg)+"_d"+to_string((int)scan.distance)+".pcd", cloud);
                ofstream ofs(path+"/"+to_string(i)+"_r"+to_string(deg)+"_d"+to_string((int)scan.distance)+".ply");
                ofs<<"ply"<<endl;
                ofs<<"format ascii 1.0"<<endl;
                ofs<<"comment exported in vitual_scanner"<<endl;
                ofs<<"element vertex "<<xx.size()<<endl;
                ofs<<"property float x"<<endl;
                ofs<<"property float y"<<endl;
                ofs<<"property float z"<<endl;
                ofs<<"property uchar red"<<endl;
                ofs<<"property uchar green"<<endl;
                ofs<<"property uchar blue"<<endl;
                ofs<<"end_header"<<endl;
                for(int idx=0;idx<xx.size();idx++)
                    ofs<<xx[idx]<<" "<<yy[idx]<<" "<<zz[idx]<<" "<<rr[idx]<<" "<<gg[idx]<<" "<<bb[idx]<<endl;
                ofs.close();
            }
        }
    }log.close();
    return 0;
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
